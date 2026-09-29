"""Query (phase 3 prototype): naive hit counts and containment per unit.

Reads are streamed through the Rust kernel at the index's ``t_max``; every sampled hash
is looked up in tier 2 and counts for each unit in its set whose own threshold it passes
(``hash <= max_hash_g``), which also discards most fingerprint false hits. No query
sketch is stored: only per-(unit, hash) counts are kept. Indexes imported from sourmash
signatures (``meta["hash"] == "sourmash"``) hash reads with sourmash instead; ``frames``
and ``genetic_code`` then do not apply.

Shared k-mers count for every unit that holds them in ``kmers_hit``; ``kmers_unique`` is
what a gather-style greedy assignment leaves each unit (the phase-4 detection baseline).
``coverage_em`` re-splits the hits of the units gather keeps by EM (phase-4 quantification);
``coverage_zi`` fits coverage to the k-mers present only (zero-inflated EM).
Two more baselines give each hit k-mer to one unit (:func:`assign_best`): winner-take-all
(``kmers_wta``) and uniqueness-first (``kmers_ufirst``), each with the coverage its k-mers'
hits give.
"""

import heapq
from collections.abc import Iterable
from pathlib import Path
from typing import Final

import numpy as np
import polars as pl
from scipy.sparse import coo_array, csr_array
from scipy.sparse.csgraph import connected_components

from kmer_functional_profiler import _core
from kmer_functional_profiler.compat import sourmash_hits
from kmer_functional_profiler.index import PIN_BITS, Index, IndexParams, PackedTable


def unit_hits(
    table: PackedTable, max_hash_g: np.ndarray, hashes: np.ndarray, reads: np.ndarray
) -> pl.DataFrame:
    """Expand sampled query hashes to (unit, hash, read, pin_q, holders) hits through ``table``.

    A hit counts for a unit only if the hash passes that unit's ``max_hash_g``; ``pin_q`` is
    the k-mer's quantised ``p_in`` in that unit and ``holders`` the number of index units
    it counts for.
    """
    set_ids = table.lookup(hashes)
    found = set_ids >= 0
    hashes, reads, set_ids = hashes[found], reads[found], set_ids[found]
    starts = table.set_offsets[set_ids].astype(np.int64)
    lengths = table.set_offsets[set_ids + 1].astype(np.int64) - starts
    # Positions of every set member, hit by hit.
    within = np.arange(lengths.sum()) - np.repeat(np.cumsum(lengths) - lengths, lengths)
    values = table.set_values[np.repeat(starts, lengths) + within]
    units = (values >> np.uint64(PIN_BITS)).astype(np.uint32)
    pin_q = (values & np.uint64(2**PIN_BITS - 1)).astype(np.uint8)
    hit = np.repeat(np.arange(len(hashes)), lengths)
    keep = hashes[hit] <= max_hash_g[units]
    holders = np.bincount(hit[keep], minlength=len(hashes)).astype(np.uint32)
    hit = hit[keep]
    return pl.DataFrame(
        {
            "unit": units[keep],
            "hash": hashes[hit],
            "read": reads[hit],
            "pin_q": pin_q[keep],
            "holders": holders[hit],
        }
    )


def gather(kmers: pl.DataFrame, t_g: np.ndarray) -> pl.DataFrame:
    """Assign each hit k-mer to one unit, greedily, as ``sourmash gather`` does.

    Repeatedly take the unit with the most unassigned hit k-mers, scaled by ``1 / t_g`` so
    units sampled at different rates compare as estimated k-mer counts, and give it those
    k-mers. Units left with none are explained away and get no row. Returns ``unit``,
    ``kmers_unique`` (k-mers assigned) and ``gather_rank`` (0 = taken first).
    """
    remaining = {u: set(h) for u, h in kmers.group_by("unit").agg("hash").iter_rows()}
    # Scores only fall as k-mers are taken, so a lazy heap needs only stale-top rechecks.
    heap = [(-len(h) / t_g[u], u) for u, h in remaining.items()]
    heapq.heapify(heap)
    taken: set[int] = set()
    rows: list[tuple[int, int, int]] = []
    while heap:
        _, unit = heapq.heappop(heap)
        mine = remaining[unit] = remaining[unit] - taken  # iterates the unit's set, not taken
        if not mine:
            continue
        score = len(mine) / t_g[unit]
        if heap and score < -heap[0][0]:
            heapq.heappush(heap, (-score, unit))
            continue
        taken |= mine
        rows.append((unit, len(mine), len(rows)))
    return pl.DataFrame(
        rows,
        schema={"unit": pl.UInt32, "kmers_unique": pl.UInt32, "gather_rank": pl.UInt32},
        orient="row",
    )


def assign_best(kmers: pl.DataFrame, score: pl.Expr) -> pl.DataFrame:
    """Give each hit k-mer to the unit holding it with the highest ``score`` (ties: lowest).

    ``score`` is evaluated per row of ``kmers`` (``unit``, ``hash`` and any other columns)
    and must be constant within a unit. Static unit scores make this one pass, unlike
    :func:`gather`, whose scores fall as k-mers are taken. The rules it serves:

    - winner-take-all, as sylph: score = containment (k-mers hit / ``m_g``);
    - uniqueness-first: score = sum of 1 / (units holding the k-mer) over the unit's hit
      k-mers, divided by ``t_g`` so units sampled at different rates compare fairly.

    Returns the rows of ``kmers`` that won.
    """
    return (
        kmers.with_columns(_score=score)
        .sort(["hash", "_score", "unit"], descending=[False, True, False])
        .unique("hash", keep="first", maintain_order=True)
        .drop("_score")
    )


WTA_SCORE: Final = pl.col("hash").n_unique().over("unit") / pl.col("m_g")
UFIRST_SCORE: Final = (1 / pl.len().over("hash")).sum().over("unit") / pl.col("t_g")


EXPLAINED_AWAY: Final = 1e-3  # expected hits below which EM reports coverage 0


def _div(num: np.ndarray, den: np.ndarray) -> np.ndarray:
    """num / den, 0 where den is 0 (units or k-mers whose weight has vanished)."""
    out = np.zeros_like(num, dtype=np.float64)
    np.divide(num, den, out=out, where=den > 0)
    return out


def em(
    kmers: pl.DataFrame,
    m_g: np.ndarray,
    *,
    zero_inflated: bool = False,
    prior: tuple[float, float] | None = None,
    tol: float = 1e-6,
    max_iter: int = 1000,
) -> pl.DataFrame:
    """Per-unit k-mer ``coverage`` (and ``present`` fraction) by EM over k-mer hit counts.

    ``kmers`` has one row per (``unit``, ``hash``) with the k-mer's ``hits``. Each k-mer's
    hits are Poisson with mean the sum of ``coverage`` over the units holding it, and each
    unit's ``m_g`` kept k-mers (hit or not) all count in its expectation; EM finds the
    maximum-likelihood coverages. Components share no k-mers, so their updates are
    independent even though they run in one sparse product.

    ``zero_inflated`` lets only a fraction ``present`` of a unit's kept k-mers occur in the
    sample (the rest are structural zeros: the sample's strain differs there), so coverage
    is fitted to the k-mers present. Hits are shared in proportion to coverage x present,
    and ``present`` is the unit's share of hit k-mers over the kept k-mers expected to be
    hit at that coverage, capped at 1. Without shared k-mers the fixed point is the
    zero-truncated Poisson MLE sylph uses: mean hits per hit k-mer = c / (1 - exp(-c)).
    With no excess zeros ``present`` stays 1 and the result equals plain EM. Units left
    with fewer than ``EXPLAINED_AWAY`` expected hits (their k-mers explained by other units)
    get coverage 0.

    ``prior`` (a, b >= 1) puts a Beta prior on ``present`` (see :func:`fit_present_prior`),
    which pulls it towards the prior mean when few k-mers could be hit, i.e. at low coverage
    or small ``m_g``; ``present`` is then updated by EM over the unhit k-mers' presence.
    """
    units, col = np.unique(kmers["unit"].to_numpy(), return_inverse=True)
    hashes, row = np.unique(kmers["hash"].to_numpy(), return_inverse=True)
    hits = np.zeros(len(hashes))
    hits[row] = kmers["hits"].to_numpy()
    a = csr_array((np.ones(len(row)), (row, col)), shape=(len(hashes), len(units)))
    m = m_g[units].astype(np.float64)
    lam, pi = (a.T @ hits) / m, np.ones(len(units))
    for _ in range(max_iter):
        w = lam * pi
        mu = a @ w
        attributed = w * (a.T @ _div(hits, mu))  # expected hits from each unit
        if zero_inflated:
            kmers_hit = w * (a.T @ _div(np.ones_like(mu), mu))  # expected hit k-mers per unit
            seen = -np.expm1(-lam)  # chance a present k-mer is hit
            if prior is None:
                new_pi = np.minimum(1.0, _div(kmers_hit, m * seen))
            else:
                # Unhit k-mers are present with odds pi (1 - seen) : (1 - pi).
                odds = pi * (1 - seen) / np.maximum(1 - pi * seen, 1e-300)
                unhit_present = np.maximum(m - kmers_hit, 0) * odds
                new_pi = (kmers_hit + unhit_present + prior[0] - 1) / (m + sum(prior) - 2)
        else:
            new_pi = pi
        new = _div(attributed, m * new_pi)
        done = np.abs(new - lam).max(initial=0) <= tol * new.max(initial=0) and np.allclose(
            new_pi, pi, rtol=0, atol=tol
        )
        lam, pi = new, new_pi
        if done:
            break
    # Units explained away by others converge towards 0 without reaching it.
    lam[attributed < EXPLAINED_AWAY] = 0.0
    return pl.DataFrame(
        {"unit": units, "coverage": lam, "present": pi},
        schema={"unit": pl.UInt32, "coverage": pl.Float64, "present": pl.Float64},
    )


# p_in of each quantised level (pin_q = round(15 p_in)), kept inside (0, 1): a k-mer certain to
# be present (p_in = 1 at s_g = 1) would stay so forever, since EM then reads every unhit copy
# of it as present; one certain to be absent could never be hit.
PIN_P: Final = np.clip(
    np.arange(2**PIN_BITS) / (2**PIN_BITS - 1), 0.5 / (2**PIN_BITS - 1), 1 - 0.5 / (2**PIN_BITS - 1)
)


def em_pin(
    kmers: pl.DataFrame, pin_hist: np.ndarray, *, tol: float = 1e-6, max_iter: int = 1000
) -> pl.DataFrame:
    """Zero-inflated EM where each k-mer's presence follows its ``p_in``.

    k-mer x of unit g is present in the sample with probability s_g * p_in(x): a strain
    like a random member carries each k-mer with probability p_in (s_g = 1), and a more
    divergent one loses k-mers uniformly (s_g < 1). So core k-mers are expected present
    more often than private ones, which changes both how shared hits are split (by
    coverage x presence of that k-mer in each unit) and how much an unhit k-mer argues for
    low coverage rather than absence. ``p_in`` is read from its quantised level (``PIN_P``);
    s_g is capped so presence stays <= 1.

    ``kmers`` has ``unit``, ``hash``, ``hits`` and ``pin_q``; ``pin_hist[unit]`` counts the
    unit's kept k-mers per level (its sum is ``m_g``). Coverage is attributed hits over
    expected present k-mers (hit ones by their share, unhit ones by their posterior
    presence), and s_g is expected present k-mers over the sum of p_in. Like :func:`em`, a
    hit k-mer counts as present for a unit by its share of the hits. Returns ``unit``,
    ``coverage`` and ``present`` (expected fraction of kept k-mers present); units
    explained away get coverage 0.
    """
    units, col = np.unique(kmers["unit"].to_numpy(), return_inverse=True)
    hashes, row = np.unique(kmers["hash"].to_numpy(), return_inverse=True)
    level = kmers["pin_q"].to_numpy().astype(np.intp)
    hits = np.zeros(len(hashes))
    hits[row] = kmers["hits"].to_numpy()
    n = len(units)
    hist = pin_hist[units].astype(np.float64)
    m, expected = hist.sum(axis=1), hist @ PIN_P  # kept k-mers; present ones at s_g = 1
    lam, scale = np.bincount(col, weights=hits[row], minlength=n) / m, np.ones(n)
    for _ in range(max_iter):
        weight = lam[col] * scale[col] * PIN_P[level]
        share = _div(weight, np.bincount(row, weights=weight, minlength=len(hashes))[row])
        attributed = np.bincount(col, weights=share * hits[row], minlength=n)
        hit = np.zeros_like(hist)
        np.add.at(hit, (col, level), share)
        seen = -np.expm1(-lam)[:, None]  # chance a present k-mer is hit
        pi = scale[:, None] * PIN_P
        odds = pi * (1 - seen) / np.maximum(1 - pi * seen, 1e-300)
        present = hit.sum(axis=1) + (np.maximum(hist - hit, 0) * odds).sum(axis=1)
        new = _div(attributed, present)
        new_scale = np.clip(_div(present, expected), 0.0, 1 / PIN_P[-1])
        done = np.abs(new - lam).max(initial=0) <= tol * new.max(initial=0) and np.allclose(
            new_scale, scale, rtol=0, atol=tol
        )
        lam, scale = new, new_scale
        if done:
            break
    lam[attributed < EXPLAINED_AWAY] = 0.0
    return pl.DataFrame(
        {"unit": units, "coverage": lam, "present": scale * expected / m},
        schema={"unit": pl.UInt32, "coverage": pl.Float64, "present": pl.Float64},
    )


def fit_present_prior(
    fit: pl.DataFrame, min_coverage: float = 3.0, min_units: int = 10
) -> tuple[float, float] | None:
    """Beta (a, b) matching the mean and variance of ``present`` over well-covered units.

    ``fit`` is a zero-inflated :func:`em` result; at ``coverage >= min_coverage`` nearly all
    present k-mers are hit, so ``present`` is measured, not guessed. None if fewer than
    ``min_units`` qualify. a, b are floored at 1 (a unimodal prior).
    """
    well = fit.filter(pl.col("coverage") >= min_coverage)["present"]
    if well.len() < min_units:
        return None
    mean, var = float(well.mean()), float(well.var())  # type: ignore[arg-type]
    strength = mean * (1 - mean) / max(var, 1e-9) - 1
    return max(mean * strength, 1.0), max((1 - mean) * strength, 1.0)


# Background detections per unit per read (pair) per unit of sampling rate t_g, and the
# geometric ratio of own k-mers per background detection (a homolog gives a few at once);
# calibrated on the fmh benchmark (step 18): P(true | own k-mers) within 0.05 on all indexes.
BACKGROUND: Final = 1e-6
CLUMP: Final = 0.3
H_CAP: Final = 20  # hit counts pooled at and above this in the present-unit distribution


def presence(
    own: pl.DataFrame,
    t_g: np.ndarray,
    n_reads: int,
    n_index_units: int,
    *,
    background: float = BACKGROUND,
    clump: float = CLUMP,
    max_iter: int = 500,
    tol: float = 1e-9,
) -> pl.DataFrame:
    """Probability that each unit gather keeps is present rather than hit by background.

    ``own`` has one row per k-mer gather gave a unit (``unit``, ``hash``, ``holders``), so
    hits other units explain are not evidence. An absent unit is still detected by
    background (off-target homologs, error k-mers) with probability 1 - exp(-mu_g), mu_g =
    ``background`` x ``n_reads`` x ``t_g``, then with h own k-mers ~ Geometric (ratio
    ``clump``), each hit k-mer's odds scaled by its index ``holders`` (k-mers many units
    share are conserved motifs, which unindexed genes carry too). A present unit has h
    drawn from a distribution f fitted to the sample.
    With A absent units (``n_index_units`` minus the expected present), P(present | h) =
    w_h / (w_h + A x P(h | absent)), w_h the expected present units with h own k-mers; w
    and A are fitted by fixed-point iteration, at which w_h is about the units seen with h
    own k-mers less the background expected to give h.

    The background rate is calibrated, not fitted: fitted per sample it is not
    identifiable (a free f, or a free per-unit hit rate, absorbs any unit with few hits,
    so fits collapse to no background or to everything being background). Scaling by
    t_g rather than m_g matches the fmh benchmark, where false positives are small,
    often floored, KOs.
    """
    per_unit = (
        own.group_by("unit")
        .agg(h=pl.len(), log_holders=pl.col("holders").cast(pl.Float64).log().sum())
        .sort("unit")
    )
    units = per_unit["unit"].to_numpy()
    h = per_unit["h"].to_numpy().astype(np.float64)
    mu = background * n_reads * t_g[units]
    capped = np.minimum(h, H_CAP).astype(np.intp)
    # P(own hits | absent); units at the cap are never taken for background.
    absent = np.where(
        capped < H_CAP,
        -np.expm1(-mu)
        * (1 - clump)
        * clump ** (h - 1)
        * np.exp(per_unit["log_holders"].to_numpy()),
        0.0,
    )
    prob = np.ones(len(units))
    for _ in range(max_iter):
        w = np.bincount(capped, weights=prob, minlength=H_CAP + 1)
        new = w[capped] / (w[capped] + (n_index_units - prob.sum()) * absent)
        done = np.abs(new - prob).max(initial=0) <= tol
        prob = new
        if done:
            break
    return pl.DataFrame(
        {"unit": units, "present_prob": prob},
        schema={"unit": pl.UInt32, "present_prob": pl.Float64},
    )


# Floor of the copies error (sd on the log scale), calibrated on the fmh benchmark (step 18):
# brings abundance_zi intervals to 95% at kfp_s100.
COPIES_ERROR: Final = 0.2
MH_STEPS: Final = 5  # Metropolis steps on coverage per sweep
PI_GRID: Final = (np.arange(512) + 0.5) / 512  # present-fraction grid for exact draws


def _ztp_log_post(lam: np.ndarray, hits: np.ndarray, hit: np.ndarray) -> np.ndarray:
    """Log posterior of coverage from ``hits`` over ``hit`` k-mers, zero-truncated Poisson."""
    out: np.ndarray = hits * np.log(lam) - hit * (lam + np.log(-np.expm1(-lam))) - 0.01 * lam
    return out


def posterior_zi(
    hit_reads: pl.DataFrame,
    m_g: np.ndarray,
    pin_sum: np.ndarray,
    draws: int,
    *,
    present_prob: np.ndarray | None = None,
    len_cv: np.ndarray | None = None,
    copies_error: float = COPIES_ERROR,
    sweeps: int = 10,
    level: float = 0.95,
    shared_evidence: float = 0.5,
    seed: int = 0,
) -> pl.DataFrame:
    """Intervals and ambiguity groups for the zero-inflated model: resampled reads + Gibbs.

    ``hit_reads`` has one row per (``unit``, ``hash``, ``read``) with ``n`` hits. Each of
    ``draws`` replicates

    1. reweights every read (pair) by a Poisson(1) draw and rebuilds the k-mer counts (one
       read hits many neighbouring k-mers, so counts are not independent), and refits
       zero-inflated :func:`em`: this carries the read-sampling uncertainty;
    2. runs ``sweeps`` Gibbs sweeps from that fit on those counts and keeps the last state,
       which adds the uncertainty of splitting shared k-mers between similar units. A sweep
       splits each hit k-mer's count among its holders multinomially, in proportion to
       coverage x present fraction (EM's shares, approximating the holders' joint presence);
       updates coverage from the counts on the k-mers each unit was given, a zero-truncated
       Poisson likelihood free of the present fraction (Metropolis on log coverage,
       Gamma(1, 0.01) prior); and draws the present fraction exactly on a grid from hit
       k-mers ~ Binomial(m, pi (1 - e^-lambda)), plus the present unhit k-mers for copies.

    A unit given no hits in a draw is absent in it (coverage and abundance 0), and so is one
    dropped with probability 1 - ``present_prob[unit]`` (:func:`presence`), so weak or
    background-like evidence widens intervals down to 0. Copies (present k-mers over
    ``pin_sum``) assume an average member; each draw scales them by exp(N(0, s^2)), s^2 =
    log(1 + ``len_cv[unit]``^2 / copies) + ``copies_error``^2, the spread of the mean
    kept k-mers of that many members plus a calibrated floor. Where units
    cannot be told apart, their shared counts move between them across draws, so their
    intervals widen. Ambiguity groups come from shared evidence: over all draws, a unit is
    linked to another when at least ``shared_evidence`` of the hits allocated to it lie on
    k-mers both hold, and linked units form a group. This catches lopsided pairs, e.g. a
    small unit living on hits it shares with a well-supported one, whose draws barely
    correlate. ``own_evidence`` is a unit's share of allocated hits on k-mers no other
    detected unit holds. A group's total is usually far better determined than its members.
    Returns per unit ``coverage_zi_lo``/``_hi``, ``abundance_zi_lo``/``_hi`` at ``level``,
    ``own_evidence``, and ``ambiguity_group`` (the group's smallest unit id; null when
    alone), ``group_size`` and the group totals ``group_coverage_zi_lo``/``_hi`` and
    ``group_abundance_zi_lo``/``_hi``.
    """
    rng = np.random.default_rng(seed)
    # Sorted, so the seeded draws do not depend on row order (Polars group_by does not fix it).
    keyed = (
        hit_reads.sort("unit", "hash", "read")
        .group_by("unit", "hash", maintain_order=True)
        .agg("read", "n")
    )
    units, col = np.unique(keyed["unit"].to_numpy(), return_inverse=True)
    hashes, row = np.unique(keyed["hash"].to_numpy(), return_inverse=True)
    n_units, n_rows = len(units), len(hashes)
    key_of_hit = np.repeat(np.arange(keyed.height), keyed["read"].list.len().to_numpy())
    reads, read_of = np.unique(
        keyed["read"].explode(empty_as_null=False).to_numpy(), return_inverse=True
    )
    n = keyed["n"].explode(empty_as_null=False).to_numpy().astype(np.int64)
    # Entries sorted by k-mer, with their position within it, for sequential binomial splits.
    order = np.argsort(row, kind="stable")
    row_s, col_s = row[order], col[order]
    starts = np.searchsorted(row_s, np.arange(n_rows))
    position = np.arange(len(row_s)) - starts[row_s]
    holders = np.bincount(row_s, minlength=n_rows)
    m = m_g[units].astype(np.float64)
    keep_prob = np.ones(n_units) if present_prob is None else present_prob[units]
    cv2 = np.zeros(n_units) if len_cv is None else len_cv[units] ** 2
    coverage = np.zeros((draws, n_units))
    abundance, group_coverage, group_abundance = (np.zeros_like(coverage) for _ in range(3))
    evidence = np.zeros(len(row_s))  # hits allocated per (k-mer, unit) entry, all draws
    for b in range(draws):
        weight = rng.poisson(1.0, len(reads))[read_of]
        per_key = np.rint(np.bincount(key_of_hit, weights=n * weight, minlength=keyed.height))
        table = keyed.select("unit", "hash").with_columns(hits=per_key).filter(pl.col("hits") > 0)
        if table.height == 0:
            continue
        fit = em(table, m_g, zero_inflated=True)
        at = np.searchsorted(units, fit["unit"].to_numpy())
        lam, pi = np.full(n_units, 1e-6), np.full(n_units, 1e-6)
        lam[at] = np.maximum(fit["coverage"].to_numpy(), 1e-6)
        pi[at] = np.clip(fit["present"].to_numpy(), 1e-6, 1 - 1e-6)
        # Holders of a k-mer see the same reads, so any holder's count is the k-mer's.
        count = np.zeros(n_rows, dtype=np.int64)
        count[row] = per_key.astype(np.int64)
        for _ in range(sweeps):
            rate = (lam * pi)[col_s]
            rest = np.bincount(row_s, weights=rate, minlength=n_rows)  # rate not yet visited
            left = count.copy()
            given = np.zeros(len(row_s), dtype=np.int64)
            for j in range(int(holders.max(initial=0))):
                here = position == j
                r = row_s[here]
                last = holders[r] == j + 1
                share = np.where(last, 1.0, _div(rate[here], rest[r]))
                given[here] = rng.binomial(left[r], np.clip(share, 0.0, 1.0))
                left[r] -= given[here]
                rest[r] -= rate[here]
            hits = np.bincount(col_s, weights=given, minlength=n_units)
            hit = np.bincount(col_s, weights=given > 0, minlength=n_units)  # k-mers given
            alive = hit > 0
            for _ in range(MH_STEPS):
                step = rng.normal(0.0, 1.5 / np.sqrt(hits + 1))
                new = lam * np.exp(step)
                gain = _ztp_log_post(new, hits, hit) - _ztp_log_post(lam, hits, hit) + step
                lam = np.where(alive & (np.log(rng.random(n_units)) < gain), new, lam)
            seen = -np.expm1(-lam)
            log_post = hit[:, None] * np.log(PI_GRID) + (m - hit)[:, None] * np.log1p(
                -PI_GRID * seen[:, None]
            )
            cdf = np.cumsum(np.exp(log_post - log_post.max(axis=1, keepdims=True)), axis=1)
            pick = (cdf < rng.random(n_units)[:, None] * cdf[:, -1:]).sum(axis=1)
            pi = PI_GRID[np.minimum(pick, len(PI_GRID) - 1)]
        unhit_odds = pi * (1 - seen) / (1 - pi * seen)
        present = hit + rng.binomial(np.maximum(m - hit, 0).astype(np.int64), unhit_odds)
        copies = present / pin_sum[units]
        sd = np.sqrt(np.log1p(cv2 / np.maximum(copies, 1.0)) + copies_error**2)
        scaled = lam * copies * np.exp(rng.normal(0.0, sd))
        # Group totals keep a dropped member's share: under absence its hits go to others.
        group_coverage[b] = np.where(alive, lam, 0.0)
        group_abundance[b] = np.where(alive, scaled, 0.0)
        alive &= rng.random(n_units) < keep_prob
        coverage[b] = np.where(alive, lam, 0.0)
        abundance[b] = np.where(alive, scaled, 0.0)
        evidence += given

    # Ambiguity groups: link a unit to another holder of its hit k-mers when most of its
    # allocated hits lie on k-mers the two share.
    total_evidence = np.bincount(col_s, weights=evidence, minlength=n_units)
    own = _div(np.bincount(col_s, weights=evidence * (holders[row_s] == 1), minlength=n_units),
               total_evidence)  # fmt: skip
    entries = pl.DataFrame({"row": row_s, "unit": col_s, "evidence": evidence})
    shared = (
        entries.join(entries.select("row", other="unit"), on="row")
        .filter(pl.col("unit") != pl.col("other"))
        .group_by("unit", "other")
        .agg(pl.col("evidence").sum())
    )
    unit_a, unit_b = shared["unit"].to_numpy(), shared["other"].to_numpy()
    # A unit never given hits (explained away in every draw) has total 0, so it links to
    # every co-holder of its hit k-mers.
    link = shared["evidence"].to_numpy() >= shared_evidence * total_evidence[unit_a]
    graph = coo_array(
        (np.ones(int(link.sum())), (unit_a[link], unit_b[link])),
        shape=(n_units, n_units),
    )
    _, group = connected_components(graph, directed=False)
    size = np.bincount(group)[group]
    group_id = np.full(n_units, len(units), dtype=np.int64)
    np.minimum.at(group_id, group, np.arange(n_units))  # smallest member index per group
    group_unit = units[group_id[group]]
    total = [np.stack([np.bincount(group, w, minlength=group.max() + 1) for w in x])[:, group]
             for x in (group_coverage, group_abundance)]  # fmt: skip
    total = [np.where(size > 1, t, x) for t, x in zip(total, (coverage, abundance), strict=True)]
    tails = [(1 - level) / 2, (1 + level) / 2]
    columns = {"unit": units.astype(np.uint32)}
    for name, x in (("coverage_zi", coverage), ("abundance_zi", abundance),
                    ("group_coverage_zi", total[0]), ("group_abundance_zi", total[1])):  # fmt: skip
        columns[f"{name}_lo"], columns[f"{name}_hi"] = np.quantile(x, tails, axis=0)
    return pl.DataFrame(columns).with_columns(
        ambiguity_group=pl.when(pl.Series(size) > 1).then(pl.Series(group_unit.astype(np.uint32))),
        group_size=pl.Series(size.astype(np.uint32)),
        own_evidence=pl.Series(own),
    )


def profile(
    index: Index,
    r1: str | Path,
    r2: str | Path | None = None,
    *,
    genetic_code: int = 11,
    frames: str = "stopfree",
    batch_reads: int = 100_000,
    draws: int = 0,
    kmers_out: str | Path | None = None,
) -> pl.DataFrame:
    """Per-unit hits, distinct k-mers hit, reads hit, containment and mean coverage.

    ``containment`` is the fraction of the unit's kept k-mers seen at least once;
    ``coverage`` is hits per kept k-mer; ``kmers_unique`` and ``gather_rank`` come from
    :func:`gather` (0 and null for units explained away) and ``coverage_em`` from :func:`em`
    over the units gather keeps (0 for the rest), with ``coverage_zi`` and ``present_zi``
    from its zero-inflated form, and ``coverage_zip``/``present_zip`` from :func:`em_pin`.
    ``copies_zi`` is present k-mers over an average member's kept k-mers (``pin_sum``), the
    member-equivalents present, and ``abundance_zi`` = ``coverage_zi`` x ``copies_zi``: for a
    unit whose members come from many genomes (a KO), total depth over its gene copies.
    ``kmers_wta``/``coverage_wta`` and ``kmers_ufirst``/``coverage_ufirst`` are the k-mers
    :func:`assign_best` gives each unit and their hits per kept k-mer. Units without hits
    are omitted.

    With a dense tier, the reads are streamed a second time at its rate and the EM
    estimates (``coverage_em``, ``_zi``, ``_zib``) are fitted on the dense hits of the units
    gather keeps, over their ``m_dense`` k-mers; ``kmers_dense`` counts those hit.

    ``present_prob`` (:func:`presence`) is the probability that a unit gather keeps is
    present rather than hit by background, from the tier-2 k-mers gather gave it; units
    gather drops get 0.

    ``draws`` > 0 adds 95% posterior intervals for ``coverage_zi`` and ``abundance_zi`` and
    ambiguity groups from that many Gibbs sweeps (:func:`posterior_zi`). ``kmers_out``
    writes the tier-2 hits per (``unit``, ``hash``) with ``hits`` and ``holders`` to Parquet.
    """
    params = IndexParams(**index.meta["params"])

    def stream(max_hash: int) -> Iterable[dict[str, np.ndarray]]:
        if index.meta.get("hash") == "sourmash":
            return sourmash_hits(r1, r2, params.k, max_hash, batch_reads)
        return _core.FastxHits(
            r1,
            r2,
            k=params.k,
            alphabet=params.alphabet,
            genetic_code=genetic_code,
            frames=frames,
            max_hash=max_hash,
            batch_reads=batch_reads,
        )

    max_hash_g = index.units["max_hash_g"].to_numpy()

    def by_read(hits: pl.DataFrame) -> pl.DataFrame:
        return hits.group_by("unit", "hash", "read").agg(
            n=pl.len(), pin_q=pl.col("pin_q").first(), holders=pl.col("holders").first()
        )

    def per_kmer(per_read: pl.DataFrame) -> pl.DataFrame:
        return per_read.group_by("unit", "hash").agg(
            hits=pl.col("n").sum(), pin_q=pl.col("pin_q").first(), holders=pl.col("holders").first()
        )

    empty = np.empty(0, dtype=np.uint64)
    batches, n_reads = [], 0
    for b in [{"hash": empty, "read": empty}, *stream(index.tier2.max_hash)]:
        batches.append(by_read(unit_hits(index.tier2, max_hash_g, b["hash"], b["read"])))
        if len(b["read"]):  # reads are numbered in input order; the last has sampled hashes
            n_reads = max(n_reads, int(b["read"].max()) + 1)
    per_read = pl.concat(batches)
    kmer_hits = per_kmer(per_read)
    if kmers_out is not None:
        kmer_hits.write_parquet(kmers_out)
    assigned = gather(kmer_hits.select("unit", "hash"), index.units["t_g"].to_numpy())
    detected_reads = per_read.join(assigned.select("unit"), on="unit", how="semi")
    detected = kmer_hits.join(assigned.select("unit"), on="unit")
    m_g = index.units["m_g"].to_numpy()
    pin_hist = index.units["pin_hist"].to_numpy()
    pin_sum = index.units["pin_sum"]
    len_cv = index.units["len_cv"].to_numpy()
    # Each hit k-mer is gather's: the first unit in gather order holding it.
    own = (
        kmer_hits.join(assigned.select("unit", "gather_rank"), on="unit")
        .sort("gather_rank")
        .unique("hash", keep="first")
        .select("unit", "hash", "holders")
    )
    present_prob = presence(own, index.units["t_g"].to_numpy(), n_reads, index.units.height)
    dense = index.dense
    if dense is not None:
        # Second pass: every k-mer at the dense rate, for the detected units only.
        max_hash_dense = index.units["max_hash_dense"].to_numpy()
        keep = assigned.select("unit")
        detected_reads = pl.concat(
            [
                by_read(
                    unit_hits(dense, max_hash_dense, b["hash"], b["read"]).join(
                        keep, on="unit", how="semi"
                    )
                )
                for b in [{"hash": empty, "read": empty}, *stream(dense.max_hash)]
            ]
        )
        detected = per_kmer(detected_reads)
        m_g = index.units["m_dense"].to_numpy()
        pin_hist = index.units["pin_hist_dense"].to_numpy()
        pin_sum = index.units["pin_sum_dense"]
        len_cv = index.units["len_cv_dense"].to_numpy()
    plain = em(detected, m_g).select("unit", coverage_em="coverage")
    inflated = em(detected, m_g, zero_inflated=True)
    prior = fit_present_prior(inflated)
    shrunk = (em(detected, m_g, zero_inflated=True, prior=prior) if prior else inflated).select(
        "unit", coverage_zib="coverage", present_zib="present"
    )
    # Present k-mers over an average member's kept k-mers: member-equivalents present.
    copies = pl.col("present") * pl.col("m") / pl.col("pin_sum")
    inflated = inflated.join(
        pl.DataFrame({"unit": np.arange(len(m_g), dtype=np.uint32), "m": m_g, "pin_sum": pin_sum}),
        on="unit",
    ).select(
        "unit",
        coverage_zi="coverage",
        present_zi="present",
        copies_zi=copies,
        abundance_zi=pl.col("coverage") * copies,
    )
    weighted = em_pin(detected, pin_hist).select(
        "unit", coverage_zip="coverage", present_zip="present"
    )
    rated = kmer_hits.join(index.units.select("unit", "m_g", "t_g"), on="unit")
    result = (
        per_read.group_by("unit")
        .agg(
            hits=pl.col("n").sum().cast(pl.UInt64),
            kmers_hit=pl.col("hash").n_unique().cast(pl.UInt32),
            reads=pl.col("read").n_unique().cast(pl.UInt64),
        )
        .join(index.units.select(pl.exclude("^(pin_(hist|sum)|len_cv).*$")), on="unit")
        .join(assigned, on="unit", how="left")
        .join(plain, on="unit", how="left")
        .join(inflated, on="unit", how="left")
        .join(shrunk, on="unit", how="left")
        .join(weighted, on="unit", how="left")
        .join(present_prob, on="unit", how="left")
    )
    if draws > 0:
        prob = np.zeros(len(m_g))
        prob[present_prob["unit"].to_numpy()] = present_prob["present_prob"].to_numpy()
        intervals = posterior_zi(
            detected_reads, m_g, pin_sum.to_numpy(), draws, present_prob=prob, len_cv=len_cv
        )
        result = result.join(intervals, on="unit", how="left")
    if dense is not None:
        result = result.join(
            detected.group_by("unit").agg(kmers_dense=pl.len().cast(pl.UInt32)),
            on="unit",
            how="left",
        ).with_columns(pl.col("kmers_dense").fill_null(0))
    for rule, score in (("wta", WTA_SCORE), ("ufirst", UFIRST_SCORE)):
        won = (
            assign_best(rated, score)
            .group_by("unit")
            .agg(
                pl.len().cast(pl.UInt32).alias(f"kmers_{rule}"),
                pl.col("hits").sum().alias(f"_{rule}"),
            )
        )
        result = result.join(won, on="unit", how="left")
    return (
        result.with_columns(
            pl.col("^(coverage|present|copies|abundance)_(em|zi|zib|zip)(_lo|_hi)?$").fill_null(
                0.0
            ),
            pl.col("present_prob").fill_null(0.0),
            containment=pl.col("kmers_hit") / pl.col("m_g"),
            coverage=pl.col("hits") / pl.col("m_g"),
            kmers_unique=pl.col("kmers_unique").fill_null(0),
            kmers_wta=pl.col("kmers_wta").fill_null(0),
            coverage_wta=pl.col("_wta").fill_null(0) / pl.col("m_g"),
            kmers_ufirst=pl.col("kmers_ufirst").fill_null(0),
            coverage_ufirst=pl.col("_ufirst").fill_null(0) / pl.col("m_g"),
        )
        .drop("_wta", "_ufirst")
        .sort("unit")
    )
