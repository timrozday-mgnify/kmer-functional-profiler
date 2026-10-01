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

import ctypes
import json
import os
import resource
import sys
import threading
import time
from collections.abc import Callable, Iterable, Iterator
from contextlib import contextmanager
from dataclasses import fields
from pathlib import Path
from typing import Final, NamedTuple

import numpy as np
import polars as pl
from scipy.sparse import coo_array
from scipy.sparse.csgraph import connected_components

from kmer_functional_profiler import _core
from kmer_functional_profiler.compat import sourmash_hits
from kmer_functional_profiler.index import PIN_BITS, Index, IndexParams, PackedTable

DISTINCT_SAMPLE: Final = 256  # distinct sampled k-mers are counted on 1 in this of hash space


def peak_rss() -> int:
    """Peak resident set size of this process so far, in bytes."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak if sys.platform == "darwin" else peak * 1024  # Linux reports KiB


def rss() -> int:
    """Current resident set size in bytes (Linux); the peak elsewhere."""
    try:
        return int(Path("/proc/self/statm").read_text().split()[1]) * resource.getpagesize()
    except OSError:
        return peak_rss()


def _footprint() -> int | None:
    """macOS ``phys_footprint`` (what Activity Monitor calls Memory): anonymous and
    compressed memory, without clean pages of memory-mapped files."""
    buf = ctypes.create_string_buffer(512)  # rusage_info_v2; phys_footprint at byte 72
    if _libc is None or _libc.proc_pid_rusage(os.getpid(), 2, buf) != 0:
        return None
    return int.from_bytes(buf.raw[72:80], "little")


_libc = ctypes.CDLL(None) if sys.platform == "darwin" else None


def anon_rss() -> int | None:
    """Current anonymous RSS in bytes (Linux ``RssAnon``; macOS ``phys_footprint``): memory
    the process owns, without the page cache of memory-mapped files (the index), which the
    kernel can reclaim. None elsewhere."""
    if _libc is not None:
        return _footprint()
    try:
        with open("/proc/self/status") as f:
            return next(int(line.split()[1]) * 1024 for line in f if line.startswith("RssAnon:"))
    except (OSError, StopIteration):
        return None


ANON_SAMPLE_S: Final = 0.05  # anonymous RSS sampling interval while a stage runs


class Timer:
    """Wall time, CPU time and peak RSS per query stage, and counts (the query cost study).

    A stage entered once per read batch accumulates its times. ``peak_rss`` is the
    process's peak at the stage's last exit, so the stage that raises it shows as a step.
    It counts resident pages of memory-mapped files too, so where the OS reports it (Linux,
    macOS) ``peak_anon`` is the same high-water mark for anonymous memory only, sampled every
    ``ANON_SAMPLE_S`` while a stage runs (the kernel keeps no anonymous peak).

    With ``path``, the stats are rewritten there after every stage, with ``running`` set to
    the stage in progress, so a process killed mid-query (out of memory) leaves the stages
    it finished and the one it died in. With ``log``, each stage's start and end are printed
    to stderr with the current and peak RSS.
    """

    def __init__(self, path: Path | None = None, *, log: bool = False) -> None:
        self.stages: dict[str, dict[str, float]] = {}
        self.counts: dict[str, int] = {}
        self.path, self.log, self.running = path, log, ""
        self.start = time.perf_counter()
        self.peak_anon = anon_rss()

    def _sample_anon(self, stop: threading.Event) -> None:
        while not stop.wait(ANON_SAMPLE_S):
            self._update_anon()

    def _update_anon(self) -> None:
        now = anon_rss()
        if now is not None and self.peak_anon is not None:
            self.peak_anon = max(self.peak_anon, now)

    def _log(self, event: str, stage: str) -> None:
        if self.log:
            gib = 2**30
            print(
                f"[{time.perf_counter() - self.start:9.1f}s] {event:5} {stage:<11}"
                f" rss {rss() / gib:7.2f} GiB  peak {peak_rss() / gib:7.2f} GiB"
                + (
                    "" if self.peak_anon is None else f"  anon peak {self.peak_anon / gib:7.2f} GiB"
                ),
                file=sys.stderr,
                flush=True,
            )

    def write(self) -> None:
        if self.path is not None:
            self.path.write_text(json.dumps(self.as_dict(), indent=2))

    @contextmanager
    def __call__(self, stage: str) -> Iterator[None]:
        outer, self.running = self.running, stage
        self._log("start", stage)
        self.write()
        stop = threading.Event()
        if self.peak_anon is not None and not outer:  # one sampler, for the outermost stage
            threading.Thread(target=self._sample_anon, args=(stop,), daemon=True).start()
        wall, cpu = time.perf_counter(), time.process_time()
        try:
            yield
        finally:
            stop.set()
            self._update_anon()
            s = self.stages.setdefault(stage, {"wall_s": 0.0, "cpu_s": 0.0})
            s["wall_s"] += time.perf_counter() - wall
            s["cpu_s"] += time.process_time() - cpu
            s["peak_rss"] = peak_rss()
            if self.peak_anon is not None:
                s["peak_anon"] = self.peak_anon
            self.running = outer
            self._log("end", stage)
            self.write()

    def as_dict(self) -> dict[str, object]:
        return {"running": self.running, "stages": self.stages, "counts": self.counts}


def _ids(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Distinct values of ``x`` and each element's index into them; the index is int32 when
    it fits (half the memory of ``np.unique``'s int64), which changes no result."""
    distinct, inverse = np.unique(x, return_inverse=True)
    return distinct, inverse.astype(np.int32 if len(x) < 2**30 else np.int64, copy=False)


def components(kmers: pl.DataFrame) -> tuple[int, int, int]:
    """Connected components of hit units linked by shared k-mers.

    Returns the number of components, the largest's units and its (unit, hash) pairs.
    """
    units, col = _ids(kmers["unit"].to_numpy())
    _, row = _ids(kmers["hash"].to_numpy())
    n = len(units)
    if n == 0:
        return 0, 0, 0
    label = _unit_components(col, row, n)
    per_unit = np.bincount(label)
    largest = int(per_unit.argmax())
    n_units = int(per_unit[largest])
    return len(per_unit), n_units, int((label[col] == largest).sum())


def _unit_components(col: np.ndarray, row: np.ndarray, n: int) -> np.ndarray:
    """Component label of each of ``n`` units, linked by pairs (unit ``col``, k-mer ``row``);
    components are numbered in order of their smallest unit."""
    # Units only, each linked to the next holder of the same k-mer: a unit graph with fewer
    # edges than pairs, instead of the bipartite unit-k-mer graph.
    order = np.argsort(row, kind="stable")
    by_kmer = col[order]
    same = row[order][1:] == row[order][:-1]
    graph = coo_array(
        (np.ones(int(same.sum()), dtype=np.int8), (by_kmer[:-1][same], by_kmer[1:][same])),
        shape=(n, n),
    )
    label: np.ndarray = connected_components(graph, directed=False)[1]
    return label


MAX_BATCH_PAIRS: Final = 2_000_000  # (unit, hash) rows per batch of components fitted at once


def component_batches(kmers: pl.DataFrame, max_pairs: int | None = None) -> list[pl.DataFrame]:
    """``kmers`` (one row per ``unit``, ``hash``) split into batches of whole components of
    units linked by shared k-mers, about ``max_pairs`` (default ``MAX_BATCH_PAIRS``) rows
    each (a larger component is a batch of its own). Fits that stop per component
    (:func:`em`, :func:`em_pin`) give the same result per batch as on the whole, with the
    working memory of one batch."""
    if kmers.height == 0:
        return [kmers]
    units, col = _ids(kmers["unit"].to_numpy())
    _, row = _ids(kmers["hash"].to_numpy())
    label = _unit_components(col, row, len(units))
    size = np.bincount(label[col])  # rows per component
    batch = ((np.cumsum(size) - size) // (max_pairs or MAX_BATCH_PAIRS))[
        label[col]
    ]  # by each component's start
    del units, row, label
    order = np.argsort(batch, kind="stable")
    rows = kmers[order]
    ends = np.cumsum(np.unique(batch, return_counts=True)[1])
    return [rows.slice(lo, hi - lo) for lo, hi in zip(np.r_[0, ends[:-1]], ends, strict=True)]


def _ranges(starts: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    """Concatenated ``arange(start, start + length)`` for each pair."""
    ends = np.cumsum(lengths)
    total = int(ends[-1]) if len(ends) else 0
    out: np.ndarray = np.arange(total) + np.repeat(starts - ends + lengths, lengths)
    return out


class _Block(NamedTuple):
    """The components still being fitted: their units (grouped by component), the pairs of
    those units, each pair's local unit ``c`` and local k-mer ``r``, and the k-mers."""

    unit: np.ndarray
    pairs: np.ndarray
    c: np.ndarray
    r: np.ndarray
    kmer: np.ndarray


# (coverage, second parameter) -> (new coverage, new second parameter, expected hits)
Step = Callable[[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray, np.ndarray]]


def _fit_components(
    col: np.ndarray,
    row: np.ndarray,
    lam: np.ndarray,
    other: np.ndarray,
    make_step: Callable[[_Block], Step],
    tol: float,
    max_iter: int,
    report: dict[str, int] | None = None,
) -> np.ndarray:
    """Iterate an EM per component of units linked by shared k-mers (pairs of unit ``col``
    and k-mer ``row``); return each unit's expected hits at its last step.

    Components share no k-mer, so each is its own fit and stops at its own convergence:
    its coverages ``lam`` change by at most ``tol`` x its largest, and ``other`` (presence
    or scale) by at most ``tol``. An iteration only touches components still moving (most
    converge in a few dozen, a few run to ``max_iter``); they are compacted into a new
    ``_Block``, and ``make_step`` called on it, once half of the current block has
    converged. ``lam`` and ``other`` are updated in place.
    """
    n = len(lam)
    attributed = np.zeros(n)
    label = _unit_components(col, row, n)
    active = np.argsort(label, kind="stable")  # units grouped by component
    seg = np.diff(np.flatnonzero(np.r_[True, np.diff(label[active]) != 0, True]))
    by_unit = np.argsort(col, kind="stable").astype(col.dtype)
    pair_start = np.r_[0, np.cumsum(np.bincount(col, minlength=n))]
    local = np.empty(n, dtype=col.dtype)
    it = 0
    while len(active) and it < max_iter:
        pairs = by_unit[_ranges(pair_start[active], np.diff(pair_start)[active])]
        local[active] = np.arange(len(active))
        kmer, r = _ids(row[pairs])
        step = make_step(_Block(active, pairs, local[col[pairs]], r, kmer))
        starts = np.r_[0, np.cumsum(seg)[:-1]]
        la, se, done = lam[active], other[active], np.zeros(len(active), dtype=bool)
        while it < max_iter and done.sum() * 2 <= len(active):
            it += 1
            new, new_se, att = step(la, se)
            converged = (
                np.maximum.reduceat(np.abs(new - la), starts)
                <= tol * np.maximum.reduceat(new, starts)
            ) & (np.maximum.reduceat(np.abs(new_se - se), starts) <= tol)
            la, se = new, new_se
            now = np.repeat(converged, seg) & ~done
            lam[active[now]], other[active[now]], attributed[active[now]] = (
                la[now],
                se[now],
                att[now],
            )
            done |= now
        rest = ~done  # at max_iter, units still moving keep their last values
        lam[active[rest]], other[active[rest]], attributed[active[rest]] = (
            la[rest],
            se[rest],
            att[rest],
        )
        moving = ~done[starts]
        active, seg = active[np.repeat(moving, seg)], seg[moving]
    if report is not None:  # units still moving here stopped at max_iter, not converged
        report["em_iterations"] = max(report.get("em_iterations", 0), it)
        report["em_unconverged_units"] = report.get("em_unconverged_units", 0) + len(active)
    return attributed


def unit_hits(
    table: PackedTable, max_hash_g: np.ndarray, hashes: np.ndarray, reads: np.ndarray
) -> pl.DataFrame:
    """Expand sampled query hashes to (unit, hash, read, pin_q, holders) hits through ``table``.

    A hit counts for a unit only if the hash passes that unit's ``max_hash_g``; ``pin_q`` is
    the k-mer's quantised ``p_in`` in that unit and ``holders`` the number of index units
    it counts for. Rows follow the input order, then each set's order. The table's arrays
    are read in place (``_core.unit_hits``): one bucket and one set per hash, no keys.
    """
    u64 = np.uint64
    return pl.DataFrame(
        _core.unit_hits(
            table,
            np.asarray(max_hash_g, dtype=u64),
            np.asarray(hashes, dtype=u64),
            np.asarray(reads, dtype=u64),
        )
    )


def gather(kmers: pl.DataFrame, t_g: np.ndarray) -> pl.DataFrame:
    """Assign each hit k-mer to one unit, greedily, as ``sourmash gather`` does.

    Repeatedly take the unit with the most unassigned hit k-mers, scaled by ``1 / t_g`` so
    units sampled at different rates compare as estimated k-mer counts, and give it those
    k-mers. Units left with none are explained away and get no row. Returns ``unit``,
    ``kmers_unique`` (k-mers assigned) and ``gather_rank`` (0 = taken first).
    """
    got = _core.gather(
        kmers["unit"].cast(pl.UInt32).to_numpy(),
        kmers["hash"].cast(pl.UInt64).to_numpy(),
        np.asarray(t_g, dtype=np.float64),
    )
    return pl.DataFrame(got).with_columns(
        gather_rank=pl.int_range(len(got["unit"]), dtype=pl.UInt32)
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
    tol: float = 1e-8,
    max_iter: int = 1000,
    report: dict[str, int] | None = None,
) -> pl.DataFrame:
    """Per-unit k-mer ``coverage`` (and ``present`` fraction) by EM over k-mer hit counts.

    ``kmers`` has one row per (``unit``, ``hash``) with the k-mer's ``hits``. Each k-mer's
    hits are Poisson with mean the sum of ``coverage`` over the units holding it, and each
    unit's ``m_g`` kept k-mers (hit or not) all count in its expectation; EM finds the
    maximum-likelihood coverages. Components share no k-mers, so each is fitted and
    stops on its own (:func:`_fit_components`); ``tol`` 1e-8 per component matches the
    accuracy the joint fit had at 1e-6.

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

    ``report`` collects ``em_iterations`` (the most any component took) and
    ``em_unconverged_units`` (units of components stopped at ``max_iter``), summed over calls.
    """
    kmers = kmers.sort("unit", "hash")  # sums in a fixed order: results independent of input order
    units, col = _ids(kmers["unit"].to_numpy())
    hashes, row = _ids(kmers["hash"].to_numpy())
    hits = np.zeros(len(hashes))
    hits[row] = kmers["hits"].to_numpy()
    del kmers, hashes  # the sorted copy is not needed through the fit
    m = m_g[units].astype(np.float64)
    lam = np.bincount(col, weights=hits[row], minlength=len(units)) / m
    pi = np.ones(len(units))

    def make_step(b: _Block) -> Step:
        h, mb, per_unit = hits[b.kmer], m[b.unit], len(b.unit)

        def over_units(x: np.ndarray) -> np.ndarray:  # per unit, sum of x over its k-mers
            return np.bincount(b.c, weights=x[b.r], minlength=per_unit)

        def step(la: np.ndarray, p: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            w = la * p
            mu = np.bincount(b.r, weights=w[b.c], minlength=len(b.kmer))
            attributed = w * over_units(_div(h, mu))  # expected hits from each unit
            if zero_inflated:
                kmers_hit = w * over_units(_div(np.ones_like(mu), mu))  # expected hit k-mers
                seen = -np.expm1(-la)  # chance a present k-mer is hit
                if prior is None:
                    new_p = np.minimum(1.0, _div(kmers_hit, mb * seen))
                else:
                    # Unhit k-mers are present with odds p (1 - seen) : (1 - p).
                    odds = p * (1 - seen) / np.maximum(1 - p * seen, 1e-300)
                    unhit_present = np.maximum(mb - kmers_hit, 0) * odds
                    new_p = (kmers_hit + unhit_present + prior[0] - 1) / (mb + sum(prior) - 2)
            else:
                new_p = p
            return _div(attributed, mb * new_p), new_p, attributed

        return step

    attributed = _fit_components(col, row, lam, pi, make_step, tol, max_iter, report)
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
    kmers: pl.DataFrame,
    pin_hist: np.ndarray,
    *,
    tol: float = 1e-8,
    max_iter: int = 1000,
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
    unit's kept k-mers per level (its sum is ``m_g``). Coverage is attributed hits
    over expected present k-mers (hit ones by their share, unhit ones by their posterior
    presence), and s_g is expected present k-mers over the sum of p_in. Like :func:`em`, a
    hit k-mer counts as present for a unit by its share of the hits. Returns ``unit``,
    ``coverage`` and ``present`` (expected fraction of kept k-mers present); units explained
    away get coverage 0.
    """
    kmers = kmers.sort("unit", "hash")  # sums in a fixed order: results independent of input order
    units, col = _ids(kmers["unit"].to_numpy())
    hashes, row = _ids(kmers["hash"].to_numpy())
    level = kmers["pin_q"].to_numpy().astype(np.intp)
    hits = np.zeros(len(hashes))
    hits[row] = kmers["hits"].to_numpy()
    del kmers, hashes  # the sorted copy is not needed through the fit
    n = len(units)
    hist = pin_hist[units].astype(np.float64)
    levels = hist.shape[1]
    m, expected = hist.sum(axis=1), hist @ PIN_P  # kept k-mers; present ones at s_g = 1
    lam, scale = np.bincount(col, weights=hits[row], minlength=n) / m, np.ones(n)

    def make_step(b: _Block) -> Step:
        h, lv, hs, ex = hits[b.kmer], level[b.pairs], hist[b.unit], expected[b.unit]
        cell = b.c * levels + lv

        def step(la: np.ndarray, sc: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            weight = la[b.c] * sc[b.c] * PIN_P[lv]
            share = _div(weight, np.bincount(b.r, weights=weight, minlength=len(b.kmer))[b.r])
            att = np.bincount(b.c, weights=share * h[b.r], minlength=len(b.unit))
            hit = np.bincount(cell, weights=share, minlength=hs.size).reshape(hs.shape)
            seen = -np.expm1(-la)[:, None]  # chance a present k-mer is hit
            pi = sc[:, None] * PIN_P
            odds = pi * (1 - seen) / np.maximum(1 - pi * seen, 1e-300)
            present = hit.sum(axis=1) + (np.maximum(hs - hit, 0) * odds).sum(axis=1)
            return _div(att, present), np.clip(_div(present, ex), 0.0, 1 / PIN_P[-1]), att

        return step

    # tol 1e-8 per component matches the joint fit's accuracy at 1e-6 (phase 6, step 15).
    attributed = _fit_components(col, row, lam, scale, make_step, tol, max_iter)
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
    for i in range(max_iter):
        w = np.bincount(capped, weights=prob, minlength=H_CAP + 1)
        new = w[capped] / (w[capped] + (n_index_units - prob.sum()) * absent)
        done = np.abs(new - prob).max(initial=0) <= tol
        prob = new
        if done:
            print(f"presence converged in {i + 1} iterations", file=sys.stderr, flush=True)
            break
    else:
        print(
            f"presence failed to converge in {max_iter} iterations, "
            f"max change {np.abs(new - prob).max(initial=0):.2e}",
            file=sys.stderr,
            flush=True,
        )
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
    pairs: pl.DataFrame,
    hash_reads: pl.DataFrame,
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

    ``pairs`` has one row per hit (``unit``, ``hash``) and ``hash_reads`` one per (``hash``,
    ``read``) with the ``n`` times the read hit the k-mer (every unit holding a k-mer sees the
    same reads, so they are stored once per k-mer). Each of ``draws`` replicates

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

    Units in different components (linked by shared k-mers) never interact, so the draws run
    per batch of components (:func:`component_batches`) sized to ``POSTERIOR_BATCH_BYTES``,
    each with its own generator seeded by ``seed`` and the batch's smallest unit. A read's
    Poisson weight is a function of ``seed``, the draw and its id (:func:`_poisson1`), so
    reads hitting several batches weigh the same in each.
    """
    per_pair = 1600 + 16 * draws  # bytes per (unit, hash) row at peak, measured (step 28)
    parts = component_batches(
        pairs.select("unit", "hash"), max(1, POSTERIOR_BATCH_BYTES // per_pair)
    )
    return pl.concat(
        [
            _posterior_batch(
                part,
                # This batch's k-mers' rows: a scan per batch, no copy of all rows.
                hash_reads.join(part.select("hash").unique(), on="hash", how="semi"),
                m_g,
                pin_sum,
                draws,
                present_prob=present_prob,
                len_cv=len_cv,
                copies_error=copies_error,
                sweeps=sweeps,
                level=level,
                shared_evidence=shared_evidence,
                seed=seed,
            )
            for part in parts
        ]
    ).sort("unit")


POSTERIOR_BATCH_BYTES: Final = 2**29  # working memory of one batch of the posterior's draws

# Poisson(1) CDF at 0..24; beyond, a probability below 1e-24.
POISSON1_CDF: Final = np.cumsum(np.exp(-1.0) / np.cumprod(np.r_[1.0, np.arange(1.0, 25.0)]))


def _poisson1(seed: int, draw: int, reads: np.ndarray) -> np.ndarray:
    """Poisson(1) weight of each read (id) in a draw: a SplitMix64 hash of (``seed``,
    ``draw``, read) to a uniform, then the inverse CDF, so the same read has the same weight
    wherever it is drawn."""
    with np.errstate(over="ignore"):
        x = np.asarray(reads, dtype=np.uint64) + np.uint64(
            (seed * 1_000_003 + draw) % 2**64
        ) * np.uint64(0x9E3779B97F4A7C15)
        for shift, mult in ((30, 0xBF58476D1CE4E5B9), (27, 0x94D049BB133111EB)):
            x = (x ^ (x >> np.uint64(shift))) * np.uint64(mult)
        x ^= x >> np.uint64(31)
    uniform = (x >> np.uint64(11)).astype(np.float64) * 2.0**-53
    return np.searchsorted(POISSON1_CDF, uniform, side="right")


def _posterior_batch(
    pairs: pl.DataFrame,
    hash_reads: pl.DataFrame,
    m_g: np.ndarray,
    pin_sum: np.ndarray,
    draws: int,
    *,
    present_prob: np.ndarray | None,
    len_cv: np.ndarray | None,
    copies_error: float,
    sweeps: int,
    level: float,
    shared_evidence: float,
    seed: int,
) -> pl.DataFrame:
    """:func:`posterior_zi` on one batch of whole components."""
    # Sorted, so the seeded draws do not depend on row order (Polars group_by does not fix it).
    keyed = pairs.select("unit", "hash").sort("unit", "hash")
    units, col = np.unique(keyed["unit"].to_numpy(), return_inverse=True)
    rng = np.random.default_rng([seed, int(units[0])])
    hashes, row = np.unique(keyed["hash"].to_numpy(), return_inverse=True)
    n_units, n_rows = len(units), len(hashes)
    by_hash = hash_reads.join(pl.DataFrame({"hash": hashes}), on="hash", how="semi").sort(
        "hash", "read"
    )
    hash_of_hit: np.ndarray = np.searchsorted(hashes, by_hash["hash"].to_numpy()).astype(
        np.int32 if n_rows < 2**31 else np.int64
    )
    reads, read_of = _ids(by_hash["read"].to_numpy())
    n = by_hash["n"].to_numpy().astype(np.int64)
    del by_hash
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
        weight = _poisson1(seed, b, reads)[read_of]
        per_hash = np.rint(np.bincount(hash_of_hit, weights=n * weight, minlength=n_rows))
        per_key = per_hash[row]
        table = keyed.select("unit", "hash").with_columns(hits=per_key).filter(pl.col("hits") > 0)
        if table.height == 0:
            continue
        fit = em(table, m_g, zero_inflated=True)
        at = np.searchsorted(units, fit["unit"].to_numpy())
        lam, pi = np.full(n_units, 1e-6), np.full(n_units, 1e-6)
        lam[at] = np.maximum(fit["coverage"].to_numpy(), 1e-6)
        pi[at] = np.clip(fit["present"].to_numpy(), 1e-6, 1 - 1e-6)
        # Holders of a k-mer see the same reads, so any holder's count is the k-mer's.
        count = per_hash.astype(np.int64)
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
            u = rng.random(n_units)
            pick = np.empty(n_units, dtype=np.intp)
            # Units x grid arrays, a chunk of units at a time (8 GB per 2 M units at once).
            chunk = max(1, CHUNK_BYTES // (8 * len(PI_GRID)))
            for lo in range(0, n_units, chunk):
                part = slice(lo, lo + chunk)
                log_post = hit[part, None] * np.log(PI_GRID) + (m - hit)[part, None] * np.log1p(
                    -PI_GRID * seen[part, None]
                )
                cdf = np.cumsum(np.exp(log_post - log_post.max(axis=1, keepdims=True)), axis=1)
                pick[part] = (cdf < u[part, None] * cdf[:, -1:]).sum(axis=1)
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
    total = [np.stack([np.bincount(group, w, minlength=n_units) for w in x])[:, group]
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


CHUNK_BYTES: Final = 2**25  # units x grid arrays of the posterior are built in chunks of this


class _Summed:
    """Per-batch aggregates summed into one table, re-aggregated each time the batches
    waiting outgrow the running total, so memory stays about twice the total's.

    The table is split into ``parts`` by the first key (mod ``parts``), each merged on its
    own, so a merge's transient memory (polars' group-by over the concatenation) is that of
    one part, not the whole table."""

    min_rows = 1_000_000  # rows waiting before the first merge, over all parts

    def __init__(self, keys: list[str], aggs: dict[str, pl.Expr], parts: int = 16) -> None:
        self.keys, self.aggs, self.parts = keys, aggs, parts
        self.totals: list[pl.DataFrame | None] = [None] * parts
        self.pending: list[list[pl.DataFrame]] = [[] for _ in range(parts)]
        self.rows = [0] * parts

    def add(self, df: pl.DataFrame) -> None:
        self.empty = df.clear()
        split = df.with_columns(_part=pl.col(self.keys[0]) % self.parts).partition_by(
            "_part", as_dict=True, include_key=False
        )
        for key, part in split.items():
            p = int(key[0])
            self.pending[p].append(part)
            self.rows[p] += part.height
            total = self.totals[p]
            if self.rows[p] > max(
                self.min_rows // self.parts, 0 if total is None else total.height
            ):
                self.totals[p], self.pending[p], self.rows[p] = self._merged(p), [], 0

    def _merged(self, p: int) -> pl.DataFrame:
        total = self.totals[p]
        parts = [total, *self.pending[p]] if total is not None else self.pending[p]
        return pl.concat([self.empty, *parts]).group_by(self.keys).agg(**self.aggs)

    def total(self) -> pl.DataFrame:
        """The summed table. It releases the parts, so it is called once, at the end."""
        merged = [self._merged(p) for p in range(self.parts)]
        self.totals, self.pending = [None] * self.parts, [[] for _ in range(self.parts)]
        return pl.concat(merged, rechunk=False)


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
    timer: Timer | None = None,
    all_estimators: bool = False,
    low_memory: bool = False,
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
    are omitted. ``coverage_em`` is the shipped estimate; the others (``_zi`` and its
    ``copies``/``abundance``, ``_zib``, ``_zip``, ``_wta``, ``_ufirst``) are fitted only
    with ``all_estimators`` (benchmarks), except that ``draws`` > 0 fits ``_zi``, whose
    intervals the posterior gives.

    With a dense tier, the reads are streamed a second time at its rate and the EM
    estimates (``coverage_em``, ``_zi``, ``_zib``) are fitted on the dense hits of the units
    gather keeps, over their ``m_dense`` k-mers; ``kmers_dense`` counts those hit.

    ``present_prob`` (:func:`presence`) is the probability that a unit gather keeps is
    present rather than hit by background, from the tier-2 k-mers gather gave it; units
    gather drops get 0.

    ``draws`` > 0 adds 95% posterior intervals for ``coverage_zi`` and ``abundance_zi`` and
    ambiguity groups from that many Gibbs sweeps (:func:`posterior_zi`). Its per-read rows
    are kept from the first pass over the reads, or with ``low_memory`` rebuilt by a second
    pass for the detected units' k-mers only (same result; a read pass more, the rows of
    all hit units never held). ``kmers_out`` writes
    the tier-2 hits per (``unit``, ``hash``) with ``hits`` and ``holders`` to Parquet.

    ``timer`` records each stage's time and peak RSS and the counts the query's cost
    hinges on (sampled k-mers, distinct ones estimated on 1 in ``DISTINCT_SAMPLE`` of hash
    space, hit k-mers, hit rows, (unit, hash) pairs, detected units, component sizes).
    """
    record = timer is not None
    timer = timer or Timer()
    counts = timer.counts
    # Format-1 indexes also record tier1_per_unit.
    params = IndexParams(**{f.name: index.meta["params"][f.name] for f in fields(IndexParams)})

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

    max_hash_g = index.units["max_hash_g"]

    def by_read(hits: pl.DataFrame) -> pl.DataFrame:
        return hits.group_by("unit", "hash", "read").agg(
            n=pl.len(), pin_q=pl.col("pin_q").first(), holders=pl.col("holders").first()
        )

    def reread(kmers: np.ndarray) -> pl.DataFrame:
        """(hash, read, n) rows for the sampled ``kmers``: times each read has each, which is
        also a hit's rows over its holders, so these are the first pass's rows for them."""
        wanted = np.sort(kmers).astype(np.uint64)
        parts = [pl.DataFrame(schema={"hash": pl.UInt64, "read": pl.UInt64, "n": pl.UInt32})]
        for b in stream(index.tier2.max_hash):
            at = np.minimum(np.searchsorted(wanted, b["hash"]), max(len(wanted) - 1, 0))
            keep = wanted[at] == b["hash"] if len(wanted) else np.zeros(len(at), dtype=bool)
            found = pl.DataFrame({"hash": b["hash"][keep], "read": b["read"][keep]})
            parts.append(found.group_by("hash", "read").agg(n=pl.len()))
        return pl.concat(parts)

    def per_kmer(per_read: pl.DataFrame) -> pl.DataFrame:
        return per_read.group_by("unit", "hash").agg(
            hits=pl.col("n").sum(), pin_q=pl.col("pin_q").first(), holders=pl.col("holders").first()
        )

    # Per-(unit, hash) hits and per-unit reads are summed batch by batch (read ids never
    # span batches), so per-read rows are kept only for the posterior.
    first = {"pin_q": pl.col("pin_q").first(), "holders": pl.col("holders").first()}
    pairs = _Summed(["unit", "hash"], {"hits": pl.col("hits").sum(), **first})
    reads_summed = _Summed(["unit"], {"reads": pl.col("reads").sum()})
    empty = np.empty(0, dtype=np.uint64)
    batches, n_reads, subsample = [], 0, []
    # The posterior's per-read rows come from this pass, unless a dense tier's pass gives
    # them or ``low_memory`` re-reads the reads for the detected k-mers only.
    keep_reads = draws > 0 and index.dense is None and not low_memory
    sampled = hit_kmers = hit_rows = 0
    reads = iter(stream(index.tier2.max_hash))
    b: dict[str, np.ndarray] | None = {"hash": empty, "read": empty}
    while b is not None:
        with timer("lookup"):
            hits = unit_hits(index.tier2, max_hash_g, b["hash"], b["read"])
        with timer("aggregate"):
            pairs.add(hits.group_by("unit", "hash").agg(hits=pl.len(), **first))
            reads_summed.add(
                hits.select("unit", "read").unique().group_by("unit").agg(reads=pl.len())
            )
            if keep_reads:
                # Once per (hash, read): each hit has a row per unit holding the k-mer.
                batches.append(
                    hits.group_by("hash", "read").agg(n=pl.len() // pl.col("holders").first())
                )
        if len(b["read"]):  # reads are numbered in input order; the last has sampled hashes
            n_reads = max(n_reads, int(b["read"].max()) + 1)
        sampled += len(b["hash"])
        subsample.append(b["hash"][b["hash"] <= index.tier2.max_hash // DISTINCT_SAMPLE])
        hit_rows += hits.height
        hit_kmers += round((1 / hits["holders"]).sum()) if hits.height else 0  # rows per hit
        with timer("hash"):
            b = next(reads, None)
    with timer("aggregate"):
        per_read = pl.concat(batches) if batches else None
        kmer_hits, unit_reads = pairs.total(), reads_summed.total()
    counts |= {
        "reads": n_reads,
        "sampled_kmers": sampled,
        "distinct_sampled_kmers_est": len(np.unique(np.concatenate(subsample))) * DISTINCT_SAMPLE,
        "hit_kmers": hit_kmers,
        "hit_rows": hit_rows,
        "read_rows": 0 if per_read is None else per_read.height,
        "unit_kmer_pairs": kmer_hits.height,
        "hit_units": kmer_hits["unit"].n_unique(),
    }
    with timer("hit_units"):  # rows gathered from the unit columns, never the whole table
        hit_info = index.units.rows(np.sort(kmer_hits["unit"].unique().to_numpy()))
    if record:
        with timer("components"):
            (
                counts["components"],
                counts["largest_component_units"],
                counts["largest_component_pairs"],
            ) = components(kmer_hits)
    if kmers_out is not None:
        kmer_hits.write_parquet(kmers_out)
    # From here units are numbered 0.. in the hit units' order (``index`` keeps the index's
    # ids), so per-unit arrays come from their rows only: a full column is ~4 GB per 8 bytes.
    ids = hit_info.select(
        pl.col("unit").cast(kmer_hits.schema["unit"]),
        index=pl.int_range(pl.len(), dtype=pl.UInt32),
    )

    def renumber(df: pl.DataFrame) -> pl.DataFrame:
        return df.join(ids, on="unit").drop("unit").rename({"index": "unit"})

    kmer_hits, unit_reads = renumber(kmer_hits), renumber(unit_reads)
    hit_info = hit_info.with_columns(index=pl.col("unit"), unit=ids["index"])
    t_g = hit_info["t_g"].to_numpy()
    with timer("gather"):
        assigned = gather(kmer_hits.select("unit", "hash"), t_g)
    counts["detected_units"] = assigned.height
    detected = kmer_hits.join(assigned.select("unit"), on="unit")
    m_g = hit_info["m_g"].to_numpy()
    pin_hist = hit_info["pin_hist"].to_numpy()
    pin_sum = hit_info["pin_sum"]
    len_cv = hit_info["len_cv"].to_numpy()
    # Each hit k-mer is gather's: the first unit in gather order holding it.
    own = (
        kmer_hits.join(assigned.select("unit", "gather_rank"), on="unit")
        .sort("gather_rank")
        .unique("hash", keep="first")
        .select("unit", "hash", "holders")
    )
    with timer("presence"):
        present_prob = presence(own, t_g, n_reads, index.units.height)
    del own
    dense = index.dense
    if dense is not None:
        # Second pass: every k-mer at the dense rate, for the detected units only.
        max_hash_dense = index.units["max_hash_dense"]
        keep = assigned.select("unit")  # detected units are hit units, so renumber keeps all
        with timer("dense"):
            detected_reads = pl.concat(
                [
                    by_read(
                        renumber(unit_hits(dense, max_hash_dense, b["hash"], b["read"])).join(
                            keep, on="unit", how="semi"
                        )
                    )
                    for b in [{"hash": empty, "read": empty}, *stream(dense.max_hash)]
                ]
            )
            detected = per_kmer(detected_reads)
            if draws > 0:
                per_read = detected_reads.group_by("hash", "read").agg(pl.col("n").first())
        m_g = hit_info["m_dense"].to_numpy()
        pin_hist = hit_info["pin_hist_dense"].to_numpy()
        pin_sum = hit_info["pin_sum_dense"]
        len_cv = hit_info["len_cv_dense"].to_numpy()
    with timer("fit_em"):
        # Fits stop per component, so batches of components give the same units' results.
        parts = component_batches(detected)
        detected = pl.concat(parts, rechunk=False)  # the batches' rows, not a second copy

        def per_batch(fit: Callable[[pl.DataFrame], pl.DataFrame]) -> pl.DataFrame:
            return pl.concat([fit(part) for part in parts]).sort("unit")

        counts["fit_batches"] = len(parts)
        counts["fit_largest_batch_pairs"] = max(part.height for part in parts)
        plain = per_batch(lambda part: em(part, m_g, report=counts)).select(
            "unit", coverage_em="coverage"
        )
    fits = [plain]
    if all_estimators or draws > 0:
        with timer("fit_zi"):
            inflated = per_batch(lambda part: em(part, m_g, zero_inflated=True))
        # Present k-mers over an average member's kept k-mers: member-equivalents present.
        copies = pl.col("present") * pl.col("m") / pl.col("pin_sum")
        fits.append(
            inflated.join(
                pl.DataFrame(
                    {"unit": np.arange(len(m_g), dtype=np.uint32), "m": m_g, "pin_sum": pin_sum}
                ),
                on="unit",
            ).select(
                "unit",
                coverage_zi="coverage",
                present_zi="present",
                copies_zi=copies,
                abundance_zi=pl.col("coverage") * copies,
            )
        )
    if all_estimators:
        with timer("fit_zib"):
            prior = fit_present_prior(inflated)
            fitted = (
                per_batch(lambda part: em(part, m_g, zero_inflated=True, prior=prior))
                if prior
                else inflated
            )
            fits.append(fitted.select("unit", coverage_zib="coverage", present_zib="present"))
        with timer("fit_zip"):
            fits.append(
                per_batch(lambda part: em_pin(part, pin_hist)).select(
                    "unit", coverage_zip="coverage", present_zip="present"
                )
            )
    unit_info = hit_info.select(pl.exclude("^(pin_(hist|sum)|len_cv).*$"))
    if "name" not in unit_info.columns:  # sourmash imports name units, builds by cluster_rep
        unit_info = unit_info.with_columns(name=pl.col("cluster_rep").cast(pl.String))
    with timer("result"):
        result = (
            kmer_hits.group_by("unit")
            .agg(hits=pl.col("hits").sum().cast(pl.UInt64), kmers_hit=pl.len().cast(pl.UInt32))
            .join(unit_reads.with_columns(pl.col("reads").cast(pl.UInt64)), on="unit")
            .join(unit_info, on="unit")
            .join(assigned, on="unit", how="left")
            .join(present_prob, on="unit", how="left")
        )
        for fit in fits:
            result = result.join(fit, on="unit", how="left")
    if draws > 0 and per_read is None:
        with timer("reread"):
            per_read = reread(detected["hash"].unique().to_numpy())
    if draws > 0:
        prob = np.zeros(len(m_g))
        prob[present_prob["unit"].to_numpy()] = present_prob["present_prob"].to_numpy()
        assert per_read is not None  # per-read rows are kept when draws > 0
        with timer("posterior"):
            intervals = posterior_zi(
                detected.select("unit", "hash"),
                per_read,
                m_g,
                pin_sum.to_numpy(),
                draws,
                present_prob=prob,
                len_cv=len_cv,
            )
        result = result.join(intervals, on="unit", how="left")
    if dense is not None:
        result = result.join(
            detected.group_by("unit").agg(kmers_dense=pl.len().cast(pl.UInt32)),
            on="unit",
            how="left",
        ).with_columns(pl.col("kmers_dense").fill_null(0))
    if all_estimators:
        rated = kmer_hits.join(hit_info.select("unit", "m_g", "t_g"), on="unit")
        with timer("baselines"):
            for rule, score in (("wta", WTA_SCORE), ("ufirst", UFIRST_SCORE)):
                won = (
                    assign_best(rated, score)
                    .group_by("unit")
                    .agg(
                        pl.len().cast(pl.UInt32).alias(f"kmers_{rule}"),
                        pl.col("hits").sum().alias(f"_{rule}"),
                    )
                )
                result = result.join(won, on="unit", how="left").with_columns(
                    pl.col(f"kmers_{rule}").fill_null(0),
                    (pl.col(f"_{rule}").fill_null(0) / pl.col("m_g")).alias(f"coverage_{rule}"),
                )
            result = result.drop("_wta", "_ufirst")
    return (
        result.with_columns(
            pl.col("^(coverage|present|copies|abundance)_(em|zi|zib|zip)(_lo|_hi)?$").fill_null(
                0.0
            ),
            pl.col("present_prob").fill_null(0.0),
            containment=pl.col("kmers_hit") / pl.col("m_g"),
            coverage=pl.col("hits") / pl.col("m_g"),
            kmers_unique=pl.col("kmers_unique").fill_null(0),
        )
        .with_columns(unit=pl.col("index"))  # back to the index's unit ids
        .drop("index")
        .sort("unit")
    )
