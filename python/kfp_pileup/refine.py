"""Post-processing of the pileups (plan, phase 13, steps 7-9), run by ``call`` after the
crude pass (:mod:`kfp_pileup.pileup`) on the alignments it keeps.

- *Step 7, base qualities and a learned error model:* qualities re-read from the FASTQ for
  recruited reads only; a codon's error is 1 - prod(1 - e_base) over its three bases. At deep
  target columns with an unambiguous consensus, mismatches against it by quality and read
  cycle (thirds of the read), fitted as a non-increasing error(quality) per cycle bin (pool
  adjacent violators, shrunk to the nominal error at sparse bins): ``error_model.tsv``. The
  error is carried through translation by :func:`substitution_shares` (the chance a base
  error turns a codon of a into b). Variants are tested LoFreq-style: the alt residue's
  weight against a Poisson with mean sum(w e) x share(null -> alt) (the Poisson-binomial's
  approximation), BH-adjusted, with strand (Fisher) and read-position (Mann-Whitney on the
  distance to the read's nearer end) bias: ``calls.tsv``. Linkage is step 9's haplotypes.
- *Step 8, the joint EM over read sources:* each pair's source among the units it aligns to
  (targets and neighbours), its component's reversed decoys and its component's novel
  relative. A mate's likelihood under a source is binomial in the differences of its best
  alignment there (``mism``, plus read residues the alignment leaves out), at the source's
  divergence d_s: d^m (1 - d)^(n - m) / (A - 1)^m. Units' d lies in [error, ``novel``],
  the novel relative's in [``novel``, 0.5] on the pair's best alignment (a relative missing
  from the reference), decoys' in [0.3, 0.95]. A mate with no alignment to a source takes its
  worst alignment's differences (an optimistic bound, as DIAMOND's ``--top`` cut it). Pairs
  whose mates sit on one sequence get an insert-size term, flat within 2 sigma. pi is
  started from ``coverage_zi`` x effective length; the prior fit adds ``PRIOR_PAIRS``
  pseudo-pairs per component spread the same way, the likelihood-only fit none. Pairs whose
  best source differs between the two are flagged; pileups use the likelihood-only fit.
  ``sources.tsv``: pairs, abundance (pairs over effective length, positions a fragment can
  start from while overlapping the target by ``OVERLAP`` nt) and divergence per source.
- *Step 9, consensus updates and haplotypes:* where a target column's weighted majority
  differs from the current sequence, the sequence is updated and its alignments re-scored
  along their paths (substitutions only), then the EM re-run, up to ``ROUNDS`` times. Then
  per target sequence, a mixture of K <= ``K_MAX`` haplotypes over its variable columns, K
  by BIC; a haplotype diverging from the original sequence by more than ``novel`` over the
  covered columns is a novel relative, and its reads' weight moves to the novel class:
  ``haplotypes.tsv``.
"""

import gzip
from collections.abc import Sequence
from pathlib import Path
from typing import Final

import numpy as np
import polars as pl
from scipy import special, stats

from kmer_functional_profiler.reference import BASES, CODE_11

NOVEL: Final = 0.1  # least divergence of a novel relative: MGnify90's cluster radius
PRIOR_PAIRS: Final = 2.0  # weak Dirichlet: pseudo-pairs per component from coverage_zi
OVERLAP: Final = 60  # nt a mate shares with a target to align (~20 aa)
DEEP: Final = 10.0  # error model: least column depth ...
UNAMBIGUOUS: Final = 0.9  # ... and consensus frequency
PSEUDO: Final = 10.0  # error model: pseudo-observations at the nominal error
DEFAULT_ERROR: Final = 0.01  # residue error without qualities or deep columns
CYCLE_BINS: Final = 3
MAX_Q: Final = 45
ROUNDS: Final = 5
K_MAX: Final = 3
MIN_COVERED: Final = 10  # novel relatives: least covered columns
MIN_ALT_READS: Final = 2.0  # haplotypes: least weight of a variable column's second residue
EM_ITER, EM_TOL = 2000, 1e-8
DECOY_D: Final = (0.3, 0.95)
GAP: Final = ord("-")


def substitution_shares(protein: bool) -> np.ndarray:
    """256 x 256: the share of a's errors observed as b. Amino acids: single-base changes of
    a's codons (synonymous codons equally likely) that give b; bases: 1/3 each."""
    m = np.zeros((256, 256))
    if protein:
        for i, a in enumerate(CODE_11):
            codon = [BASES[i // 16], BASES[i // 4 % 4], BASES[i % 4]]
            for pos in range(3):
                for b in BASES.replace(codon[pos], ""):
                    c = [*codon[:pos], b, *codon[pos + 1 :]]
                    t = CODE_11[BASES.index(c[0]) * 16 + BASES.index(c[1]) * 4 + BASES.index(c[2])]
                    if t != a:
                        m[ord(a), ord(t)] += 1
    else:
        for a in "ACGT":
            for b in "ACGT".replace(a, ""):
                m[ord(a), ord(b)] = 1
    total = m.sum(1, keepdims=True)
    return np.divide(m, total, out=np.zeros_like(m), where=total > 0)


def read_qualities(paths: Sequence[str | Path], names: set[str]) -> dict[tuple[str, int], bytes]:
    """Phred+33 quality strings of the named reads, keyed (name, mate); one pass per file."""
    out: dict[tuple[str, int], bytes] = {}
    for mate, path in enumerate(paths):
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as f:
            name = ""
            for i, line in enumerate(f):
                if i % 4 == 0:
                    name = line[1:].split()[0]
                    if name.endswith(("/1", "/2")):
                        name = name[:-2]
                elif i % 4 == 3 and name in names:
                    out[(name, mate)] = line.rstrip("\n").encode()
    return out


def _decreasing(y: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Weighted least-squares non-increasing fit (pool adjacent violators)."""
    blocks: list[list[float]] = []  # value, weight, size
    for v, wt in zip(y, w, strict=True):
        blocks.append([v, wt, 1])
        while len(blocks) > 1 and blocks[-2][0] < blocks[-1][0]:
            v2, w2, n2 = blocks.pop()
            v1, w1, n1 = blocks.pop()
            blocks.append([(v1 * w1 + v2 * w2) / (w1 + w2), w1 + w2, n1 + n2])
    return np.repeat([b[0] for b in blocks], [int(b[2]) for b in blocks])


def _site_errors(
    cols: pl.DataFrame, hits: pl.DataFrame, quals: dict[tuple[str, int], bytes], protein: bool
) -> pl.DataFrame:
    """Per target column row: nominal error (null without qualities), its Phred bin (-1
    without), cycle bin, distance to the read's nearer end and the read's strand."""
    keys = hits.select("hit", "read_name", "mate", "reverse", "qlen")
    rows = cols.join(keys, on="hit")
    k = rows.filter(pl.col("qpos") >= 0)
    qpos = k["qpos"].to_numpy()
    p = qpos >> 1 if protein else qpos
    width = 3 if protein else 1
    lo = p - 2 * (qpos & 1) if protein else p
    if quals:
        ids = {key: i for i, key in enumerate(quals)}
        offsets = np.cumsum([0, *(len(q) for q in quals.values())])
        qual = np.frombuffer(b"".join(quals.values()), np.uint8).astype(float) - 33
        idx = np.array(
            [ids.get(key, -1) for key in zip(k["read_name"], k["mate"], strict=True)], np.int64
        )
        length = np.where(idx >= 0, offsets[idx + 1] - offsets[idx], 0)
        e_ok = 1.0 - np.ones(len(p))
        known = idx >= 0
        for j in range(width):
            at = np.clip(lo + j, 0, np.maximum(length - 1, 0))
            e = 10 ** (-qual[np.where(known, offsets[idx] + at, 0)] / 10)
            e_ok = np.where(known, 1 - (1 - e_ok) * (1 - e), np.nan)
        nominal = e_ok
    else:
        length = (k["qlen"].to_numpy() * width).astype(np.int64)
        known = np.zeros(len(p), bool)
        nominal = np.full(len(p), np.nan)
    length = np.maximum(length, p + width)
    qbin = np.where(
        known, np.clip(np.round(-10 * np.log10(np.maximum(nominal, 1e-10))), 0, MAX_Q), -1
    )
    return k.select("hit", "col", "qpos", "residue", "reverse").with_columns(
        nominal=pl.Series(nominal),
        qbin=pl.Series(qbin.astype(np.int16)),
        cycle=pl.Series(np.minimum(CYCLE_BINS * p // length, CYCLE_BINS - 1).astype(np.int8)),
        from_end=pl.Series(np.minimum(p, length - 1 - p)),
    )


def error_model(sites: pl.DataFrame, weights: pl.DataFrame) -> pl.DataFrame:
    """Error rate per (cycle, qbin) from mismatches at deep unambiguous columns (weights:
    hit, seq, weight), shrunk to the nominal error and made non-increasing in quality."""
    w = sites.filter(pl.col("residue") != GAP).join(weights, on="hit")
    cons = (
        w.group_by("seq", "col", "residue")
        .agg(pl.col("weight").sum())
        .group_by("seq", "col")
        .agg(
            consensus=pl.col("residue").sort_by("weight").last(),
            freq=pl.col("weight").max() / pl.col("weight").sum(),
            depth=pl.col("weight").sum(),
        )
        .filter(pl.col("depth") >= DEEP, pl.col("freq") >= UNAMBIGUOUS)
    )
    seen = (
        w.join(cons, on=["seq", "col"])
        .group_by("cycle", "qbin")
        .agg(
            weight=pl.col("weight").sum(),
            mismatches=(pl.col("weight") * (pl.col("residue") != pl.col("consensus"))).sum(),
        )
    )
    total = seen.select(pl.col("mismatches").sum() / pl.col("weight").sum()).item()
    overall = total if total is not None and np.isfinite(total) else DEFAULT_ERROR
    grid = (
        sites.group_by("cycle", "qbin")
        .agg(nominal=pl.col("nominal").mean())
        .join(seen, on=["cycle", "qbin"], how="left")
        .with_columns(pl.col("weight", "mismatches").fill_null(0.0))
        .with_columns(prior=pl.col("nominal").fill_nan(None).fill_null(overall))
        .with_columns(
            raw=(pl.col("mismatches") + PSEUDO * pl.col("prior")) / (pl.col("weight") + PSEUDO)
        )
        .sort("cycle", "qbin")
    )
    fitted = [
        _decreasing(g["raw"].to_numpy(), (g["weight"] + PSEUDO).to_numpy())
        for _, g in grid.group_by("cycle", maintain_order=True)
    ]
    return grid.with_columns(error=pl.Series(np.concatenate(fitted))).drop("prior")


def _em(
    pair: np.ndarray,
    src: np.ndarray,
    diffs: np.ndarray,
    match: np.ndarray,
    ins: np.ndarray,
    init: np.ndarray,
    prior: np.ndarray,
    lo: np.ndarray,
    hi: np.ndarray,
    alpha: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Responsibilities (per row; rows sorted by pair), pairs per source and divergences.
    ponytail: plain EM; SQUAREM (as kfp's) if iterations dominate at scale."""
    starts = np.flatnonzero(np.r_[True, pair[1:] != pair[:-1]])
    sizes = np.diff(np.r_[starts, len(pair)])
    n_src = len(init)
    pi, d = init / init.sum(), np.clip(np.sqrt(lo * hi), lo, hi)
    r = np.zeros(len(pair))
    for _ in range(EM_ITER):
        ll = (
            np.log(np.maximum(pi[src], 1e-300))
            + diffs * np.log(d[src] / (alpha - 1))
            + match * np.log1p(-d[src])
            + ins
        )
        e = np.exp(ll - np.repeat(np.maximum.reduceat(ll, starts), sizes))
        r = e / np.repeat(np.add.reduceat(e, starts), sizes)
        n = np.bincount(src, r, n_src)
        new = (n + prior) / (n + prior).sum()
        dm, dn = np.bincount(src, r * diffs, n_src), np.bincount(src, r * (diffs + match), n_src)
        d = np.clip(np.divide(dm, dn, out=d.copy(), where=dn > 0), lo, hi)
        done = np.abs(new - pi).max() < EM_TOL
        pi = new
        if done:
            break
    return r, np.bincount(src, r, n_src), d


class _Model:
    """Step 8's EM inputs from the alignments: sources, pairs and their candidate rows."""

    def __init__(
        self,
        hits: pl.DataFrame,
        refs: pl.DataFrame,
        prof: pl.DataFrame | None,
        protein: bool,
        error: float,
        novel: float,
    ) -> None:
        self.alpha = 20 if protein else 4
        nt = 3 if protein else 1
        h = hits.join(refs.select("seq", "component"), on="seq").with_columns(
            pair=pl.col("read_name").rank("dense").cast(pl.Int64) - 1
        )
        # sources: units (targets, neighbours), then a decoy and a novel class per component
        units = (
            refs.filter(pl.col("role") != 2)
            .select("unit", "component", "role")
            .unique("unit", maintain_order=True)
        )
        comps = refs["component"].unique().sort()
        n_u, n_c = units.height, comps.len()
        src_unit = dict(zip(units["unit"], range(n_u), strict=True))
        comp_i = dict(zip(comps, range(n_c), strict=True))
        h = h.with_columns(
            src=pl.when(pl.col("role") == "decoy")
            .then(n_u + pl.col("component").replace_strict(comp_i, return_dtype=pl.Int64))
            .otherwise(pl.col("unit").replace_strict(src_unit, default=0, return_dtype=pl.Int64))
        )
        n_read = pl.col("qlen").max().over("pair", "mate")
        h = h.with_columns(m=pl.col("mism") + n_read - pl.col("qlen"), n=n_read)
        self.hits = h
        # fragment length from pairs whose mates' best alignments share a sequence
        best = h.sort("score", descending=True).group_by("pair", "mate", "src").first()
        same = (
            best.group_by("pair", "src")
            .agg(
                seqs=pl.col("seq").n_unique(),
                mates=pl.len(),
                span=(pl.col("end").max() - pl.col("start").min()) * nt,
            )
            .filter(pl.col("mates") == 2, pl.col("seqs") == 1)
        )
        spans = same["span"].to_numpy().astype(float)
        self.fragment = float(np.median(spans)) if len(spans) else 0.0
        sigma = 1.4826 * float(np.median(np.abs(spans - self.fragment))) if len(spans) else 0.0
        # novel class: each mate's best non-decoy alignment; the pair's component from its best
        non_decoy = h.filter(pl.col("role") != "decoy").sort("score", descending=True)
        pair_comp = (
            non_decoy.group_by("pair")
            .first()
            .select("pair", novel=n_u + n_c + pl.col("component").replace_strict(comp_i))
        )
        novel_rows = (
            non_decoy.group_by("pair", "mate")
            .first()
            .join(pair_comp, on="pair")
            .select("pair", "mate", src=pl.col("novel").cast(pl.Int64), m="m")
        )
        per_src = pl.concat([best.select("pair", "mate", "src", "m"), novel_rows])
        mates = h.group_by("pair", "mate").agg(n=pl.col("n").first(), worst=pl.col("m").max())
        cand = (
            per_src.select("pair", "src")
            .unique()
            .join(mates, on="pair")
            .join(per_src, on=["pair", "mate", "src"], how="left")
            .with_columns(m=pl.col("m").fill_null(pl.col("worst")))
            .group_by("pair", "src")
            .agg(diffs=pl.col("m").sum(), match=(pl.col("n") - pl.col("m")).sum())
            .join(same.select("pair", "src", "span"), on=["pair", "src"], how="left")
        )
        if sigma > 0:
            z = ((pl.col("span") - self.fragment).abs() / sigma - 2).clip(lower_bound=0)
            cand = cand.with_columns(ins=(-0.5 * z**2).fill_null(0.0))
        else:
            cand = cand.with_columns(ins=pl.lit(0.0))
        self.cand = cand.sort("pair", "src")
        # per-source bounds, effective lengths, priors and starting points
        lengths = refs.filter(pl.col("rep"), pl.col("role") != 2).select("unit", "length")
        meta = units.join(lengths, on="unit", how="left").with_columns(
            eff=(pl.col("length") * nt + self.fragment - 2 * OVERLAP).clip(lower_bound=1.0)
        )
        if prof is not None and "coverage_zi" in prof.columns:
            meta = meta.join(
                prof.select(pl.col("unit").cast(pl.UInt32), "coverage_zi"), on="unit", how="left"
            )
        else:
            meta = meta.with_columns(coverage_zi=pl.lit(1.0))
        meta = meta.with_columns(
            mass=(pl.col("coverage_zi").fill_null(0.0) * pl.col("eff")).clip(lower_bound=0.0)
        ).with_columns(share=pl.col("mass") / pl.col("mass").sum().over("component"))
        share = meta["share"].fill_nan(0.0).fill_null(0.0).to_numpy()
        comp_of = meta["component"].replace_strict(comp_i, return_dtype=pl.Int64).to_numpy()
        pairs_c = np.bincount(pair_comp["novel"].to_numpy() - n_u - n_c, minlength=n_c).astype(
            float
        )
        unit_init = (share + 1e-3) * pairs_c[comp_of]
        self.units = meta
        self.n_u, self.n_c = n_u, n_c
        self.comps = comps
        self.init = np.r_[unit_init, 0.05 * pairs_c, 0.05 * pairs_c] + 1e-9
        self.prior = np.r_[PRIOR_PAIRS * share, np.zeros(2 * n_c)]
        floor = min(max(error, 1e-4), novel)
        self.lo = np.r_[np.full(n_u, floor), np.full(n_c, DECOY_D[0]), np.full(n_c, novel)]
        self.hi = np.r_[np.full(n_u, novel), np.full(n_c, DECOY_D[1]), np.full(n_c, 0.5)]

    def fit(self, prior: bool, novel: bool = True) -> tuple[pl.DataFrame, np.ndarray, np.ndarray]:
        """``novel`` False starts the novel classes at zero, which EM keeps them at."""
        c = self.cand
        init = self.init.copy()
        if not novel:
            init[self.n_u + self.n_c :] = 0.0
        r, n, d = _em(
            c["pair"].to_numpy(),
            c["src"].to_numpy(),
            c["diffs"].to_numpy().astype(float),
            c["match"].to_numpy().astype(float),
            c["ins"].to_numpy(),
            init,
            self.prior if prior else np.zeros_like(self.prior),
            self.lo,
            self.hi,
            self.alpha,
        )
        return c.select("pair", "src").with_columns(r=pl.Series(r)), n, d

    def rescore(self, hits: pl.DataFrame) -> None:
        """New differences after consensus updates (hit, mism): rebuild the candidates."""
        h = self.hits.drop("mism", "m").join(hits, on="hit")
        self.hits = h.with_columns(m=pl.col("mism") + pl.col("n") - pl.col("qlen"))
        best = self.hits.sort("score", descending=True).group_by("pair", "mate", "src").first()
        non_decoy = self.hits.filter(pl.col("role") != "decoy").sort("score", descending=True)
        novel_src = self.cand.filter(pl.col("src") >= self.n_u + self.n_c).select(
            "pair", novel="src"
        )
        novel_rows = (
            non_decoy.group_by("pair", "mate")
            .first()
            .join(novel_src, on="pair")
            .select("pair", "mate", src="novel", m="m")
        )
        per_src = pl.concat([best.select("pair", "mate", "src", "m"), novel_rows])
        mates = self.hits.group_by("pair", "mate").agg(
            n=pl.col("n").first(), worst=pl.col("m").max()
        )
        counts = (
            self.cand.select("pair", "src")
            .join(mates, on="pair")
            .join(per_src, on=["pair", "mate", "src"], how="left")
            .with_columns(m=pl.col("m").fill_null(pl.col("worst")))
            .group_by("pair", "src")
            .agg(diffs=pl.col("m").sum(), match=(pl.col("n") - pl.col("m")).sum())
        )
        self.cand = (
            self.cand.drop("diffs", "match").join(counts, on=["pair", "src"]).sort("pair", "src")
        )

    def hit_weights(self, fit: pl.DataFrame) -> pl.DataFrame:
        """Per alignment: the pair's responsibility for its source, spread over the mate's
        alignments to that source by score; and per mate the novel class's."""
        h = self.hits.join(fit, on=["pair", "src"], how="left").with_columns(
            share=2.0 ** (pl.col("score") - pl.col("score").max().over("pair", "mate", "src"))
        )
        h = h.with_columns(
            weight=pl.col("r").fill_null(0.0)
            * pl.col("share")
            / pl.col("share").sum().over("pair", "mate", "src")
        )
        novel = fit.filter(pl.col("src") >= self.n_u + self.n_c).select("pair", novel="r")
        return h.join(novel, on="pair", how="left").select(
            "hit", "pair", "mate", "seq", "src", "weight", novel=pl.col("novel").fill_null(0.0)
        )


def _pile(sites: pl.DataFrame, weights: pl.DataFrame) -> pl.DataFrame:
    return (
        sites.join(weights.select("hit", "weight"), on="hit")
        .group_by("seq", "col", "residue")
        .agg(pl.col("weight").sum())
    )


def _haplotypes(
    x: np.ndarray, w: np.ndarray, cons: np.ndarray, orig: np.ndarray, eps: float, alpha: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(haplotypes, mixing weights, responsibilities n x haplotypes) of reads x (n x V
    residue bytes, 0 where a read does not cover the column). Haplotype 0 is the original
    sequence, fixed (reads that match the reference always have a place); K <= ``K_MAX``
    free haplotypes beside it, seeded from the consensus and then the read furthest from
    the haplotypes so far, K by BIC."""
    if x.shape[1] == 0 or x.shape[0] == 0:
        return orig[None, :], np.ones(1), np.ones((x.shape[0], 1))
    covered = x > 0
    log_match, log_miss = np.log1p(-eps), np.log(eps / (alpha - 1))
    values = np.unique(x[covered])
    best: tuple[float, np.ndarray, np.ndarray, np.ndarray] | None = None
    n_obs = max(w.sum(), 2.0)
    for k in range(K_MAX + 1):
        haps = orig[None, :].copy()
        if k and (cons != orig).any():
            haps = np.vstack([haps, cons])
        while len(haps) < k + 1:  # seed: the read furthest from the haplotypes so far
            miss = ((x[:, None, :] != haps[None]) & covered[:, None, :]).sum(2).min(1)
            far = int(np.argmax(miss))
            if miss[far] == 0:
                break
            haps = np.vstack([haps, np.where(covered[far], x[far], cons)])
        if len(haps) < k + 1:
            break
        eta = np.full(k + 1, 1 / (k + 1))
        for _ in range(50):
            ll = np.where(
                covered[:, None, :],
                np.where(x[:, None, :] == haps[None], log_match, log_miss),
                0.0,
            ).sum(2) + np.log(eta)
            resp = np.exp(ll - ll.max(1, keepdims=True))
            resp /= resp.sum(1, keepdims=True)
            wr = w[:, None] * resp
            eta = np.maximum(wr.sum(0) / w.sum(), 1e-12)
            votes = np.stack([wr.T @ (x == v) for v in values])  # values x haplotypes x V
            new = np.where(votes.max(0) > 0, values[votes.argmax(0)], haps)
            new[0] = orig
            if (new == haps).all():
                break
            haps = new
        log_likelihood = float(
            w @ (ll.max(1) + np.log(np.exp(ll - ll.max(1, keepdims=True)).sum(1)))
        )
        bic = -2 * log_likelihood + (k + k * x.shape[1]) * np.log(n_obs)
        if best is None or bic < best[0]:
            best = (bic, haps, eta, resp)
    assert best is not None
    return best[1], best[2], best[3]


def refine(
    hits: pl.DataFrame,
    cols: pl.DataFrame,
    refs: pl.DataFrame,
    sequences: dict[str, str],
    out: Path,
    protein: bool,
    profile: str | Path | None,
    reads: Sequence[str | Path],
    min_depth: float,
    min_freq: float,
    novel: float = NOVEL,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Steps 7-9 on the crude pass's alignments (``hits``: assignments with hit ids, spans,
    ``qlen`` and ``mism``; ``cols``: target columns). Writes ``error_model.tsv``,
    ``sources.tsv``, ``haplotypes.tsv`` and ``calls.tsv``; returns per alignment
    (hit, em_posterior, em_novel, em_flag) and per unit (em_pairs, em_abundance,
    em_divergence)."""
    alpha = 20 if protein else 4
    # step 7: qualities and the error model, from the crude posteriors
    quals = read_qualities(reads, set(hits["read_name"].unique())) if reads else {}
    sites = _site_errors(cols, hits, quals, protein)
    crude = hits.filter("kept", pl.col("role") == "target").select("hit", "seq", weight="posterior")
    model = error_model(sites, crude)
    model.write_csv(out / "error_model.tsv", separator="\t")
    sites = sites.join(model.select("cycle", "qbin", err="error"), on=["cycle", "qbin"])
    error = float(sites["err"].to_numpy().mean()) if sites.height else DEFAULT_ERROR

    # step 8: the joint EM; step 9: consensus updates and re-scoring between fits
    prof = pl.read_csv(profile, separator="\t", infer_schema_length=None) if profile else None
    em = _Model(hits, refs, prof, protein, error, novel)
    offsets = np.cumsum([0, *refs["length"].to_list()])
    original = np.frombuffer("".join(sequences[s] for s in refs["seq_id"]).encode(), np.uint8)
    current = original.copy()
    start_of = dict(enumerate(offsets[:-1].tolist()))
    position = pl.col("col").cast(pl.Int64) + pl.col("seq").replace_strict(
        start_of, return_dtype=pl.Int64
    )
    target_cols = cols.join(hits.select("hit", "seq"), on="hit").with_columns(g=position)
    g = target_cols["g"].to_numpy()
    res = target_cols["residue"].to_numpy()
    hit_of = target_cols["hit"].to_numpy()
    mism = hits.select("hit", "mism").sort("hit")
    sites = sites.join(hits.select("hit", "seq"), on="hit")
    updated = 0
    for _ in range(ROUNDS):
        # without the novel class: a strain's reads must reach its consensus before the
        # final fit and the haplotypes judge novelty against the original sequence
        fit, _, _ = em.fit(prior=False, novel=False)
        weights = em.hit_weights(fit).filter(pl.col("src") < em.n_u)
        pile = _pile(sites.filter(pl.col("residue") != GAP), weights).with_columns(g=position)
        pile = pile.with_columns(orig=pl.Series(original[pile["g"].to_numpy()]))
        # a new majority, where the original residue has lost its support (not where a
        # misplaced relative outnumbers the unit's own strain)
        top = (
            pile.group_by("seq", "col", "g")
            .agg(
                consensus=pl.col("residue").sort_by("weight").last(),
                freq=pl.col("weight").max() / pl.col("weight").sum(),
                depth=pl.col("weight").sum(),
                orig_weight=pl.col("weight").filter(pl.col("residue") == pl.col("orig")).sum(),
            )
            .filter(
                pl.col("depth") >= min_depth,
                pl.col("freq") > 0.5,
                pl.col("orig_weight")
                < pl.max_horizontal(MIN_ALT_READS, min_freq * pl.col("depth")),
            )
        )
        at = top["g"].to_numpy()
        change = top["consensus"].to_numpy() != current[at]
        if not change.any():
            break
        before = current.copy()
        current[at[change]] = top["consensus"].to_numpy()[change]
        updated += int(change.sum())
        delta = (res != current[g]).astype(np.int64) - (res != before[g])
        moved = np.bincount(hit_of, delta, minlength=hits.height)
        mism = mism.with_columns(
            mism=pl.col("mism") + pl.Series(moved[mism["hit"].to_numpy()].astype(np.int32))
        )
        em.rescore(mism)
    fit, n_free, d = em.fit(prior=False)
    fit_prior, n_prior, _ = em.fit(prior=True)
    flagged = (
        _best(fit)
        .join(_best(fit_prior), on="pair", suffix="_prior")
        .select("pair", em_flag=pl.col("src") != pl.col("src_prior"))
    )
    weights = em.hit_weights(fit)

    # step 9: haplotypes per target sequence; novel ones move their reads' weight
    on_targets = (
        weights.filter(pl.col("src") < em.n_u)
        .join(hits.select("hit", "role"), on="hit")
        .filter(pl.col("role") == "target", pl.col("weight") > 1e-3)
    )
    pile = _pile(sites.filter(pl.col("residue") != GAP), on_targets)
    columns = pile.group_by("seq", "col").agg(
        depth=pl.col("weight").sum(),
        second=pl.col("weight").sort(descending=True).slice(1, 1).first().fill_null(0.0),
    )
    covered = columns.filter(pl.col("depth") >= min_depth).with_columns(g=position)
    at = covered["g"].to_numpy()
    covered = covered.with_columns(updated=pl.Series(current[at] != original[at]))
    variable = covered.filter(
        (pl.col("second") >= pl.max_horizontal(min_freq * pl.col("depth"), MIN_ALT_READS))
        | pl.col("updated")
    )
    eps = max(error, 0.01)
    # each mate's divergence from the nearest of the unit's original sequences
    nearest = (
        hits.filter(pl.col("role") != "decoy")
        .with_columns(n=pl.col("qlen").max().over("read_name", "mate"))
        .with_columns(d=(pl.col("mism") + pl.col("n") - pl.col("qlen")) / pl.col("n"))
        .with_columns(d_unit=pl.col("d").min().over("read_name", "mate", "unit"))
        .select("hit", "d_unit")
    )
    hap_rows, novel_hits = [], []
    site_res = sites.filter(pl.col("residue") != GAP).select("hit", "seq", "col", "residue")
    for (seq,), cov in covered.group_by("seq"):
        off = int(offsets[seq])
        cov_cols = np.sort(cov["col"].to_numpy())
        var_cols = np.sort(variable.filter(pl.col("seq") == seq)["col"].to_numpy())
        reads_w = (
            on_targets.filter(pl.col("seq") == seq)
            .select("hit", "weight")
            .join(nearest, on="hit", how="left", maintain_order="left")
        )
        x = np.zeros((reads_w.height, len(var_cols)), np.uint8)
        if len(var_cols):
            obs = site_res.filter(
                pl.col("seq") == seq, pl.col("col").is_in(var_cols.tolist())
            ).join(reads_w.with_row_index("i"), on="hit")
            x[obs["i"].to_numpy(), np.searchsorted(var_cols, obs["col"].to_numpy())] = obs[
                "residue"
            ].to_numpy()
        cons = current[off + var_cols]
        w = reads_w["weight"].to_numpy()
        haps, eta, resp = _haplotypes(x, w, cons, original[off + var_cols], eps, alpha)
        is_novel = np.zeros(len(haps), bool)
        d_unit = reads_w["d_unit"].fill_null(0.0).to_numpy()
        for k, hap in enumerate(haps):
            full = current[off + cov_cols].copy()
            full[np.searchsorted(cov_cols, var_cols)] = hap
            div = float((full != original[off + cov_cols]).mean())
            mass = float(w @ resp[:, k])
            # its reads, against the nearest of the unit's sequences, sequencing error removed
            read_div = float((w * resp[:, k]) @ d_unit / mass) - error if mass > 0 else 0.0
            is_novel[k] = len(cov_cols) >= MIN_COVERED and min(div, read_div) > novel
            hap_rows.append(
                (int(seq), k, float(eta[k]), mass, div, read_div, len(cov_cols), is_novel[k])
            )
        if is_novel.any():
            novel_hits.append(reads_w.select("hit", p_novel=pl.Series(resp[:, is_novel].sum(1))))
    haplotypes = pl.DataFrame(
        hap_rows,
        schema={
            "seq": pl.UInt32,
            "haplotype": pl.Int64,
            "share": pl.Float64,
            "reads": pl.Float64,
            "divergence": pl.Float64,
            "read_divergence": pl.Float64,
            "covered": pl.Int64,
            "novel": pl.Boolean,
        },
        orient="row",
    )
    refs.select("seq", "seq_id", "unit", "cluster_rep").join(haplotypes, on="seq").drop("seq").sort(
        "seq_id", "haplotype"
    ).write_csv(out / "haplotypes.tsv", separator="\t")
    to_novel = (
        pl.concat(novel_hits)
        if novel_hits
        else pl.DataFrame(schema={"hit": pl.UInt32, "p_novel": pl.Float64})
    )
    weights = (
        weights.join(to_novel, on="hit", how="left")
        .with_columns(p_novel=pl.col("p_novel").fill_null(0.0))
        .with_columns(
            gone=pl.col("weight") * pl.col("p_novel"),
            weight=pl.col("weight") * (1 - pl.col("p_novel")),
        )
        .with_columns(novel=pl.col("novel") + pl.col("gone").sum().over("pair", "mate"))
    )

    # step 7: variant tests on the final pileup
    on_targets = weights.join(hits.select("hit", "role"), on="hit").filter(
        pl.col("role") == "target"
    )
    _calls(sites, on_targets, refs, original, position, protein, min_depth, min_freq).write_csv(
        out / "calls.tsv", separator="\t"
    )

    # outputs
    kinds = ["unit"] * em.n_u + ["decoy"] * em.n_c + ["novel"] * em.n_c
    comps = em.units["component"].to_list() + em.comps.to_list() * 2
    src_units = em.units["unit"].to_list() + [None] * (2 * em.n_c)
    eff = np.r_[em.units["eff"].to_numpy(), np.full(2 * em.n_c, np.nan)]
    sources = pl.DataFrame(
        {
            "kind": kinds,
            "unit": pl.Series(src_units, dtype=pl.UInt32),
            "component": comps,
            "pairs": n_free,
            "pairs_prior": n_prior,
            "divergence": d,
            "eff_len": eff,
        }
    ).with_columns(abundance=pl.col("pairs") / pl.col("eff_len"))
    sources.write_csv(out / "sources.tsv", separator="\t")
    per_hit = weights.join(flagged, on="pair", how="left").select(
        "hit", em_posterior="weight", em_novel="novel", em_flag="em_flag"
    )
    per_unit = sources.filter(pl.col("kind") == "unit").select(
        "unit", em_pairs="pairs", em_abundance="abundance", em_divergence="divergence"
    )
    print(
        f"refine: fragment {em.fragment:.0f} nt, error {error:.4f}, {updated} consensus "
        f"updates, {int(haplotypes['novel'].sum())} novel haplotypes"
    )
    return per_hit, per_unit


def _best(fit: pl.DataFrame) -> pl.DataFrame:
    """Each pair's most responsible source."""
    return fit.sort("r", descending=True).group_by("pair").first()


def _as_char(column: str) -> pl.Expr:
    return pl.col(column).map_elements(chr, return_dtype=pl.String)


def _calls(
    sites: pl.DataFrame,
    weights: pl.DataFrame,
    refs: pl.DataFrame,
    reference: np.ndarray,
    position: pl.Expr,
    protein: bool,
    min_depth: float,
    min_freq: float,
) -> pl.DataFrame:
    """Per piled target column with an alt residue at ``min_freq``: depth, consensus, alt and
    its frequency, the error test (p, BH q), strand and read-position bias, PASS."""
    shares = substitution_shares(protein)
    rows = (
        sites.filter(pl.col("residue") != GAP)
        .join(weights.select("hit", "weight"), on="hit")
        .filter(pl.col("weight") > 0)
    )
    rows = rows.with_columns(g=position)
    rows = rows.with_columns(ref=pl.Series(reference[rows["g"].to_numpy()]))
    per = rows.group_by("seq", "col", "residue").agg(
        weight=pl.col("weight").sum(), ref=pl.col("ref").first()
    )
    cols = (
        per.group_by("seq", "col")
        .agg(
            depth=pl.col("weight").sum(),
            consensus=pl.col("residue").sort_by("weight").last(),
            ref=pl.col("ref").first(),
            alt=pl.col("residue")
            .filter(pl.col("residue") != pl.col("ref"))
            .sort_by(pl.col("weight").filter(pl.col("residue") != pl.col("ref")))
            .last(),
            alt_weight=pl.col("weight").filter(pl.col("residue") != pl.col("ref")).max(),
        )
        .join(
            rows.group_by("seq", "col").agg(expected=(pl.col("weight") * pl.col("err")).sum()),
            on=["seq", "col"],
        )
        .filter(pl.col("depth") >= min_depth)
    )
    n_tested = cols.height  # every covered column is a test (p = 1 without an alt)
    cols = cols.filter(pl.col("alt").is_not_null()).with_columns(
        alt_freq=pl.col("alt_weight") / pl.col("depth")
    )
    null = np.where(
        cols["alt"].to_numpy() == cols["consensus"].to_numpy(),
        cols["ref"].to_numpy(),
        cols["consensus"].to_numpy(),
    )
    lam = np.maximum(cols["expected"].to_numpy() * shares[null, cols["alt"].to_numpy()], 1e-12)
    p = np.clip(special.gammainc(cols["alt_weight"].to_numpy(), lam), 0, 1)  # P(Poisson >= k)
    q = np.asarray(stats.false_discovery_control(np.r_[p, np.ones(n_tested - len(p))]))[: len(p)]
    cols = cols.with_columns(pvalue=pl.Series(p), qvalue=pl.Series(q)).filter(
        pl.col("alt_freq") >= min_freq
    )
    # strand and read-position bias, per reported column
    detail = rows.join(cols.select("seq", "col", "alt"), on=["seq", "col"]).with_columns(
        is_alt=pl.col("residue") == pl.col("alt")
    )
    strand_p, pos_p = [], []
    groups = {k: g for k, g in detail.group_by("seq", "col")}
    for seq, col in cols.select("seq", "col").iter_rows():
        g = groups[(seq, col)]
        alt, other = g.filter("is_alt"), g.filter(~pl.col("is_alt"))
        table = [
            [round(part.filter(pl.col("reverse") == r)["weight"].sum()) for r in (False, True)]
            for part in (alt, other)
        ]
        strand_p.append(float(stats.fisher_exact(table).pvalue))
        pos_p.append(
            float(
                stats.mannwhitneyu(alt["from_end"].to_numpy(), other["from_end"].to_numpy()).pvalue
            )
            if alt.height >= 3 and other.height >= 3
            else None
        )
    return (
        cols.with_columns(
            strand_p=pl.Series(strand_p, dtype=pl.Float64),
            position_p=pl.Series(pos_p, dtype=pl.Float64),
        )
        .with_columns(
            passed=(pl.col("qvalue") <= 0.01)
            & (pl.col("strand_p") >= 1e-3)
            & (pl.col("position_p").fill_null(1.0) >= 1e-3)
        )
        .join(refs.select("seq", "seq_id", "unit", "cluster_rep", "rep"), on="seq")
        .select(
            "seq_id",
            "unit",
            "cluster_rep",
            "rep",
            pos=pl.col("col") + 1,
            ref=_as_char("ref"),
            depth="depth",
            consensus=_as_char("consensus"),
            alt=_as_char("alt"),
            alt_freq="alt_freq",
            pvalue="pvalue",
            qvalue="qvalue",
            strand_p="strand_p",
            position_p="position_p",
            passed="passed",
        )
        .sort("seq_id", "pos")
    )
