"""Read pairs to MGnify90 units through their HMM columns (docs/pfam-anchor-plan.md).

``oracle`` (phase 1): no index. Each read pair's hits (every family, ``-E 10``) are scored
against every GA domain of the same family in the reference (``prepare.py``'s
``domains.parquet``), by identity over the HMM columns both have:

- per hit and unit, the unit's best domain, if it shares at least ``MIN_SHARE`` of the hit's
  columns;
- per pair, family and unit, matches summed over the pair's hits, over the shared columns
  (``--score shared``) or over all the pair's columns in that family (``read``), if the shared
  columns are at least ``MIN_SHARE`` of the pair's;
- the pair's call: the units with the best identity over its families, ties kept.

Beside it, DIAMOND's call (``prepare.py diamond``, ``--top 10``): per unit the sum over mates
of the best bitscore, best units tied. *Restricted* keeps only units with a GA domain in a
family the pair hits, the oracle's universe; *all* keeps every unit.

Alignment agreement (phase 1, step 3): where a mate's ORF hit and DIAMOND both align it to
the truth protein (strains: the source; relatives: the nearest remaining relative), the
read-to-protein residue pairs the two HMM alignments imply are checked against DIAMOND's.

Writes per read set ``oracle.parquet`` (one row per pair with a hit), and ``oracle.tsv`` and
``alignment.tsv`` under ``--data``, by read set and identity bin.
"""

import argparse
import gzip
import itertools
import json
import math
import re
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import polars as pl
from Bio.Align import substitution_matrices
from numpy.typing import NDArray

sys.path.insert(0, str(Path(__file__).resolve().parent))
import columns  # noqa: E402
import prepare  # noqa: E402

MIN_SHARE = 0.6  # phase 1b's MIN_OVERLAP
SCALE = 10_000  # matches and shared columns per hit stay below this: packed into one int64
EPS = 1e-9
# identity over shared columns (phase 1b), or over all the read's columns
SCORES = ("shared", "read")
OP = re.compile(r"(\d+)([MID])")


class Family:
    """One family's reference GA domains on its HMM columns, rows grouped by unit."""

    def __init__(self, length: int, domains: pl.DataFrame) -> None:
        domains = domains.sort("unit")
        self.matrix = np.stack(
            [
                columns.column_residues(*r, length)[0]
                for r in domains.select("hmm_from", "ali_from", "cigar", "sequence").rows()
            ]
        )
        unit = domains["unit"].to_numpy()
        self.starts = np.flatnonzero(np.r_[True, unit[1:] != unit[:-1]])
        self.units = unit[self.starts]  # sorted
        self.proteins = (
            domains["protein_id"].to_numpy() if "protein_id" in domains.columns else None
        )

    def score(self, codes: NDArray[np.uint8], by: str) -> tuple[NDArray, NDArray, int]:
        """Per unit, matches and shared columns of its best domain against one hit's columns
        (0, 0 when no domain shares ``MIN_SHARE`` of them), and the hit's column count. Best is
        by identity over shared columns (``by="shared"``) or by matches (``"read"``)."""
        idx = np.flatnonzero(codes)
        c = codes[idx]
        sub = self.matrix[:, idx]
        shared = (sub > 0).sum(1, dtype=np.int64)
        match = ((sub == c) & (c != columns.OTHER)).sum(1, dtype=np.int64)
        ident = match * 1_000_000 // np.maximum(shared, 1) if by == "shared" else match
        key = np.where(
            shared >= MIN_SHARE * len(idx), ident * SCALE**2 + shared * SCALE + match, -1
        )
        best = np.maximum.reduceat(key, self.starts)
        ok = best >= 0
        return np.where(ok, best % SCALE, 0), np.where(ok, best // SCALE % SCALE, 0), len(idx)


def lengths(hmm: Path) -> dict[int, int]:
    return {acc: e.shape[0] - 1 for acc, e in columns.read_hmms(str(hmm)).items()}


def weak_bar(data: Path, fdr: float = 0.05) -> float:
    """The lowest domain score at which reversed reference proteins' domains are at most ``fdr``
    of the real domains below GA scoring as high (``prepare.py nullref``; pfam-reannotate's
    null arm), counted over all mini-slice proteins."""
    real = pl.read_parquet(data / "hits.parquet").filter(~pl.col("above_ga"))["score"].to_numpy()
    null = np.sort(pl.read_parquet(data / "null_hits.parquet")["score"].to_numpy())
    for s in np.unique(np.round(real, 1)):
        n_real = (real >= s).sum()
        if n_real and (len(null) - np.searchsorted(null, s)) / n_real <= fdr:
            return float(s)
    return float("inf")


def load_families(data: Path, ref_arm: str = "ga") -> tuple[dict[int, Family], pl.DataFrame]:
    """Families, and the reference domains with unit and sequence: ``ga`` (GA domains), ``weak``
    (also those below GA scoring at least :func:`weak_bar`) or ``reps`` (GA domains of
    representatives only)."""
    ref = pl.read_parquet(
        data / "holdout/members.parquet", columns=["protein_id", "cluster_rep", "sequence"]
    )
    keep = pl.col("above_ga")
    if ref_arm == "weak":
        keep = keep | (pl.col("score") >= weak_bar(data))
    domains = (
        pl.read_parquet(data / "domains.parquet")
        .filter(keep)
        .join(ref.rename({"cluster_rep": "unit"}), on="protein_id")
    )
    if ref_arm == "reps":
        domains = domains.filter(pl.col("protein_id") == pl.col("unit"))
    size = lengths(data / "hmms.hmm")
    fams = {acc: Family(size[acc], g) for (acc,), g in domains.group_by(["pfam_accession"])}
    return fams, domains


def truth_of(data: Path, name: str) -> pl.DataFrame:
    """Per read pair: identity bin, truth unit, the source's GA families and the truth protein."""
    truth = pl.read_parquet(data / name / "truth.parquet")
    if name == "relatives":
        rel = pl.read_csv(data / "holdout/relatives.tsv", separator="\t").select(
            "protein_id", "relative"
        )
        truth = truth.join(rel, on="protein_id")
    else:
        truth = truth.with_columns(relative=pl.col("protein_id"))
    return (
        pl.read_parquet(data / name / "reads.parquet", columns=["read", "strain", "pair_bases"])
        .join(truth, on="strain")
        .select(
            "read",
            "identity_bin",
            "unit",
            "protein_id",
            "relative",
            read_len=pl.col("pair_bases") // 2,
        )
    )


def read_hits(data: Path, name: str) -> pl.DataFrame:
    orfs = pl.read_parquet(data / name / "orfs.parquet")
    return (
        pl.read_parquet(data / name / "hits.parquet")
        .rename({"protein_id": "orf_id"})
        .join(orfs, on="orf_id")
        .sort("read", "pfam_accession", "orf_id", "domain")
    )


def oracle_calls(
    hits: pl.DataFrame, fams: dict[int, Family], truth: dict[int, int], by: str
) -> pl.DataFrame:
    rows = []
    for (read,), pair in hits.group_by(["read"], maintain_order=True):
        target = truth[read]
        best, tied, truth_ident, families = -1.0, set(), -1.0, 0
        for (acc,), g in pair.group_by(["pfam_accession"], maintain_order=True):
            fam = fams.get(acc)
            if fam is None:  # no reference GA domain in this family
                continue
            families += 1
            M = np.zeros(len(fam.units), dtype=np.int64)
            S = np.zeros(len(fam.units), dtype=np.int64)
            N = 0
            for r in g.select("hmm_from", "ali_from", "cigar", "sequence").rows():
                codes = columns.column_residues(*r, fam.matrix.shape[1] - 1)[0]
                m, s, n = fam.score(codes, by)
                M, S, N = M + m, S + s, N + n
            ok = (S >= MIN_SHARE * N) & (S > 0)
            if not ok.any():
                continue
            ident = np.where(ok, M / (np.maximum(S, 1) if by == "shared" else N), -1.0)
            top = ident.max()
            if top > best + EPS:
                best, tied = top, set()
            if top >= best - EPS:
                tied |= set(fam.units[ident >= top - EPS].tolist())
            i = np.searchsorted(fam.units, target)
            if i < len(fam.units) and fam.units[i] == target:
                truth_ident = max(truth_ident, ident[i])
        rows.append(
            (read, families, best, len(tied), target in tied, truth_ident >= 0, truth_ident)
        )
    schema = ["read", "families", "best_identity", "ties", "right", "candidate", "truth_identity"]
    return pl.DataFrame(rows, schema=schema, orient="row").with_columns(
        pl.col("best_identity", "truth_identity").replace(-1.0, None)
    )


def diamond_calls(
    dia: pl.DataFrame, truth: pl.DataFrame, keep: pl.DataFrame | None
) -> pl.DataFrame:
    """Per pair: whether the truth unit is among the units of best summed mate bitscore.
    ``keep``: (read, unit) allowed, else every unit."""
    if keep is not None:
        dia = dia.join(keep, on=["read", "unit"], how="semi")
    per_unit = (
        dia.group_by("read", "mate", "unit").agg(pl.col("bitscore").max())
        .group_by("read", "unit").agg(pl.col("bitscore").sum())
    )  # fmt: skip
    best = per_unit.filter(pl.col("bitscore") == pl.col("bitscore").max().over("read"))
    return (
        best.join(truth.select("read", truth_unit="unit"), on="read")
        .group_by("read")
        .agg(right=(pl.col("unit") == pl.col("truth_unit")).any(), ties=pl.len())
    )


def diamond_residues(qstart: int, qframe: int, sstart: int, cigar: str) -> dict[int, int]:
    """Each aligned codon's first base on its strand (1-based, forward coordinates) -> the
    subject residue (0-based) DIAMOND pairs it with."""
    sign = 1 if qframe > 0 else -1
    out, q, j = {}, 0, sstart - 1
    for n, op in OP.findall(cigar):
        n = int(n)
        if op == "M":
            for k in range(n):
                out[qstart + sign * 3 * (q + k)] = j + k
            q, j = q + n, j + n
        elif op == "I":
            q += n
        else:
            j += n
    return out


def agreement(
    hits: pl.DataFrame, truth: pl.DataFrame, dia: pl.DataFrame, domains: pl.DataFrame, size: dict
) -> pl.DataFrame:
    """Per (ORF hit, truth protein) both methods align: residues compared and residues placed
    differently, against the truth protein's domain in that family that agrees best."""
    best = (
        dia.sort("bitscore", descending=True)
        .unique(["read", "mate", "protein_id"], keep="first")
        .select(
            "read",
            "mate",
            relative="protein_id",
            qstart="qstart",
            qframe="qframe",
            sstart="sstart",
            dcigar="cigar",
        )
    )
    prot = domains.select(
        relative="protein_id", pfam_accession="pfam_accession", p_from="hmm_from", p_ali="ali_from",
        p_cigar="cigar", p_seq="sequence",
    )  # fmt: skip
    pairs = (
        hits.join(truth.select("read", "relative", "identity_bin", "read_len"), on="read")
        .join(best, on=["read", "mate", "relative"])
        .join(prot, on=["relative", "pfam_accession"])
    )
    rows: dict[tuple[int, int], tuple] = {}
    for r in pairs.iter_rows(named=True):
        length = size[r["pfam_accession"]]
        _, orf_pos = columns.column_residues(
            r["hmm_from"], r["ali_from"], r["cigar"], r["sequence"], length
        )
        _, p_pos = columns.column_residues(
            r["p_from"], r["p_ali"], r["p_cigar"], r["p_seq"], length
        )
        both = (orf_pos >= 0) & (p_pos >= 0)
        i, j = orf_pos[both], p_pos[both]
        # each ORF residue's codon: first base on its strand, 1-based forward coordinates
        plus = r["strand"] == "+"
        key = r["start"] + 3 * i + 1 if plus else r["read_len"] - r["start"] - 3 * i
        d = diamond_residues(r["qstart"], r["qframe"], r["sstart"], r["dcigar"])
        compared = differ = 0
        for k, jj in zip(key.tolist(), j.tolist(), strict=True):
            if k in d:
                compared += 1
                differ += d[k] != jj
        k = (r["orf_id"], r["pfam_accession"])
        if k not in rows or (compared - differ) > (rows[k][2] - rows[k][3]):
            rows[k] = (r["read"], r["identity_bin"], compared, differ)
    return pl.DataFrame(
        list(rows.values()), schema=["read", "identity_bin", "compared", "differ"], orient="row"
    )


def oracle(a: argparse.Namespace) -> None:
    fams, domains = load_families(a.data)
    size = lengths(a.data / "hmms.hmm")
    unit_fams = domains.select("unit", "pfam_accession").unique()
    ga_all = (
        pl.read_parquet(a.data / "hits.parquet")
        .filter("above_ga")
        .select("protein_id", "pfam_accession")
        .unique()
    )
    unit_of = domains.select("protein_id", "unit").unique()
    summary, align = [], []
    for name in prepare.SETS:
        truth = truth_of(a.data, name)
        hits = read_hits(a.data, name)
        target = dict(truth.select("read", "unit").iter_rows())
        calls = pl.concat(
            oracle_calls(hits, fams, target, by).with_columns(score=pl.lit(by)) for by in a.score
        )
        # the pair hits a family of its source protein (any score; at GA)
        src = (
            hits.select("read", "pfam_accession", "above_ga")
            .join(truth.select("read", "protein_id"), on="read")
            .join(ga_all, on=["protein_id", "pfam_accession"])
            .group_by("read")
            .agg(hit_source=pl.lit(True), ga_source=pl.col("above_ga").any())
        )
        ref = pl.read_parquet(
            a.data / "holdout/members.parquet", columns=["protein_id", "cluster_rep"]
        )
        dia = pl.read_parquet(a.data / name / "diamond.parquet").join(
            ref.rename({"cluster_rep": "unit"}), on="protein_id"
        )
        pair_fams = hits.select("read", "pfam_accession").unique()
        keep = pair_fams.join(unit_fams, on="pfam_accession").select("read", "unit").unique()
        restricted = diamond_calls(dia, truth, keep).rename(
            {"right": "dia_right", "ties": "dia_ties"}
        )
        every = diamond_calls(dia, truth, None).rename(
            {"right": "dia_all_right", "ties": "dia_all_ties"}
        )
        calls = (
            calls.join(truth.select("read", "identity_bin"), on="read")
            .join(src, on="read", how="left")
            .join(restricted, on="read", how="left")
            .join(every, on="read", how="left")
            .with_columns(
                pl.col("hit_source", "ga_source", "dia_right", "dia_all_right").fill_null(False)
            )
        )
        calls.write_parquet(a.data / name / "oracle.parquet")
        on_source = calls.filter("hit_source")
        summary.append(
            on_source.group_by("score", "identity_bin")
            .agg(
                pairs=pl.len(),
                candidate=pl.col("candidate").mean(),
                right=pl.col("right").mean(),
                credit=(pl.col("right") / pl.col("ties")).mean(),
                ties_median=pl.col("ties").median(),
                dia_right=pl.col("dia_right").mean(),
                dia_ties_median=pl.col("dia_ties").median(),
                dia_credit=(pl.col("dia_right") / pl.col("dia_ties")).fill_null(0).mean(),
                ties_mean=pl.col("ties").mean(),
                dia_ties_mean=pl.col("dia_ties").mean(),
                dia_all_right=pl.col("dia_all_right").mean(),
                right_ga=pl.col("right").filter("ga_source").mean(),
                dia_right_ga=pl.col("dia_right").filter("ga_source").mean(),
            )
            .with_columns(set=pl.lit(name), gap=pl.col("dia_right") - pl.col("right"))
        )
        agree = agreement(hits, truth, dia, domains.join(unit_of, on=["protein_id", "unit"]), size)
        align.append(
            agree.group_by("identity_bin")
            .agg(
                hits=pl.len(),
                no_overlap=(pl.col("compared") == 0).mean(),
                residues=pl.col("compared").sum(),
                differ=pl.col("differ").sum() / pl.col("compared").sum(),
                hits_all_agree=((pl.col("differ") == 0) & (pl.col("compared") > 0)).mean(),
            )
            .with_columns(set=pl.lit(name))
        )
    for rows, out in ((summary, "oracle.tsv"), (align, "alignment.tsv")):
        df = pl.concat(rows)
        keys = [k for k in ("set", "score", "identity_bin") if k in df.columns]
        df = df.select(*keys, pl.exclude(keys)).sort(
            keys, descending=[*[False] * (len(keys) - 1), True]
        )
        df.write_csv(a.data / out, separator="\t", float_precision=3)
        with pl.Config(tbl_rows=-1, tbl_cols=-1, float_precision=3, tbl_width_chars=200):
            print(df)


# Phase 2: column seeds. --------------------------------------------------------------------------

MURPHY10 = "LVIM C A G ST P FYW EDNQ KR H"  # tree-index alphabets.py's murphy10
SAMPLING = (1, 2, 4)  # index keeps seeds whose first column is a multiple of s
MIN_SEEDS = (1, 2)  # matching seeds a unit needs to be a candidate


def alphabet(groups: str | None) -> NDArray[np.int64]:
    """Residue code (0-21) -> letter of the seed alphabet; -1 for no residue and OTHER."""
    out = np.full(22, -1, dtype=np.int64)
    if groups is None:
        out[1:21] = np.arange(20)
    else:
        for g, letters in enumerate(groups.split()):
            for x in letters:
                out[columns.AMINO.index(x) + 1] = g
    return out


def shapes(length: int, ref: NDArray[np.uint8]) -> dict[str, NDArray[np.int64]]:
    """Per shape, its seeds: (seeds x width) HMM columns. Anchored shapes have one seed per
    first column; greedy ones are phase 1b's (``columns.greedy_seeds``)."""
    out = {}
    for name, pattern in (("k5", "11111"), ("k6", "111111"), ("k8", "11111111"),
                          ("s1101011", "1101011"), ("s11010010011", "11010010011")):  # fmt: skip
        offsets = np.flatnonzero(np.frombuffer(pattern.encode(), dtype=np.uint8) == 49)
        out[name] = np.arange(1, length - len(pattern) + 2)[:, None] + offsets
    unit_view = SimpleNamespace(ref=ref, length=length)
    unit_view.pairs = columns.Family.nearest_pairs(unit_view)
    for k in (3, 4, 5):
        seeds = columns.greedy_seeds(unit_view, 12, k)
        out[f"greedy_k{k}"] = np.array(seeds, dtype=np.int64).reshape(-1, k)
    return out


ARMS = (  # (arm, shape, alphabet); greedy seeds are not anchored, so only s = 1
    *((s, s, None) for s in ("k5", "k6", "k8", "s1101011", "s11010010011")),
    *((f"greedy_k{k}", f"greedy_k{k}", None) for k in (3, 4, 5)),
    ("murphy_k6", "k6", MURPHY10),
    ("murphy_k8", "k8", MURPHY10),
)


class SeedIndex:
    """One family's postings for one shape and alphabet: sorted ``code << bits | unit``, with
    ``code = seed * A^width + residues``; ``s`` keeps seeds whose first column divides by it.
    ``anchored=False`` (the ablation): ``code = residues``, so a seed matches at any column of
    the family; sampling by ``s`` is then fixed at build, and queries use every seed."""

    def __init__(
        self,
        fam: Family,
        seeds: NDArray[np.int64],
        letters: NDArray[np.int64],
        anchored: bool = True,
        s: int = 1,
    ) -> None:
        self.seeds, self.letters, self.anchored = seeds, letters, anchored
        self.base = int(letters.max()) + 1
        self.bits = max(1, int(len(fam.units) - 1).bit_length())
        unit_of_row = np.repeat(
            np.arange(len(fam.units)), np.diff(np.r_[fam.starts, len(fam.matrix)])
        )
        parts = []
        for lo in range(
            0, len(fam.matrix), 2000
        ):  # ponytail: row chunks bound the seeds x rows block
            code = self.codes(fam.matrix[lo : lo + 2000])
            ok = (code >= 0) & (seeds[:, 0] % s == 0)
            parts.append((code[ok] << self.bits) | unit_of_row[lo : lo + 2000][np.nonzero(ok)[0]])
        self.post = np.unique(np.concatenate(parts)) if parts else np.zeros(0, dtype=np.int64)

    def codes(self, rows: NDArray[np.uint8]) -> NDArray[np.int64]:
        """(rows x seeds) codes, -1 where a seed column has no residue."""
        x = self.letters[rows[:, self.seeds]]  # rows x seeds x width
        packed = (x * self.base ** np.arange(self.seeds.shape[1])).sum(2)
        if self.anchored:
            packed += np.arange(len(self.seeds)) * self.base ** self.seeds.shape[1]
        return np.where((x >= 0).all(2), packed, -1)

    def sampled(self, s: int) -> NDArray[np.int64]:
        if s == 1 or not self.anchored:
            return self.post
        seed = (self.post >> self.bits) // self.base ** self.seeds.shape[1]
        return self.post[self.seeds[seed, 0] % s == 0]

    def counts(self, post: NDArray[np.int64], query: NDArray[np.uint8], n_units: int) -> NDArray:
        """Per unit, how many of the query's distinct seed codes it has."""
        q = np.unique(self.codes(query))
        q = q[q >= 0]
        lo = np.searchsorted(post, q << self.bits)
        hi = np.searchsorted(post, (q + 1) << self.bits)
        n = hi - lo
        if not n.sum():
            return np.zeros(n_units, dtype=np.int64)
        idx = np.repeat(lo - np.r_[0, np.cumsum(n)[:-1]], n) + np.arange(n.sum())
        return np.bincount(post[idx] & ((1 << self.bits) - 1), minlength=n_units)


def seeds(a: argparse.Namespace) -> None:
    fams, _ = load_families(a.data)
    sets = []
    for name in prepare.SETS:
        oracle_right = pl.read_parquet(a.data / name / "oracle.parquet").filter(
            "hit_source", pl.col("score") == "read"
        )
        hits = read_hits(a.data, name).join(oracle_right.select("read"), on="read")
        truth = truth_of(a.data, name)
        sets.append(
            hits.join(truth.select("read", "unit"), on="read").with_columns(set=pl.lit(name))
        )
    hits = pl.concat(sets)
    rows, speed, sizes = [], [], []
    for (acc,), fam_hits in hits.group_by(["pfam_accession"], maintain_order=True):
        fam = fams.get(acc)
        if fam is None:
            continue
        t0 = time.perf_counter()
        length = fam.matrix.shape[1] - 1
        shape = shapes(length, fam.matrix[fam.starts])
        letters = {None: alphabet(None), MURPHY10: alphabet(MURPHY10)}
        indexes = {arm: SeedIndex(fam, shape[sh], letters[al], a.anchored) for arm, sh, al in ARMS}
        build = time.perf_counter() - t0
        rows_per_unit = np.diff(np.r_[fam.starts, len(fam.matrix)])
        posts = {
            (arm, s): indexes[arm].sampled(s)
            if a.anchored
            else SeedIndex(fam, shape[sh], letters[al], False, s).post
            for arm, sh, al in ARMS
            for s in (SAMPLING if not arm.startswith("greedy") else (1,))
        }
        sizes += [(acc, arm, s, len(p), len(fam.units)) for (arm, s), p in posts.items()]
        brute = query = 0.0
        n_hits = 0
        cand_rows: dict[tuple, int] = {}
        for (st, read), pair in fam_hits.group_by(["set", "read"], maintain_order=True):
            target = pair["unit"][0]
            hits_codes = [
                columns.column_residues(*r, length)[0]
                for r in pair.select("hmm_from", "ali_from", "cigar", "sequence").rows()
            ]
            t1 = time.perf_counter()
            M = np.zeros(len(fam.units), dtype=np.int64)
            S = np.zeros(len(fam.units), dtype=np.int64)
            N = 0
            for codes in hits_codes:
                m, s_, n = fam.score(codes, "read")
                M, S, N = M + m, S + s_, N + n
            brute += time.perf_counter() - t1
            n_hits += len(hits_codes)
            ok = (S >= MIN_SHARE * N) & (S > 0)
            ident = np.where(ok, M / max(N, 1), -1.0)
            i = np.searchsorted(fam.units, target)
            has_truth = i < len(fam.units) and fam.units[i] == target
            arms = [("oracle", 1, 1, np.ones(len(fam.units), dtype=bool))]
            query_rows = np.stack(hits_codes)  # distinct seeds over all the pair's hits
            for (arm, s), post in posts.items():
                t2 = time.perf_counter()
                count = indexes[arm].counts(post, query_rows, len(fam.units))
                query += time.perf_counter() - t2
                for h in MIN_SEEDS:
                    cand = count >= h
                    arms.append((arm, s, h, cand))
                    key = (arm, s, h)
                    cand_rows[key] = cand_rows.get(key, 0) + int(rows_per_unit[cand].sum())
            for arm, s, h, cand in arms:
                okc = ok & cand
                best = ident[okc].max() if okc.any() else -1.0
                tied = okc & (ident >= best - EPS)
                rows.append((st, read, acc, arm, s, h, best, bool(has_truth and tied[i]),
                             int(tied.sum()), int(cand.sum()), len(fam.units)))  # fmt: skip
        for (arm, s, h), n_rows in cand_rows.items():
            speed.append(
                (
                    acc,
                    len(fam.units),
                    len(fam.matrix),
                    n_hits,
                    brute,
                    query,
                    build,
                    arm,
                    s,
                    h,
                    n_rows,
                )
            )
        print(
            f"PF{acc:05d}: {len(fam.units)} units, {n_hits} hits, build {build:.1f} s", flush=True
        )
    cols = ["set", "read", "pfam_accession", "arm", "s", "h", "best", "truth_tied", "ties",
            "candidates", "units"]  # fmt: skip
    per_fam = pl.DataFrame(rows, schema=cols, orient="row")
    per_pair = (
        per_fam.with_columns(top=pl.col("best").max().over("set", "read", "arm", "s", "h"))
        .group_by("set", "read", "arm", "s", "h")
        .agg(
            right=(pl.col("truth_tied") & (pl.col("best") >= pl.col("top") - EPS)).any(),
            ties=pl.col("ties").filter(pl.col("best") >= pl.col("top") - EPS).sum(),
            share=pl.col("candidates").sum() / pl.col("units").sum(),
            units=pl.col("units").sum(),
        )
    )
    tag = "" if a.anchored else "_unanchored"
    per_pair.write_parquet(a.data / f"seeds{tag}.parquet")
    oracle = per_pair.filter(pl.col("arm") == "oracle").select("set", "read", oracle_right="right")
    index = (
        pl.DataFrame(
            sizes, schema=["pfam_accession", "arm", "s", "postings", "units"], orient="row"
        )
        .group_by("arm", "s")
        .agg(postings_per_unit=pl.col("postings").sum() / pl.col("units").sum())
    )
    big = pl.col("units") >= a.large
    table = (
        per_pair.filter(pl.col("arm") != "oracle")
        .join(oracle, on=["set", "read"])
        .group_by("arm", "s", "h")
        .agg(
            *(
                expr.alias(f"{st}_{label}")
                for st in prepare.SETS
                for label, expr in (
                    ("kept", pl.col("right").filter(pl.col("set") == st, "oracle_right").mean()),
                    ("share", pl.col("share").filter(pl.col("set") == st).median()),
                    (
                        "kept_large",
                        pl.col("right").filter(pl.col("set") == st, "oracle_right", big).mean(),
                    ),
                    ("share_large", pl.col("share").filter(pl.col("set") == st, big).median()),
                )
            )
        )
        .join(index, on=["arm", "s"], how="left")
        .sort("arm", "s", "h")
    )
    table.write_csv(a.data / f"seeds{tag}.tsv", separator="\t", float_precision=4)
    sp = pl.DataFrame(
        speed,
        schema=[
            "pfam_accession",
            "units",
            "rows",
            "hits",
            "brute_s",
            "query_s",
            "build_s",
            "arm",
            "s",
            "h",
            "candidate_rows",
        ],  # fmt: skip
        orient="row",
    )
    sp.write_csv(a.data / f"seed_speed{tag}.tsv", separator="\t", float_precision=6)
    with pl.Config(tbl_rows=-1, tbl_cols=-1, float_precision=3, tbl_width_chars=250):
        print(table)


# Anchor residues: the n most informative columns a read covers, as a first check. ----------

ANCHOR_N = (8, 12, 16, 24)  # columns per read pair and family
ANCHOR_SHARES = (0.005, 0.01, 0.02, 0.05, 0.1)  # units kept, by matches at the anchors


def column_ranks(fam: Family, rng) -> dict[str, NDArray[np.float64]]:
    """Per column, how informative it is, from one domain per unit: entropy across units, that
    per bit of retention lost between nearest units (phase 1b's "efficient"), or random."""
    ref = fam.matrix[fam.starts]
    counts = np.stack([(ref == x).sum(0) for x in range(1, 21)], 1).astype(float) + 0.1
    p = counts / counts.sum(1, keepdims=True)
    entropy = -(p * np.log2(p)).sum(1)
    view = SimpleNamespace(ref=ref, length=ref.shape[1] - 1)
    same, seen = columns.Family.nearest_pairs(view)
    retention = (same.sum(0) + 1) / (seen.sum(0) + 2)
    lost = np.maximum(-np.log2(np.clip(retention, 1e-3, 0.999)), 0.02)
    return {"entropy": entropy, "efficient": entropy / lost, "random": rng.random(len(entropy))}


def anchors(a: argparse.Namespace) -> None:
    rng = np.random.default_rng(1)
    fams, _ = load_families(a.data)
    hits = pl.concat(
        read_hits(a.data, name)
        .join(
            pl.read_parquet(a.data / name / "oracle.parquet")
            .filter("hit_source", pl.col("score") == "read")
            .select("read"),
            on="read",
        )
        .join(truth_of(a.data, name).select("read", "unit"), on="read")
        .with_columns(set=pl.lit(name))
        for name in prepare.SETS
    )
    rows = []
    for (acc,), fam_hits in hits.group_by(["pfam_accession"], maintain_order=True):
        fam = fams.get(acc)
        if fam is None:
            continue
        length, n_units = fam.matrix.shape[1] - 1, len(fam.units)
        ranks = column_ranks(fam, rng)
        for (st, read), pair in fam_hits.group_by(["set", "read"], maintain_order=True):
            target = pair["unit"][0]
            hits_codes = [
                columns.column_residues(*r, length)[0]
                for r in pair.select("hmm_from", "ali_from", "cigar", "sequence").rows()
            ]
            M = np.zeros(n_units, dtype=np.int64)
            S = np.zeros(n_units, dtype=np.int64)
            N = 0
            for codes in hits_codes:
                m, s_, n = fam.score(codes, "read")
                M, S, N = M + m, S + s_, N + n
            ok = (S >= MIN_SHARE * N) & (S > 0)
            ident = np.where(ok, M / max(N, 1), -1.0)
            top = ident.max()
            i = np.searchsorted(fam.units, target)
            if top < 0 or i >= n_units or fam.units[i] != target:
                rows.append((st, read, acc, None, None, top, False, None, n_units))
                continue
            truth_tied = bool(ok[i] and ident[i] >= top - EPS)
            merged = np.max(np.stack(hits_codes), 0)  # the pair's residue per column
            covered = np.flatnonzero(merged)
            for rank_name, rank in ranks.items():
                order = covered[np.argsort(-rank[covered], kind="stable")]
                for n in ANCHOR_N:
                    cols = order[:n]
                    match = (fam.matrix[:, cols] == merged[cols]).sum(1)
                    per_unit = np.maximum.reduceat(match, fam.starts)
                    above = np.sort(per_unit)[::-1]
                    for q in ANCHOR_SHARES:
                        cut = above[min(n_units, max(1, int(np.ceil(q * n_units)))) - 1]
                        kept = per_unit >= cut  # ties at the cut are kept
                        hit = truth_tied and bool(kept[i])
                        arm = f"{rank_name}_n{n}"
                        rows.append((st, read, acc, arm, q, top, hit, int(kept.sum()), n_units))
        print(f"PF{acc:05d}: {n_units} units", flush=True)
    per_fam = pl.DataFrame(
        rows,
        schema=[
            "set",
            "read",
            "pfam_accession",
            "arm",
            "q",
            "top",
            "truth_kept",
            "candidates",
            "units",
        ],
        orient="row",
    )
    pair_top = per_fam.group_by("set", "read").agg(best=pl.col("top").max())
    arms = per_fam.filter(pl.col("arm").is_not_null()).join(pair_top, on=["set", "read"])
    table = (
        arms.group_by("set", "read", "arm", "q")
        .agg(
            kept=(pl.col("truth_kept") & (pl.col("top") >= pl.col("best") - EPS)).any(),
            share=pl.col("candidates").sum() / pl.col("units").sum(),
        )
        .join(
            pl.concat(
                pl.read_parquet(a.data / name / "oracle.parquet")
                .filter("hit_source", "right", pl.col("score") == "read")
                .select("read", set=pl.lit(name))
                for name in prepare.SETS
            ),
            on=["set", "read"],
        )
        .group_by("arm", "q")
        .agg(
            *(
                expr.alias(f"{st}_{label}")
                for st in prepare.SETS
                for label, expr in (
                    ("kept", pl.col("kept").filter(pl.col("set") == st).mean()),
                    ("share", pl.col("share").filter(pl.col("set") == st).median()),
                )
            )
        )
        .sort("arm", "q")
    )
    table.write_csv(a.data / "anchors.tsv", separator="\t", float_precision=4)
    with pl.Config(tbl_rows=-1, tbl_cols=-1, float_precision=3, tbl_width_chars=200):
        print(table)


# The anchored route (phases 3-4): anchor residues, then column identity on the units kept. --

AMBIGUOUS = 0.02  # identity below the best that still makes a unit a neighbour
STATUS_ORDER = {"ga": 0, "one": 1, "set": 2, "weak": 3}  # a call's status: its strongest hit's
HIT_SOURCES = ("hmmsearch", "ga", "recruited", "recruited_top")
REF_ARMS = ("ga", "weak", "reps")


class Route:
    """One family, laid out for both stages: column-major rows for the anchor stage, row ranges
    per unit for the identity stage, and the columns ranked by entropy across units."""

    def __init__(self, fam: Family) -> None:
        self.fam = fam
        self.cols = np.ascontiguousarray(fam.matrix.T)  # columns x rows
        self.ends = np.r_[fam.starts[1:], len(fam.matrix)]
        ref = fam.matrix[fam.starts]
        counts = np.stack([(ref == x).sum(0) for x in range(1, 21)], 1) + 0.1
        p = counts / counts.sum(1, keepdims=True)
        self.entropy = -(p * np.log2(p)).sum(1)

    def candidates(self, merged: NDArray[np.uint8], n: int, keep: float) -> NDArray[np.int64]:
        """Stage 1: units whose best domain matches the read at its ``n`` highest-entropy
        columns at least as often as the top ``keep`` of units (ties at the cut kept)."""
        covered = np.flatnonzero(merged)
        anchors = covered[np.argsort(-self.entropy[covered], kind="stable")[:n]]
        count = np.zeros(self.cols.shape[1], dtype=np.uint8)
        for c in anchors:
            count += self.cols[c] == merged[c]
        per_unit = np.maximum.reduceat(count, self.fam.starts)
        k = max(1, int(np.ceil(keep * len(per_unit))))
        cut = np.partition(per_unit, len(per_unit) - k)[len(per_unit) - k]
        return np.flatnonzero(per_unit >= cut)

    def identity(self, cand: NDArray[np.int64], hits_codes: list) -> tuple:
        """Stage 2: per candidate unit, matches over the read's columns summed over its hits
        (each hit on the unit's best domain), whether it shares ``MIN_SHARE`` of them, the
        read's column count, and the row of the unit's domain with the most matches."""
        lengths = self.ends[cand] - self.fam.starts[cand]
        offsets = np.r_[0, np.cumsum(lengths)[:-1]]
        rows = np.repeat(self.fam.starts[cand] - offsets, lengths) + np.arange(lengths.sum())
        block = self.fam.matrix[rows]
        M = np.zeros(len(cand), dtype=np.int64)
        S = np.zeros(len(cand), dtype=np.int64)
        per_row = np.zeros(len(rows), dtype=np.int64)
        N = 0
        for codes in hits_codes:
            idx = np.flatnonzero(codes)
            c = codes[idx]
            sub = block[:, idx]
            shared = (sub > 0).sum(1, dtype=np.int64)
            match = ((sub == c) & (c != columns.OTHER)).sum(1, dtype=np.int64)
            valid = shared >= MIN_SHARE * len(idx)
            per_row += np.where(valid, match, 0)
            key = np.where(valid, match * SCALE + shared, -1)
            best = np.maximum.reduceat(key, offsets)
            ok = best >= 0
            M += np.where(ok, best // SCALE, 0)
            S += np.where(ok, best % SCALE, 0)
            N += len(idx)
        local = np.arange(len(rows))
        win = np.maximum.reduceat(per_row * 2**24 + (2**24 - 1 - local), offsets)
        return M, (S >= MIN_SHARE * N) & (S > 0), N, rows[2**24 - 1 - win % 2**24]


def route_calls(hits: pl.DataFrame, routes: dict[int, Route], n: int, keep: float) -> pl.DataFrame:
    """Per read pair: the best units over its families, ties kept, with its matches (the score),
    columns, identity, each tied unit's winning protein, the units within ``AMBIGUOUS`` of the
    best, and the call's status (its family's strongest hit: ``ga``, ``one``, ``set``, ``weak``)."""
    rows = []
    for (read,), pair in hits.group_by(["read"], maintain_order=True):
        best, call = -1.0, None
        for (acc,), g in pair.group_by(["pfam_accession"], maintain_order=True):
            route = routes.get(acc)
            if route is None:
                continue
            length = route.fam.matrix.shape[1] - 1
            hits_codes = [
                columns.column_residues(*r, length)[0]
                for r in g.select("hmm_from", "ali_from", "cigar", "sequence").rows()
            ]
            cand = route.candidates(np.max(np.stack(hits_codes), 0), n, keep)
            M, ok, N, win = route.identity(cand, hits_codes)
            if not ok.any():
                continue
            ident = np.where(ok, M / N, -1.0)
            top = ident.max()
            if top > best + EPS:
                units = route.fam.units[cand]
                is_top = ident >= top - EPS
                near = units[(ident >= top - AMBIGUOUS) & ~is_top]
                status = min(g["status"].to_list(), key=STATUS_ORDER.__getitem__)
                best = top
                call = (acc, int(M[is_top][0]), N, top, units[is_top].tolist(),
                        route.fam.proteins[win[is_top]].tolist(), near.tolist(), status,
                        len(cand))  # fmt: skip
        if call is not None:
            rows.append((read, *call))
    schema = {
        "read": pl.Int64, "pfam_accession": pl.Int32, "matches": pl.Int64, "columns": pl.Int64,
        "identity": pl.Float64, "units": pl.List(pl.Int64), "proteins": pl.List(pl.Int64),
        "near": pl.List(pl.Int64), "status": pl.String, "candidates": pl.Int64,
    }  # fmt: skip
    return pl.DataFrame(rows, schema=schema, orient="row")


def write_handoff(dst: Path, calls: pl.DataFrame, ref: pl.DataFrame, reads_dir: Path) -> None:
    """``kfp-pileup reference``'s format, as tree-index tier 2 writes it (``tier2.handoff``).
    Targets: each called unit's winning proteins and its representative. Neighbours: up to
    ``NEIGHBOURS`` near units per family that no pair called, by representative (else their
    first member). ``assignments.parquet`` splits a pair equally over its tied units, each on its
    winning protein."""
    neighbours = 5  # tier2.NEIGHBOURS
    dst.mkdir(parents=True, exist_ok=True)
    won = (
        calls.with_columns(weight=1 / pl.col("units").list.len())
        .explode(["units", "proteins"], empty_as_null=True)
        .select("read", "matches", "pfam_accession", "weight", unit="units", protein_id="proteins")
    )
    called = won.select("pfam_accession", "unit").unique("unit")
    members = ref.rename({"cluster_rep": "unit"})
    rep_or_first = members.sort(pl.col("protein_id") != pl.col("unit")).unique(
        "unit", keep="first", maintain_order=True
    )
    near = (
        calls.explode("near", empty_as_null=True).select("pfam_accession", unit="near")
        .drop_nulls().join(called.select("unit"), on="unit", how="anti").unique("unit")
        .filter(pl.int_range(pl.len()).over("pfam_accession") < neighbours)
    )  # fmt: skip
    targets = (
        pl.concat(
            [
                members.join(won.select("protein_id").unique(), on="protein_id"),
                rep_or_first.join(called.select("unit"), on="unit"),
            ]
        )
        .unique("protein_id")
        .join(called, on="unit")
    )
    refs = pl.concat(
        [
            targets.with_columns(role=pl.lit("target")),
            rep_or_first.join(near, on="unit").with_columns(role=pl.lit("neighbour")),
        ]
    ).sort("unit", "protein_id").with_columns(
        cluster_rep=pl.col("unit"),
        unit=pl.col("unit").rank("dense").cast(pl.Int64) - 1,
        rep=pl.col("protein_id") == pl.col("unit"),
        component=pl.col("pfam_accession"),
        seq_id=pl.format("s{}", pl.int_range(pl.len())),
        length=pl.col("sequence").str.len_bytes(),
        pos=pl.lit(-1),
    )  # fmt: skip
    with open(dst / "ref.fa", "w") as f:
        f.writelines(f">{i}\n{s}\n" for i, s in refs.select("seq_id", "sequence").iter_rows())
    refs.select(
        "unit", "cluster_rep", "protein_id", "rep", "component", "role", "seq_id", "length", "pos"
    ).write_csv(dst / "refs.tsv", separator="\t")
    (dst / "molecule.txt").write_text("protein\n")
    won.join(refs.select("protein_id", "seq_id", "unit"), on="protein_id").select(
        pair="read", seq_id="seq_id", unit="unit_right", score=pl.col("matches").cast(pl.Float64),
        weight="weight",
    ).sort("pair", "seq_id").write_parquet(dst / "assignments.parquet")  # fmt: skip
    want = set(calls["read"].to_list())
    for m in (1, 2):
        with (
            gzip.open(reads_dir / f"reads_{m}.fq.gz", "rt") as f,
            gzip.open(dst / f"reads_{m}.fq.gz", "wt", compresslevel=1) as out,
        ):
            for i, rec in enumerate(zip(*[f] * 4, strict=True)):
                if i in want:
                    out.writelines(rec)


def route_hits(d: Path, prefix: str, source: str) -> pl.DataFrame | None:
    """A read set's hits for the route, with ``status``: every ``hmmsearch`` hit (``ga`` or
    ``weak``), GA ones only, or kfpfam's recruited loci (every target in a set, or ``top`` only).
    None for the reversed reads under recruitment, which has no null run of its own here."""
    orfs = pl.read_parquet(d / f"{prefix}orfs.parquet")
    if source in ("hmmsearch", "ga"):
        hits = pl.read_parquet(d / f"{prefix}hits.parquet").rename({"protein_id": "orf_id"})
        if source == "ga":
            hits = hits.filter("above_ga")
        hits = hits.with_columns(
            status=pl.when("above_ga").then(pl.lit("ga")).otherwise(pl.lit("weak"))
        )
        return hits.join(orfs, on="orf_id").sort("read", "pfam_accession", "orf_id", "domain")
    if prefix:
        return None
    hits = pl.read_parquet(d / "recruited.parquet").filter(pl.col("status") != "none")
    if source == "recruited_top":
        hits = hits.filter(pl.col("top") == 1)
    return (
        hits.drop("read", "mate").join(orfs, on="orf_id")
        .select("read", "orf_id", "pfam_accession", "hmm_from", "ali_from", "cigar", "sequence",
                "above_ga", "status", "locus")
        .sort("read", "pfam_accession", "orf_id", "locus")
    )  # fmt: skip


def tag_of(a: argparse.Namespace) -> str:
    return "" if (a.hits, a.ref) == ("hmmsearch", "ga") else f"_{a.hits}_{a.ref}"


def query(a: argparse.Namespace) -> None:
    t0 = time.perf_counter()
    fams, _ = load_families(a.data, a.ref)
    routes = {acc: Route(f) for acc, f in fams.items()}
    timing = {"build_s": time.perf_counter() - t0, "hits": a.hits, "ref": a.ref}
    ref = pl.read_parquet(
        a.data / "holdout/members.parquet", columns=["protein_id", "cluster_rep", "sequence"]
    )
    tag = tag_of(a)
    for name, prefix in itertools.product(prepare.SETS, ("", "null_")):
        d = a.data / name
        hits = route_hits(d, prefix, a.hits)
        key = f"{name}{'_null' if prefix else ''}"
        t1 = time.perf_counter()
        calls = (
            route_calls(hits, routes, a.anchors, a.keep) if hits is not None
            else route_calls(pl.DataFrame(schema={"read": pl.Int64}), routes, a.anchors, a.keep)
        )  # fmt: skip
        timing[f"{key}_s"] = time.perf_counter() - t1
        timing[f"{key}_pairs"] = pl.read_parquet(d / "reads.parquet").height
        calls.write_parquet(d / f"{prefix}anchored{tag}.parquet")
        if not prefix:
            write_handoff(d / f"handoff{tag}", calls, ref, d)
        print(f"{key}: {calls.height} pairs called in {timing[f'{key}_s']:.1f} s", flush=True)
    (a.data / f"anchored{tag}.json").write_text(json.dumps(timing, indent=1))


def compare(a: argparse.Namespace) -> None:
    """Phases 3-4: the anchored route against DIAMOND (all units), by read set and identity bin,
    over every pair and over pairs with a hit on a family of their source; and by the call's
    status."""
    tag = tag_of(a)
    timing = json.loads((a.data / f"anchored{tag}.json").read_text())
    ref = pl.read_parquet(a.data / "holdout/members.parquet", columns=["protein_id", "cluster_rep"])
    out, by_status = [], []
    for name in prepare.SETS:
        d = a.data / name
        truth = truth_of(a.data, name)
        calls = pl.read_parquet(d / f"anchored{tag}.parquet")
        null = pl.read_parquet(d / f"null_anchored{tag}.parquet")
        n_pairs = truth.height
        # the bar: fewest matches such that reversed pairs are called at most --null-rate
        bar = 0
        for m in sorted(set(null["matches"].to_list())):
            if (null["matches"] >= bar).sum() / n_pairs <= a.null_rate:
                break
            bar = m + 1
        right = pl.col("units").list.contains(pl.col("unit"))
        anchored = (
            calls.filter(pl.col("matches") >= bar)
            .join(truth.select("read", "unit"), on="read")
            .select(
                "read",
                "status",
                right=right,
                credit=right / pl.col("units").list.len(),
                called=pl.lit(True),
            )  # fmt: skip
        )
        dia = pl.read_parquet(d / "diamond.parquet").join(
            ref.rename({"cluster_rep": "unit"}), on="protein_id"
        )
        dia_calls = diamond_calls(dia, truth, None).select(
            "read", dia_called=pl.lit(True), dia_right="right",
            dia_credit=pl.col("right") / pl.col("ties"),
        )  # fmt: skip
        src = (
            pl.read_parquet(d / "oracle.parquet")
            .filter(pl.col("score") == "read")
            .select("read", "hit_source")
        )
        df = (
            truth.select("read", "identity_bin")
            .join(anchored, on="read", how="left")
            .join(dia_calls, on="read", how="left")
            .join(src, on="read", how="left")
            .with_columns(
                pl.col("called", "right", "dia_called", "dia_right", "hit_source").fill_null(False),
                pl.col("credit", "dia_credit").fill_null(0.0),
            )
        )
        for scope, sub in (("all", df), ("hit_source", df.filter("hit_source"))):
            out.append(
                sub.group_by("identity_bin")
                .agg(
                    pairs=pl.len(),
                    called=pl.col("called").mean(),
                    right=pl.col("right").mean(),
                    dia_right=pl.col("dia_right").mean(),
                    credit=pl.col("credit").mean(),
                    dia_credit=pl.col("dia_credit").mean(),
                    precision=pl.col("right").filter("called").mean(),
                    dia_precision=pl.col("dia_right").filter("dia_called").mean(),
                    dia_right_no_hit=(pl.col("dia_right") & ~pl.col("hit_source")).sum()
                    / pl.col("dia_right").sum(),
                )
                .with_columns(
                    set=pl.lit(name),
                    scope=pl.lit(scope),
                    bar=pl.lit(bar),
                    null_called=pl.lit((null["matches"] >= bar).sum() / n_pairs),
                )
            )
        by_status.append(
            df.filter("called").group_by("status")
            .agg(pairs=pl.len(), precision=pl.col("right").mean(), credit=pl.col("credit").mean(),
                 dia_right=pl.col("dia_right").mean())
            .with_columns(set=pl.lit(name), share_of_pairs=pl.col("pairs") / n_pairs)
        )  # fmt: skip
        timing[f"{name}_ms_per_pair"] = timing[f"{name}_s"] / timing[f"{name}_pairs"] * 1e3
    table = pl.concat(out).select(
        "set", "scope", "identity_bin", pl.exclude("set", "scope", "identity_bin")
    )
    table = table.sort("set", "scope", "identity_bin", descending=[False, False, True])
    table.write_csv(a.data / f"phase3{tag}.tsv", separator="\t", float_precision=4)
    status = (
        pl.concat(by_status)
        .select("set", "status", pl.exclude("set", "status"))
        .sort("set", "status")
    )
    status.write_csv(a.data / f"status{tag}.tsv", separator="\t", float_precision=4)
    (a.data / f"anchored{tag}.json").write_text(json.dumps(timing, indent=1))
    with pl.Config(tbl_rows=-1, tbl_cols=-1, float_precision=3, tbl_width_chars=250):
        print(table)
        print(status)
    print(json.dumps(timing, indent=1))


# From pseudo-alignment to alignment (phase 4b): fill between the HMM's pairs, then extend. ---

X_DROP = 15  # BLOSUM62 points below the best before an ungapped extension stops
GAP_OPEN, GAP_EXTEND = -11, -1


def blosum() -> NDArray[np.int64]:
    """BLOSUM62 on byte codes, for scoring residue pairs."""
    m = substitution_matrices.load("BLOSUM62")
    out = np.full((256, 256), -4, dtype=np.int64)
    for a in m.alphabet:
        for b in m.alphabet:
            out[ord(a), ord(b)] = int(m[a][b])
    return out


def build_alignment(orf: str, prot: str, i: NDArray, j: NDArray, sub: NDArray, aligner) -> dict:
    """ORF residue -> protein residue, per stage: ``implied`` (the HMM's match columns both
    use), ``filled`` (gaps between consecutive pairs: along the diagonal when both sides skip
    the same number, else a global DP of the two stretches), ``extended`` (ungapped from the
    first and last pair to the ORF's ends, X-drop). Also the DP cells used."""
    implied = dict(zip(i.tolist(), j.tolist(), strict=True))
    filled = dict(implied)
    cells = 0
    for (i1, j1), (i2, j2) in itertools.pairwise(zip(i.tolist(), j.tolist(), strict=True)):
        di, dj = i2 - i1, j2 - j1
        if di == dj:
            filled.update((i1 + k, j1 + k) for k in range(1, di))
        elif di > 1 and dj > 1:
            cells += (di - 1) * (dj - 1)
            aln = aligner.align(orf[i1 + 1 : i2], prot[j1 + 1 : j2])[0]
            for (qa, qb), (sa, _) in zip(*aln.aligned, strict=True):
                filled.update((i1 + 1 + qa + k, j1 + 1 + sa + k) for k in range(qb - qa))
    extended = dict(filled)
    o, p = np.frombuffer(orf.encode(), np.uint8), np.frombuffer(prot.encode(), np.uint8)
    for start, step in ((min(filled), -1), (max(filled), 1)):
        s0, t0 = start, filled[start]
        score = best = 0
        keep = 0
        k = 1
        while 0 <= s0 + step * k < len(o) and 0 <= t0 + step * k < len(p):
            score += sub[o[s0 + step * k], p[t0 + step * k]]
            if score > best:
                best, keep = score, k
            elif best - score > X_DROP:
                break
            k += 1
        extended.update((s0 + step * m, t0 + step * m) for m in range(1, keep + 1))
    return {"implied": implied, "filled": filled, "extended": extended, "cells": cells}


def align(a: argparse.Namespace) -> None:
    """For each ORF hit that DIAMOND also aligns to the truth protein (the settled route's hits:
    kfpfam's recruited, top only), the alignment built from the HMM's pairs at each stage,
    against DIAMOND's: residues of the ORF aligned, the share placed as DIAMOND places them
    (precision), and DIAMOND's pairs reproduced (recall)."""
    from Bio.Align import PairwiseAligner

    aligner = PairwiseAligner(mode="global", open_gap_score=GAP_OPEN, extend_gap_score=GAP_EXTEND)
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    sub = blosum()
    _, domains = load_families(a.data, "weak")
    size = lengths(a.data / "hmms.hmm")
    ref = pl.read_parquet(a.data / "holdout/members.parquet", columns=["protein_id", "sequence"])
    prot_seq = dict(ref.iter_rows())
    rows = []
    t0 = time.perf_counter()
    n_built = 0
    for name in prepare.SETS:
        d = a.data / name
        truth = truth_of(a.data, name)
        hits = route_hits(d, "", "recruited_top")
        dia = (
            pl.read_parquet(d / "diamond.parquet")
            .sort("bitscore", descending=True)
            .unique(["read", "mate", "protein_id"], keep="first")
            .select(
                "read",
                "mate",
                relative="protein_id",
                qstart="qstart",
                qframe="qframe",
                sstart="sstart",
                dcigar="cigar",
                dlen="length",
            )  # fmt: skip
        )
        dom = domains.select(
            relative="protein_id", pfam_accession="pfam_accession", p_from="hmm_from",
            p_ali="ali_from", p_cigar="cigar",
        )  # fmt: skip
        pairs = (
            hits.join(
                pl.read_parquet(d / "orfs.parquet").select("orf_id", "mate", "strand", "start"),
                on="orf_id",
            )
            .join(truth.select("read", "relative", "identity_bin", "read_len"), on="read")
            .join(dia, on=["read", "mate", "relative"])
            .join(dom, on=["relative", "pfam_accession"])
        )  # fmt: skip
        best: dict[tuple, tuple] = {}
        for r in pairs.iter_rows(named=True):
            length = size[r["pfam_accession"]]
            pseq = prot_seq[r["relative"]]
            _, orf_pos = columns.column_residues(
                r["hmm_from"], r["ali_from"], r["cigar"], r["sequence"], length
            )
            _, p_pos = columns.column_residues(r["p_from"], r["p_ali"], r["p_cigar"], pseq, length)
            both = (orf_pos >= 0) & (p_pos >= 0)
            if both.sum() < 2:
                continue
            t1 = time.perf_counter()
            built = build_alignment(r["sequence"], pseq, orf_pos[both], p_pos[both], sub, aligner)
            n_built += 1
            elapsed = time.perf_counter() - t1
            dmap = diamond_residues(r["qstart"], r["qframe"], r["sstart"], r["dcigar"])
            n_orf = len(r["sequence"])
            s0, rl = r["start"], r["read_len"]
            keys = [s0 + 3 * x + 1 if r["strand"] == "+" else rl - s0 - 3 * x for x in range(n_orf)]
            in_dia = {x: dmap[k] for x, k in enumerate(keys) if k in dmap}
            out = [
                r["read"],
                r["mate"],
                r["orf_id"],
                r["identity_bin"],
                n_orf,
                len(in_dia),
                built["cells"],
                elapsed,
            ]
            for stage in ("implied", "filled", "extended"):
                m = built[stage]
                shared = [x for x in m if x in in_dia]
                agree = sum(m[x] == in_dia[x] for x in shared)
                out += [len(m), len(shared), agree]
            k = (name, r["orf_id"], r["pfam_accession"])
            # repeats: several of the protein's domains are in the family; keep the best agreeing
            if k not in best or out[-1] > best[k][-1]:
                best[k] = (name, *out)
        rows += best.values()
    cols = [
        "set",
        "read",
        "mate",
        "orf_id",
        "identity_bin",
        "orf_len",
        "diamond_aligned",
        "dp_cells",
        "seconds",
    ]
    for stage in ("implied", "filled", "extended"):
        cols += [f"{stage}_aligned", f"{stage}_shared", f"{stage}_agree"]
    df = pl.DataFrame(rows, schema=cols, orient="row")
    df.write_parquet(a.data / "align.parquet")
    summary = (
        df.group_by("set", "identity_bin")
        .agg(
            *(
                expr
                for stage in ("implied", "filled", "extended")
                for expr in (
                    (pl.col(f"{stage}_aligned").sum() / pl.col("orf_len").sum()).alias(
                        f"{stage}_cover"
                    ),
                    (pl.col(f"{stage}_agree").sum() / pl.col(f"{stage}_shared").sum()).alias(
                        f"{stage}_precision"
                    ),
                    (pl.col(f"{stage}_agree").sum() / pl.col("diamond_aligned").sum()).alias(
                        f"{stage}_recall"
                    ),
                )
            ),
            hits=pl.len(),
            orf_len=pl.col("orf_len").mean(),
            diamond_cover=pl.col("diamond_aligned").sum() / pl.col("orf_len").sum(),
            dp_hits=(pl.col("dp_cells") > 0).mean(),
            dp_cells=pl.col("dp_cells").mean(),
            us_per_hit=pl.col("seconds").mean() * 1e6,
        )
        .sort("set", "identity_bin", descending=[False, True])
    )
    summary.write_csv(a.data / "align.tsv", separator="\t", float_precision=4)
    with pl.Config(tbl_rows=-1, tbl_cols=-1, float_precision=3, tbl_width_chars=260):
        print(summary)
    print(f"{n_built} alignments built in {time.perf_counter() - t0:.1f} s")


# Phase 4b: the route's alignments as kfp-pileup call's DIAMOND lines. ------------------------

LAMBDA, K = 0.267, 0.041  # BLOSUM62 with gaps 11/1, as DIAMOND scores (gapped Karlin-Altschul)
ROUTE_TAG = "_recruited_top_weak"  # the settled route (phase 4)
STAGE_TAG = {"extended": "", "filled": "_filled", "implied": "_implied"}


def lines_tag(a: argparse.Namespace) -> str:
    """The lines' file suffix: alignment stage, then ``_c<N>`` with competitors."""
    return STAGE_TAG[a.stage] + (f"_c{a.competitors}" if a.competitors else "")


MIN_DEPTH_PILE = 3.0  # kfp-pileup's MIN_DEPTH: pileup columns compared


def diamond_line(name: str, sseqid: str, orf: str, prot: str, pairs: dict, strand: str,
                 start: int, read_len: int, sub: NDArray) -> str:  # fmt: skip
    """One alignment in ``kfp-pileup call``'s DIAMOND fields: qseqid sseqid bitscore qstart qend
    sstart send qframe btop. Read bases are 1-based; on the reverse strand qstart is the higher
    one (DIAMOND's convention, checked against its own lines)."""
    xs = sorted(pairs)
    btop, run, raw = [], 0, 0
    prev_x = prev_j = None
    for x in xs:
        j = pairs[x]
        if prev_x is not None:
            gq, gs = x - prev_x - 1, j - prev_j - 1  # residues in the query, subject only
            for g, side in ((gq, "q"), (gs, "s")):
                if g > 0:
                    if run:
                        btop.append(str(run))
                        run = 0
                    raw += GAP_OPEN + GAP_EXTEND * g
                    skipped = orf[prev_x + 1 : x] if side == "q" else prot[prev_j + 1 : j]
                    btop += [f"{r}-" if side == "q" else f"-{r}" for r in skipped]
        a, b = orf[x], prot[j]
        raw += int(sub[ord(a), ord(b)])
        if a == b:
            run += 1
        else:
            if run:
                btop.append(str(run))
                run = 0
            btop.append(a + b)
        prev_x, prev_j = x, j
    if run:
        btop.append(str(run))
    bits = (LAMBDA * raw - math.log(K)) / math.log(2)
    if strand == "+":
        qstart, qend = start + 3 * xs[0] + 1, start + 3 * xs[-1] + 3
        frame = (qstart - 1) % 3 + 1
    else:
        qstart, qend = read_len - start - 3 * xs[0], read_len - start - 3 * xs[-1] - 2
        frame = -((read_len - qstart) % 3 + 1)
    sstart, send = pairs[xs[0]] + 1, pairs[xs[-1]] + 1
    fields = (name, sseqid, f"{bits:.1f}", qstart, qend, sstart, send, frame, "".join(btop))
    return "\t".join(map(str, fields)) + "\n"


def read_names(path: Path) -> list[str]:
    with gzip.open(path, "rt") as f:
        return [line[1:].split()[0] for i, line in enumerate(f) if i % 4 == 0]


def lines(a: argparse.Namespace) -> None:
    """Per read set, ``lines_{1,2}.tsv`` beside the settled route's handoff: each called pair's
    mates aligned (``build_alignment``) to every tied unit's winning protein and to the near
    units' neighbour sequences, through the mate's hits in the call's family."""
    from Bio.Align import PairwiseAligner

    aligner = PairwiseAligner(mode="global", open_gap_score=GAP_OPEN, extend_gap_score=GAP_EXTEND)
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    sub = blosum()
    _, domains = load_families(a.data, "weak")
    size = lengths(a.data / "hmms.hmm")
    dom = {}
    for k, g in domains.group_by("protein_id", "pfam_accession"):
        dom[k] = g.select("hmm_from", "ali_from", "cigar", "sequence").rows()
    for name in prepare.SETS:
        d = a.data / name
        handoff = d / f"handoff{ROUTE_TAG}"
        t0 = time.perf_counter()
        refs = pl.read_csv(handoff / "refs.tsv", separator="\t")
        seq_of = dict(refs.select("protein_id", "seq_id").iter_rows())
        rep_of = dict(
            refs.filter(pl.col("role") == "neighbour")
            .select("cluster_rep", "protein_id")
            .iter_rows()
        )
        calls = pl.read_parquet(d / f"anchored{ROUTE_TAG}.parquet")
        hits = route_hits(d, "", "recruited_top").join(
            pl.read_parquet(d / "orfs.parquet").select("orf_id", "mate", "strand", "start"),
            on="orf_id",
        )
        by_pair = {k[0]: g for k, g in hits.group_by(["read"])}
        # competitors: every reference sequence of a family, by its domains' column residues
        family_rows = {}
        for (acc,), g in refs.filter(pl.col("role").is_in(["target", "neighbour"])).group_by(
            ["component"]
        ):
            pids, rows = [], []
            for pid in g["protein_id"].to_list():
                for r in dom.get((pid, acc), []):
                    pids.append(pid)
                    rows.append(columns.column_residues(*r, size[acc])[0])
            if rows:
                family_rows[acc] = (np.array(pids), np.stack(rows))
        names = read_names(d / "reads_1.fq.gz")
        read_len = pl.read_parquet(d / "reads.parquet")["pair_bases"][0] // 2
        out = {1: [], 2: []}
        n_lines = 0
        for c in calls.sort("read").iter_rows(named=True):
            targets = list(c["proteins"]) + [rep_of[u] for u in c["near"] if u in rep_of]
            mate_hits = by_pair[c["read"]].filter(pl.col("pfam_accession") == c["pfam_accession"])
            length = size[c["pfam_accession"]]
            for h in mate_hits.iter_rows(named=True):
                _, orf_pos = columns.column_residues(
                    h["hmm_from"], h["ali_from"], h["cigar"], h["sequence"], length
                )
                hit_targets = list(targets)
                if a.competitors and c["pfam_accession"] in family_rows:
                    pids, rows = family_rows[c["pfam_accession"]]
                    codes = columns.column_residues(
                        h["hmm_from"], h["ali_from"], h["cigar"], h["sequence"], length
                    )[0]
                    match = ((rows == codes) & (codes > 0)).sum(1)
                    order = np.argsort(-match, kind="stable")
                    seen = set(hit_targets)
                    for pid in pids[order].tolist():
                        if len(hit_targets) >= len(targets) + a.competitors:
                            break
                        if pid not in seen:
                            seen.add(pid)
                            hit_targets.append(pid)
                for pid in hit_targets:
                    best = None
                    for r in dom.get((pid, c["pfam_accession"]), []):
                        _, p_pos = columns.column_residues(*r, length)
                        both = (orf_pos >= 0) & (p_pos >= 0)
                        if both.sum() >= 2 and (best is None or both.sum() > best[0]):
                            best = (both.sum(), orf_pos[both], p_pos[both], r[3])
                    if best is None:
                        continue
                    built = build_alignment(h["sequence"], best[3], best[1], best[2], sub, aligner)
                    if len(built[a.stage]) < 2:
                        continue
                    out[h["mate"]].append(diamond_line(
                        names[c["read"]], seq_of[pid], h["sequence"], best[3], built[a.stage],
                        h["strand"], h["start"], read_len, sub,
                    ))  # fmt: skip
                    n_lines += 1
        for m in (1, 2):
            (handoff / f"lines{lines_tag(a)}_{m}.tsv").write_text("".join(out[m]))
        print(
            f"{name}: {n_lines} lines for {calls.height} pairs in {time.perf_counter() - t0:.1f} s",
            flush=True,
        )


def pileup(a: argparse.Namespace) -> None:
    """Per read set, on the settled route's handoff: DIAMOND as the pileup workflow runs it
    (``blastx --top 10 -F 15 --sensitive``, its fields) and the route's own lines, each through
    ``kfp-pileup call`` -> ``call_diamond/`` and ``call_anchored/``, with times in
    ``pileup_times.json``."""
    kfp = Path(sys.executable).with_name("kfp-pileup")
    times = {}
    for name in prepare.SETS:
        h = a.data / name / f"handoff{ROUTE_TAG}"
        cap = f"{a.max_gb:g}g"
        docker = ["docker", "run", "--rm", f"--memory={cap}", f"--memory-swap={cap}", "--cpus=4",
                  "-v", f"{h.resolve()}:/w", "-w", "/w", prepare.DIAMOND, "diamond"]  # fmt: skip
        if not a.skip_diamond:
            prepare.run(*docker, "makedb", "--in", "ref.fa", "--db", "ref", "--threads", 4)
            t0 = time.perf_counter()
            for m in (1, 2):
                prepare.run(*docker, "blastx", "--db", "ref", "--query", f"reads_{m}.fq.gz",
                            "--out", f"diamond_{m}.tsv", "--outfmt", 6, "qseqid", "sseqid",
                            "bitscore", "qstart", "qend", "sstart", "send", "qframe", "btop",
                            "--top", 10, "-F", 15, "--sensitive", "--threads", 4,
                            "-b", 0.5)  # fmt: skip
            times[f"{name}_diamond_s"] = time.perf_counter() - t0
        arms = [("diamond", "diamond")] * (not a.skip_diamond) + [
            (f"anchored{lines_tag(a)}", f"lines{lines_tag(a)}")
        ]
        for arm, prefix in arms:
            t0 = time.perf_counter()
            tables = [h / f"{prefix}_{m}.tsv" for m in (1, 2)]
            prepare.run(
                kfp,
                "call",
                h,
                h / f"call_{arm}",
                *tables,
                "--reads",
                h / "reads_1.fq.gz",
                "--reads",
                h / "reads_2.fq.gz",
                "--stats",
                h / f"call_{arm}.json",
            )
            times[f"{name}_call_{arm}_s"] = time.perf_counter() - t0
    out = a.data / f"pileup_times{lines_tag(a)}.json"
    out.write_text(json.dumps(times, indent=1))
    print(json.dumps(times, indent=1))


def variant_truth(data: Path, name: str, h: Path, variants: pl.DataFrame) -> float | None:
    """Strains only: of the variants on a strain's own source protein, the share whose alt is
    the strain's residue there (``truth.tsv``'s back-translated CDS, translated; ``pos``
    1-based). None for relatives, whose reads are the held-out protein itself."""
    if name != "strains":
        return None
    refs = pl.read_csv(h / "refs.tsv", separator="\t").select("seq_id", "protein_id")
    truth = pl.read_csv(data / name / "truth.tsv", separator="\t").select(
        protein_id=pl.col("member").cast(pl.Int64), cds="cds"
    )
    strain = {
        pid: "".join(prepare.CODON.get(cds[i : i + 3], "X") for i in range(0, len(cds) - 2, 3))
        for pid, cds in truth.iter_rows()
    }
    x = variants.join(refs, on="seq_id").filter(pl.col("protein_id").is_in(list(strain)))
    good = [strain[pid][pos - 1 : pos] == alt for _, pos, alt, pid in x.iter_rows()]
    return sum(good) / len(good) if good else None


def pileupcmp(a: argparse.Namespace) -> None:
    """``call_diamond`` against ``call_anchored`` (``pileup``; ``--stage``), per read set ->
    ``pileup_compare.tsv``: mates kept; each mate's best unit by ``call``'s EM posterior,
    against each other and the truth unit; pileup columns (depth >= ``MIN_DEPTH`` in both)
    by depth and consensus; per-unit abundance and identity; variants, and on strains their
    truth (:func:`variant_truth`)."""
    rows = []
    for name in prepare.SETS:
        h = a.data / name / f"handoff{ROUTE_TAG}"
        refs = pl.read_csv(h / "refs.tsv", separator="\t").select("unit", "cluster_rep").unique()
        strain_unit = pl.read_parquet(a.data / name / "truth.parquet").select(
            "strain", truth="unit"
        )
        best, piles, units, variants, true_share = {}, {}, {}, {}, {}
        anchored = h / f"call_anchored{lines_tag(a)}"
        for arm, c in (("diamond", h / "call_diamond"), ("anchored", anchored)):
            asg = (
                pl.read_parquet(c / "assignments.parquet")
                .filter("kept")
                .sort("em_posterior", descending=True)
                .unique(["read_name", "mate"], keep="first")
                .join(refs, on="unit")
                .with_columns(strain=pl.col("read_name").str.extract(r"^s(\d+)_").cast(pl.Int64))
                .join(strain_unit, on="strain")
            )
            best[arm] = asg.select(
                "read_name", "mate", "cluster_rep", right=pl.col("cluster_rep") == pl.col("truth")
            )
            p = (
                pl.read_parquet(c / "pileup.parquet")
                .group_by("seq_id", "col")
                .agg(
                    depth=pl.col("weight").sum(),
                    consensus=pl.col("residue").sort_by("weight").last(),
                )
            )
            piles[arm] = p
            units[arm] = pl.read_csv(c / "units.tsv", separator="\t").select(
                "cluster_rep", "em_abundance", "identity"
            )
            variants[arm] = pl.read_csv(c / "variants.tsv", separator="\t").select(
                "seq_id", "pos", "alt"
            )
            true_share[arm] = variant_truth(a.data, name, h, variants[arm])
        b = best["diamond"].join(
            best["anchored"], on=["read_name", "mate"], how="full", suffix="_a"
        )
        both = b.drop_nulls(["cluster_rep", "cluster_rep_a"])
        p = piles["diamond"].join(piles["anchored"], on=["seq_id", "col"], suffix="_a")
        deep = p.filter((pl.col("depth") >= MIN_DEPTH_PILE) & (pl.col("depth_a") >= MIN_DEPTH_PILE))
        u = (
            units["diamond"]
            .join(units["anchored"], on="cluster_rep", how="full", suffix="_a")
            .with_columns(pl.col("em_abundance", "em_abundance_a").fill_nan(0.0).fill_null(0.0))
        )
        v = variants["diamond"].join(variants["anchored"], on=["seq_id", "pos", "alt"], how="inner")
        rows.append({
            "set": name,
            "mates_diamond": best["diamond"].height,
            "mates_anchored": best["anchored"].height,
            "mates_both": both.height,
            "same_best_unit": (both["cluster_rep"] == both["cluster_rep_a"]).mean(),
            "right_diamond": best["diamond"]["right"].mean(),
            "right_anchored": best["anchored"]["right"].mean(),
            "columns_diamond": piles["diamond"].height,
            "columns_anchored": piles["anchored"].height,
            "columns_deep_both": deep.height,
            "depth_corr": deep.select(pl.corr("depth", "depth_a")).item(),
            "same_consensus": (deep["consensus"] == deep["consensus_a"]).mean(),
            "abundance_corr": u.select(pl.corr("em_abundance", "em_abundance_a")).item(),
            "identity_mean_diff": u.drop_nulls(["identity", "identity_a"]).select(
                (pl.col("identity_a") - pl.col("identity")).mean()
            ).item(),
            "variants_diamond": variants["diamond"].height,
            "variants_anchored": variants["anchored"].height,
            "variants_shared": v.height,
            "variants_true_diamond": true_share["diamond"],
            "variants_true_anchored": true_share["anchored"],
        })  # fmt: skip
    df = pl.DataFrame(rows)
    df.write_csv(a.data / f"pileup_compare{lines_tag(a)}.tsv", separator="\t", float_precision=4)
    with pl.Config(tbl_cols=-1, tbl_width_chars=250, float_precision=3):
        print(df.transpose(include_header=True, column_names="set"))


# Confirm, then map every read. ---------------------------------------------------------------

MEMBERS_PER_UNIT = 3  # the pileup workflow's --members: sequences per target unit


def confirmed_dir(a: argparse.Namespace) -> str:
    """A read set's confirmed reference, by the confirmation threshold."""
    return f"confirmed_w{a.min_weight:g}"


NEIGHBOURS_PER_FAMILY = 5  # tier2.NEIGHBOURS


def confirm(a: argparse.Namespace) -> None:
    """Per read set: units confirmed by the settled route's calls (each called pair's weight
    split over its tied units; confirmed at ``--min-weight``) -> ``confirmed.parquet`` (unit,
    family, weight, truth) and ``confirmed/``, a ``kfp-pileup reference`` directory: up to
    ``MEMBERS_PER_UNIT`` members per confirmed unit (representative first) as targets, and up to
    ``NEIGHBOURS_PER_FAMILY`` near units per family that were not confirmed, by representative."""
    ref = pl.read_parquet(
        a.data / "holdout/members.parquet", columns=["protein_id", "cluster_rep", "sequence"]
    ).rename({"cluster_rep": "unit"})
    for name in prepare.SETS:
        d = a.data / name
        calls = pl.read_parquet(d / f"anchored{ROUTE_TAG}.parquet")
        weight = (
            calls.with_columns(w=1 / pl.col("units").list.len())
            .explode("units", empty_as_null=True)
            .group_by(unit="units")
            .agg(pl.col("w").sum().alias("weight"), pl.col("pfam_accession").mode().first())
        )
        present = pl.read_parquet(d / "truth.parquet").select("unit").unique()
        confirmed = weight.filter(pl.col("weight") >= a.min_weight).with_columns(
            truth=pl.col("unit").is_in(present["unit"].implode())
        )
        confirmed.write_parquet(d / f"{confirmed_dir(a)}.parquet")
        near = (
            calls.explode("near", empty_as_null=True).select("pfam_accession", unit="near")
            .drop_nulls().join(confirmed.select("unit"), on="unit", how="anti").unique("unit")
            .filter(pl.int_range(pl.len()).over("pfam_accession") < NEIGHBOURS_PER_FAMILY)
        )  # fmt: skip
        ordered = ref.sort(pl.col("protein_id") != pl.col("unit"), "protein_id")
        targets = ordered.join(confirmed.select("unit", "pfam_accession"), on="unit").filter(
            pl.int_range(pl.len()).over("unit") < MEMBERS_PER_UNIT
        )
        neighbours = ordered.join(near, on="unit").filter(pl.int_range(pl.len()).over("unit") < 1)
        refs = pl.concat(
            [
                targets.with_columns(role=pl.lit("target")),
                neighbours.with_columns(role=pl.lit("neighbour")),
            ]  # fmt: skip
        ).with_columns(
            cluster_rep=pl.col("unit"),
            unit=pl.col("unit").rank("dense").cast(pl.Int64) - 1,
            rep=pl.col("protein_id") == pl.col("unit"),
            component=pl.col("pfam_accession"),
            seq_id=pl.format("s{}", pl.int_range(pl.len())),
            length=pl.col("sequence").str.len_bytes(),
            pos=pl.lit(-1),
        )
        out = d / confirmed_dir(a)
        out.mkdir(exist_ok=True)
        with open(out / "ref.fa", "w") as f:
            f.writelines(f">{i}\n{s}\n" for i, s in refs.select("seq_id", "sequence").iter_rows())
        refs.select(
            "unit", "cluster_rep", "protein_id", "rep", "component", "role", "seq_id", "length",
            "pos",
        ).write_csv(out / "refs.tsv", separator="\t")  # fmt: skip
        (out / "molecule.txt").write_text("protein\n")
        n_present = present.height
        found = confirmed.filter("truth").height
        print(
            f"{name}: {confirmed.height} units confirmed ({found} present: recall"
            f" {found / n_present:.3f} of {n_present}, precision {found / confirmed.height:.3f});"
            f" {refs.filter(pl.col('role') == 'target').height} targets,"
            f" {neighbours.height} neighbours",
            flush=True,
        )


def mapall(a: argparse.Namespace) -> None:
    """Per read set: every read pair aligned to the confirmed reference (``confirm``) by
    DIAMOND, as the pileup workflow runs it, then ``kfp-pileup call`` -> ``confirmed/call_diamond``;
    times in ``mapall_times.json``."""
    kfp = Path(sys.executable).with_name("kfp-pileup")
    times = {}
    for name in prepare.SETS:
        d = a.data / name
        c = d / confirmed_dir(a)
        cap = f"{a.max_gb:g}g"
        docker = ["docker", "run", "--rm", f"--memory={cap}", f"--memory-swap={cap}", "--cpus=4",
                  "-v", f"{d.resolve()}:/w", "-w", "/w", prepare.DIAMOND, "diamond"]  # fmt: skip
        prepare.run(
            *docker,
            "makedb",
            "--in",
            f"{confirmed_dir(a)}/ref.fa",
            "--db",
            f"{confirmed_dir(a)}/ref",
        )
        t0 = time.perf_counter()
        cd = confirmed_dir(a)
        for m in (1, 2):
            prepare.run(*docker, "blastx", "--db", f"{cd}/ref", "--query", f"reads_{m}.fq.gz",
                        "--out", f"{cd}/diamond_{m}.tsv", "--outfmt", 6, "qseqid", "sseqid",
                        "bitscore", "qstart", "qend", "sstart", "send", "qframe", "btop",
                        "--top", 10, "-F", 15, "--sensitive", "--threads", 4,
                        "-b", 0.5)  # fmt: skip
        times[f"{name}_diamond_s"] = time.perf_counter() - t0
        t0 = time.perf_counter()
        tables = [c / f"diamond_{m}.tsv" for m in (1, 2)]
        prepare.run(kfp, "call", c, c / "call_diamond", *tables, "--reads", d / "reads_1.fq.gz",
                    "--reads", d / "reads_2.fq.gz", "--stats", c / "call_diamond.json")  # fmt: skip
        times[f"{name}_call_s"] = time.perf_counter() - t0
    (a.data / f"mapall_times_{confirmed_dir(a)}.json").write_text(json.dumps(times, indent=1))
    print(json.dumps(times, indent=1))


def mapcmp(a: argparse.Namespace) -> None:
    """Per read set and arm (``confirmed/call_*``): over every read pair, a mate is *mapped* when
    ``call`` keeps an alignment for it and *right* when its best unit by EM posterior is the truth
    unit; split by whether the pair has a hit on a family of its source; variants and their truth
    on strains (:func:`variant_truth`) -> ``mapall.tsv``."""
    rows = []
    for name in prepare.SETS:
        d = a.data / name
        c = d / confirmed_dir(a)
        refs = pl.read_csv(c / "refs.tsv", separator="\t").select("unit", "cluster_rep").unique()
        reads = pl.read_parquet(d / "reads.parquet").select("read", "strain")
        names = pl.DataFrame(
            {"read": range(len(n := read_names(d / "reads_1.fq.gz"))), "read_name": n}
        )
        truth = pl.read_parquet(d / "truth.parquet").select("strain", truth="unit")
        src = (
            pl.read_parquet(d / "oracle.parquet").filter(pl.col("score") == "read")
            .select("read", "hit_source")
        )  # fmt: skip
        mates = (
            reads.join(names, on="read").join(truth, on="strain").join(src, on="read", how="left")
            .with_columns(pl.col("hit_source").fill_null(False))
            .join(pl.DataFrame({"mate": [0, 1]}), how="cross")
        )  # fmt: skip
        for arm_dir in sorted(c.glob("call_*/")):
            arm = arm_dir.name.removeprefix("call_")
            best = (
                pl.read_parquet(arm_dir / "assignments.parquet").filter("kept")
                .sort("em_posterior", descending=True)
                .unique(["read_name", "mate"], keep="first")
                .join(refs, on="unit")
                .select("read_name", pl.col("mate").cast(pl.Int64), "cluster_rep")
            )  # fmt: skip
            df = mates.join(best, on=["read_name", "mate"], how="left").with_columns(
                mapped=pl.col("cluster_rep").is_not_null(),
                right=(pl.col("cluster_rep") == pl.col("truth")).fill_null(False),
            )
            variants = pl.read_csv(arm_dir / "variants.tsv", separator="\t").select(
                "seq_id", "pos", "alt"
            )
            for scope, sub in (
                ("all", df),
                ("hit", df.filter("hit_source")),
                ("no_hit", df.filter(~pl.col("hit_source"))),
            ):
                rows.append({
                    "set": name, "arm": arm, "scope": scope, "mates": sub.height,
                    "mapped": sub["mapped"].mean(), "right": sub["right"].mean(),
                    "precision": sub.filter("mapped")["right"].mean(),
                    "variants": variants.height,
                    "variants_true": variant_truth(a.data, name, c, variants),
                })  # fmt: skip
    df = pl.DataFrame(rows)
    df.write_csv(a.data / f"mapall_{confirmed_dir(a)}.tsv", separator="\t", float_precision=4)
    with pl.Config(tbl_rows=-1, tbl_cols=-1, float_precision=3, tbl_width_chars=200):
        print(df)


TREE = prepare.ROOT / "workflows/tree-index"


def treeindex(a: argparse.Namespace) -> None:
    """Per read set: the tree index's standard tier on the confirmed reference's sequences (the
    same targets and neighbours DIAMOND searches), as one part (build.py, merge.py, store.py),
    then every read pair through tier 1 (query.py --decoy, verify.py), tier 1 then tier 2
    (tier2.py), and tier 2 alone (tier2.py --all-leaves) -> ``<confirmed>/tree/``; each step's
    time in ``tree_times.json``."""
    times = {}
    for name in prepare.SETS:
        d = a.data / name
        c = d / confirmed_dir(a)
        tree = c / "tree"
        tree.mkdir(exist_ok=True)
        refs = pl.read_csv(c / "refs.tsv", separator="\t").select("protein_id")
        members = pl.read_parquet(a.data / "holdout/members.parquet").join(refs, on="protein_id")
        members.write_parquet(tree / "members.parquet")
        reads = [d / "reads_1.fq.gz", d / "reads_2.fq.gz"]
        py = [sys.executable]
        steps = [
            ("build", [*py, TREE / "build.py", "--members", tree / "members.parquet",
                       "--out", tree / "part0000", "--database-proteins", members.height]),
            ("merge", [*py, TREE / "merge.py", "--parts", tree / "part0000",
                       "--out", tree / "index"]),
            ("store", [*py, TREE / "store.py", "build",
                       "--proteins", tree / "part0000/proteins.parquet",
                       "--members", tree / "members.parquet", "--parts", tree / "index/parts.json",
                       "--name", "part0000", "--format", "flat", "--out", tree / "store"]),
            ("tier1_query", [*py, TREE / "query.py", "--index", tree / "index", "--reads", *reads,
                             "--decoy", "--out", tree / "q"]),
            ("tier1_verify", [*py, TREE / "verify.py", "--index", tree / "index", "--store",
                              tree / "store", "--reads", *reads, "--query", tree / "q"]),
            ("tier2_after_tier1", [*py, TREE / "tier2.py", "--index", tree / "index", "--store",
                                   tree / "store", "--reads", *reads, "--query", tree / "q",
                                   "--memory-gb", 3, "--no-handoff", "--out", tree / "t2"]),
            ("tier2_alone", [*py, TREE / "tier2.py", "--index", tree / "index", "--store",
                             tree / "store", "--reads", *reads, "--query", tree / "q",
                             "--all-leaves", "--memory-gb", 3, "--no-handoff", "--out",
                             tree / "t2_alone"]),
        ]  # fmt: skip
        for step, cmd in steps:
            t0 = time.perf_counter()
            prepare.run(*cmd)
            times[f"{name}_{step}_s"] = time.perf_counter() - t0
    (a.data / f"tree_times_{confirmed_dir(a)}.json").write_text(json.dumps(times, indent=1))
    print(json.dumps(times, indent=1))


def mapscore(a: argparse.Namespace) -> None:
    """Every mapper on the same confirmed reference, before the gather, per read pair: DIAMOND
    (each unit's best bit score summed over the mates, best units tied), tier 1 (``verify.py``),
    tier 1 then tier 2, and tier 2 alone (``tier2.py``; best-scoring leaves, ties kept). A pair is
    *called* at the mapper's bar (the tree index's decoy bar at 1%; any line for DIAMOND) and
    *right* when the truth unit is among its best -> ``mapscore_<confirmed>.tsv``, by read set,
    identity bin and whether the pair has a Pfam hit on its source's family."""
    rows = []
    for name in prepare.SETS:
        d = a.data / name
        c = d / confirmed_dir(a)
        tree = c / "tree"
        truth = truth_of(a.data, name).select("read", "identity_bin", "unit")
        src = (
            pl.read_parquet(d / "oracle.parquet").filter(pl.col("score") == "read")
            .select("read", "hit_source")
        )  # fmt: skip
        base = truth.join(src, on="read", how="left").with_columns(
            pl.col("hit_source").fill_null(False)
        )
        names = read_names(d / "reads_1.fq.gz")
        pair_of = pl.DataFrame({"read_name": names, "read": range(len(names))})
        refs = pl.read_csv(c / "refs.tsv", separator="\t").select("seq_id", unit_of="cluster_rep")
        arms = {}
        lines = pl.concat(
            pl.read_csv(c / f"diamond_{m}.tsv", separator="\t", has_header=False,
                        new_columns=["read_name", "seq_id", "bitscore"], columns=[0, 1, 2])
            .with_columns(mate=pl.lit(m))
            for m in (1, 2)
        )  # fmt: skip
        per_unit = (
            lines.join(pair_of, on="read_name").join(refs, on="seq_id")
            .group_by("read", "mate", "unit_of").agg(pl.col("bitscore").max())
            .group_by("read", "unit_of").agg(pl.col("bitscore").sum())
        )  # fmt: skip
        arms["diamond"] = (
            per_unit.filter(pl.col("bitscore") == pl.col("bitscore").max().over("read"))
            .group_by("read")
            .agg(units=pl.col("unit_of"))
        )
        pos_units = (
            pl.read_parquet(tree / "index/proteins.parquet").filter(pl.col("pos") >= 0)
            .join(pl.read_parquet(tree / "members.parquet").select("protein_id", "cluster_rep"),
                  on="protein_id")
            .select("pos", unit_of="cluster_rep").unique()
        )  # fmt: skip
        for arm, qdir in (("tier1", tree / "q"), ("tier1_tier2", tree / "t2"),
                          ("tier2_alone", tree / "t2_alone")):  # fmt: skip
            bar = json.loads((qdir / "verify.json").read_text())["read_bar_0.01"]
            hits = pl.read_parquet(qdir / "reads.parquet", columns=["read", "pos", "score"])
            top = hits.filter(pl.col("score") == pl.col("score").max().over("read"))
            arms[arm] = (
                top.filter(pl.col("score") >= bar).join(pos_units, on="pos")
                .group_by("read").agg(units=pl.col("unit_of").unique())
            )  # fmt: skip
        for arm, called in arms.items():
            df = base.join(called, on="read", how="left").with_columns(
                called=pl.col("units").is_not_null(),
                right=pl.col("units").list.contains(pl.col("unit")).fill_null(False),
                ties=pl.col("units").list.len(),
            )
            for scope, sub in (("all", df), ("hit", df.filter("hit_source")),
                               ("no_hit", df.filter(~pl.col("hit_source")))):  # fmt: skip
                rows.append({"set": name, "arm": arm, "scope": scope, "pairs": sub.height,
                             "called": sub["called"].mean(), "right": sub["right"].mean(),
                             "precision": sub.filter("called")["right"].mean(),
                             "median_ties": sub.filter("called")["ties"].median()})  # fmt: skip
    df = pl.DataFrame(rows)
    df.write_csv(a.data / f"mapscore_{confirmed_dir(a)}.tsv", separator="\t", float_precision=4)
    with pl.Config(tbl_rows=-1, tbl_cols=-1, float_precision=3, tbl_width_chars=200):
        print(df)


STEPS = {
    "oracle": oracle,
    "seeds": seeds,
    "anchors": anchors,
    "query": query,
    "align": align,
    "lines": lines,
    "pileup": pileup,
    "pileupcmp": pileupcmp,
    "compare": compare,
    "confirm": confirm,
    "mapall": mapall,
    "mapcmp": mapcmp,
    "treeindex": treeindex,
    "mapscore": mapscore,
}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("step", choices=list(STEPS))
    p.add_argument("--data", type=Path, default=prepare.ROOT / "data/pfam-anchor")
    p.add_argument("--score", nargs="+", choices=SCORES, default=list(SCORES))
    p.add_argument("--large", type=int, default=1500, help="units in a large family (seeds)")
    p.add_argument("--anchors", type=int, default=16, help="query: anchor columns per read")
    p.add_argument("--keep", type=float, default=0.01, help="query: share of units kept by anchors")
    p.add_argument("--null-rate", type=float, default=0.001, help="compare: reversed pairs called")
    p.add_argument("--stage", choices=list(STAGE_TAG), default="extended", help="lines, pileup")
    p.add_argument("--min-weight", type=float, default=1.0, help="confirm: called pairs per unit")
    p.add_argument("--skip-diamond", action="store_true", help="pileup: anchored lines only")
    p.add_argument(
        "--competitors",
        type=int,
        default=0,
        help="lines: also align each mate to this many more sequences of its family",
    )
    p.add_argument("--hits", choices=HIT_SOURCES, default="hmmsearch", help="query, compare")
    p.add_argument(
        "--ref", choices=REF_ARMS, default="ga", help="query, compare: reference domains"
    )
    p.add_argument(
        "--unanchored",
        dest="anchored",
        action="store_false",
        help="seeds: keys without their HMM column (the ablation)",
    )
    p.add_argument("--max-gb", type=float, default=prepare.MAX_GB, help="memory watchdog's cap")
    a = p.parse_args()
    prepare.watchdog(a.max_gb)
    STEPS[a.step](a)


if __name__ == "__main__":
    main()
