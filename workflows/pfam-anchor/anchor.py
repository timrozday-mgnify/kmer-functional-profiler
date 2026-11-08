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
import re
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import polars as pl
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


def load_families(data: Path) -> tuple[dict[int, Family], pl.DataFrame]:
    """Families, and the reference GA domains with unit and sequence."""
    ref = pl.read_parquet(
        data / "holdout/members.parquet", columns=["protein_id", "cluster_rep", "sequence"]
    )
    domains = (
        pl.read_parquet(data / "domains.parquet")
        .filter("above_ga")
        .join(ref.rename({"cluster_rep": "unit"}), on="protein_id")
    )
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
    ``code = seed * A^width + residues``; ``s`` keeps seeds whose first column divides by it."""

    def __init__(self, fam: Family, seeds: NDArray[np.int64], letters: NDArray[np.int64]) -> None:
        self.seeds, self.letters = seeds, letters
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
            ok = code >= 0
            parts.append((code[ok] << self.bits) | unit_of_row[lo : lo + 2000][np.nonzero(ok)[0]])
        self.post = np.unique(np.concatenate(parts)) if parts else np.zeros(0, dtype=np.int64)

    def codes(self, rows: NDArray[np.uint8]) -> NDArray[np.int64]:
        """(rows x seeds) codes, -1 where a seed column has no residue."""
        x = self.letters[rows[:, self.seeds]]  # rows x seeds x width
        packed = (x * self.base ** np.arange(self.seeds.shape[1])).sum(2)
        packed += np.arange(len(self.seeds)) * self.base ** self.seeds.shape[1]
        return np.where((x >= 0).all(2), packed, -1)

    def sampled(self, s: int) -> NDArray[np.int64]:
        if s == 1:
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
        indexes = {arm: SeedIndex(fam, shape[sh], letters[al]) for arm, sh, al in ARMS}
        build = time.perf_counter() - t0
        rows_per_unit = np.diff(np.r_[fam.starts, len(fam.matrix)])
        posts = {
            (arm, s): indexes[arm].sampled(s)
            for arm, _, _ in ARMS
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
    per_pair.write_parquet(a.data / "seeds.parquet")
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
    table.write_csv(a.data / "seeds.tsv", separator="\t", float_precision=4)
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
    sp.write_csv(a.data / "seed_speed.tsv", separator="\t", float_precision=6)
    with pl.Config(tbl_rows=-1, tbl_cols=-1, float_precision=3, tbl_width_chars=250):
        print(table)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("step", choices=["oracle", "seeds"])
    p.add_argument("--data", type=Path, default=prepare.ROOT / "data/pfam-anchor")
    p.add_argument("--score", nargs="+", choices=SCORES, default=list(SCORES))
    p.add_argument("--large", type=int, default=1500, help="units in a large family (seeds)")
    p.add_argument("--max-gb", type=float, default=prepare.MAX_GB, help="memory watchdog's cap")
    a = p.parse_args()
    prepare.watchdog(a.max_gb)
    {"oracle": oracle, "seeds": seeds}[a.step](a)


if __name__ == "__main__":
    main()
