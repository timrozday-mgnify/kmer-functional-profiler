"""Truth for the Zymo mocks' reads: read -> gene (by sequence) -> the gene's MGnify90 cluster -> was
it called?

A read's gene is found among the standard's genes holding a domain of the read's family, by its
aligned
residues (at most 5% mismatches, no gaps). A gene's clusters are the representatives its whole
protein
scores best against, over all 630,402 (Smith-Waterman, exhaustive). Calls come from `follow`'s
calls.tsv.
"""

import sys
from pathlib import Path

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
from analyze import REP, blosum  # noqa: E402
from prep import AA, CODE, OUT, ZYMO  # noqa: E402


def family_units() -> list[np.ndarray]:
    """As Rust's `family_units` over both arms: each family's units, sorted."""
    fams = pl.read_csv(
        OUT / "bin/families.tsv",
        separator="\t",
        has_header=False,
        new_columns=["acc", "len", "reps", "members"],
    )
    out = [[] for _ in range(fams.height)]
    for arm in ("reps", "members"):
        unit = np.fromfile(OUT / f"bin/{arm}.unit", "<u4")
        at = np.concatenate([[0], np.cumsum(fams[arm].to_numpy())])
        for i in range(fams.height):
            out[i].append(unit[at[i] : at[i + 1]])
    return [np.unique(np.concatenate(x)) for x in out], fams["acc"].to_list()


def scores(d: str):
    """truth.fam of one bin directory: a function query -> (units, scores)."""
    units, accs = family_units()
    fam = np.fromfile(OUT / d / "q.fam", "<u4")
    off = np.concatenate([[0], np.cumsum([len(units[f]) for f in fam])])
    s = np.memmap(OUT / d / "truth.fam", "<i2", "r")
    return lambda i: (units[fam[i]], s[off[i] : off[i + 1]])


def find(residues: np.ndarray, gene: np.ndarray) -> float:
    """Best ungapped identity of the read's residues within the gene."""
    n = len(residues)
    if len(gene) < n:
        return 0.0
    w = np.lib.stride_tricks.sliding_window_view(gene, n)
    return float((w == residues).mean(1).max())


def main() -> None:
    best = np.dtype([("s", "<i2"), ("n", "<u4")])
    genes = pl.read_parquet(OUT / "zymo_queries.parquet")
    gbest = np.fromfile(OUT / "binz/truth.best", best)
    d = np.diag(blosum())
    code = {a: i for i, a in enumerate(AA)}
    gene_scores, read_scores = scores("binz"), scores("bin")
    q = (
        pl.read_parquet(OUT / "queries.parquet")
        .sort("qid")
        .filter((pl.col("set") == "real") & pl.col("run").is_in(list(ZYMO)))
    )
    coded = [CODE[np.frombuffer(s.encode(), np.uint8)] for s in genes["sequence"]]
    by = {}
    for i, (b, f) in enumerate(genes.select("bundle", "pfam_accession").iter_rows()):
        by.setdefault((b, f), []).append(i)
    g_rows = []
    for i, (g, n, s) in enumerate(genes.select("genome", "gene", "sequence").iter_rows()):
        self_scores = sum(d[code[c]] for c in s if c in code)
        g_rows.append(
            (
                i,
                g,
                n,
                len(s),
                int(gbest["s"][i]) / self_scores,
                int(gene_scores(i)[1].max()) == int(gbest["s"][i]),
            )
        )
    gt = pl.DataFrame(
        g_rows, schema=["gid", "genome", "gene", "length", "similarity", "in_family"], orient="row"
    )
    gt.write_csv(REP / "zymo_genes.tsv", separator="\t")
    print(
        gt.select(
            pl.len(),
            (pl.col("similarity") > 0.9).mean().alias("over .9"),
            pl.col("similarity").median(),
            pl.col("in_family").mean(),
        )
    )

    truth = {}
    for qid, run, f, s in q.select("qid", "run", "pfam_accession", "sequence").iter_rows():
        r = CODE[np.frombuffer(s.encode(), np.uint8)]
        c = [(find(r, coded[i]), i) for i in by.get((ZYMO[run], f), [])]
        top = max(c, default=(0, None))[0]
        if top >= 0.95:
            right = set()
            for x, i in c:
                if x == top:
                    u, sc = gene_scores(i)
                    right |= set(u[sc == sc.max()].tolist())
            truth[qid] = (right, gt["similarity"][[i for x, i in c if x == top][0]], top)
    print(f"{len(truth)} of {q.height} reads have a gene")
    out = []
    for qid, (right, sim, _top) in truth.items():
        u, sc = read_scores(qid)
        ties = u[sc == sc.max()]
        out.append((qid, "exhaustive, family", len(right & set(ties.tolist())) > 0, len(ties), sim))
    p = OUT / "follow/calls.tsv"
    if p.exists():
        c = pl.read_csv(p, separator="\t", schema_overrides={"units": pl.String})
        for opt, qid, _score, ties, units in c.select(
            "option", "qid", "score", "ties", "units"
        ).iter_rows():
            if qid in truth:
                called = set(map(int, units.split(","))) if units else set()
                out.append((qid, opt, len(truth[qid][0] & called) > 0, ties, truth[qid][1]))
    o = pl.DataFrame(
        out, schema=["qid", "method", "right", "ties", "gene_similarity"], orient="row"
    ).join(q.select(pl.col("qid").cast(pl.Int64), "run"), on="qid")
    o.write_csv(REP / "zymo_calls.tsv", separator="\t")
    o = o.with_columns(
        gene=pl.when(pl.col("gene_similarity") >= 0.9).then(pl.lit("≥.9")).otherwise(pl.lit("<.9")),
        credit=pl.col("right") / pl.col("ties").clip(1),
    )
    for k in (["method"], ["method", "gene"], ["method", "run"]):
        t = (
            o.group_by(k)
            .agg(pl.len(), pl.col("right").mean(), pl.col("credit").mean(), pl.col("ties").median())
            .sort(k)
        )
        t.write_csv(REP / f"zymo_{'_'.join(k)}.tsv", separator="\t")
        print(t)


if __name__ == "__main__":
    main()
