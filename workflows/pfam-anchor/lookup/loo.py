"""Members, fairly: simulated reads were cut from member proteins, so a read's own source protein
must not be in the
members index. Rows that only the source protein gave are left out for that read; rows another
member shares stay.

-> report/loo.tsv: per simulated read, is the source cluster among the top ties, and how many ties,
    for
   reps (GA representatives), members (leaky: source protein kept), members_loo, and the two-step
   route
   (top 500 clusters by bits on representatives, then members of those, leave-one-out).
"""

import sys
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze import REP, blosum  # noqa: E402
from prep import OUT, SLICE  # noqa: E402

B = np.full((22, 22), -1, np.int16)
B[:20, :20] = blosum()
B[:, 21] = 0  # a column the target lacks adds nothing, as in the prototype's column score


def sole() -> tuple[np.ndarray, np.ndarray]:
    """Per members row (pack's order): the one protein it comes from, or -1; and its cluster."""
    members = pl.read_parquet(
        SLICE / "members.parquet", columns=["protein_id", "cluster_rep", "sequence"]
    )
    d = (
        pl.read_parquet(OUT / "domains.parquet")
        .join(members, on="protein_id")
        .sort("pfam_accession", "cluster_rep", "protein_id", "domain")
    )
    d = d.filter("above_ga").with_columns(
        key=pl.col("sequence").str.slice(
            pl.col("ali_from") - 1, pl.col("ali_to") - pl.col("ali_from") + 1
        )
    )
    k = ["pfam_accession", "cluster_rep", "key", "hmm_from"]
    g = d.group_by(k, maintain_order=True).agg(
        n=pl.col("protein_id").n_unique(), first=pl.col("protein_id").first()
    )
    return g.select(pl.when(pl.col("n") == 1).then(pl.col("first")).otherwise(-1))[
        "first"
    ].to_numpy(), g["cluster_rep"].to_numpy()


def main() -> None:
    only, cluster = sole()
    units = pl.read_parquet(OUT / "bin/units.parquet")["unit"].to_numpy()
    mu = np.fromfile(OUT / "bin/members.unit", "<u4")
    assert len(only) == len(mu) and (units[mu] == cluster).all(), "row order differs from pack's"
    ru, ga = np.fromfile(OUT / "bin/reps.unit", "<u4"), np.fromfile(OUT / "bin/reps.ga", "u1")
    fams = pl.read_csv(
        OUT / "bin/families.tsv",
        separator="\t",
        has_header=False,
        new_columns=["acc", "len", "reps", "members"],
    )
    rr, mm = (
        np.memmap(OUT / "bin/reps.rows", "u1", "r"),
        np.memmap(OUT / "bin/members.rows", "u1", "r"),
    )
    q = pl.read_parquet(OUT / "queries.parquet").sort("qid")
    qrows, qat = np.fromfile(OUT / "bin/q.rows", "u1"), 0
    starts, length = [], dict(fams.select("acc", "len").iter_rows())
    for a in q["pfam_accession"]:
        starts.append(qat)
        qat += length[a]
    unit_at = {int(u): i for i, u in enumerate(units)}
    out = []
    r0 = m0 = rb = mb = 0
    for acc, n, a, b in fams.iter_rows():
        g = ga[r0 : r0 + a] == 1
        R = np.asarray(rr[rb : rb + a * n]).reshape(a, n)[g].T.copy()
        RU = ru[r0 : r0 + a][g]
        M = np.asarray(mm[mb : mb + b * n]).reshape(b, n).T.copy()
        MU, MO = mu[m0 : m0 + b], only[m0 : m0 + b]
        r0, m0, rb, mb = r0 + a, m0 + b, rb + a * n, mb + b * n
        count = np.stack([(x == R).sum(1) for x in range(20)], 1).astype(float)
        bits = np.zeros((n, 22))
        bits[:, :20] = -np.log2((count + 0.05) / (count.sum(1, keepdims=True) + 1))
        x = q.filter((pl.col("pfam_accession") == acc) & (pl.col("set") == "sim"))
        for qid, pid, su, ident in x.select("qid", "source", "source_unit", "identity").iter_rows():
            row = qrows[starts[qid] : starts[qid] + n]
            cols = np.flatnonzero(row < 20)
            src = unit_at[su]
            keep = pid != MO

            def top(rows_u, s, src=src):
                t = np.unique(rows_u[s == s.max()])
                return src in t, len(t)

            def score(X, kind, cols=cols, row=row):
                if kind == "identity":
                    return (X[cols] == row[cols, None]).sum(0)
                return B[row[cols][:, None], X[cols]].sum(0)

            rb_ = ((R[cols] == row[cols, None]) * bits[cols, row[cols]][:, None]).sum(0)
            best = np.zeros(len(units), np.float32)
            np.maximum.at(best, RU, rb_.astype(np.float32))
            have = np.unique(RU)
            top500 = have[np.argsort(-best[have])[:500]]
            two = np.isin(MU, top500) & keep
            for kind in ("identity", "blosum"):
                rs, ms = score(R, kind), score(M, kind)
                res = {
                    "reps": top(RU, rs),
                    "members": top(MU, ms),
                    "members_loo": top(MU[keep], ms[keep]),
                }
                U2, S2 = (
                    np.concatenate([MU[two], RU[np.isin(RU, top500)]]),
                    np.concatenate([ms[two], rs[np.isin(RU, top500)]]),
                )
                res["two_step_loo"] = top(U2, S2)
                for arm, (right, ties) in res.items():
                    out.append((qid, ident, kind, arm, right, ties, bool((pid == MO).any())))
        print(acc, x.height, flush=True)
    o = pl.DataFrame(
        out, schema=["qid", "identity", "score", "arm", "right", "ties", "sole_row"], orient="row"
    )
    o.write_csv(REP / "loo.tsv", separator="\t")
    o = o.with_columns(credit=pl.col("right") / pl.col("ties"))
    pl.Config.set_tbl_rows(60)
    print(
        o.group_by("score", "arm")
        .agg(
            pl.len(),
            pl.col("right").mean(),
            pl.col("credit").mean(),
            pl.col("ties").median(),
            pl.col("sole_row").mean(),
        )
        .sort("score", "arm")
    )
    print(
        o.filter(pl.col("score") == "identity")
        .pivot(on="arm", index="identity", values="right", aggregate_function="mean")
        .sort("identity")
    )


if __name__ == "__main__":
    main()
