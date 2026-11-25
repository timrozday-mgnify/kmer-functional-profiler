"""What members cost if each is held as its differences from its cluster's representative row.

For each GA member row (members arm) with a representative GA row in the same family and cluster:
columns that
differ from the nearest of the representative's rows, over the columns the member covers.
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import polars as pl  # noqa: E402
from prep import NONE, OUT  # noqa: E402

fams = pl.read_csv(
    OUT / "bin/families.tsv",
    separator="\t",
    has_header=False,
    new_columns=["acc", "len", "reps", "members"],
)
ga = np.fromfile(OUT / "bin/reps.ga", "u1")
ru, mu = np.fromfile(OUT / "bin/reps.unit", "<u4"), np.fromfile(OUT / "bin/members.unit", "<u4")
rr, mm = np.memmap(OUT / "bin/reps.rows", "u1", "r"), np.memmap(OUT / "bin/members.rows", "u1", "r")
n = cov = diffs = runs = norep = same = 0
res_m = res_r = 0
r0 = m0 = rb = mb = 0
for _acc, length, a, b in fams.iter_rows():
    R = np.asarray(rr[rb : rb + a * length]).reshape(a, length)
    g = ga[r0 : r0 + a] == 1
    R, U = R[g], ru[r0 : r0 + a][g]
    M, V = np.asarray(mm[mb : mb + b * length]).reshape(b, length), mu[m0 : m0 + b]
    res_m += int((M != NONE).sum())
    res_r += int((R != NONE).sum())
    order = np.argsort(U, kind="stable")
    U, R = U[order], R[order]
    lo, hi = np.searchsorted(U, V), np.searchsorted(U, V, "right")
    for k in range(0, b, 20000):
        sl = slice(k, k + 20000)
        has = lo[sl] < hi[sl]
        norep += int((~has).sum())
        X, st, en = M[sl][has], lo[sl][has], hi[sl][has]
        best = None
        # the representative's rows tried one by one: nearly always one or two
        for j in range(int((en - st).max(initial=0))):
            ok = st + j < en
            d = np.where(ok, ((X != NONE) & (R[np.minimum(st + j, len(R) - 1)] != X)).sum(1), 10**6)
            best = d if best is None else np.minimum(best, d)
        if best is not None:
            diffs += int(best.sum())
            same += int((best == 0).sum())
            n += len(best)
            cov += int((X != NONE).sum())
    r0, m0, rb, mb = r0 + a, m0 + b, rb + a * length, mb + b * length
out = {
    "rep_ga_rows": int(ga.sum()),
    "rep_ga_residues": res_r,
    "member_rows": int(len(mu)),
    "member_residues": res_m,
    "member_rows_without_rep_row": norep,
    "members_compared": n,
    "identical_to_rep_row": same / n,
    "mean_columns": cov / n,
    "mean_differences": diffs / n,
    "differing_share": diffs / cov,
}
(OUT / "report/members.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
