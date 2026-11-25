"""Size of an annotated MGnify90 index, by how far targets are collapsed -> report/annot_sizes.tsv.
Estimates."""

import json
import sys
from pathlib import Path

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[1]))
from analyze import OUT, REP, queries  # noqa: E402
from zymo import family_units, scores  # noqa: E402

E = json.loads((REP / "extrapolation.json").read_text())
G = json.loads((REP / "annot_genomes.json").read_text())
rows, res = E["release_rows_GA_reps"], E["release_row_residues"]
c = pl.read_csv(REP / "annot_collapse.tsv", separator="\t").filter(pl.col("reads") == "real")
tr = pl.read_csv(REP / "annot_transfer.tsv", separator="\t")
edges = [0.4, 0.6, 0.8, 0.95]
labels = ["<.4", ".4-.6", ".6-.8", ".8-.95", "≥.95"]
by = {
    r["bin"]: (r["recall"], r["precision"])
    for r in tr.filter((pl.col("label") == "ko") & pl.col("in_pfam"))
    .group_by("bin")
    .agg(
        recall=pl.col("right").sum() / pl.col("has").sum(),
        precision=pl.col("right").sum() / pl.col("called").sum(),
    )
    .iter_rows(named=True)
}
print({k: (round(a, 3), round(b, 3)) for k, (a, b) in sorted(by.items())})
q = (
    queries()
    .select("qid", "set", "self_score", "overall", "pfam_accession")
    .filter(pl.col("set") == "real")
)
units, accs = family_units()
sc = scores("bin")


def expected(t: float | None) -> tuple[float, float]:
    """Expected KO recall and precision of the real reads, from their similarity to the best kept
    target."""
    masks = None
    if t:
        raw, at, masks = np.fromfile(OUT / f"bin/collapse.{t:g}.units", "<u4"), 0, []
        for u in units:
            n = int(raw[at])
            masks.append(np.isin(u, raw[at + 1 : at + 1 + n]))
            at += 1 + n
    f = []
    for qid, _, ss, best, acc in q.iter_rows():
        x = sc(qid)[1]
        f.append(float(x[masks[accs.index(acc)]].max() if t else best) / ss)
    k = np.digitize(f, edges)
    return tuple(float(np.mean([by[labels[i]][j] for i in k])) for j in (0, 1))


out = []
for t, rl, sl in [(None, 1.0, 1.0)] + [
    (r["t"], r["rows_left"], r["residues_left"])
    for r in c.sort("t", descending=True).iter_rows(named=True)
]:
    rec, pre = expected(t)
    n, s = rows * rl, res * sl
    out.append(
        {
            "collapse": t or 0.0,
            "rows_M": n / 1e6,
            "GB_5bit": (s * 5 / 8 + 8 * n) / 1e9,
            "GB_1byte": (s + 8 * n) / 1e9,
            "KO_rows_only_5bit": (s * 5 / 8 + 8 * n) / 1e9 * G["p_ko_given_pfam"],
            "CAZy_rows_only_5bit": (s * 5 / 8 + 8 * n) / 1e9 * G["p_cazy_given_pfam"],
            "label_GB_per_scheme": 2 * n / 1e9,
            "scan_ms": E["scan_ms_per_read_domain"] * rl,
            "ko_recall": rec,
            "ko_precision": pre,
        }
    )
o = pl.DataFrame(out)
o.write_csv(REP / "annot_sizes.tsv", separator="\t")
pl.Config.set_tbl_cols(12)
pl.Config.set_tbl_width_chars(220)
print(o.with_columns(pl.all().round(3)))
