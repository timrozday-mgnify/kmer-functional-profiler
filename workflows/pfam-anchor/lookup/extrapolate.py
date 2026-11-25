"""From the slice to the whole of MGnify proteins 90 -> report/extrapolation.tsv and .json.

Measured inputs: shares of a family's units each option lets through (constant over 1-16x
subsampling), postings
per residue, stage costs (scalar Rust, M4 Pro, one core), scan cost per row (bench).
Estimated: rows per family over the release, hence candidates, CPU and bytes. Assumptions are in
NOTES.
"""

import json
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze import OUT, REP  # noqa: E402

SZ = json.loads((OUT / "sizing.json").read_text())
REPS = 1_661_142_385
REP_RESIDUES = 2.5e11  # STATS.md
FIX = REP_RESIDUES / (
    SZ["sample_mean_length"] * REPS
)  # the sample's reps are longer than the release's
NOTES = (
    "rows = release reps x sample's domains per rep x length correction; "
    "every read domain is scanned against its family"
)
SCAN_NS = {8: 0.385, 16: 0.6, 40: 1.55}  # ns a row, bench, one M4 core (27 GB/s)
SCAN_REST = 0.04  # ns a row per column beyond 40
STAGE_NS = {}  # filled from follow/stage_seconds.tsv: ns a row
POSTING_BYTES = (4.0, 6.0, 8.0)  # low, used, high
HIT_NS = 2.0  # a posting read in order
LOOKUP_NS = 200.0  # a key looked up in a structure too big for cache


def main() -> None:
    for r in pl.read_csv(OUT / "follow/stage_seconds.tsv", separator="\t").iter_rows(named=True):
        STAGE_NS[r["stage"]] = r["seconds"] / r["rows"] * 1e9
    rows = SZ["release_rows"] * FIX
    residues = SZ["release_row_residues"] * FIX
    mean_rows = SZ["rows_per_read_hit_mean"] * FIX
    q = pl.read_csv(OUT / "full_reps/queries.tsv", separator="\t")
    local = q["family_units"].mean()  # real and simulated together; per-set below
    cols = float(q["columns"].mean())
    s = pl.read_csv(REP / "cand_full_reps.tsv", separator="\t").filter(
        (pl.col("bin") == "all") & (pl.col("set") == "real")
    )
    ix = pl.read_csv(REP / "index_full_reps.tsv", separator="\t").with_columns(
        ppr=pl.col("postings_per_residue")
    )
    n_local = 16330  # units per read domain on the slice, real reads
    out = []
    for r in s.iter_rows(named=True):
        o = r["option"]
        share = r["cands_mean"] / n_local
        kind = "scan" if o.startswith("scan") or o.startswith("rare") else "seed"
        cand = (
            share * mean_rows
            if kind == "seed" or "/T" in o or "/N+" in o or "/self" in o or "share" in o
            else r["cands_mean"]
        )
        p = ix.filter(
            pl.col("option") == ("/".join(o.split("/")[:2]) if kind == "seed" else "")
        ).select("ppr")
        ppr = p["ppr"][0] if len(p) else 0.0
        if kind == "seed":
            hits = r["hits_mean"] / n_local * mean_rows
            us = (cols * LOOKUP_NS + hits * HIT_NS + cand * STAGE_NS["columns"]) / 1e3
        else:
            k = 16 if "16" in o else 8 if "bits8" in o else cols
            ns = (
                SCAN_NS[16]
                if k == 16
                else SCAN_NS[8]
                if k == 8
                else SCAN_NS[40] + SCAN_REST * max(0, cols - 40)
            )
            us = mean_rows * ns / 1e3
        out.append(
            {
                "option": o,
                "real_recall": r["best"],
                "cands_local": r["cands_mean"],
                "share_of_family": share,
                "cands_release": cand,
                "index_GB_low": residues * ppr * POSTING_BYTES[0] / 1e9,
                "index_GB": residues * ppr * POSTING_BYTES[1] / 1e9,
                "index_GB_high": residues * ppr * POSTING_BYTES[2] / 1e9,
                "us_first_stage_per_read_domain": us,
            }
        )
    t = pl.DataFrame(out)
    t.write_csv(REP / "extrapolation.tsv", separator="\t")
    info = {
        "length_correction": FIX,
        "release_rows_GA_reps": rows,
        "release_row_residues": residues,
        "mean_rows_per_read_domain": mean_rows,
        "rows_per_read_domain_quantiles": {
            k: v * FIX for k, v in SZ["rows_per_read_hit_quantiles"].items()
        },
        "matrix_GB_1_byte": (residues + 8 * rows) / 1e9,
        "matrix_GB_5_bit": (residues * 5 / 8 + 8 * rows) / 1e9,
        "stage_ns_per_row": STAGE_NS,
        "columns_per_query_mean": cols,
        "scan_ms_per_read_domain": mean_rows * (SCAN_NS[40] + SCAN_REST * max(0, cols - 40)) / 1e6,
        "scan16_ms_per_read_domain": mean_rows * SCAN_NS[16] / 1e6,
        "local_family_units_mean": local,
    }
    (REP / "extrapolation.json").write_text(json.dumps(info, indent=1))
    print(json.dumps(info, indent=1))
    pl.Config.set_tbl_rows(40)
    pl.Config.set_tbl_cols(12)
    want = [
        "golomb4b/s1/h1",
        "golomb5/s1/h1",
        "k5/s1/h1",
        "bits8/s1/h1",
        "bits10/s1/h1",
        "bits14/s1/h1",
        "bits14/s1/h2",
        "golomb4b/s4/h1",
        "bits8/s4/h1",
        "target16/s1/h1",
        "target12/s1/h1",
        "rare3/T8",
        "scan_bits/top100",
        "scan_bits/top500",
        "scan_bits/best-20",
        "scan_bits/selft0.3",
        "scan_ent16/top1000",
        "scan_bits16/top1000",
    ]
    print(t.filter(pl.col("option").is_in(want)).with_columns(pl.exclude("option").round(3)))


if __name__ == "__main__":
    main()
