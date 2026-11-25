"""Release-wide sizing: what a MGnify90 representative carries, and how big the family is that a
read lands in.

Measured: 148,471 release representatives (the head of mgy_clusters.fa.gz) against full Pfam 38.2 at
GA
(kfpfam data/truth/rest.*), and the families of the 11 B11 libraries' read domains (Pfam 38.2, nail
arm).
Estimated: rows per family over the release, as release representatives x share of rows; families of
fewer
than 1 row in the sample get the mean of the unseen mass.
"""

import glob
import json
from pathlib import Path

import numpy as np
import polars as pl

HOME = Path.home() / "Documents"
OUT = HOME / "kmer-functional-profiler/data/pfam-anchor/lookup"
REPS = 1_661_142_385  # release 2026_07
TRUTH = HOME / "kfpfam/data/truth"


def main() -> None:
    n = res = 0
    with open(TRUTH / "rest.faa") as f:
        for line in f:
            if line[0] == ">":
                n += 1
            else:
                res += len(line) - 1
    rows = []
    with open(TRUTH / "rest.domtbl") as f:
        for line in f:
            if line[0] != "#":
                c = line.split(maxsplit=22)
                rows.append(
                    (
                        c[0],
                        int(c[4][2:].split(".")[0]),
                        int(c[15]),
                        int(c[16]),
                        int(c[17]),
                        int(c[18]),
                    )
                )
    d = pl.DataFrame(
        rows, schema=["protein", "pfam", "hmm_from", "hmm_to", "ali_from", "ali_to"], orient="row"
    )
    cols = (d["hmm_to"] - d["hmm_from"] + 1).sum()
    ali = (d["ali_to"] - d["ali_from"] + 1).sum()
    per_rep = d.height / n
    out = {
        "sample_reps": n,
        "sample_mean_length": res / n,
        "sample_reps_with_domain": d["protein"].n_unique() / n,
        "domains_per_rep": per_rep,
        "domains_per_rep_with_domain": d.height / d["protein"].n_unique(),
        "columns_per_domain": cols / d.height,
        "residues_per_domain": ali / d.height,
        "share_of_residues_in_domains": ali / res,
        "families_seen": d["pfam"].n_unique(),
        "release_rows": per_rep * REPS,
        "release_row_residues": ali / n * REPS,
    }
    size = d.group_by("pfam").agg(pl.len().alias("k")).with_columns(rows=pl.col("k") / n * REPS)
    reads = []
    for p in sorted(
        glob.glob(
            str(
                HOME
                / "nar-review-response/data/B11/results_main_nail/*/function-summary/pfam"
                / "*_alignments.parquet"
            )
        )
    ):
        x = pl.scan_parquet(p).group_by("pfam").agg(pl.len().alias("reads")).collect()
        reads.append(
            x.with_columns(
                run=pl.lit(Path(p).name.split("_")[0]),
                pfam=pl.col("pfam").str.slice(2).cast(pl.Int64),
            )
        )
    r = pl.concat(reads)
    unseen = 0.5 / n * REPS  # a family with no domain in the sample: at most about this
    j = r.join(size, on="pfam", how="left").with_columns(pl.col("rows").fill_null(unseen))
    out["read_domains"] = int(j["reads"].sum())
    out["read_domains_per_library"] = {k[0]: int(v["reads"].sum()) for k, v in r.group_by("run")}
    s = j.group_by("pfam").agg(pl.col("reads").sum(), pl.col("rows").first()).sort("rows")
    w = s["reads"].to_numpy() / s["reads"].sum()
    v = s["rows"].to_numpy()
    out["rows_per_read_hit_mean"] = float((w * v).sum())
    cw = np.cumsum(w)
    out["rows_per_read_hit_quantiles"] = {
        q: float(v[np.searchsorted(cw, q)]) for q in (0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
    }
    out["rows_max_family"] = float(v.max())
    out["read_domains_in_unseen_families"] = float(w[v <= unseen].sum())
    # the same without rows weighting by reads, for per-family streaming: total rows of families
    # that reads hit
    for run, g in j.group_by("run"):
        out.setdefault("rows_of_hit_families_by_library", {})[run[0]] = float(g["rows"].sum())
    s.write_parquet(OUT / "family_sizes.parquet")
    (OUT / "sizing.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
