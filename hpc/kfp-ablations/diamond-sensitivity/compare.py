"""Step 3b's arms against the baseline arm (plan, phase 13): compare.tsv and decision.tsv.

Per benchmark (``name=dir``, each arm's outputs in ``dir/results-dmnd-<arm>``) and sample
(``<sample>/pileup``, ``<sample>/variants``): mates aligned, best unit per mate, per-unit EM
abundance, freebayes records and Floria's multi-strain contigs against the baseline; truth
metrics from the arm's ``summary.tsv`` where there is one (the spike benchmark); cost from
``trace.tsv``. The rule (fixed before the run) picks the cheapest passing ``_F15`` arm.
"""

import argparse
import gzip
from pathlib import Path
from typing import Any

import polars as pl

TRUTH = [
    "read_precision_em", "sibling_caught_em", "sibling_caught_em_parent_target",
    "strain_novel_em", "l1_targets_em", "f1", "consensus_identity_0.9", "consensus_identity_0.95",
]  # fmt: skip
# the rule: (column, baseline-relative floor or absolute floor)
PRECISION_DROP, SIBLING_DROP, MIN_SPEARMAN, MIN_JACCARD = 0.005, 0.02, 0.99, 0.98


def best_unit(pileup: Path) -> pl.DataFrame:
    """Each aligned mate's unit with the largest EM posterior."""
    return (
        pl.scan_parquet(pileup / "assignments.parquet")
        .select("read_name", "mate", "unit", "em_posterior")
        .sort("em_posterior", descending=True)
        .group_by("read_name", "mate")
        .first()
        .select("read_name", "mate", "unit")
        .collect()
    )


def vcf_keys(path: Path) -> set[tuple[str, ...]]:
    if not path.exists():
        return set()
    with gzip.open(path, "rt") as f:
        return {tuple(line.split("\t", 5)[:5]) for line in f if not line.startswith("#")}


def multi_strain(variants: Path) -> int | None:
    info = variants / "floria" / "contig_ploidy_info.tsv"
    if not info.exists():
        return None
    table = pl.read_csv(info, separator="\t")
    return table.filter(pl.col("average_straincount") >= 2).height


def cost(trace: Path) -> dict[str, Any]:
    if not trace.exists():
        return {}
    t = pl.read_csv(trace, separator="\t", infer_schema_length=None, null_values=["-"])
    t = t.filter(pl.col("status").is_in(["COMPLETED", "CACHED"])).with_columns(
        cpu_h=pl.col("realtime").cast(pl.Float64) / 3.6e6 * pl.col("%cpu").cast(pl.Float64) / 100
    )
    step = pl.col("name").str.extract(r"PILEUP:(\w+)")
    after_reference = ~step.is_in(["QUERY", "REFERENCE"]) & step.is_not_null()
    return {
        "diamond_cpu_h": t.filter(step == "DIAMOND")["cpu_h"].sum(),
        "pileup_cpu_h": t.filter(after_reference)["cpu_h"].sum(),
        "call_peak_rss_gb": (t.filter(step == "CALL")["peak_rss"].cast(pl.Float64).max() or 0)
        / 1e9,
    }


def against(arm: Path, base: Path, sample: str) -> dict[str, Any]:
    a, b = best_unit(arm / sample / "pileup"), best_unit(base / sample / "pileup")
    both = a.join(b, on=["read_name", "mate"], suffix="_base")
    ua, ub = (
        pl.read_csv(d / sample / "pileup" / "units.tsv", separator="\t").select(
            "unit", pl.col("em_abundance").fill_null(0.0)
        )
        for d in (arm, base)
    )
    u = ua.join(ub, on="unit", how="full", coalesce=True, suffix="_base").fill_null(0.0)
    x, y = pl.col("em_abundance"), pl.col("em_abundance_base")
    va, vb = (vcf_keys(d / sample / "variants" / "freebayes.vcf.gz") for d in (arm, base))
    union = va | vb
    return {
        "mates_aligned_ratio": a.height / b.height if b.height else None,
        "best_unit_agree": (both["unit"] == both["unit_base"]).mean() if both.height else None,
        "abundance_spearman": u.select(pl.corr(x, y, method="spearman")).item(),
        "abundance_l1": u.select((x / x.sum() - y / y.sum()).abs().sum()).item(),
        "freebayes_records": len(va),
        "freebayes_jaccard": len(va & vb) / len(union) if union else 1.0,
        "floria_multi_strain": multi_strain(arm / sample / "variants"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Step 3b's arms against the baseline arm")
    parser.add_argument("--base", required=True)
    parser.add_argument("--arms", nargs="+", required=True)
    parser.add_argument("--bench", action="append", required=True, help="name=directory")
    parser.add_argument("--out", type=Path, default=Path("."))
    args = parser.parse_args()
    rows = []
    for bench in args.bench:
        name, directory = bench.split("=", 1)
        base = Path(directory) / f"results-dmnd-{args.base}"
        samples = sorted(p.parent.name for p in base.glob("*/pileup"))
        for arm in args.arms:
            out = Path(directory) / f"results-dmnd-{arm}"
            if not out.is_dir():
                print(f"{name}: no {out}, skipped")
                continue
            common = {"bench": name, "arm": arm, **cost(out / "trace.tsv")}
            summary = out / "summary.tsv"
            if summary.exists():
                s = pl.read_csv(summary, separator="\t", infer_schema_length=None).row(
                    0, named=True
                )
                common |= {c: s.get(c) for c in TRUTH}
            for sample in samples:
                rows.append({**common, "sample": sample, **against(out, base, sample)})
    table = pl.DataFrame(rows, infer_schema_length=None)
    table.write_csv(args.out / "compare.tsv", separator="\t")

    # cost is per benchmark run, repeated on its samples' rows
    runs = table.group_by("arm", "bench").agg(pl.col("diamond_cpu_h", "pileup_cpu_h").first())
    per_arm = (
        table.group_by("arm")
        .agg(
            pl.col("abundance_spearman", "freebayes_jaccard").min(),
            pl.col("read_precision_em", "sibling_caught_em").drop_nulls().first(),
        )
        .join(runs.group_by("arm").agg(pl.col("diamond_cpu_h", "pileup_cpu_h").sum()), on="arm")
    )
    ref = per_arm.filter(pl.col("arm") == args.base).row(0, named=True)
    decision = per_arm.with_columns(
        passes=pl.col("arm").str.ends_with("_F15")
        & (pl.col("read_precision_em") >= ref["read_precision_em"] - PRECISION_DROP)
        & (pl.col("sibling_caught_em") >= ref["sibling_caught_em"] - SIBLING_DROP)
        & (pl.col("abundance_spearman") >= MIN_SPEARMAN)
        & (pl.col("freebayes_jaccard") >= MIN_JACCARD)
    ).sort("diamond_cpu_h")
    decision.write_csv(args.out / "decision.tsv", separator="\t")
    with pl.Config(tbl_cols=-1, tbl_rows=-1, tbl_width_chars=200):
        print(decision)
    chosen = decision.filter("passes")
    print("chosen:", chosen["arm"][0] if chosen.height else args.base)


if __name__ == "__main__":
    main()
