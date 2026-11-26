"""`--profile`'s tables and figures, from profile.py's chunks.

  summary --proteins DIR --pairs DIR --runs DIR --out PREFIX   (polars; the workflow runs it)
      PREFIX.proteins.parquet, PREFIX.pairs.parquet   every chunk's, one file each
      PREFIX.lengths.tsv   by protein length: CPU-s of each tool and kfpfam's stages, speed-up,
                           recall
      PREFIX.models.tsv    by the true model's length: true pairs, kfpfam's recall, HMMER's cost
                           per candidate
      PREFIX.run.json      the chunks' totals: per-protein sums, and the batch runs they scale to
  plot PREFIX   (matplotlib: uv run --with matplotlib python profile_plot.py plot results/profile)
      PREFIX.png           cost against protein length, speed-up and recall by length, kfpfam's
                           stages per residue, and by model length HMMER's cost per candidate and
                           recall of the true pairs

Lengths are binned; a bin's speed-up is its hmmsearch CPU over its kfpfam CPU (filter and HMMER).
"""

import argparse
import json
import sys
from pathlib import Path

import polars as pl

LENGTHS = [50, 100, 150, 200, 300, 400, 600, 800, 1200, 1600]
MODELS = [25, 50, 75, 100, 150, 200, 300, 500]
STAGES = ["kf_gate_s", "kf_sort_s", "kf_diagonal_s", "kf_bound_s", "kf_gapped_s", "kf_hmmer_s"]


def binned(column: str, edges: list[int]) -> pl.Expr:
    labels = (
        [f"<{edges[0]}"]
        + [f"{a}-{b}" for a, b in zip(edges, edges[1:], strict=False)]
        + [f">={edges[-1]}"]
    )
    return pl.col(column).cut(edges, labels=labels, left_closed=True).alias("bin")


def by_length(proteins: pl.DataFrame) -> pl.DataFrame:
    out = (
        proteins.with_columns(
            binned("length", LENGTHS),
            true_n=pl.col("true_models").list.len(),
            missed_n=pl.col("kf_missed").list.len(),
            extra_n=pl.col("kf_extra").list.len(),
        )
        .group_by("bin")
        .agg(
            pl.len().alias("proteins"),
            pl.col("length").sum().alias("residues"),
            pl.col("hmmsearch_s").sum(),
            pl.col("kf_filter_s").sum(),
            *[pl.col(s).sum() for s in STAGES],
            pl.col("kf_s").sum(),
            pl.col("kf_candidates").sum(),
            pl.col("true_n").sum().alias("true_pairs"),
            pl.col("missed_n").sum().alias("missed"),
            pl.col("extra_n").sum().alias("extra"),
        )
        .with_columns(
            speedup=pl.col("hmmsearch_s") / pl.col("kf_s"),
            recall=1 - pl.col("missed") / pl.col("true_pairs"),
            hmmsearch_us_per_residue=pl.col("hmmsearch_s") / pl.col("residues") * 1e6,
            kf_us_per_residue=pl.col("kf_s") / pl.col("residues") * 1e6,
        )
        .sort("bin")
    )
    return out


def by_model(pairs: pl.DataFrame) -> pl.DataFrame:
    return (
        pairs.with_columns(binned("model_length", MODELS))
        .group_by("bin")
        .agg(
            pl.col("true").sum().alias("true_pairs"),
            (pl.col("true") & pl.col("found")).sum().alias("found"),
            (pl.col("true") & ~pl.col("candidate")).sum().alias("not_candidates"),
            pl.col("candidate").sum().alias("candidates"),
            pl.col("hmmer_s").sum(),
            pl.col("hmmer_s").filter(pl.col("true")).median().alias("hmmer_s_true_median"),
            pl.col("hmmer_s").filter(~pl.col("true")).median().alias("hmmer_s_false_median"),
        )
        .with_columns(
            recall=pl.col("found") / pl.col("true_pairs"),
            hmmer_us_per_candidate=pl.col("hmmer_s") / pl.col("candidates") * 1e6,
        )
        .sort("bin")
    )


def summary(a: argparse.Namespace) -> None:
    proteins = pl.concat([pl.read_parquet(f) for f in sorted(a.proteins.glob("*.parquet"))])
    pairs = pl.concat([pl.read_parquet(f) for f in sorted(a.pairs.glob("*.parquet"))])
    proteins.write_parquet(f"{a.out}.proteins.parquet", compression="zstd")
    pairs.write_parquet(f"{a.out}.pairs.parquet", compression="zstd")
    by_length(proteins).write_csv(f"{a.out}.lengths.tsv", separator="\t", float_precision=4)
    by_model(pairs).write_csv(f"{a.out}.models.tsv", separator="\t", float_precision=4)
    runs = [json.loads(f.read_text()) for f in sorted(a.runs.glob("*.json"))]
    total = {
        "chunks": len(runs),
        "proteins": sum(r["proteins"] for r in runs),
        "residues": sum(r["residues"] for r in runs),
        "hmmsearch_s": sum(r["hmmsearch_s"]["batch"] for r in runs),
        "kfpfam_filter_s": sum(r["kfpfam_filter_s"]["per_protein_sum"] for r in runs),
        "kfpfam_hmmer_s": sum(r["kfpfam_hmmer_s"]["batch"] for r in runs),
        "true_pairs": sum(r["true_pairs"] for r in runs),
        "kfpfam_missed": sum(r["kfpfam_missed"] for r in runs),
        "kfpfam_extra": sum(r["kfpfam_extra"] for r in runs),
        "hmmsearch_scale": [r["hmmsearch_s"]["scale"] for r in runs],
        "kfpfam_hmmer_scale": [r["kfpfam_hmmer_s"]["scale"] for r in runs],
        "kernel": sorted({r["kfpfam_filter_s"]["kernel"] for r in runs}),
    }
    total["speedup"] = total["hmmsearch_s"] / (total["kfpfam_filter_s"] + total["kfpfam_hmmer_s"])
    Path(f"{a.out}.run.json").write_text(json.dumps(total, indent=1) + "\n")
    print(json.dumps({k: v for k, v in total.items() if not k.endswith("scale")}))


def plot(a: argparse.Namespace) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    proteins = pl.read_parquet(f"{a.prefix}.proteins.parquet")
    pairs = pl.read_parquet(f"{a.prefix}.pairs.parquet")
    lengths, models = by_length(proteins), by_model(pairs)
    blue, orange = "#2a6fb0", "#d0702a"
    fig, ax = plt.subplots(2, 2, figsize=(13, 9.5))

    p = ax[0, 0]
    sample = proteins.sample(min(20000, proteins.height), seed=0)
    p.scatter(
        sample["length"], sample["hmmsearch_s"], s=2, alpha=0.15, color=blue, label="hmmsearch"
    )
    p.scatter(sample["length"], sample["kf_s"], s=2, alpha=0.15, color=orange, label="kfpfam")
    med = (
        proteins.with_columns(binned("length", LENGTHS))
        .group_by("bin")
        .agg(pl.col("length").median(), pl.col("hmmsearch_s").median(), pl.col("kf_s").median())
        .sort("length")
    )
    p.plot(med["length"], med["hmmsearch_s"], color=blue, lw=2)
    p.plot(med["length"], med["kf_s"], color=orange, lw=2)
    p.set(
        xscale="log",
        yscale="log",
        xlabel="protein length (residues)",
        ylabel="CPU-s per protein",
        title="Cost per protein (lines: median per length bin)",
    )
    p.legend(markerscale=6)

    p = ax[0, 1]
    x = range(lengths.height)
    p.bar(x, lengths["speedup"], color=orange)
    p.set_xticks(list(x), lengths["bin"].cast(pl.Utf8).to_list(), rotation=45, ha="right")
    p.set(
        ylabel="hmmsearch CPU / kfpfam CPU",
        xlabel="protein length",
        title="Speed-up and recall by protein length",
    )
    r = p.twinx()
    r.plot(list(x), lengths["recall"], "k.-")
    r.set_ylabel("kfpfam recall of hmmsearch's pairs")
    r.set_ylim(min(0.95, lengths["recall"].drop_nulls().min() or 1), 1.001)

    p = ax[1, 0]
    bottom = [0.0] * lengths.height
    per = lambda c: (lengths[c] / lengths["residues"] * 1e6).to_list()  # noqa: E731
    for stage in STAGES:
        p.bar(x, per(stage), bottom=bottom, label=stage.removeprefix("kf_").removesuffix("_s"))
        bottom = [b + v for b, v in zip(bottom, per(stage), strict=True)]
    p.plot(list(x), per("hmmsearch_s"), "k.--", label="hmmsearch")
    p.set_xticks(list(x), lengths["bin"].cast(pl.Utf8).to_list(), rotation=45, ha="right")
    p.set(
        ylabel="CPU-µs per residue",
        xlabel="protein length",
        title="kfpfam's stages by protein length, per residue",
    )
    p.legend(fontsize=8)

    p = ax[1, 1]
    x = range(models.height)
    p.bar(x, models["hmmer_us_per_candidate"], color=orange)
    p.set_xticks(list(x), models["bin"].cast(pl.Utf8).to_list(), rotation=45, ha="right")
    p.set(
        ylabel="kfpfam's HMMER, µs per candidate",
        xlabel="model length (columns)",
        title="By model length: HMMER per candidate, and recall of true pairs",
    )
    r = p.twinx()
    r.plot(list(x), models["recall"], "k.-")
    r.set_ylabel("recall")
    r.set_ylim(min(0.95, models["recall"].drop_nulls().min() or 1), 1.001)

    fig.tight_layout()
    fig.savefig(f"{a.prefix}.png", dpi=130)
    print(f"{a.prefix}.png")


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = p.add_subparsers(required=True)
    s = sub.add_parser("summary")
    s.set_defaults(run=summary)
    for name in ("proteins", "pairs", "runs"):
        s.add_argument(f"--{name}", type=Path, required=True)
    s.add_argument("--out", required=True)
    x = sub.add_parser("plot")
    x.set_defaults(run=plot)
    x.add_argument("prefix")
    a = p.parse_args()
    a.run(a)


if __name__ == "__main__":
    sys.exit(main())
