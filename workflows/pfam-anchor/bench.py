"""Phase 6's benchmark (docs/pfam-anchor-plan.md): recall and precision by depth x AAI, every
miss attributed to the stage that lost it, and the cost per stage.

- ``score``: per read set and ``call`` arm of the confirmed reference ->
  ``bench_reads.parquet`` (one row per true read pair: depth, AAI, the stage it reached),
  ``bench_recall.tsv`` (per depth x AAI cell: read pairs and strains with Wilson intervals,
  recall against all pairs and against the resolvable ones), ``bench_misses.tsv`` (the
  waterfall: pairs lost at each stage per cell, with example reads), ``bench_units.tsv``
  (unit precision) and ``bench_background.tsv`` (real background reads against DIAMOND).
- ``cost``: Nextflow's trace (``--trace``) and the scripts' timers -> ``bench_cost.tsv``:
  per stage CPU-s, wall-s, peak RSS, core-ms per read pair and per Mbp; index sizes on disk.

A true read pair's stages, in order; it is counted at the first it fails:

=================  ==========================================================================
``no_coding``      the pair holds no coding base (flank only): no method can place it
``no_orf``         neither mate gave an ORF of ``--min-orf`` residues
``no_candidate``   neither route listed the truth unit (anchored ``units``; tier 1 or tier 2
                   leaves' units); split by whether the pair had a Pfam hit
``not_confirmed``  the truth unit missed the confirmed reference
``not_aligned``    ``call`` has no alignment of either mate to the truth unit
``not_kept``       aligned, but no alignment to the truth unit kept
``em_near``        the EM's best unit is another unit in the truth unit's component
``em_novel``       the EM gave the pair to the component's novel class
``em_other``       the EM's best unit is elsewhere
``right``          the EM's best unit is the truth unit
=================  ==========================================================================

*Resolvable*: DIAMOND (every mate against the whole reference, top 10%) puts the truth unit
among its best-scoring units for the pair. A pair that isn't is out of reach of any
reference-based method at this read length, so recall is reported against both denominators.
Wilson intervals are over strains (units) as well as pairs, as a strain's pairs aren't
independent.
"""

import argparse
import json
import math
from pathlib import Path

import anchor
import polars as pl
import prepare

DEPTH_EDGES = [0.75, 1.5, 3.5, 7.5, 20.0]  # bins around the grid 0.5, 1, 2, 5, 10, 30
DEPTH_LABELS = ["0.5", "1", "2", "5", "10", "30"]
STAGES = ["no_coding", "no_orf", "no_candidate", "not_confirmed", "not_aligned", "not_kept",
          "em_near", "em_novel", "em_other", "right"]  # fmt: skip
EXAMPLES = 5  # example read names per miss bucket


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float | None, float | None]:
    if n == 0:
        return None, None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def unit_of(data: Path) -> pl.DataFrame:
    return pl.read_parquet(data / "holdout/members.parquet", columns=["protein_id", "cluster_rep"])


def pos_units(data: Path) -> pl.DataFrame:
    """The whole-reference tree index's leaf positions -> units (``anchor.unanchored``'s)."""
    return (
        pl.read_parquet(data / "whole/index/proteins.parquet").filter(pl.col("pos") >= 0)
        .join(unit_of(data), on="protein_id").select("pos", unit="cluster_rep").unique()
    )  # fmt: skip


def candidates(data: Path, name: str, tag: str) -> pl.DataFrame:
    """(read, unit) listed by either route: anchored ``units``, tier 1 and tier 2 leaves."""
    anchored = (
        pl.read_parquet(data / name / f"anchored{tag}.parquet", columns=["read", "units"])
        .explode("units", empty_as_null=True).drop_nulls()
        .select(pl.col("read").cast(pl.Int64), unit="units")
    )  # fmt: skip
    pu = pos_units(data)
    tiers = [
        pl.read_parquet(f, columns=["read", "pos"])
        .join(pu, on="pos")
        .select(pl.col("read").cast(pl.Int64), "unit")
        for f in (data / f"whole/q-{name}/reads.parquet", data / f"whole/t2-{name}/reads.parquet")
        if f.exists()
    ]
    return pl.concat([anchored, *tiers]).unique()


def diamond_best(data: Path, name: str) -> pl.DataFrame:
    """Per read pair, the units at DIAMOND's best bitscore (either mate)."""
    dia = pl.read_parquet(
        data / name / "diamond.parquet", columns=["read", "protein_id", "bitscore"]
    )
    return (
        dia.join(unit_of(data), on="protein_id")
        .filter(pl.col("bitscore") == pl.col("bitscore").max().over("read"))
        .select("read", unit="cluster_rep").unique()
    )  # fmt: skip


def best_calls(arm_dir: Path, refs: pl.DataFrame) -> pl.DataFrame:
    """Per read name: the kept alignment with the highest EM posterior (its unit and component),
    and whether the novel class outbid it."""
    return (
        pl.read_parquet(arm_dir / "assignments.parquet").filter("kept")
        .sort("em_posterior", descending=True).unique("read_name", keep="first")
        .join(refs, on="unit")
        .select("read_name", called="cluster_rep", called_comp="component",
                novel=pl.col("em_novel") > pl.col("em_posterior"))
    )  # fmt: skip


def stage_of(df: pl.DataFrame) -> pl.DataFrame:
    s = (
        pl.when(pl.col("coding") <= 0).then(pl.lit("no_coding"))
        .when(~pl.col("has_orf")).then(pl.lit("no_orf"))
        .when(~pl.col("candidate")).then(pl.lit("no_candidate"))
        .when(~pl.col("confirmed")).then(pl.lit("not_confirmed"))
        .when(~pl.col("aligned")).then(pl.lit("not_aligned"))
        .when(~pl.col("kept")).then(pl.lit("not_kept"))
        .when(pl.col("called") == pl.col("unit")).then(pl.lit("right"))
        .when(pl.col("novel")).then(pl.lit("em_novel"))
        .when(pl.col("called_comp") == pl.col("comp")).then(pl.lit("em_near"))
        .otherwise(pl.lit("em_other"))
    )  # fmt: skip
    return df.with_columns(stage=s)


def score(a: argparse.Namespace) -> None:
    tag = anchor.tag_of(a)
    c_name = anchor.confirmed_dir(a)
    reads_all, bg_rows, unit_rows = [], [], []
    for name in prepare.SETS:
        d = a.data / name
        c = d / c_name
        refs = pl.read_csv(c / "refs.tsv", separator="\t")
        comp = refs.select("cluster_rep", "component").unique("cluster_rep")
        refs = refs.select("unit", "cluster_rep", "component").unique("unit")
        truth = pl.read_csv(d / "truth.tsv", separator="\t").select("strain", "depth")
        truth = truth.join(pl.read_parquet(d / "truth.parquet"), on="strain")
        names = anchor.read_names(d / "reads_1.fq.gz")
        pairs = pl.read_parquet(d / "reads.parquet", columns=["read", "strain", "coding"])
        pairs = pairs.with_columns(read_name=pl.Series(names))
        orf = pl.read_parquet(d / "orfs.parquet", columns=["read"]).unique()
        hit = pl.read_parquet(d / "hits.parquet", columns=["protein_id"]).unique()
        hit = (
            pl.read_parquet(d / "orfs.parquet", columns=["orf_id", "read"])
            .join(hit.rename({"protein_id": "orf_id"}), on="orf_id").select("read").unique()
        )  # fmt: skip
        cand = candidates(a.data, name, tag)
        confirmed = pl.read_parquet(d / f"{c_name}.parquet", columns=["unit"]).unique()
        best_dia = diamond_best(a.data, name)
        true = (
            pairs.join(truth.select("strain", "unit", "identity_bin", "depth"), on="strain")
            .join(comp.rename({"cluster_rep": "unit", "component": "comp"}), on="unit", how="left")
            .with_columns(
                has_orf=pl.col("read").is_in(orf["read"].implode()),
                pfam_hit=pl.col("read").is_in(hit["read"].implode()),
                confirmed=pl.col("unit").is_in(confirmed["unit"].implode()),
            )
            .join(cand.with_columns(candidate=pl.lit(True)), on=["read", "unit"], how="left")
            .join(best_dia.with_columns(resolvable=pl.lit(True)), on=["read", "unit"], how="left")
            .with_columns(pl.col("candidate", "resolvable").fill_null(False))
            .with_columns(
                depth_bin=pl.col("depth").cut(DEPTH_EDGES, labels=DEPTH_LABELS),
                coding=pl.col("coding").fill_null(0),
            )
        )  # fmt: skip
        for arm_dir in sorted(c.glob("call_*/")):
            arm = arm_dir.name.removeprefix("call_")
            asg = pl.read_parquet(
                arm_dir / "assignments.parquet", columns=["read_name", "unit", "kept"]
            )
            asg = asg.join(refs.select("unit", "cluster_rep"), on="unit")
            to_truth = (
                true.select("read_name", "unit")
                .join(
                    asg.rename({"unit": "u"}),
                    left_on=["read_name", "unit"],
                    right_on=["read_name", "cluster_rep"],
                )  # fmt: skip
                .group_by("read_name")
                .agg(aligned=pl.lit(True), kept=pl.col("kept").any())
            )
            df = (
                true.join(to_truth, on="read_name", how="left")
                .join(best_calls(arm_dir, refs), on="read_name", how="left")
                .with_columns(pl.col("aligned", "kept", "novel").fill_null(False))
            )
            df = stage_of(df).with_columns(set=pl.lit(name), arm=pl.lit(arm))
            reads_all.append(df.drop("comp", "called_comp"))
            # units: called at >= half a pair, true or false
            units = pl.read_csv(arm_dir / "units.tsv", separator="\t")
            called = set(units.filter(pl.col("em_pairs") >= 0.5)["cluster_rep"].to_list())
            truth_units = set(truth["unit"].to_list())
            bg_units = set()
            if name == "strains":
                bg = pairs.filter(pl.col("strain") < 0)
                if bg.height:
                    bgc = bg.join(best_calls(arm_dir, refs), on="read_name", how="left").join(
                        best_dia.group_by("read").agg(dia=pl.col("unit")), on="read", how="left"
                    )
                    bgc = bgc.with_columns(
                        mapped=pl.col("called").is_not_null(),
                        dia_hit=pl.col("dia").is_not_null(),
                        agree=pl.col("dia").list.contains(pl.col("called")).fill_null(False),
                    )
                    bg_units = set(bgc.filter("agree")["called"].to_list())
                    m, h = int(bgc["mapped"].sum()), int(bgc["dia_hit"].sum())
                    ag = int(bgc["agree"].sum())
                    only_mapped = bgc.filter(pl.col("mapped") & ~pl.col("dia_hit")).height
                    only_dia = bgc.filter(~pl.col("mapped") & pl.col("dia_hit")).height
                    bg_rows.append({
                        "arm": arm, "background_pairs": bg.height, "mapped": m,
                        "diamond_hit": h, "mapped_and_diamond_agree": ag,
                        "mapped_without_diamond_hit": only_mapped,
                        "diamond_hit_not_mapped": only_dia,
                        "agreement_precision": ag / m if m else None,
                        "agreement_recall": ag / h if h else None,
                    })  # fmt: skip
            tp = len(called & truth_units)
            fp_supported = len((called - truth_units) & bg_units)
            fp = len(called - truth_units - bg_units)
            unit_rows.append({
                "set": name, "arm": arm, "called": len(called), "true_called": tp,
                "true_units": len(truth_units), "background_supported": fp_supported,
                "false": fp, "precision": tp / max(tp + fp, 1),
                "precision_lo": wilson(tp, tp + fp)[0], "precision_hi": wilson(tp, tp + fp)[1],
            })  # fmt: skip

    reads = pl.concat(reads_all, how="diagonal_relaxed")
    reads.write_parquet(a.data / "bench_reads.parquet")
    keys = ["set", "arm", "depth_bin", "identity_bin"]
    rows = []
    for k, g in reads.group_by(keys, maintain_order=True):
        n, right = g.height, int((g["stage"] == "right").sum())
        res = g.filter("resolvable")
        nr, rr = res.height, int((res["stage"] == "right").sum())
        strains = g.group_by("strain").agg(found=(pl.col("stage") == "right").any())
        ns, fs = strains.height, int(strains["found"].sum())
        rows.append({
            **dict(zip(keys, k, strict=True)), "pairs": n, "recall": right / n,
            "recall_lo": wilson(right, n)[0], "recall_hi": wilson(right, n)[1],
            "resolvable": nr, "recall_resolvable": rr / nr if nr else None,
            "ceiling": nr / n, "strains": ns, "strains_found": fs / ns,
            "strains_lo": wilson(fs, ns)[0], "strains_hi": wilson(fs, ns)[1],
        })  # fmt: skip
    recall = pl.DataFrame(rows).sort(keys)
    recall.write_csv(a.data / "bench_recall.tsv", separator="\t", float_precision=4)
    misses = (
        reads.group_by([*keys, "stage", "pfam_hit"]).agg(
            pairs=pl.len(), examples=pl.col("read_name").head(EXAMPLES).str.join(",")
        )
        .with_columns(share=pl.col("pairs") / pl.col("pairs").sum().over(keys))
        .with_columns(order=pl.col("stage").replace_strict(STAGES, list(range(len(STAGES)))))
        .sort([*keys, "order", "pfam_hit"]).drop("order")
    )  # fmt: skip
    misses.write_csv(a.data / "bench_misses.tsv", separator="\t", float_precision=4)
    pl.DataFrame(unit_rows).write_csv(a.data / "bench_units.tsv", separator="\t", float_precision=4)
    if bg_rows:
        pl.DataFrame(bg_rows).write_csv(
            a.data / "bench_background.tsv", separator="\t", float_precision=4
        )
    overall = (
        reads.group_by("set", "arm", "stage").agg(pairs=pl.len())
        .with_columns(share=pl.col("pairs") / pl.col("pairs").sum().over("set", "arm"))
        .pivot("stage", index=["set", "arm"], values="share").fill_null(0)
    )  # fmt: skip
    with pl.Config(tbl_rows=-1, tbl_cols=-1, float_precision=3, tbl_width_chars=250):
        print(overall.select("set", "arm", *[s for s in STAGES if s in overall.columns]))
        print(pl.DataFrame(unit_rows))


def cost(a: argparse.Namespace) -> None:
    """Per Nextflow task: CPU-s, wall-s, peak RSS, normalised per read pair and per Mbp of
    reads; the scripts' own per-step timers beside them; index sizes on disk."""
    pairs = bases = 0
    for name in prepare.SETS:
        r = pl.read_parquet(a.data / name / "reads.parquet", columns=["pair_bases"])
        pairs, bases = pairs + r.height, bases + int(r["pair_bases"].sum())
    rows = []
    if a.trace and a.trace.exists():
        tr = pl.read_csv(a.trace, separator="\t", null_values=["-"], infer_schema_length=0)
        for r in tr.filter(pl.col("status").is_in(["COMPLETED", "CACHED"])).iter_rows(named=True):
            wall = float(r["realtime"] or 0) / 1000
            cpu = wall * float(r["%cpu"] or 0) / 100 if r["%cpu"] else None
            rss = float(r["peak_rss"]) / 2**30 if r["peak_rss"] else None
            rows.append({"source": "trace", "stage": r["name"], "wall_s": wall, "cpu_s": cpu,
                         "peak_rss_gb": rss})  # fmt: skip
    for f in sorted(a.data.glob("**/*times*.json")) + sorted(a.data.glob("anchored*.json")):
        for k, v in json.loads(f.read_text()).items():
            if isinstance(v, float) and (k.endswith("_s") or "times" in f.name):
                rows.append({"source": str(f.relative_to(a.data)), "stage": k, "wall_s": v,
                             "cpu_s": None, "peak_rss_gb": None})  # fmt: skip
    df = pl.DataFrame(rows, schema={"source": pl.String, "stage": pl.String, "wall_s": pl.Float64,
                                    "cpu_s": pl.Float64, "peak_rss_gb": pl.Float64})  # fmt: skip
    df = df.with_columns(
        core_ms_per_pair=pl.col("cpu_s") * 1000 / max(pairs, 1),
        core_s_per_mbp=pl.col("cpu_s") / max(bases / 1e6, 1e-9),
        wall_ms_per_pair=pl.col("wall_s") * 1000 / max(pairs, 1),
    )
    sizes = []
    for rel in ("whole/index", "whole/store", "reference.dmnd", "hmms.hmm"):
        p = a.data / rel
        if p.exists():
            n = (
                sum(f.stat().st_size for f in p.rglob("*") if f.is_file())
                if p.is_dir()
                else p.stat().st_size
            )
            sizes.append({"source": "disk", "stage": rel, "bytes": n})
    df = pl.concat([df, pl.DataFrame(sizes)], how="diagonal_relaxed") if sizes else df
    df.write_csv(a.data / "bench_cost.tsv", separator="\t", float_precision=4)
    with pl.Config(tbl_rows=-1, tbl_cols=-1, float_precision=3, tbl_width_chars=250):
        print(f"{pairs} read pairs, {bases / 1e6:.1f} Mbp")
        print(df)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("step", choices=["score", "cost"])
    p.add_argument("--data", type=Path, required=True)
    p.add_argument("--hits", choices=anchor.HIT_SOURCES, default="hmmsearch")
    p.add_argument("--ref", choices=anchor.REF_ARMS, default="weak")
    p.add_argument("--min-weight", type=float, default=0.1)
    p.add_argument("--no-unanchored", dest="with_unanchored", action="store_false")
    p.add_argument("--trace", type=Path, help="cost: Nextflow's trace.tsv")
    a = p.parse_args()
    {"score": score, "cost": cost}[a.step](a)


if __name__ == "__main__":
    main()
