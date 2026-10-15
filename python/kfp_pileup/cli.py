"""kfp-pileup CLI: a per-sample reference from a kfp profile, then pileups from alignments."""

import json
from pathlib import Path
from typing import Annotated

import polars as pl
import typer

from kfp_pileup.check import MIN_READS, aai_check
from kfp_pileup.pileup import BATCH_READS, MIN_DEPTH, MIN_FREQ, SCORE_NATS, call
from kfp_pileup.reference import MIN_PRESENT, build_reference
from kfp_pileup.refine import NOVEL
from kmer_functional_profiler.store import MEMBERS

app = typer.Typer(no_args_is_help=True, add_completion=False)


@app.command()
def reference(
    profile: Annotated[Path, typer.Argument(help="kfp profile.tsv")],
    index_dir: Path,
    out_dir: Path,
    min_present: Annotated[
        float, typer.Option(help="Targets: units with present_prob at least this")
    ] = MIN_PRESENT,
    members: Annotated[int, typer.Option(help="Sequences per unit")] = MEMBERS,
    store: Annotated[
        list[Path] | None,
        typer.Option(help="Members tables (default the index's targets.parquet); repeatable"),
    ] = None,
    withhold: Annotated[
        Path | None,
        typer.Option(help="Units (cluster_rep, one per line) left out of the reference"),
    ] = None,
) -> None:
    """Write ref.fa (targets, component neighbours, reversed decoys), refs.tsv and
    molecule.txt; print the sequences per role."""
    held = withhold.read_text().split() if withhold else ()
    counts = build_reference(profile, index_dir, out_dir, min_present, members, store or (), held)
    typer.echo(json.dumps(counts))


@app.command(name="call")
def call_command(
    ref_dir: Annotated[Path, typer.Argument(help="Output of `reference`")],
    out_dir: Path,
    alignments: Annotated[
        list[Path],
        typer.Argument(help="DIAMOND tabular output per mate (in mate order), or one SAM"),
    ],
    aligner: Annotated[str, typer.Option(help="diamond or bwa")] = "diamond",
    profile: Annotated[
        Path | None, typer.Option(help="The profile, whose aai/ani columns join units.tsv")
    ] = None,
    min_depth: Annotated[float, typer.Option(help="variants.tsv: least column depth")] = MIN_DEPTH,
    min_freq: Annotated[float, typer.Option(help="variants.tsv: least alt frequency")] = MIN_FREQ,
    score_nats: Annotated[
        float, typer.Option(help="bwa: nats per alignment score point")
    ] = SCORE_NATS,
    batch_reads: int = BATCH_READS,
    reads: Annotated[
        list[Path] | None,
        typer.Option(help="FASTQ of each mate (in mate order), for base qualities; repeatable"),
    ] = None,
    novel: Annotated[
        float, typer.Option(help="Least divergence of a novel relative from the reference")
    ] = NOVEL,
) -> None:
    """Read posteriors (assignments.parquet), pileups (pileup.parquet, variants.tsv) and
    per-unit identity (units.tsv) from alignments to the reference; then the error model,
    the joint EM, consensus updates and haplotypes (error_model.tsv, sources.tsv,
    haplotypes.tsv, calls.tsv); print counts."""
    counts = call(
        alignments, ref_dir, out_dir, aligner, profile, min_depth, min_freq, score_nats,
        batch_reads, reads or (), novel,
    )  # fmt: skip
    typer.echo(json.dumps(counts))


@app.command(name="aai-check")
def aai_check_command(
    units: Annotated[Path, typer.Argument(help="units.tsv of `call`")],
    out: Annotated[Path, typer.Option(help="TSV to write")] = Path("aai_check.tsv"),
    min_reads: Annotated[
        float, typer.Option(help="Units with fewer reads are left out")
    ] = MIN_READS,
) -> None:
    """Step 6's aai check: the profile's aai against the reads' alignment identity, by
    identity, depth, length and cluster size; prints the table (flagged: step 6's rule)."""
    table = aai_check(units, min_reads)
    table.write_csv(out, separator="\t")
    with pl.Config(tbl_rows=-1, tbl_width_chars=150):
        typer.echo(table)
