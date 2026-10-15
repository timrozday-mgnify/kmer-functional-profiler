"""kfp-pileup CLI: a per-sample reference from a kfp profile, then pileups from alignments."""

import json
from pathlib import Path
from typing import Annotated

import typer

from kfp_pileup.pileup import BATCH_READS, MIN_DEPTH, MIN_FREQ, SCORE_NATS, call
from kfp_pileup.reference import MIN_PRESENT, build_reference
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
) -> None:
    """Write ref.fa (targets, component neighbours, reversed decoys), refs.tsv and
    molecule.txt; print the sequences per role."""
    counts = build_reference(profile, index_dir, out_dir, min_present, members, store or ())
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
) -> None:
    """Read posteriors (assignments.parquet), pileups (pileup.parquet, variants.tsv) and
    per-unit identity (units.tsv) from alignments to the reference; print counts."""
    counts = call(
        alignments, ref_dir, out_dir, aligner, profile, min_depth, min_freq, score_nats,
        batch_reads,
    )  # fmt: skip
    typer.echo(json.dumps(counts))
