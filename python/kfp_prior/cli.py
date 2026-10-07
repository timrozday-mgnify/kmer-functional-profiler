"""kfp-prior: genome-informed unit presence, a companion of kmer-functional-profiler."""

from pathlib import Path
from typing import Annotated

import polars as pl
import typer

from kfp_prior.model import Carriage, update

app = typer.Typer(no_args_is_help=True)


@app.callback()
def main() -> None:
    """Genome-informed unit presence; the carriage table is a species index (the
    profiler's `species-index --genomes`)."""


@app.command(name="update")
def update_command(
    profile: Annotated[Path, typer.Argument(help="profile.tsv of `query`")],
    genomes: Annotated[Path, typer.Argument(help="genomes.tsv of `genomes`")],
    genome_index: Annotated[Path, typer.Argument(help="Output of `annotate-genomes`")],
    carriage: Annotated[
        Path, typer.Argument(help="Output of the profiler's `species-index --genomes`")
    ],
    out: Annotated[Path, typer.Option(help="Per-unit presence TSV")] = Path("presence.tsv"),
    pfam_out: Annotated[Path, typer.Option(help="Per-Pfam presence TSV")] = Path(
        "pfam_presence.tsv"
    ),
) -> None:
    """Unit presence updated with a prior from the detected genomes (profile unchanged)."""
    presence, by_pfam = update(
        pl.read_csv(profile, separator="\t"),
        pl.read_csv(genomes, separator="\t"),
        genome_index,
        Carriage(carriage),
    )
    presence.write_csv(out, separator="\t")
    if by_pfam is not None:
        by_pfam.write_csv(pfam_out, separator="\t")
    typer.echo(f"{int(presence['imputed'].sum())} units imputed -> {out}", err=True)
