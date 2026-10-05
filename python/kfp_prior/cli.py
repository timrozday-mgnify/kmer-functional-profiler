"""kfp-prior: genome-informed unit presence, a companion of kmer-functional-profiler."""

import json
from pathlib import Path
from typing import Annotated

import polars as pl
import typer

from kfp_prior.model import Carriage, build_carriage, update

app = typer.Typer(no_args_is_help=True)


@app.command()
def build(
    genome_index: Annotated[Path, typer.Argument(help="Output of `annotate-genomes`")],
    out_dir: Path,
) -> None:
    """Carriage table: units carried per family, genus and species, α fitted per rank."""
    typer.echo(json.dumps(build_carriage(genome_index, out_dir), indent=2))


@app.command(name="update")
def update_command(
    profile: Annotated[Path, typer.Argument(help="profile.tsv of `query`")],
    genomes: Annotated[Path, typer.Argument(help="genomes.tsv of `genomes`")],
    genome_index: Annotated[Path, typer.Argument(help="Output of `annotate-genomes`")],
    carriage: Annotated[Path, typer.Argument(help="Output of `kfp-prior build`")],
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
