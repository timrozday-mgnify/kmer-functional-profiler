"""kfp-genomes: genome profile and functional update, a companion of kmer-functional-profiler
(plan: Genome profile and functional update, phase 12)."""

import json
from pathlib import Path
from typing import Annotated

import typer

from kfp_genomes.panel import (
    MAX_PER_SPECIES,
    NEIGHBOURS,
    panel_from_catalogue,
    panel_from_genomes,
)

app = typer.Typer(no_args_is_help=True)


@app.callback()
def main() -> None:
    """Place the sample's species and strains among reference genomes from the functional
    profile (`place`, on a `panel`), then update the functional profile from them (`update`).
    The profile is never modified."""


@app.command()
def panel(
    index_dir: Annotated[Path, typer.Argument(help="The index profiles are made with")],
    out_dir: Annotated[Path, typer.Argument(help="Panel directory to write")],
    catalogue: Annotated[
        Path | None,
        typer.Option(
            help="MGnify genome catalogue directory (genomes-all_metadata.tsv and "
            "species_catalogue/ as on the FTP site): pangenome families annotated"
        ),
    ] = None,
    genomes: Annotated[
        Path | None,
        typer.Option(help="Output of `annotate-genomes` with taxonomy (species from s__), instead"),
    ] = None,
    exclude: Annotated[
        Path | None,
        typer.Option(
            help="Genome names (one per line) left out of the panel and every count; their "
            "carriage goes to held_out.parquet (benchmark truth)"
        ),
    ] = None,
    species: Annotated[
        Path | None,
        typer.Option(help="With --catalogue: species representatives (one per line) to keep"),
    ] = None,
    max_per_species: Annotated[
        int, typer.Option(help="Panel genomes per species, by farthest-point sampling")
    ] = MAX_PER_SPECIES,
    neighbours: Annotated[
        int, typer.Option(help="Nearest neighbours per panel genome (the placement graph)")
    ] = NEIGHBOURS,
    alpha: Annotated[
        float | None,
        typer.Option(help="Fix the shrinkage α at every rank (ablation; ~0: no shrinkage)"),
    ] = None,
    completeness: Annotated[
        bool, typer.Option(help="Weight genomes by completeness (off: an ablation)")
    ] = True,
) -> None:
    """Build a genome panel: per species, reference genomes' carriage, copies and k-mer
    survival per unit, their nearest neighbours, and the rank-shrunk carriage frequency."""
    if (catalogue is None) == (genomes is None):
        raise typer.BadParameter("give one of --catalogue and --genomes")
    options = {
        "exclude": frozenset(exclude.read_text().split()) if exclude else frozenset(),
        "max_per_species": max_per_species,
        "neighbours": neighbours,
        "alpha": alpha,
        "completeness": completeness,
    }
    if catalogue is not None:
        keep = set(species.read_text().split()) if species is not None else None
        meta = panel_from_catalogue(index_dir, catalogue, out_dir, keep, **options)
    else:
        assert genomes is not None
        meta = panel_from_genomes(index_dir, genomes, out_dir, **options)
    typer.echo(json.dumps(meta, indent=2))
