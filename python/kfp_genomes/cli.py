"""kfp-genomes: genome profile and functional update, a companion of kmer-functional-profiler
(plan: Genome profile and functional update, phase 12)."""

import json
from pathlib import Path
from typing import Annotated

import polars as pl
import typer

from kfp_genomes.evidence import Histogram
from kfp_genomes.panel import (
    MAX_PER_SPECIES,
    NEIGHBOURS,
    Panel,
    check_profile,
    panel_from_catalogue,
    panel_from_genomes,
)
from kfp_genomes.report import genome_profile

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


@app.command()
def place(
    profile_tsv: Annotated[Path, typer.Argument(help="profile.tsv of `query` (phase 12)")],
    panel_dir: Annotated[Path, typer.Argument(help="Output of `kfp-genomes panel`")],
    out_dir: Annotated[Path, typer.Argument(help="Directory for the outputs")] = Path("place"),
    own_hist: Annotated[
        Path | None,
        typer.Option(
            help="The profile's own-k-mer histogram (`query --own-hist`): strain mixtures' evidence"
        ),  # fmt: skip
    ] = None,
    max_strains: Annotated[
        int, typer.Option(help="Strains per species considered (K ≥ 2 by MCMC; slower)")
    ] = 1,
    index: Annotated[
        Path | None, typer.Option(help="The profile's index, checked against the panel's")
    ] = None,
    seed: int = 0,
    beta: Annotated[
        float | None,
        typer.Option(help="Fix the background prior (ablation) instead of fitting it"),
    ] = None,
) -> None:
    """The genome profile: species present, their strains' depths, and where each strain
    sits among the panel's reference genomes, from the functional profile (unchanged).
    Writes genome_profile.tsv, candidates.tsv, placements.parquet, strain_units.parquet and
    place_summary.json."""
    panel_ = Panel(panel_dir)
    prof = pl.read_csv(profile_tsv, separator="\t")
    try:
        if index is not None:
            panel_.check_index(index)
        check_profile(prof, panel_)
    except ValueError as e:
        raise typer.BadParameter(str(e)) from e
    hist = Histogram.read(own_hist) if own_hist is not None else None
    result = genome_profile(prof, panel_, max_strains=max_strains, hist=hist, seed=seed,
                            beta=beta)  # fmt: skip
    out_dir.mkdir(parents=True, exist_ok=True)
    for key in ("genome_profile", "candidates"):
        table = result[key]
        assert isinstance(table, pl.DataFrame)
        table.write_csv(out_dir / f"{key}.tsv", separator="\t")
    for key in ("placements", "strain_units"):
        table = result[key]
        if isinstance(table, pl.DataFrame):
            table.write_parquet(out_dir / f"{key}.parquet")
    (out_dir / "place_summary.json").write_text(json.dumps(result["summary"], indent=2) + "\n")
    typer.echo(f"{result['summary']['strains']} strains -> {out_dir}", err=True)  # type: ignore[index]
