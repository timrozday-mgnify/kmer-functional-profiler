"""Prototype CLI; mirrors the planned Rust one."""

import json
import sys
from pathlib import Path
from typing import Annotated

import typer

from kmer_functional_profiler import __version__
from kmer_functional_profiler.index import IndexParams, build_index

app = typer.Typer(no_args_is_help=True)
DEFAULTS = IndexParams()


@app.command()
def version() -> None:
    """Print the version."""
    typer.echo(__version__)


@app.command()
def index(
    members: Annotated[Path, typer.Argument(help="Parquet: protein_id, cluster_rep, ...")],
    out_dir: Path,
    pfam: Annotated[Path | None, typer.Option(help="Parquet: protein_id, pfam_accession")] = None,
    k: int = DEFAULTS.k,
    alphabet: str = DEFAULTS.alphabet,
    t_base: float = DEFAULTS.t_base,
    n_min: int = DEFAULTS.n_min,
    max_groups: int = DEFAULTS.max_groups,
    tier1_per_unit: int = DEFAULTS.tier1_per_unit,
    fp_bits: int = DEFAULTS.fp_bits,
    batch_residues: int = DEFAULTS.batch_residues,
) -> None:
    """Build an index from a members table; print its stats as JSON."""
    params = IndexParams(
        k, alphabet, t_base, n_min, max_groups, tier1_per_unit, fp_bits, batch_residues
    )
    typer.echo(json.dumps(build_index(members, out_dir, params, pfam), indent=2))


@app.command()
def query() -> None:
    """Profile reads against an index (phase 3)."""
    sys.exit("not implemented yet")
