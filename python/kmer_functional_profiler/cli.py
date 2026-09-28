"""Prototype CLI; mirrors the planned Rust one."""

import sys

import typer

from kmer_functional_profiler import __version__

app = typer.Typer(no_args_is_help=True)


@app.command()
def version() -> None:
    """Print the version."""
    typer.echo(__version__)


@app.command()
def index() -> None:
    """Build an index (phase 2)."""
    sys.exit("not implemented yet")


@app.command()
def query() -> None:
    """Profile reads against an index (phase 3)."""
    sys.exit("not implemented yet")
