"""Prototype CLI; mirrors the planned Rust one."""

import json
from pathlib import Path
from typing import Annotated

import typer

from kmer_functional_profiler import __version__
from kmer_functional_profiler.compat import import_signatures
from kmer_functional_profiler.index import Index, IndexParams, build_index
from kmer_functional_profiler.query import profile

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
def import_sourmash(
    signatures: Annotated[Path, typer.Argument(help="sourmash signatures, e.g. .sig.zip")],
    out_dir: Path,
    ksize: int = 11,
) -> None:
    """Build an index from protein FracMinHash signatures (queried with sourmash hashing)."""
    typer.echo(json.dumps(import_signatures(signatures, out_dir, ksize), indent=2))


@app.command()
def query(
    index_dir: Path,
    r1: Path,
    r2: Annotated[Path | None, typer.Argument()] = None,
    out: Annotated[Path, typer.Option(help="TSV of per-unit counts")] = Path("profile.tsv"),
    frames: str = "stopfree",
    genetic_code: int = 11,
) -> None:
    """Profile reads (FASTA/FASTQ, optionally paired) against an index."""
    result = profile(Index.load(index_dir), r1, r2, genetic_code=genetic_code, frames=frames)
    result.write_csv(out, separator="\t")
    typer.echo(f"{result.height} units hit -> {out}", err=True)
