"""Prototype CLI; mirrors the planned Rust one."""

import json
from pathlib import Path
from typing import Annotated

import typer

from kmer_functional_profiler import __version__
from kmer_functional_profiler.compat import import_signatures
from kmer_functional_profiler.index import Index, IndexParams, build_index
from kmer_functional_profiler.query import Timer, profile

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
    t_cap: float = DEFAULTS.t_cap,
    oversample: float = DEFAULTS.oversample,
    mask_adapters: bool = DEFAULTS.mask_adapters,
    max_groups: int = DEFAULTS.max_groups,
    fp_bits: int = DEFAULTS.fp_bits,
    t_dense: float = DEFAULTS.t_dense,
    batch_residues: int = DEFAULTS.batch_residues,
    postings: Annotated[bool, typer.Option(help="Also write postings.parquet")] = False,
) -> None:
    """Build an index from a members table; print its stats as JSON."""
    params = IndexParams(
        k=k,
        alphabet=alphabet,
        t_base=t_base,
        n_min=n_min,
        t_cap=t_cap,
        oversample=oversample,
        mask_adapters=mask_adapters,
        max_groups=max_groups,
        fp_bits=fp_bits,
        t_dense=t_dense,
        batch_residues=batch_residues,
    )
    typer.echo(json.dumps(build_index(members, out_dir, params, pfam, postings), indent=2))


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
    draws: Annotated[
        int, typer.Option(help="Posterior draws for 95% intervals and ambiguity groups")
    ] = 0,
    kmers: Annotated[
        Path | None, typer.Option(help="Parquet of tier-2 hits per unit and k-mer (diagnostics)")
    ] = None,
    stats: Annotated[
        Path | None,
        typer.Option(
            help="JSON of time and peak RSS per stage, and hit counts; rewritten after each "
            "stage, and each stage's start and end are logged to stderr with its RSS"
        ),
    ] = None,
) -> None:
    """Profile reads (FASTA/FASTQ, optionally paired) against an index."""
    timer = Timer(stats, log=stats is not None)
    with timer("total"):
        with timer("load"):
            loaded = Index.load(index_dir)
        result = profile(
            loaded,
            r1,
            r2,
            genetic_code=genetic_code,
            frames=frames,
            draws=draws,
            kmers_out=kmers,
            timer=timer if stats else None,
        )
        result.write_csv(out, separator="\t")
    timer.write()
    typer.echo(f"{result.height} units hit -> {out}", err=True)
