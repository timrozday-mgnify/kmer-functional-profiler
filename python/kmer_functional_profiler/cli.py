"""Prototype CLI; mirrors the planned Rust one."""

import json
from pathlib import Path
from typing import Annotated

import typer

from kmer_functional_profiler import __version__
from kmer_functional_profiler.compat import import_signatures
from kmer_functional_profiler.index import (
    AAI_CALIBRATION,
    Index,
    IndexParams,
    build_index,
    write_unit_columns,
)
from kmer_functional_profiler.mask import Mask, build_mask
from kmer_functional_profiler.query import Timer, check_aai_calibration, profile

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
def unit_columns(index_dir: Path) -> None:
    """Write the unit table's numeric columns as .npy files, for an index built without them."""
    write_unit_columns(index_dir)


@app.command()
def calibrate_aai(
    index_dir: Path,
    calibration: Annotated[
        Path | None,
        typer.Argument(help="JSON from the benchmark's aai-calibrate; omit to remove the map"),
    ] = None,
) -> None:
    """Attach an aai calibration to an index (``aai_calibration.json``): its queries then
    report ``aai`` as alignment identity and keep the raw estimate as ``aai_raw``."""
    target = index_dir / AAI_CALIBRATION
    if calibration is None:
        target.unlink(missing_ok=True)
        return
    cal = json.loads(calibration.read_text())
    meta = json.loads((index_dir / "meta.json").read_text())
    try:
        check_aai_calibration(cal, meta["params"])
    except ValueError as e:
        raise typer.BadParameter(str(e)) from e
    target.write_text(json.dumps(cal, indent=1) + "\n")


@app.command(name="mask")
def mask_command(
    genome: Annotated[Path, typer.Argument(help="Host genome FASTA (plain or gzip)")],
    index_dir: Path,
    out_dir: Path,
) -> None:
    """Build a mask sidecar: the index's k-mers in the six-frame translated genome."""
    typer.echo(json.dumps(build_mask(genome, index_dir, out_dir), indent=2))


@app.command()
def query(
    index_dir: Path,
    r1: Path,
    r2: Annotated[Path | None, typer.Argument()] = None,
    out: Annotated[Path, typer.Option(help="TSV of per-unit counts")] = Path("profile.tsv"),
    frames: Annotated[
        str,
        typer.Option(
            help="Frames hashed: stopfree; edges (also the terminal segments, >= 20 aa, of "
            "frames with stops; edges:M sets the length); all"
        ),
    ] = "stopfree",
    min_qual: Annotated[
        int, typer.Option(help="Mask bases below this Phred quality as N (0 = off)")
    ] = 0,
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
    in_memory: Annotated[
        bool,
        typer.Option(
            help="Read the index tiers into memory instead of memory-mapping them: no page "
            "faults during lookup, at the cost of their full size in RSS"
        ),
    ] = False,
    low_memory: Annotated[
        bool,
        typer.Option(
            help="Trade time for memory where results are unchanged: with --draws, re-read "
            "the reads for the detected units' k-mers instead of keeping per-read rows of "
            "every hit unit through the query"
        ),
    ] = False,
    aai: Annotated[
        bool,
        typer.Option(
            help="Containment AAI of the detected units (aai, with a closed-form 95% interval "
            "aai_lo/aai_hi; no --draws needed); fits the zero-inflated EM. aai_naive is "
            "reported for every hit unit regardless"
        ),
    ] = False,
    min_aai: Annotated[
        float, typer.Option(help="Drop units with aai_naive below this (0 = keep all)")
    ] = 0.0,
    extra_index: Annotated[
        list[Path] | None,
        typer.Option(
            help="Another index (same k, alphabet, hash) queried jointly; repeatable. Its "
            "units compete with the first index's; output ids are offset, with a source column"
        ),
    ] = None,
    mask: Annotated[
        Path | None,
        typer.Option(help="Mask sidecar (from `mask`) of the first index: host k-mers dropped"),
    ] = None,
    summary: Annotated[
        Path | None,
        typer.Option(
            help="JSON of sample-level explained and unknown fractions: model-based "
            "(explained_fraction) and a census of known k-mers (census_containment); fits "
            "the zero-inflated EM"
        ),
    ] = None,
    all_estimators: Annotated[
        bool,
        typer.Option(
            help="Also fit the benchmark estimators (zero-inflated, empirical-Bayes, p_in-"
            "weighted EM; winner-take-all, uniqueness-first); coverage_em is the default"
        ),
    ] = False,
) -> None:
    """Profile reads (FASTA/FASTQ, optionally paired) against an index."""
    timer = Timer(stats, log=stats is not None)
    sample: dict[str, float | int | None] | None = None if summary is None else {}
    with timer("total"):
        with timer("load"):
            loaded = Index.load(index_dir, mmap=not in_memory)
            extra = [Index.load(d, mmap=not in_memory) for d in extra_index or []]
            masked = None if mask is None else Mask(mask)
        result = profile(
            loaded,
            r1,
            r2,
            genetic_code=genetic_code,
            frames=frames,
            min_qual=min_qual,
            draws=draws,
            kmers_out=kmers,
            timer=timer if stats else None,
            all_estimators=all_estimators,
            low_memory=low_memory,
            with_aai=aai,
            min_aai=min_aai,
            extra=extra,
            mask=masked,
            summary=sample,
        )
        result.write_csv(out, separator="\t")
        if summary is not None:
            summary.write_text(json.dumps(sample, indent=2) + "\n")
    timer.write()
    typer.echo(f"{result.height} units hit -> {out}", err=True)
