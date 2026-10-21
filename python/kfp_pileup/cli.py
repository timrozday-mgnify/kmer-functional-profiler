"""kfp-pileup CLI: a per-sample reference from a kfp profile, then pileups from alignments."""

import json
from pathlib import Path
from typing import Annotated

import polars as pl
import typer

from kfp_pileup.check import MIN_READS, aai_check
from kfp_pileup.handoff import MIN_POSTERIOR, handoff
from kfp_pileup.pileup import BATCH_READS, MIN_DEPTH, MIN_FREQ, SCORE_NATS, call
from kfp_pileup.reference import MAX_FDR, MAX_TARGETS, NEIGHBOURS, build_reference
from kfp_pileup.refine import NOVEL
from kfp_pileup.route import BLOCK_LETTERS, MIN_KMERS, ROUTE_K, route
from kmer_functional_profiler.store import MEMBERS

app = typer.Typer(no_args_is_help=True, add_completion=False)


@app.command()
def reference(
    profile: Annotated[Path, typer.Argument(help="kfp profile.tsv")],
    index_dir: Path,
    out_dir: Path,
    max_fdr: Annotated[
        float,
        typer.Option(
            help="Targets: the top units by present_prob whose expected false share is at "
            "most this (present_q <= max_fdr)"
        ),
    ] = MAX_FDR,
    members: Annotated[int, typer.Option(help="Sequences per unit")] = MEMBERS,
    store: Annotated[
        list[Path] | None,
        typer.Option(help="Members tables (default the index's targets.parquet); repeatable"),
    ] = None,
    withhold: Annotated[
        Path | None,
        typer.Option(help="Units (cluster_rep, one per line) left out of the reference"),
    ] = None,
    neighbours: Annotated[
        int, typer.Option(help="Neighbours per target in a component, the most-hit first")
    ] = NEIGHBOURS,
    max_targets: Annotated[
        int,
        typer.Option(help="At most this many targets, the most hit first (0 = no cap)"),
    ] = MAX_TARGETS,
) -> None:
    """Write ref.fa (targets and component neighbours), refs.tsv and
    molecule.txt; print the sequences per role."""
    held = withhold.read_text().split() if withhold else ()
    counts = build_reference(
        profile, index_dir, out_dir, max_fdr, members, store or (), held, neighbours, max_targets
    )
    typer.echo(json.dumps(counts))


@app.command(name="route")
def route_command(
    ref_dir: Annotated[Path, typer.Argument(help="Output of `reference` (amino acids)")],
    index_dir: Annotated[Path, typer.Argument(help="The kfp index (its k and alphabet)")],
    r1: Path,
    r2: Path,
    out_dir: Path,
    block_letters: Annotated[
        int, typer.Option(help="Residues per block of whole components")
    ] = BLOCK_LETTERS,
    k: Annotated[int, typer.Option(help="Routing k-mer length (0: the index's k)")] = ROUTE_K,
    min_kmers: Annotated[
        int, typer.Option(help="Distinct k-mers a pair shares with a block to be sent there")
    ] = MIN_KMERS,
    genetic_code: int = 11,
) -> None:
    """Split the reference into blocks of components and each read pair to the blocks it
    shares at least --min-kmers distinct k-mers with (other pairs dropped): b<i>.fa,
    b<i>_1/2.fq.gz, letters.txt (DIAMOND's --dbsize), route.json; print counts."""
    typer.echo(
        json.dumps(
            route(ref_dir, index_dir, (r1, r2), out_dir, block_letters, k, min_kmers, genetic_code)
        )
    )


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
    novel: Annotated[
        float, typer.Option(help="Least divergence of a novel relative from the reference")
    ] = NOVEL,
    reads: Annotated[
        list[Path] | None,
        typer.Option(help="diamond: the FASTQs (mate order), for the nucleotide pileup"),
    ] = None,
    stats: Annotated[
        Path | None, typer.Option(help="JSON of each stage's time and memory (as query's)")
    ] = None,
) -> None:
    """Read posteriors (assignments.parquet), pileups (pileup.parquet, variants.tsv) and
    per-unit identity (units.tsv) from alignments to the reference; then the joint EM over
    read sources, consensus updates and haplotypes (sources.tsv, haplotypes.tsv); with
    DIAMOND alignments and --reads, the nucleotide pileup and consensus (pileup_nt.parquet,
    consensus.*, recruited_*.fq.gz); print counts."""
    counts = call(
        alignments, ref_dir, out_dir, aligner, profile, min_depth, min_freq, score_nats,
        batch_reads, novel, reads or (), stats,
    )  # fmt: skip
    typer.echo(json.dumps(counts))


@app.command(name="handoff")
def handoff_command(
    pileup_dir: Annotated[Path, typer.Argument(help="Output of `call`")],
    ref_dir: Annotated[Path, typer.Argument(help="Output of `reference`")],
    sam: Annotated[Path, typer.Argument(help="bwa-mem2 SAM (consensus.fna or ref.fa)")],
    out: Annotated[Path, typer.Option(help="SAM to write")] = Path("handoff.sam"),
    min_posterior: Annotated[
        float, typer.Option(help="Least responsibility of the pair's unit")
    ] = MIN_POSTERIOR,
) -> None:
    """Keep each mate's best alignment to its pair's assigned unit (step 8), MAPQ from the
    posterior: the BAM for step 6c's callers; print counts."""
    typer.echo(json.dumps(handoff(pileup_dir, ref_dir, sam, out, min_posterior)))


@app.command(name="aai-check")
def aai_check_command(
    units: Annotated[Path, typer.Argument(help="units.tsv of `call`")],
    out: Annotated[Path, typer.Option(help="TSV to write")] = Path("aai_check.tsv"),
    min_reads: Annotated[
        float, typer.Option(help="Units with fewer reads are left out")
    ] = MIN_READS,
) -> None:
    """Step 6's aai check: the profile's aai against the reads' alignment identity, by
    identity, depth, length and cluster size; prints the table (flagged: step 6's rule)."""
    table = aai_check(units, min_reads)
    table.write_csv(out, separator="\t")
    with pl.Config(tbl_rows=-1, tbl_width_chars=150):
        typer.echo(table)
