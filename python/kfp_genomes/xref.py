"""Read cross-reference: function x genome from two independent read classifications
(plan: Genome profile and functional update, Read cross-reference pipeline; phase 12).

Inputs, per sample (``workflows/read-xref``):

- *read -> function:* kfp-pileup's ``assignments.parquet`` and ``refs.tsv`` (DIAMOND or
  bwa-mem2 on the recruited reference). A mate whose best alignment is a target is spread
  over the target units by its posteriors, renormalised over its target alignments;
- *read -> genome:* bowtie2 ``-k`` hits as TSV (``HITS_COLUMNS``: one row per aligned mate
  line) against the genomes sylph picked, with ``contigs.tsv`` (``contig``, ``genome``,
  ``length``). Per pair and genome, the score is the sum over mates of the mate's best
  ``AS`` there; genomes where more mates align rank first, then by score; the pair is split
  evenly over the genomes tied at the best. Alignments below ``min_identity`` (1 - NM /
  aligned length) or ``min_mapq`` are dropped first.

Join on the pair (read name without ``/1``, ``/2``): a pair counts once, each side's
weights the mean over its assigned mates, so mates that disagree are split. Its function
weight times its genome weight goes to (unit, genome); pairs with a function but no genome
go to ``unclassified``.

Outputs:

- ``function_genome.tsv``: ``unit``, ``rank`` (``genome``, and ``species``, ``genus`` from
  the GTDB-style ``taxonomy`` of ``--genomes`` when given), ``taxon``, ``reads`` (pairs),
  ``share`` (of the unit's pairs), ``coverage`` (pairs x mean aligned bases per pair /
  unit length in bases) and ``copies`` (coverage / the taxon's depth; null for
  ``unclassified``). Rows of a unit and rank sum to the unit's pairs;
- ``genomes.tsv``: per genome (and taxon), ``reads`` (every pair bowtie2 gave it, with a
  function or not), ``function_reads``, ``length`` and ``depth`` (aligned bases /
  length).

Crude by design (``ponytail:``): a pair with one mate in the gene counts as a whole pair;
ties are exact ``AS`` ties, so near-identical genomes split only when they score the same.
Not done here: Pfam and KO roll-ups (the scoring code's job, from the index's unit tables).
"""

from pathlib import Path
from typing import Annotated, Final

import polars as pl
import typer

HITS_COLUMNS: Final = ["read_name", "flag", "contig", "mapq", "score", "mismatches", "aligned"]
MIN_IDENTITY: Final = 0.95
MIN_MAPQ: Final = 0
UNCLASSIFIED: Final = "unclassified"
PAIR: Final = pl.col("read_name").str.replace(r"/[12]$", "")


def function_weights(assignments: pl.DataFrame) -> pl.DataFrame:
    """(pair, unit, w): each pair's function weights, summing to 1 per assigned pair."""
    mates = (
        assignments.filter(pl.col("kept"), pl.col("role") == "target")
        .group_by(pair=PAIR, mate="mate", unit="unit")
        .agg(pl.col("posterior").sum())
        .with_columns(w=pl.col("posterior") / pl.col("posterior").sum().over("pair", "mate"))
    )
    return (
        mates.group_by("pair", "unit")
        .agg(pl.col("w").sum())
        .with_columns(pl.col("w") / pl.col("w").sum().over("pair"))
    )


def genome_weights(
    hits: pl.DataFrame, contigs: pl.DataFrame, min_identity: float, min_mapq: int
) -> pl.DataFrame:
    """(pair, genome, w, bases): each pair's genome weights (summing to 1) and the pair's
    aligned bases on that genome."""
    best = (
        hits.filter(
            pl.col("flag") & 4 == 0,
            pl.col("mapq") >= min_mapq,
            1 - pl.col("mismatches") / pl.col("aligned") >= min_identity,
        )
        .join(contigs.select("contig", "genome"), on="contig")
        .with_columns(pair=PAIR, mate=(pl.col("flag") & 128 > 0))
        # a mate's best alignment on each genome
        .group_by("pair", "mate", "genome")
        .agg(pl.col("score", "aligned").sort_by("score").last())
        .group_by("pair", "genome")
        .agg(pl.len().alias("mates"), pl.col("score", "aligned").sum())
    )
    top = best.filter(pl.col("mates") == pl.col("mates").max().over("pair"))
    top = top.filter(pl.col("score") == pl.col("score").max().over("pair"))
    return top.select("pair", "genome", w=1 / pl.len().over("pair"), bases=pl.col("aligned") * 1.0)


def unit_lengths(refs: pl.DataFrame, nucleotide: bool) -> pl.DataFrame:
    """(unit, length in bases) of each target unit's ``rep`` sequence (kfp's ``fetch``: the
    representative, else the longest member)."""
    return refs.filter(pl.col("role") == "target", pl.col("rep")).select(
        "unit", unit_length=pl.col("length") * (1 if nucleotide else 3)
    )


def cross_reference(
    assignments: pl.DataFrame,
    refs: pl.DataFrame,
    nucleotide: bool,
    hits: pl.DataFrame,
    contigs: pl.DataFrame,
    taxonomy: pl.DataFrame | None = None,
    min_identity: float = MIN_IDENTITY,
    min_mapq: int = MIN_MAPQ,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """(function_genome, genomes) tables; see the module docstring."""
    fw = function_weights(assignments)
    gw = genome_weights(hits, contigs, min_identity, min_mapq)
    pairs = gw.select("pair").unique().height
    bases_per_pair = gw.select((pl.col("w") * pl.col("bases")).sum()).item() / max(pairs, 1)

    joined = fw.join(gw, on="pair", how="left").select(
        "pair",
        "unit",
        genome=pl.col("genome").fill_null(UNCLASSIFIED),
        w=pl.col("w") * pl.col("w_right").fill_null(1.0),
    )
    genomes = (
        contigs.group_by("genome")
        .agg(pl.col("length").sum())
        .join(
            gw.group_by("genome").agg(
                reads=pl.col("w").sum(), bases=(pl.col("w") * pl.col("bases")).sum()
            ),
            on="genome",
            how="left",
        )
        .join(
            joined.group_by("genome").agg(function_reads=pl.col("w").sum()),
            on="genome",
            how="left",
        )
        .select(
            "genome",
            reads=pl.col("reads").fill_null(0.0),
            function_reads=pl.col("function_reads").fill_null(0.0),
            length="length",
            bases=pl.col("bases").fill_null(0.0),
        )
    )
    ranks = [("genome", pl.col("genome"))]
    if taxonomy is not None:
        genomes = genomes.join(taxonomy.select("genome", "taxonomy"), on="genome", how="left")
        joined = joined.join(taxonomy.select("genome", "taxonomy"), on="genome", how="left")
        ranks += [
            (rank, pl.col("taxonomy").str.extract(rf"{rank[0]}__([^;]+)", 1))
            for rank in ("species", "genus")
        ]

    tables, depth_tables = [], []
    for rank, taxon in ranks:
        # unclassified pairs (and genomes without that rank) stay unclassified at every rank
        def at_rank(df: pl.DataFrame, taxon: pl.Expr = taxon) -> pl.DataFrame:
            return df.with_columns(
                taxon=pl.when(pl.col("genome") == UNCLASSIFIED)
                .then(pl.lit(UNCLASSIFIED))
                .otherwise(taxon.fill_null(UNCLASSIFIED))
            )

        depth = (
            at_rank(genomes)
            .group_by("taxon")
            .agg(pl.col("reads", "function_reads", "length", "bases").sum())
            .with_columns(rank=pl.lit(rank), depth=pl.col("bases") / pl.col("length"))
        )
        depth_tables.append(depth)
        tables.append(
            at_rank(joined)
            .group_by("unit", "taxon")
            .agg(reads=pl.col("w").sum())
            .with_columns(rank=pl.lit(rank))
            .join(depth.select("taxon", "depth"), on="taxon", how="left")
        )
    lengths = unit_lengths(refs, nucleotide)
    function_genome = (
        pl.concat(tables)
        .join(lengths, on="unit", how="left")
        .with_columns(
            share=pl.col("reads") / pl.col("reads").sum().over("unit", "rank"),
            coverage=pl.col("reads") * bases_per_pair / pl.col("unit_length"),
        )
        .select(
            "unit",
            "rank",
            "taxon",
            "reads",
            "share",
            "coverage",
            copies=pl.when(pl.col("taxon") != UNCLASSIFIED).then(
                pl.col("coverage") / pl.col("depth")
            ),
        )
        .sort("unit", "rank", "reads", descending=[False, False, True])
    )
    genome_table = (
        pl.concat(depth_tables)
        .select("rank", "taxon", "reads", "function_reads", "length", "depth")
        .sort("rank", "reads", descending=[False, True])
    )
    return function_genome, genome_table


def read_hits(path: Path) -> pl.DataFrame:
    """bowtie2 hits TSV; empty when no genome was picked."""
    schema = dict.fromkeys(HITS_COLUMNS, pl.Int64) | {"read_name": pl.String, "contig": pl.String}
    if path.stat().st_size == 0:
        return pl.DataFrame(schema=schema)
    return pl.read_csv(path, separator="\t", has_header=False, new_columns=HITS_COLUMNS,
                       schema_overrides=schema)  # fmt: skip


app = typer.Typer(no_args_is_help=True, add_completion=False)


@app.command()
def main(
    pileup_dir: Annotated[Path, typer.Argument(help="kfp-pileup call's output")],
    ref_dir: Annotated[Path, typer.Argument(help="kfp-pileup reference's output")],
    hits: Annotated[Path, typer.Argument(help="bowtie2 hits TSV (HITS_COLUMNS, no header)")],
    contigs: Annotated[Path, typer.Argument(help="contigs.tsv: contig, genome, length")],
    out_dir: Path,
    genomes: Annotated[
        Path | None,
        typer.Option(help="genomes.tsv with genome and taxonomy (GTDB-style) for roll-ups"),
    ] = None,
    min_identity: float = MIN_IDENTITY,
    min_mapq: int = MIN_MAPQ,
) -> None:
    """Write function_genome.tsv and genomes.tsv to OUT_DIR."""
    taxonomy = None
    if genomes is not None:
        table = pl.read_csv(genomes, separator="\t", infer_schema=False)
        if "taxonomy" in table.columns:
            taxonomy = table
    function_genome, genome_table = cross_reference(
        pl.read_parquet(pileup_dir / "assignments.parquet"),
        pl.read_csv(ref_dir / "refs.tsv", separator="\t", infer_schema=False).with_columns(
            pl.col("unit").cast(pl.UInt32),
            pl.col("length").cast(pl.Int64),
            pl.col("rep") == "true",
        ),
        (ref_dir / "molecule.txt").read_text().strip() == "nucleotide",
        read_hits(hits),
        pl.read_csv(
            contigs,
            separator="\t",
            schema={"contig": pl.String, "genome": pl.String, "length": pl.Int64},
        ),
        taxonomy,
        min_identity,
        min_mapq,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    function_genome.write_csv(out_dir / "function_genome.tsv", separator="\t")
    genome_table.write_csv(out_dir / "genomes.tsv", separator="\t")
    typer.echo(f"{function_genome.height} function x taxon rows, {genome_table.height} taxa")


if __name__ == "__main__":
    app()
