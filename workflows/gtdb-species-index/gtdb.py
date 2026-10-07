"""Genomes for a GTDB species index (plan: phase 11, step 11): the release's metadata to a
genome set, with each genome's NCBI DNA URL and the TSV ``annotate-genomes`` reads."""

import argparse
import random
import re
from pathlib import Path

import polars as pl

NCBI: str = "https://ftp.ncbi.nlm.nih.gov/genomes/all"


def ncbi_url(accession: str, assembly_name: str, base: str = NCBI) -> str:
    """The genomic FASTA of an NCBI assembly (``GCA_000005845.2``): its directory is the
    accession and the assembly name with characters other than letters, digits, ``.``,
    ``_`` and ``-`` made ``_``."""
    stem = f"{accession}_{re.sub(r'[^A-Za-z0-9._-]', '_', assembly_name)}"
    d = accession[4:13]
    return f"{base}/{accession[:3]}/{d[:3]}/{d[3:6]}/{d[6:]}/{stem}/{stem}_genomic.fna.gz"


def pick(
    metadata: list[str],
    max_per_species: int,
    min_completeness: float,
    max_contamination: float,
    species: set[str] | None = None,
    max_species: int | None = None,
    seed: int = 0,
    base: str = NCBI,
) -> pl.DataFrame:
    """Genomes of GTDB metadata tables (``bac120``/``ar53`` ``_metadata``) of at least
    ``min_completeness`` and at most ``max_contamination`` (CheckM2, else CheckM), at most
    ``max_per_species`` per species: its representative, then the best by completeness - 5
    contamination. ``species`` (GTDB names, ``s__`` included) or ``max_species`` (sampled)
    keep a subset. Columns ``genome`` (accession without ``RS_``/``GB_``), ``url``,
    ``taxonomy``, ``completeness`` (percent), ``representative``."""
    meta = pl.concat(
        [pl.read_csv(p, separator="\t", infer_schema=False) for p in metadata], how="diagonal"
    )
    quality = {
        m: pl.coalesce([c for c in (f"checkm2_{m}", f"checkm_{m}") if c in meta.columns])
        for m in ("completeness", "contamination")
    }
    meta = meta.select(
        genome=pl.col("accession").str.slice(3),
        assembly="ncbi_assembly_name",
        taxonomy="gtdb_taxonomy",
        species=pl.col("gtdb_taxonomy").str.extract(r"(s__[^;]+)$", 1),
        representative=pl.col("gtdb_representative") == "t",
        completeness=quality["completeness"].cast(pl.Float64),
        contamination=quality["contamination"].cast(pl.Float64),
    ).filter(
        pl.col("completeness") >= min_completeness, pl.col("contamination") <= max_contamination
    )
    if species is not None:
        meta = meta.filter(pl.col("species").is_in(list(species)))
    if max_species is not None:
        names = sorted(meta["species"].unique())
        chosen = random.Random(seed).sample(names, min(max_species, len(names)))
        meta = meta.filter(pl.col("species").is_in(chosen))
    picked = (
        meta.sort(
            "species",
            "representative",
            pl.col("completeness") - 5 * pl.col("contamination"),
            "genome",
            descending=[False, True, True, False],
        )  # fmt: skip
        .group_by("species", maintain_order=True)
        .head(max_per_species)
    )
    urls = [ncbi_url(g, a, base) for g, a in zip(picked["genome"], picked["assembly"], strict=True)]
    return picked.select(
        "genome", url=pl.Series(urls), taxonomy="taxonomy", completeness="completeness",
        representative="representative",
    )  # fmt: skip


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metadata", nargs="+", help="GTDB bac120/ar53 metadata TSVs (or .gz)")
    parser.add_argument("--max-per-species", type=int, default=20)
    parser.add_argument("--min-completeness", type=float, default=50)
    parser.add_argument("--max-contamination", type=float, default=10)
    parser.add_argument("--species", help="GTDB species names (s__...), one per line, to keep")
    parser.add_argument("--max-species", type=int, help="keep this many species, sampled")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--ncbi", default=NCBI, help="NCBI genomes/all mirror")
    args = parser.parse_args()
    keep = set(Path(args.species).read_text().splitlines()) if args.species else None
    genomes = pick(
        args.metadata, args.max_per_species, args.min_completeness, args.max_contamination,
        keep, args.max_species, args.seed, args.ncbi,
    )  # fmt: skip
    genomes.write_csv("genomes.tsv", separator="\t")
    # annotate-genomes' input: proteins at faa/<genome>.faa beside it (GENES writes them)
    genomes.select(
        "genome", path=pl.format("faa/{}.faa", "genome"), taxonomy="taxonomy",
        completeness="completeness",
    ).write_csv("annotate.tsv", separator="\t")  # fmt: skip


if __name__ == "__main__":
    main()
