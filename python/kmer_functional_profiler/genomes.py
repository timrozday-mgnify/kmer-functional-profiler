"""Genome mode (phase 11): which reference genomes a sample's unit hits come from.

``annotate_genomes`` runs each genome's proteins through the query's first pass (tier 2,
each unit's own ``t_g``) with the genome as the read id. A genome's content *c_{G,u}* is its
raw hits on unit *u* (every holder of a k-mer counts it, nothing split), so a sample of
genomes at depths λ_G has expected raw unit hits Σ_G λ_G c_{G,u}. It writes, to ``out``:

- ``unit_offsets.npy``, ``genome.npy``, ``hits.npy``, ``kmers.npy``: content as a unit-major
  CSR (rows of unit *u* are ``offsets[u]:offsets[u+1]``, sorted by genome), memory-mapped
  so a fit reads the hit units' rows only;
- ``genome_offsets.npy``, ``genome_hits.npy``: the same hits genome-major, which the
  zero-inflated fit needs over every unit of a candidate (hit or not);
- ``genomes.tsv``: ``genome`` (row number), ``name``, ``path``, ``taxonomy`` (if given),
  ``proteins``, ``units``, ``hits`` (Σ_u c_{G,u});
- ``unit_pfam.parquet``: the index's Pfam labels, if it has them, so the function × taxon
  table needs no index;
- ``meta.json``: the index's ``meta.json`` SHA-256, its unit count, k and alphabet.
"""

import gzip
import json
from collections.abc import Callable, Iterator
from pathlib import Path
from typing import Final

import numpy as np
import polars as pl

from kmer_functional_profiler import _core
from kmer_functional_profiler.index import Index
from kmer_functional_profiler.mask import _sha256
from kmer_functional_profiler.query import _ranges, em, gather, unit_hits

PROTEINS_PER_BATCH: Final = 20_000
# ponytail: both defaults are guesses, to be tuned on the fmh benchmark's genome arm.
MIN_CONTAINMENT: Final = 0.1  # screen: share of a genome's content on hit units
MIN_UNITS: Final = 10  # detection: hit units gather leaves to a genome


def fasta(path: Path) -> Iterator[bytes]:
    """Sequences of a FASTA file (plain or gzip), one per record."""
    opener = gzip.open if path.suffix == ".gz" else open
    seq = bytearray()
    with opener(path, "rb") as f:
        for line in f:
            if line.startswith(b">"):
                if seq:
                    yield bytes(seq)
                seq.clear()
            else:
                seq += line.strip()
    if seq:
        yield bytes(seq)


def _batches(path: Path) -> Iterator[list[bytes]]:
    batch: list[bytes] = []
    for seq in fasta(path):
        batch.append(seq)
        if len(batch) == PROTEINS_PER_BATCH:
            yield batch
            batch = []
    if batch:
        yield batch


def genome_content(index: Index, proteins: Path) -> tuple[pl.DataFrame, int]:
    """Per unit, the raw tier-2 ``hits`` of one genome's proteins and the distinct
    ``kmers`` hit; and the number of proteins."""
    params = index.meta["params"]
    max_hash = int(index.tier2.max_hash)
    parts, n = [], 0
    for batch in _batches(proteins):
        n += len(batch)
        hashes = _core.hash_proteins(
            batch, params["k"], alphabet=params["alphabet"], max_hash=max_hash
        )["hash"]
        hits = unit_hits(
            index.tier2, index.units["max_hash_g"], hashes, np.zeros(len(hashes), np.uint64)
        )
        parts.append(hits.select("unit", "hash"))
    hits = pl.concat([pl.DataFrame(schema={"unit": pl.UInt32, "hash": pl.UInt64}), *parts])
    content = (
        hits.with_columns(pl.col("unit").cast(pl.UInt32))
        .group_by("unit")
        .agg(hits=pl.len().cast(pl.UInt32), kmers=pl.col("hash").n_unique().cast(pl.UInt32))
        .sort("unit")
    )
    return content, n


def annotate_genomes(
    index_dir: str | Path, genomes_tsv: str | Path, out: str | Path
) -> dict[str, object]:
    """Annotate the genomes of ``genomes_tsv`` (columns ``genome`` (name), ``path`` to
    protein FASTA, optional ``taxonomy``) with the index in ``index_dir``; write the genome
    index to ``out`` (see the module docstring) and return its counts."""
    index_dir, out = Path(index_dir), Path(out)
    index = Index.load(index_dir)
    if index.meta.get("hash") == "sourmash":
        raise ValueError("genome annotation needs the index's own hash (not a sourmash import)")
    table = pl.read_csv(genomes_tsv, separator="\t")
    missing = {"genome", "path"} - set(table.columns)
    if missing:
        raise ValueError(f"{genomes_tsv} lacks columns {sorted(missing)}")
    if table["genome"].n_unique() != table.height:
        raise ValueError("genome names must be unique")
    base = Path(genomes_tsv).parent
    contents, proteins = [], []
    # ponytail: all content rows are held and sorted in memory (~16 B per row, ~8 GB for
    # 10^5 genomes); write per-genome shards and merge them if a node runs short.
    for g, path in enumerate(table["path"]):
        content, n = genome_content(index, base / path)
        contents.append(content.with_columns(genome=pl.lit(g, pl.UInt32)))
        proteins.append(n)
    rows = pl.concat(contents)  # genome-major
    out.mkdir(parents=True, exist_ok=True)
    n_units, n_genomes = index.units.height, table.height
    per_genome = np.bincount(rows["genome"].to_numpy(), minlength=n_genomes)
    np.save(out / "genome_offsets.npy", np.r_[0, np.cumsum(per_genome)].astype(np.uint64))
    np.save(out / "genome_hits.npy", rows["hits"].to_numpy())
    by_unit = rows.sort("unit", "genome")
    per_unit = np.bincount(by_unit["unit"].to_numpy(), minlength=n_units)
    np.save(out / "unit_offsets.npy", np.r_[0, np.cumsum(per_unit)].astype(np.uint64))
    for column in ("genome", "hits", "kmers"):
        np.save(out / f"{column}.npy", by_unit[column].to_numpy())
    sums = rows.group_by("genome").agg(units=pl.len(), hits=pl.col("hits").sum())
    genomes = (
        table.rename({"genome": "name"})
        .with_columns(
            genome=pl.int_range(n_genomes, dtype=pl.UInt32),
            proteins=pl.Series(proteins, dtype=pl.UInt32),
        )
        .join(sums, on="genome", how="left")
        .with_columns(pl.col("units", "hits").fill_null(0).cast(pl.UInt64))
    )
    genomes.select(
        "genome", "name", "path", *(["taxonomy"] if "taxonomy" in table.columns else []),
        "proteins", "units", "hits",
    ).write_csv(out / "genomes.tsv", separator="\t")  # fmt: skip
    if (index_dir / "unit_pfam.parquet").exists():
        pl.read_parquet(index_dir / "unit_pfam.parquet").write_parquet(out / "unit_pfam.parquet")
    meta = {
        "index_meta_sha256": _sha256(index_dir / "meta.json"),
        "units": n_units,
        "genomes": n_genomes,
        "k": index.meta["params"]["k"],
        "alphabet": index.meta["params"]["alphabet"],
        "rows": rows.height,
    }
    (out / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    return meta


class GenomeIndex:
    """A genome index written by :func:`annotate_genomes`, its arrays memory-mapped."""

    def __init__(self, path: str | Path) -> None:
        path = Path(path)
        self.path = path
        self.meta = json.loads((path / "meta.json").read_text())
        self.genomes = pl.read_csv(path / "genomes.tsv", separator="\t")

        def load(name: str) -> np.ndarray:
            array: np.ndarray = np.load(path / f"{name}.npy", mmap_mode="r")
            return array

        self.unit_offsets, self.genome, self.hits, self.kmers = map(
            load, ("unit_offsets", "genome", "hits", "kmers")
        )
        self.genome_offsets, self.genome_hits = load("genome_offsets"), load("genome_hits")
        pfam = path / "unit_pfam.parquet"
        self.unit_pfam = pl.read_parquet(pfam) if pfam.exists() else None

    def check_index(self, index_dir: str | Path) -> None:
        """Raise if ``index_dir`` is not the index the genomes were annotated with."""
        if _sha256(Path(index_dir) / "meta.json") != self.meta["index_meta_sha256"]:
            raise ValueError(f"{self.path} was annotated with another index than {index_dir}")

    def content(self, units: np.ndarray) -> pl.DataFrame:
        """``unit``, ``genome``, ``hits``, ``kmers`` rows of ``units`` (index unit ids)."""
        units = np.asarray(units, dtype=np.int64)
        if len(units) and units.max() >= self.meta["units"]:
            raise ValueError("unit ids beyond the annotated index: a profile of another index?")
        starts = self.unit_offsets[units].astype(np.int64)
        lengths = self.unit_offsets[units + 1].astype(np.int64) - starts
        rows = _ranges(starts, lengths)
        return pl.DataFrame(
            {
                "unit": np.repeat(units, lengths).astype(np.uint32),
                "genome": self.genome[rows],
                "hits": self.hits[rows],
                "kmers": self.kmers[rows],
            }
        )

    def genome_values(self, genomes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Each of ``genomes``' content hits over all its units: the hits and, per hit, its
        position in ``genomes``."""
        genomes = np.asarray(genomes, dtype=np.int64)
        starts = self.genome_offsets[genomes].astype(np.int64)
        lengths = self.genome_offsets[genomes + 1].astype(np.int64) - starts
        return self.genome_hits[_ranges(starts, lengths)], np.repeat(
            np.arange(len(genomes)), lengths
        )


def _items_hit(
    gi: GenomeIndex, genomes: np.ndarray
) -> Callable[[np.ndarray, np.ndarray], np.ndarray]:
    """For :func:`em`'s zero-inflated fit: per genome (of ``genomes``, sorted), the units it
    would hit if all present at depth λ, Σ_u 1 - exp(-λ c_{G,u}) over all its units, from
    each genome's distinct content values and their counts."""
    values, at = gi.genome_values(genomes)
    pairs, count = np.unique(np.c_[at, values], axis=0, return_counts=True)
    owner, value = pairs[:, 0], pairs[:, 1].astype(np.float64)
    start = np.searchsorted(owner, np.arange(len(genomes)))
    length = np.diff(np.r_[start, len(owner)])

    def items_hit(ids: np.ndarray, la: np.ndarray) -> np.ndarray:
        pos = np.searchsorted(genomes, ids)
        rows = _ranges(start[pos], length[pos])
        rep = np.repeat(np.arange(len(ids)), length[pos])
        seen = -np.expm1(-la[rep] * value[rows])
        return np.bincount(rep, weights=count[rows] * seen, minlength=len(ids))

    return items_hit


def fit_genomes(
    prof: pl.DataFrame,
    gi: GenomeIndex,
    *,
    min_containment: float = MIN_CONTAINMENT,
    min_units: int = MIN_UNITS,
    report: dict[str, int] | None = None,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Genome depths explaining a profile's raw unit ``hits``.

    Screen: genomes whose content on hit units is at least ``min_containment`` of their
    total (Σ_{u hit} c_{G,u} / Σ_u c_{G,u}). G0: gather over (genome, hit unit) pairs, each
    unit to the genome with the most unassigned hit units; genomes left with fewer than
    ``min_units`` are dropped. G1: the weighted zero-inflated EM (:func:`em`, genome as
    holder, unit as item, weight c_{G,u}) on the rest, so h_u ~ Poisson(Σ_G λ_G c_{G,u})
    over the fraction π_G of each genome's content present.

    Returns per genome fitted: ``genome``, ``units_hit``, ``containment``, ``units_unique``,
    ``gather_rank``, ``depth`` (λ_G; 0 if explained away), ``present`` (π_G); and per hit
    unit carried by a fitted genome: ``unit``, ``genome``, ``c`` (the content), for the
    attribution.
    """
    hit = prof.filter(pl.col("hits") > 0).select(pl.col("unit").cast(pl.UInt32), "hits")
    content = gi.content(np.sort(hit["unit"].to_numpy()))
    total = gi.genomes["hits"].to_numpy().astype(np.float64)
    screened = (
        content.group_by("genome")
        .agg(units_hit=pl.len().cast(pl.UInt32), on_hit=pl.col("hits").sum())
        .with_columns(containment=pl.col("on_hit") / pl.lit(total).gather(pl.col("genome")))
        .filter(pl.col("containment") >= min_containment)
        .drop("on_hit")
    )
    content = content.join(screened.select("genome"), on="genome", how="semi")
    pairs = content.select(unit=pl.col("genome"), hash=pl.col("unit"))
    kept = (
        gather(pairs, np.ones(gi.meta["genomes"]))
        .rename({"unit": "genome", "kmers_unique": "units_unique"})
        .filter(pl.col("units_unique") >= min_units)
    )
    content = content.join(kept.select("genome"), on="genome", how="semi")
    if report is not None:
        report |= {"screened_genomes": screened.height, "gathered_genomes": kept.height}
    if kept.height == 0:
        fit = pl.DataFrame(
            schema={"unit": pl.UInt32, "coverage": pl.Float64, "present": pl.Float64}
        )
    else:
        fit = em(
            content.join(hit, on="unit").select(
                unit=pl.col("genome"),
                hash=pl.col("unit"),
                hits=pl.col("hits_right"),
                weight=pl.col("hits"),
            ),  # fmt: skip
            total,
            zero_inflated=True,
            report=report,
            items_hit=_items_hit(gi, np.sort(kept["genome"].to_numpy())),
        )
    genomes = (
        kept.join(screened, on="genome")
        .join(fit.rename({"unit": "genome", "coverage": "depth"}), on="genome")
        .sort("genome")
    )
    return genomes, content.select("unit", "genome", c=pl.col("hits"))


def explained(prof: pl.DataFrame, genomes: pl.DataFrame, carried: pl.DataFrame) -> pl.DataFrame:
    """Per hit unit: ``predicted`` raw hits Σ_G λ_G c_{G,u} (what its genomes would put
    there if present) and ``explained`` = min(``hits``, ``predicted``)."""
    predicted = (
        carried.join(genomes.select("genome", "depth"), on="genome")
        .group_by("unit")
        .agg(predicted=(pl.col("depth") * pl.col("c")).sum())
    )
    return (
        prof.filter(pl.col("hits") > 0)
        .select(pl.col("unit").cast(pl.UInt32), "hits")
        .join(predicted, on="unit", how="left")
        .with_columns(pl.col("predicted").fill_null(0.0))
        .with_columns(explained=pl.min_horizontal("hits", "predicted"))
    )


def genome_profile(
    prof: pl.DataFrame,
    gi: GenomeIndex,
    *,
    min_containment: float = MIN_CONTAINMENT,
    min_units: int = MIN_UNITS,
) -> tuple[pl.DataFrame, dict[str, float | int]]:
    """:func:`fit_genomes` as a report: per detected genome (``depth`` > 0) its name,
    taxonomy and ``relative_abundance`` (share of Σ depth); and the sample's raw ``hits``,
    ``explained_fraction`` (Σ explained / Σ hits: the rest is organisms not in the set,
    or units their genomes carry beyond the fit) and ``genome_equivalents`` (Σ depth)."""
    report: dict[str, int] = {}
    genomes, carried = fit_genomes(
        prof, gi, min_containment=min_containment, min_units=min_units, report=report
    )
    units = explained(prof, genomes, carried)
    detected = genomes.filter(pl.col("depth") > 0)
    info = gi.genomes.select(
        "genome", "name", *(["taxonomy"] if "taxonomy" in gi.genomes.columns else [])
    )
    table = (
        detected.with_columns(relative_abundance=pl.col("depth") / pl.col("depth").sum())
        .join(info.cast({"genome": pl.UInt32}), on="genome")
        .sort("depth", descending=True)
    )
    hits = float(units["hits"].sum())
    summary = report | {
        "hits": hits,
        "explained_fraction": float(units["explained"].sum()) / hits if hits else 0.0,
        "genome_equivalents": float(detected["depth"].sum()),
        "genomes_detected": detected.height,
    }
    return table, summary
