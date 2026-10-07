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
- ``genome_best.parquet``: per protein (``genome``, ``protein``: its number in the FASTA) its
  best ``unit`` by ``containment`` (k-mers hit over its k-mers sampled at the unit's rate),
  if ≥ ``BEST_MIN_CONTAINMENT``: which units each genome carries, for ``kfp-prior``;
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
from scipy import sparse
from scipy.sparse import csgraph

from kmer_functional_profiler import _core
from kmer_functional_profiler.index import Index
from kmer_functional_profiler.mask import _sha256
from kmer_functional_profiler.query import _div, _ranges, em, gather, unit_hits

PROTEINS_PER_BATCH: Final = 20_000
# ponytail: both defaults are guesses, to be tuned on the fmh benchmark's genome arm.
MIN_CONTAINMENT: Final = 0.1  # screen: share of a genome's content on hit units
# ponytail: guesses; a protein at ~90% identity keeps ~0.3 of its k = 11 k-mers.
BEST_MIN_CONTAINMENT: Final = 0.2  # a protein's best unit, for carriage (kfp-prior)
BEST_MIN_KMERS: Final = 2
MIN_UNITS: Final = 10  # detection: hit units gather leaves to a genome
AMBIGUOUS: Final = 0.95  # information cosine above which two genomes are one group
GTDB_RANKS: Final = ("domain", "phylum", "class", "order", "family", "genus", "species")
GTDB_PREFIX: Final = dict(zip("dpcofgs", GTDB_RANKS, strict=True))
REPORTED_RANKS: Final = ("family", "genus", "species", "genome")


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


def genome_content(index: Index, proteins: Path) -> tuple[pl.DataFrame, pl.DataFrame, int]:
    """Per unit, the raw tier-2 ``hits`` of one genome's proteins and the distinct
    ``kmers`` hit; per protein its best unit (:data:`BEST_MIN_CONTAINMENT`); and the number
    of proteins."""
    params = index.meta["params"]
    k, max_hash = params["k"], int(index.tier2.max_hash)
    parts, lengths, n = [], [], 0
    for batch in _batches(proteins):
        hashes = _core.hash_proteins(batch, k, alphabet=params["alphabet"], max_hash=max_hash)
        hits = unit_hits(
            index.tier2, index.units["max_hash_g"], hashes["hash"], hashes["seq"] + np.uint64(n)
        )
        parts.append(hits.select("unit", "hash", protein="read"))
        lengths += [max(len(p) - k + 1, 0) for p in batch]
        n += len(batch)
    empty = pl.DataFrame(schema={"unit": pl.UInt32, "hash": pl.UInt64, "protein": pl.UInt64})
    hits = pl.concat([empty, *parts]).with_columns(pl.col("unit").cast(pl.UInt32))
    content = (
        hits.group_by("unit")
        .agg(hits=pl.len().cast(pl.UInt32), kmers=pl.col("hash").n_unique().cast(pl.UInt32))
        .sort("unit")
    )
    # Containment of a protein in a unit: k-mers hit over the protein's k-mers sampled at the
    # unit's rate, so units at different rates compare.
    per_pair = hits.group_by("protein", "unit").agg(kmers=pl.col("hash").n_unique())
    sampled = np.asarray(lengths, dtype=np.float64)[per_pair["protein"].to_numpy()] * np.asarray(
        index.units["t_g"][per_pair["unit"].to_numpy()], dtype=np.float64
    )
    best = (
        per_pair.with_columns(containment=pl.col("kmers") / pl.Series(sampled))
        .filter(pl.col("kmers") >= BEST_MIN_KMERS, pl.col("containment") >= BEST_MIN_CONTAINMENT)
        .sort("protein", "containment", "unit", descending=[False, True, False])
        .unique("protein", keep="first", maintain_order=True)
        .select(pl.col("protein").cast(pl.UInt32), "unit", "containment")
    )
    return content, best, n


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
    contents, bests, proteins = [], [], []
    # ponytail: all content rows are held and sorted in memory (~16 B per row, ~8 GB for
    # 10^5 genomes); write per-genome shards and merge them if a node runs short.
    for g, path in enumerate(table["path"]):
        content, best, n = genome_content(index, base / path)
        contents.append(content.with_columns(genome=pl.lit(g, pl.UInt32)))
        bests.append(best.with_columns(genome=pl.lit(g, pl.UInt32)))
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
    pl.concat(bests).select("genome", "protein", "unit", "containment").write_parquet(
        out / "genome_best.parquet"
    )
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
    report: dict[str, float] | None = None,
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """Genome depths explaining a profile's raw unit ``hits``.

    Screen: genomes whose content on hit units is at least ``min_containment`` of their
    total (Σ_{u hit} c_{G,u} / Σ_u c_{G,u}). G0: gather over (genome, hit unit) pairs, each
    unit to the genome with the most unassigned hit units; genomes left with fewer than
    ``min_units`` are dropped. G1: the weighted zero-inflated EM (:func:`em`, genome as
    holder, unit as item, weight c_{G,u}) on the rest, so h_u ~ Poisson(Σ_G λ_G c_{G,u})
    over the fraction π_G of each genome's content present.

    Returns per genome fitted: ``genome``, ``units_hit``, ``containment``, ``units_unique``,
    ``gather_rank``, ``depth`` (λ_G; 0 if explained away), ``present`` (π_G); per hit unit
    and screened genome carrying it: ``unit``, ``genome``, ``c`` (the content); and the
    screened genomes' ``genome``, ``units_hit``, ``containment``.
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
    carried = content.select("unit", "genome", c=pl.col("hits"))
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
    return genomes, carried, screened


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


def uncertainty(
    detected: pl.DataFrame, carried: pl.DataFrame, screened: pl.DataFrame, units: pl.DataFrame
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Ambiguity groups and 95% ``depth_lo``/``depth_hi`` of the detected genomes, in closed
    form from the G1 fit (no draws).

    With w_{G,u} = π_G c_{G,u} and μ_u = Σ_G λ_G w_{G,u} (floored at h_u, so units the fit
    leaves unexplained still count), the Poisson information on the depths is
    I_{GH} = Σ_u h_u w_{G,u} w_{H,u} / μ_u². Intervals are λ ± 1.96 sd with covariance
    φ I⁺ (pseudo-inverse), φ = max(1, Pearson χ² / (units − genomes)): raw unit hits are
    correlated (a shared k-mer counts in every holder), so the likelihood is composite and
    φ takes up the overdispersion it leaves.

    Two genomes are ambiguous if the sample cannot tell them apart: the cosine of their
    information vectors over the hit units is ≥ ``AMBIGUOUS`` (scale-free, so it needs no
    depth: it also applies to screened genomes gather dropped) and so is the ratio of their
    ``containment`` (the share of content on hit units; content the sample never hit
    tells genomes apart too). A group is a connected set of ambiguous genomes with at least
    one detected; ``group`` is its smallest detected genome id and ``ambiguous`` the
    number of other members. Returns the detected genomes with these columns, and every
    group member's ``genome`` and ``group``.
    """
    members = pl.DataFrame(schema={"genome": pl.UInt32, "group": pl.UInt32})
    if detected.height == 0:
        empty = {"group": pl.UInt32, "ambiguous": pl.UInt32, "depth_lo": pl.Float64}
        return detected.with_columns(
            *(pl.lit(None, t).alias(c) for c, t in (empty | {"depth_hi": pl.Float64}).items())
        ), members
    mu = (
        carried.join(detected.select("genome", "depth", "present"), on="genome")
        .group_by("unit")
        .agg(mu=(pl.col("depth") * pl.col("present") * pl.col("c")).sum())
    )
    pairs = (
        carried.join(units.select("unit", "hits"), on="unit")
        .join(mu, on="unit", how="left")
        .with_columns(mu=pl.max_horizontal(pl.col("mu").fill_null(0.0), "hits"))
        .with_columns(x=pl.col("c") * pl.col("hits").sqrt() / pl.col("mu"))
    )
    candidates = screened["genome"].sort().to_numpy()
    unit_ids, row = np.unique(pairs["unit"].to_numpy(), return_inverse=True)
    col = np.searchsorted(candidates, pairs["genome"].to_numpy())
    x = sparse.csc_matrix(
        (pairs["x"].to_numpy(), (row, col)), shape=(len(unit_ids), len(candidates))
    )
    genomes = detected["genome"].to_numpy()
    at = np.searchsorted(candidates, genomes)
    # ponytail: dense D x S cosines and D x D pseudo-inverse, O(D^2 S) and O(D^3); per
    # component of genomes sharing units if thousands are detected.
    cross = (x[:, at].T @ x).toarray()
    norm = np.sqrt(np.asarray(x.multiply(x).sum(axis=0)).ravel())
    cosine = _div(cross, np.outer(norm[at], norm))
    share = screened.sort("genome")["containment"].to_numpy()
    ratio = np.minimum.outer(share[at], share) / np.maximum.outer(share[at], share)
    edges = sparse.coo_matrix((cosine >= AMBIGUOUS) & (ratio >= AMBIGUOUS))
    graph = sparse.csr_matrix(
        (np.ones(edges.nnz), (at[edges.row], edges.col)), shape=(len(candidates),) * 2
    )
    _, label = csgraph.connected_components(graph, directed=False)
    first = np.full(label.max() + 1, np.iinfo(np.int64).max)
    np.minimum.at(first, label[at], genomes.astype(np.int64))
    grouped = first[label] < np.iinfo(np.int64).max  # in a component with a detected genome
    members = pl.DataFrame(
        {"genome": candidates[grouped], "group": first[label][grouped]},
        schema={"genome": pl.UInt32, "group": pl.UInt32},
    )
    # Intervals: information over the detected genomes, w = π c.
    info = (x[:, at].T @ x[:, at]).toarray() * np.outer(
        detected["present"].to_numpy(), detected["present"].to_numpy()
    )
    with_mu = pairs.join(detected.select("genome"), on="genome", how="semi")
    per_unit = with_mu.group_by("unit").agg(pl.col("hits").first(), pl.col("mu").first())
    dof = max(per_unit.height - len(genomes), 1)
    chi2 = ((per_unit["hits"] - per_unit["mu"]) ** 2 / per_unit["mu"]).sum()
    phi = max(1.0, float(chi2) / dof)
    sd = np.sqrt(phi * np.clip(np.diag(np.linalg.pinv(info)), 0, None))
    depth = detected["depth"].to_numpy()
    sizes = members.group_by("group").agg(ambiguous=(pl.len() - 1).cast(pl.UInt32))
    out = detected.with_columns(
        group=pl.Series(first[label[at]], dtype=pl.UInt32),
        depth_lo=np.maximum(depth - 1.96 * sd, 0.0),
        depth_hi=depth + 1.96 * sd,
    ).join(sizes, on="group")
    return out, members


def lineages(table: pl.DataFrame) -> tuple[list[str], dict[int, list[str]]]:
    """Rank names and, per genome, its label at each: GTDB ranks from ``taxonomy`` (by
    prefix, ``d__`` to ``s__``; "" where missing), then its name at rank ``genome``.
    Without taxonomy, ``genome`` only."""
    names = dict(zip(table["genome"].to_list(), table["name"].to_list(), strict=True))
    if "taxonomy" not in table.columns:
        return ["genome"], {g: [n] for g, n in names.items()}
    out = {}
    for g, taxonomy in zip(table["genome"].to_list(), table["taxonomy"].to_list(), strict=True):
        parts = [x.strip() for x in (taxonomy or "").split(";")]
        by_rank = {GTDB_PREFIX[x[0]]: x for x in parts if x[1:3] == "__" and x[0] in GTDB_PREFIX}
        out[g] = [by_rank.get(r, "") for r in GTDB_RANKS] + [names[g]]
    return [*GTDB_RANKS, "genome"], out


def taxon_labels(table: pl.DataFrame) -> pl.DataFrame:
    """Per genome (``table`` has ``genome``, ``name``, ``group``, maybe ``taxonomy``) and
    reported rank (family and below, and genome), its ``taxon``. A genome in an ambiguity
    group of several takes the group's label: its members' common taxon at that rank, else
    their lowest common taxon above it, else ``unresolved`` (at genome rank, the members'
    names joined by commas). A missing rank takes the lowest taxon above it the same way."""
    ranks, lineage = lineages(table)
    members: dict[int, list[int]] = {}
    for g, grp in zip(table["genome"].to_list(), table["group"].to_list(), strict=True):
        members.setdefault(grp, []).append(g)
    rows = []
    for gs in members.values():
        paths = [lineage[g] for g in gs]
        common = 0
        while common < len(ranks) and len({p[common] for p in paths}) == 1:
            common += 1
        for i, rank in enumerate(ranks):
            if rank not in REPORTED_RANKS:
                continue
            if i < common and paths[0][i]:
                taxon = paths[0][i]
            else:
                above = [paths[0][j] for j in range(min(i, common)) if paths[0][j]]
                joined = ",".join(sorted(p[-1] for p in paths))
                taxon = above[-1] if above else joined if rank == "genome" else "unresolved"
            rows += [(g, rank, taxon) for g in gs]
    return pl.DataFrame(
        rows, schema={"genome": pl.UInt32, "rank": pl.String, "taxon": pl.String}, orient="row"
    )


def function_taxon(
    prof: pl.DataFrame,
    gi: GenomeIndex,
    detected: pl.DataFrame,
    carried: pl.DataFrame,
    units: pl.DataFrame,
    groups: pl.DataFrame,
) -> pl.DataFrame:
    """The function x taxon table: each unit's split hits (``hits_em``) attributed to taxa.

    Unit *u*'s classified share, ``explained`` / ``hits`` (Σ_G λ_G c_{G,u} capped at its raw
    hits), goes to the detected genomes carrying it by the fit's E-step responsibilities
    r_{G,u} ∝ λ_G π_G c_{G,u}; the rest, and all of units no detected genome carries, is
    ``unclassified``. Genomes roll up by :func:`taxon_labels` over their ambiguity
    ``groups`` (:func:`uncertainty`); units roll up to the genome
    index's Pfam labels (each label counts its units in full; accessions without version,
    ``PF01007``), or else to the profile's
    ``name``. Long format: ``function``, ``rank`` (``total`` and each reported rank),
    ``taxon``, ``hits_em``; per function and rank, rows sum to its ``total``.
    """
    split = prof.filter(pl.col("hits_em") > 0).select(
        pl.col("unit").cast(pl.UInt32), "hits_em", *(["name"] if "name" in prof.columns else [])
    )
    share = split.join(units.select("unit", "hits", "explained"), on="unit", how="left").select(
        "unit", classified=pl.col("hits_em") * (pl.col("explained") / pl.col("hits")).fill_null(0)
    )
    resp = (
        carried.join(detected.select("genome", "depth", "present"), on="genome")
        .with_columns(w=pl.col("depth") * pl.col("present") * pl.col("c"))
        .with_columns(r=pl.col("w") / pl.col("w").sum().over("unit"))
        .join(share, on="unit")
        .select("unit", "genome", hits_em=pl.col("classified") * pl.col("r"))
    )
    names = _names(gi)
    labels = taxon_labels(groups.join(names, on="genome"))
    if gi.unit_pfam is not None:
        accession = pl.col("pfam_accession")
        functions = gi.unit_pfam.select(
            pl.col("unit").cast(pl.UInt32),
            # MGnify stores the accession's number (1007 for PF01007); hmmsearch PF01007.23
            function=("PF" + accession.cast(pl.String).str.zfill(5))
            if gi.unit_pfam.schema["pfam_accession"].is_integer()
            else accession.str.replace(r"\.\d+$", ""),
        )
    elif "name" in split.columns:
        functions = split.select("unit", function=pl.col("name").cast(pl.String))
    else:
        functions = split.select("unit", function=pl.col("unit").cast(pl.String))
    classified = resp.join(labels, on="genome").select("unit", "rank", "taxon", "hits_em")
    ranks = [r for r in lineages(names)[0] if r in REPORTED_RANKS]
    per_rank = split.select("unit", "hits_em").join(pl.DataFrame({"rank": ranks}), how="cross")
    done = classified.group_by("unit", "rank").agg(done=pl.col("hits_em").sum())
    unclassified = (
        per_rank.join(done, on=["unit", "rank"], how="left")
        .with_columns(rest=pl.col("hits_em") - pl.col("done").fill_null(0.0))
        .filter(pl.col("rest") > 1e-9 * pl.col("hits_em"))  # rounding, not hits
        .select("unit", "rank", taxon=pl.lit("unclassified"), hits_em="rest")
    )
    total = split.select("unit", rank=pl.lit("total"), taxon=pl.lit(""), hits_em="hits_em")
    return (
        pl.concat([total, classified, unclassified])
        .join(functions, on="unit")
        .group_by("function", "rank", "taxon")
        .agg(pl.col("hits_em").sum())
        .sort("function", "rank", "taxon")
    )


def _names(gi: GenomeIndex) -> pl.DataFrame:
    columns = ["genome", "name", *(["taxonomy"] if "taxonomy" in gi.genomes.columns else [])]
    return gi.genomes.select(columns).cast({"genome": pl.UInt32})


def genome_profile(
    prof: pl.DataFrame,
    gi: GenomeIndex,
    *,
    min_containment: float = MIN_CONTAINMENT,
    min_units: int = MIN_UNITS,
) -> tuple[pl.DataFrame, dict[str, float | int], pl.DataFrame | None]:
    """:func:`fit_genomes` as a report: per detected genome (``depth`` > 0) its name,
    taxonomy, ``relative_abundance`` (share of Σ depth), ambiguity ``group`` and
    ``ambiguous``, and depth interval (:func:`uncertainty`); the sample's raw ``hits``,
    ``explained_fraction``
    (Σ explained / Σ hits: the rest is organisms not in the set, or units their genomes
    carry beyond the fit) and ``genome_equivalents`` (Σ depth); and the function x taxon
    table (:func:`function_taxon`; None for a profile without ``hits_em``)."""
    report: dict[str, float] = {}
    genomes, carried, screened = fit_genomes(
        prof, gi, min_containment=min_containment, min_units=min_units, report=report
    )
    units = explained(prof, genomes, carried)
    detected, groups = uncertainty(genomes.filter(pl.col("depth") > 0), carried, screened, units)
    detected = detected.join(_names(gi), on="genome")
    table = detected.with_columns(relative_abundance=pl.col("depth") / pl.col("depth").sum()).sort(
        "depth", descending=True
    )
    hits = float(units["hits"].sum())
    summary = report | {
        "hits": hits,
        "explained_fraction": float(units["explained"].sum()) / hits if hits else 0.0,
        "genome_equivalents": float(detected["depth"].sum()),
        "genomes_detected": detected.height,
        "ambiguity_groups": detected["group"].n_unique(),
    }
    stratified = (
        function_taxon(prof, gi, detected, carried, units, groups)
        if "hits_em" in prof.columns
        else None
    )
    return table, summary, stratified
