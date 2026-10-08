"""The genome panel (phase 12, step 2): per species, reference genomes with their carriage,
copies and k-mer survival per unit, cheap distances between them as a k-nearest-neighbour
graph, and the species' rank-shrunk carriage frequency *q*. Built once per genome set; no
trees, no DNA, no per-unit model fitting.

Per (genome, unit): carriage when one of its proteins (families) has the unit as best unit;
copies *n*, the number of such proteins; survival *f*, the share of the unit's tier-2 kept
k-mers found in them (catalogue: the best family's; a genome set: all its proteins').

Distances within a species: Jaccard on carried units, each unit one genome lacks weighted by
that genome's completeness (an incomplete MAG may have lost it), times one minus the mean
survival difference on the units both carry: d = 1 - J_c (1 - mean |f_g - f_h|). At most
``max_per_species`` genomes per species are kept, by farthest-point sampling from the
representative (catalogue) or the most complete genome, so a panel spans each species'
diversity. *q* comes from every genome of the species, not only the panel's (species.py's
shrinkage: n + α q_parent over N + α at family, genus and species, N = Σ completeness).

Files written to the panel directory:

- ``species.tsv``: ``species`` (row), ``id``, ``name``, ``taxonomy``, ``genomes`` (all kept),
  ``panel`` (in the panel);
- ``genomes.parquet``: ``genome`` (panel row), ``name``, ``species``, ``completeness``,
  ``fps_rank`` (0: the start);
- ``carriage.parquet``: ``genome``, ``unit``, ``n``, ``f`` (panel genomes' carried units);
- ``neighbours.parquet``: ``genome``, ``neighbour``, ``distance``, ``rank``;
- ``species_units.parquet``: ``species``, ``unit``, ``q``, ``f_mean``, ``n_mean`` (the
  carriers' means, else the genus', all genomes'), units with q ≥ ``Q_FLOOR`` or carried by
  a panel genome;
- ``units.parquet``: ``unit``, ``m_g`` (tier-2 kept k-mers) and ``m_dense`` if the index has
  a dense tier;
- ``held_out.parquet``: carriage of the genomes left out with ``exclude`` (benchmark truth);
- ``unit_pfam.parquet`` (if the index has it) and ``meta.json`` with the index checksum.
"""

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Final

import numpy as np
import polars as pl
from scipy import sparse

from kmer_functional_profiler.genomes import PROTEINS_PER_BATCH, best_units, protein_hits
from kmer_functional_profiler.index import Index
from kmer_functional_profiler.mask import _sha256
from kmer_functional_profiler.query import _ranges
from kmer_functional_profiler.species import (
    CATALOGUE_METADATA,
    RANKS,
    _counts,
    _load_index,
    _pangenome,
    fit_alpha,
    lineages,
    prevalence,
)

MAX_PER_SPECIES: Final = 50
NEIGHBOURS: Final = 10


# --- Carriage rows from either source ---


def _m_arrays(index: Index) -> tuple[np.ndarray, np.ndarray | None]:
    """Per index unit, tier-2 kept k-mers and dense ones (None without a dense tier)."""
    dense = index.units["m_dense"] if "m_dense" in index.units.schema else None
    return index.units["m_g"], dense


def genome_set_rows(index: Index, genome_index: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Genomes (``genome``, ``name``, ``taxonomy``, ``completeness`` in [0, 1] if given) and
    carriage rows (``genome``, ``unit``, ``n``, ``f``) of an ``annotate-genomes`` index."""
    genomes = pl.read_csv(genome_index / "genomes.tsv", separator="\t").cast({"genome": pl.UInt32})
    if "completeness" in genomes.columns:
        genomes = genomes.with_columns(pl.col("completeness").cast(pl.Float64) / 100)
    copies = (
        pl.read_parquet(genome_index / "genome_best.parquet")
        .group_by("genome", "unit")
        .agg(n=pl.len())
    )
    offsets = np.load(genome_index / "unit_offsets.npy", mmap_mode="r")
    units = np.unique(copies["unit"].to_numpy()).astype(np.int64)
    starts = offsets[units].astype(np.int64)
    lengths = offsets[units + 1].astype(np.int64) - starts
    rows = _ranges(starts, lengths)
    kmers = pl.DataFrame(
        {
            "unit": np.repeat(units, lengths).astype(np.uint32),
            "genome": np.load(genome_index / "genome.npy", mmap_mode="r")[rows],
            "kmers": np.load(genome_index / "kmers.npy", mmap_mode="r")[rows].astype(np.float64),
        }
    )
    carried = copies.with_columns(
        pl.col("unit").cast(pl.UInt32), pl.col("genome").cast(pl.UInt32)
    ).join(kmers, on=["genome", "unit"], how="left")
    m_g = np.asarray(index.units["m_g"][carried["unit"].to_numpy()], dtype=np.float64)
    carried = carried.select(
        "genome",
        "unit",
        n=pl.col("n").cast(pl.UInt16),
        f=(pl.col("kmers").fill_null(0) / pl.Series(np.maximum(m_g, 1))).clip(0.0, 1.0),
    )
    return genomes.select(
        "genome", "name", *(c for c in ("taxonomy", "completeness") if c in genomes.columns)
    ), carried


def catalogue_rows(
    index: Index, catalogue: Path, species: set[str] | None = None
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Genomes (with ``species``: the representative) and carriage rows of an MGnify genome
    catalogue: each species' Panaroo families annotated, carriage from the Rtab. Survival is
    the best family's, so a catalogue has one *f* per (species, unit) and family."""
    metadata = pl.read_csv(catalogue / CATALOGUE_METADATA, separator="\t", infer_schema=False)
    if species is not None:
        metadata = metadata.filter(pl.col("Species_rep").is_in(list(species)))
    genomes = metadata.select(
        genome=pl.int_range(metadata.height, dtype=pl.UInt32),
        name="Genome",
        species="Species_rep",
        taxonomy="Lineage",
        completeness=pl.col("Completeness").cast(pl.Float64) / 100,
    )
    genome_ids = dict(zip(genomes["name"], genomes["genome"], strict=True))
    m_g = index.units["m_g"]
    parts = []
    for rep in genomes["species"].unique(maintain_order=True):
        proteins, presence, columns = _pangenome(catalogue, rep, genome_ids)
        batches = (proteins[i : i + PROTEINS_PER_BATCH]
                   for i in range(0, len(proteins), PROTEINS_PER_BATCH))  # fmt: skip
        hits, lengths = protein_hits(index, batches)
        best = best_units(index, hits, lengths)
        if best.height == 0:
            continue
        kmers = (
            hits.with_columns(pl.col("protein").cast(pl.UInt32))
            .join(best.select("protein", "unit"), on=["protein", "unit"], how="semi")
            .group_by("protein", "unit")
            .agg(kmers=pl.col("hash").n_unique())
        )
        best = best.join(kmers, on=["protein", "unit"]).with_columns(
            f=(
                pl.col("kmers")
                / pl.Series(np.maximum(m_g[best["unit"].to_numpy()], 1)).cast(pl.Float64)
            ).clip(0.0, 1.0)
        )
        coo = presence.tocoo()
        present = pl.DataFrame(
            {"protein": coo.row.astype(np.uint32), "genome": columns[coo.col].astype(np.uint32)}
        )
        parts.append(
            best.join(present, on="protein")
            .group_by("genome", "unit")
            .agg(n=pl.len().cast(pl.UInt16), f=pl.col("f").max())
        )
    schema = pl.Schema({"genome": pl.UInt32, "unit": pl.UInt32, "n": pl.UInt16, "f": pl.Float64})
    rows = pl.concat(
        [pl.DataFrame(schema=schema), *(p.select(list(schema)).cast(schema) for p in parts)]
    )
    return genomes, rows


# --- Distances and farthest-point sampling ---


@dataclass(frozen=True)
class _Carriage:
    """One species' genomes as sparse carriage (``x``) and survival (``f``) rows."""

    x: sparse.csr_matrix
    f: sparse.csr_matrix
    size: np.ndarray
    completeness: np.ndarray

    @classmethod
    def build(cls, rows: pl.DataFrame, n_genomes: int, completeness: np.ndarray) -> "_Carriage":
        units, col = np.unique(rows["unit"].to_numpy(), return_inverse=True)
        g = rows["g"].to_numpy()
        shape = (n_genomes, len(units))
        x = sparse.csr_matrix((np.ones(len(g)), (g, col)), shape=shape)
        f = sparse.csr_matrix((rows["f"].to_numpy().astype(np.float64), (g, col)), shape=shape)
        return cls(x, f, np.asarray(x.sum(1)).ravel(), completeness)

    def distance(self, i: int, to: np.ndarray | None = None) -> np.ndarray:
        """d(genome i, each genome in ``to``, default all) (module docstring)."""
        to = np.arange(self.x.shape[0]) if to is None else to
        xi = self.x[i].toarray().ravel()
        fi = self.f[i].toarray().ravel()
        x, f = self.x[to], self.f[to]
        shared = x @ xi
        a_only, b_only = self.size[to] - shared, self.size[i] - shared
        union = shared + self.completeness[i] * a_only + self.completeness[to] * b_only
        jaccard = np.divide(shared, union, out=np.zeros(len(to)), where=union > 0)
        coo = f.tocoo()
        both = xi[coo.col] > 0
        diff = np.bincount(
            coo.row[both], np.abs(coo.data[both] - fi[coo.col[both]]), minlength=len(to)
        )
        mean_diff = np.divide(diff, shared, out=np.zeros(len(to)), where=shared > 0)
        return np.asarray(1 - jaccard * (1 - mean_diff))


def farthest_points(carriage: _Carriage, start: int, k: int) -> list[int]:
    """Up to ``k`` genomes by farthest-point sampling from ``start``."""
    n = carriage.x.shape[0]
    chosen = [start]
    nearest = carriage.distance(start)
    nearest[start] = -1
    while len(chosen) < min(k, n):
        nxt = int(np.argmax(nearest))
        chosen.append(nxt)
        nearest = np.minimum(nearest, carriage.distance(nxt))
        nearest[chosen] = -1
    return chosen


# --- Writing the panel ---


def _fallback_means(rows: pl.DataFrame, lineage: pl.DataFrame) -> dict[str, pl.DataFrame]:
    """Per rank (species, genus, root) the carriers' mean ``f`` and ``n`` per unit."""
    joined = rows.join(lineage, on="genome")
    out = {}
    for r in ("species", "genus"):
        out[r] = (
            joined.drop_nulls(r)
            .group_by(r, "unit")
            .agg(
                pl.col("f").mean().alias(f"f_{r}"),
                pl.col("n").cast(pl.Float64).mean().alias(f"n_{r}"),
            )  # fmt: skip
        )
    out["root"] = joined.group_by("unit").agg(
        f_root=pl.col("f").mean(), n_root=pl.col("n").cast(pl.Float64).mean()
    )
    return out


def write_panel(
    out: Path,
    index_dir: Path,
    m_g: np.ndarray,
    genomes: pl.DataFrame,
    rows: pl.DataFrame,
    exclude: frozenset[str] = frozenset(),
    meta: dict[str, object] | None = None,
    *,
    m_dense: np.ndarray | None = None,
    representatives: frozenset[str] = frozenset(),
    max_per_species: int = MAX_PER_SPECIES,
    neighbours: int = NEIGHBOURS,
    alpha: float | None = None,
    completeness: bool = True,
) -> dict[str, object]:
    """A panel at ``out`` (module docstring) from ``genomes`` (``genome``, ``name``,
    ``species`` or ``taxonomy``, optional ``completeness`` in [0, 1]) and carriage ``rows``
    (``genome``, ``unit``, ``n``, ``f``); ``m_g`` (``m_dense``) the index's tier-2 (dense)
    kept k-mers per unit. ``representatives`` are genome names FPS starts
    from (else the most complete genome). Ablations as ``species-index``: ``alpha`` fixed at
    every rank, ``completeness`` off."""
    out.mkdir(parents=True, exist_ok=True)
    lineage_all = lineages(genomes)
    held = genomes.filter(pl.col("name").is_in(list(exclude)))
    rows.join(held.select("genome", "name"), on="genome").join(
        lineage_all.select("genome", "species"), on="genome"
    ).select("name", "species", "unit", "n", "f").sort("name", "unit").write_parquet(
        out / "held_out.parquet"
    )
    kept = genomes.filter(~pl.col("name").is_in(list(exclude)))
    if "completeness" not in kept.columns or not completeness:
        kept = kept.with_columns(completeness=pl.lit(1.0))
    kept = kept.with_columns(pl.col("completeness").fill_null(1.0).clip(0.01, 1.0))
    lineage = lineage_all.join(kept.select("genome"), on="genome", how="semi")
    rows = rows.join(kept.select("genome"), on="genome", how="semi")

    # q from every kept genome of the species (completeness-weighted counts, α per rank)
    carried = rows.select("genome", "unit", c=pl.lit(1.0))
    weights = kept.select("genome", "completeness")
    n, sizes = _counts(carried.select("genome", "unit"), lineage, weights)
    fits = {r: fit_alpha(carried.select("genome", "unit"), lineage, r) for r in RANKS}
    alphas = {r: a if alpha is None else alpha for r, (a, _) in fits.items()}
    pairs = prevalence(carried, lineage, n, sizes, alphas).select("species", "unit", "q")

    # farthest-point sampling and neighbour lists per species
    info = kept.join(lineage.select("genome", "species"), on="genome").sort("genome")
    selected, edges = [], []
    for (sp,), group in info.group_by("species", maintain_order=True):
        group = group.sort("genome").with_row_index("g")
        local = rows.join(group.select("genome", "g"), on="genome")
        comp = group["completeness"].to_numpy()
        within = _Carriage.build(local, group.height, comp)
        names = group["name"].to_list()
        reps = [i for i, x in enumerate(names) if x in representatives]
        start = reps[0] if reps else int(np.lexsort((names, -comp))[0])
        chosen = farthest_points(within, start, max_per_species)
        selected.append(
            group[chosen].select("genome", "name", "completeness", species=pl.lit(sp))
            .with_columns(fps_rank=pl.int_range(len(chosen), dtype=pl.UInt16))
        )  # fmt: skip
        sub = np.asarray(chosen)
        for a, i in enumerate(chosen):
            d = within.distance(i, sub)
            d[a] = np.inf
            order = np.argsort(d, kind="stable")[: min(neighbours, len(chosen) - 1)]
            edges += [(group["genome"][i], group["genome"][int(sub[j])], float(d[j]), r)
                      for r, j in enumerate(order)]  # fmt: skip
    panel = pl.concat(selected)

    # species table, sorted by id
    first = lineage.unique("species", keep="first", maintain_order=True)
    taxonomy = (
        kept.select("genome", "taxonomy") if "taxonomy" in kept.columns
        else kept.select("genome", taxonomy=pl.lit(None, pl.String))
    )  # fmt: skip
    species = (
        first.join(taxonomy, on="genome", how="left")
        .join(lineage.group_by("species").agg(genomes=pl.len()), on="species")
        .join(panel.group_by("species").agg(panel=pl.len()), on="species")
        .sort("species")
        .with_row_index("row")
        .with_columns(
            name=pl.when(pl.col("species").str.starts_with("genome:"))
            .then(pl.col("species").str.slice(7))
            .otherwise(
                pl.col("taxonomy").str.extract(r"s__([^;]+)$", 1).fill_null(pl.col("species"))
            )
        )
    )
    ids = species.select("species", sp=pl.col("row").cast(pl.UInt32))
    panel = (
        panel.join(ids, on="species")
        .sort("sp", "fps_rank")
        .with_row_index("row")
        .with_columns(pl.col("row").cast(pl.UInt32))
    )
    rename = panel.select("genome", "row")
    carriage = (
        rows.join(rename, on="genome")
        .select(genome="row", unit=pl.col("unit").cast(pl.UInt32),
                n=pl.col("n").cast(pl.UInt16), f=pl.col("f").cast(pl.Float32))
        .sort("genome", "unit")
    )  # fmt: skip
    edge_table = (
        pl.DataFrame(edges, schema=["genome", "neighbour", "distance", "rank"], orient="row")
        .cast({"genome": pl.UInt32, "neighbour": pl.UInt32})
        .join(rename, on="genome")
        .join(rename.rename({"genome": "neighbour", "row": "nrow"}), on="neighbour")
        .select(genome="row", neighbour="nrow", distance=pl.col("distance").cast(pl.Float32),
                rank=pl.col("rank").cast(pl.UInt8))
        .sort("genome", "rank")
    ) if edges else pl.DataFrame(schema={"genome": pl.UInt32, "neighbour": pl.UInt32,
                                         "distance": pl.Float32, "rank": pl.UInt8})  # fmt: skip

    # species units: q >= Q_FLOOR, or carried by a panel genome (q then from the raw count)
    raw = (
        n.filter(pl.col("rank") == "species")
        .join(sizes.filter(pl.col("rank") == "species"), on=["rank", "clade"])
        .select(
            species="clade", unit=pl.col("unit").cast(pl.UInt32), q_raw=pl.col("n") / pl.col("N")
        )
    )
    in_panel = (
        carriage.join(panel.select(genome="row", sp="sp"), on="genome")
        .select("sp", "unit").unique()
        .join(ids, on="sp").select("species", "unit")
    )  # fmt: skip
    means = _fallback_means(rows, lineage)
    genus = first.select("species", "genus")
    units_table = (
        pl.concat([pairs.select("species", pl.col("unit").cast(pl.UInt32)), in_panel])
        .unique()
        .join(pairs.with_columns(pl.col("unit").cast(pl.UInt32)), on=["species", "unit"],
              how="left")
        .join(raw, on=["species", "unit"], how="left")
        .join(genus, on="species", how="left")
        .join(means["species"], on=["species", "unit"], how="left")
        .join(means["genus"], on=["genus", "unit"], how="left")
        .join(means["root"], on="unit", how="left")
        .join(ids, on="species")
        .select(
            species="sp",
            unit="unit",
            q=pl.coalesce("q", "q_raw").clip(0.0, 1.0),
            f_mean=pl.coalesce("f_species", "f_genus", "f_root", pl.lit(0.5)),
            n_mean=pl.coalesce("n_species", "n_genus", "n_root", pl.lit(1.0)),
        )
        .sort("species", "unit")
    )  # fmt: skip

    species.select(
        species=pl.col("row").cast(pl.UInt32), id="species", name="name", taxonomy="taxonomy",
        genomes="genomes", panel="panel",
    ).write_csv(out / "species.tsv", separator="\t")  # fmt: skip
    panel.select(
        genome="row", name="name", species="sp", completeness="completeness", fps_rank="fps_rank"
    ).write_parquet(out / "genomes.parquet")
    carriage.write_parquet(out / "carriage.parquet")
    edge_table.write_parquet(out / "neighbours.parquet")
    units_table.write_parquet(out / "species_units.parquet")
    used = np.unique(np.r_[units_table["unit"].to_numpy(), carriage["unit"].to_numpy()])
    pl.DataFrame(
        {"unit": used.astype(np.uint32), "m_g": np.asarray(m_g[used], dtype=np.uint32)}
        | ({} if m_dense is None else {"m_dense": np.asarray(m_dense[used], dtype=np.uint32)})
    ).write_parquet(out / "units.parquet")
    if (index_dir / "unit_pfam.parquet").exists():
        pl.read_parquet(index_dir / "unit_pfam.parquet").write_parquet(out / "unit_pfam.parquet")
    meta = (meta or {}) | {
        "index_meta_sha256": _sha256(index_dir / "meta.json"),
        "units": len(m_g),
        "species": species.height,
        "genomes": kept.height,
        "panel_genomes": panel.height,
        "excluded": held.height,
        "max_per_species": max_per_species,
        "neighbours": neighbours,
        "alpha": alphas,
        "alpha_fixed": alpha,
        "completeness_weighted": completeness,
        "held_out_log_loss": {r: losses for r, (_, losses) in fits.items()},
    }
    (out / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    return meta


def panel_from_genomes(
    index_dir: str | Path, genome_index: str | Path, out: str | Path, **kwargs: object
) -> dict[str, object]:
    """A panel from an ``annotate-genomes`` genome index (species from GTDB ``s__``)."""
    index_dir, genome_index = Path(index_dir), Path(genome_index)
    gmeta = json.loads((genome_index / "meta.json").read_text())
    if gmeta["index_meta_sha256"] != _sha256(index_dir / "meta.json"):
        raise ValueError(f"{genome_index} was annotated with another index than {index_dir}")
    index = _load_index(index_dir)
    genomes, rows = genome_set_rows(index, genome_index)
    meta: dict[str, object] = {
        "source": "genomes",
        "genome_index_meta_sha256": _sha256(genome_index / "meta.json"),
    }
    m_g, m_dense = _m_arrays(index)
    return write_panel(Path(out), index_dir, m_g, genomes, rows, meta=meta, m_dense=m_dense,
                       **kwargs)  # type: ignore[arg-type]  # fmt: skip


def panel_from_catalogue(
    index_dir: str | Path,
    catalogue: str | Path,
    out: str | Path,
    species: set[str] | None = None,
    **kwargs: object,
) -> dict[str, object]:
    """A panel from an MGnify genome catalogue (FPS starts at each species'
    representative)."""
    index_dir, catalogue = Path(index_dir), Path(catalogue)
    index = _load_index(index_dir)
    genomes, rows = catalogue_rows(index, catalogue, species)
    meta: dict[str, object] = {"source": "catalogue", "catalogue": str(catalogue)}
    reps = frozenset(genomes["species"].unique().to_list())
    m_g, m_dense = _m_arrays(index)
    return write_panel(Path(out), index_dir, m_g, genomes, rows, meta=meta, m_dense=m_dense,
                       representatives=reps, **kwargs)  # type: ignore[arg-type]  # fmt: skip


# --- Reading it ---


@dataclass(frozen=True)
class SpeciesPanel:
    """One species of a panel as dense arrays over its panel genomes (G) and units (U)."""

    species: int
    units: np.ndarray  # U index unit ids
    q: np.ndarray
    f_mean: np.ndarray
    n_mean: np.ndarray
    m_g: np.ndarray
    genomes: np.ndarray  # G panel rows
    names: list[str]
    completeness: np.ndarray
    carried: np.ndarray  # G x U bool
    xt: np.ndarray  # G x U: carriage corrected for completeness
    n: np.ndarray  # G x U copies (the carriers' mean where not carried)
    f: np.ndarray  # G x U survival (the carriers' mean where not carried)
    edges: np.ndarray  # E x 2 directed (g, h) positions in genomes; (0, 0) for one genome
    edge_prior: np.ndarray  # E: 1 / (G x neighbours of g)

    def carriage(self, g: int, h: int, t: float, ell: float) -> np.ndarray:
        """p_u at a point t between references g and h, pulled ell towards the species."""
        return np.asarray((1 - ell) * ((1 - t) * self.xt[g] + t * self.xt[h]) + ell * self.q)


def corrected_carriage(q: np.ndarray, completeness: np.ndarray) -> np.ndarray:
    """x̃ of a unit a genome of ``completeness`` c does not show: q (1 - c) / (1 - q c)."""
    return np.asarray(q * (1 - completeness) / np.maximum(1 - q * completeness, 1e-12))


class Panel:
    """A panel written by :func:`write_panel`."""

    def __init__(self, path: str | Path) -> None:
        path = Path(path)
        self.path = path
        self.meta = json.loads((path / "meta.json").read_text())
        self.species = pl.read_csv(
            path / "species.tsv", separator="\t", schema_overrides={"taxonomy": pl.String}
        )
        self.genomes = pl.read_parquet(path / "genomes.parquet")
        self.carriage_rows = pl.read_parquet(path / "carriage.parquet")
        self.neighbours = pl.read_parquet(path / "neighbours.parquet")
        self.species_units = pl.read_parquet(path / "species_units.parquet")
        self.units = pl.read_parquet(path / "units.parquet")
        pfam = path / "unit_pfam.parquet"
        self.unit_pfam = pl.read_parquet(pfam) if pfam.exists() else None

    def check_index(self, index_dir: str | Path) -> None:
        """Raise if ``index_dir`` is not the index the panel was built with."""
        if _sha256(Path(index_dir) / "meta.json") != self.meta["index_meta_sha256"]:
            raise ValueError(f"{self.path} was built with another index than {index_dir}")

    def species_panel(self, species: int) -> SpeciesPanel:
        su = self.species_units.filter(pl.col("species") == species).sort("unit")
        units = su["unit"].to_numpy().astype(np.int64)
        m_g = (
            su.select("unit").join(self.units, on="unit", how="left")["m_g"]
            .fill_null(1).to_numpy().astype(np.float64)
        )  # fmt: skip
        g = self.genomes.filter(pl.col("species") == species).sort("genome")
        rows = g["genome"].to_numpy().astype(np.int64)
        q = su["q"].to_numpy().astype(np.float64)
        f_mean = su["f_mean"].to_numpy().astype(np.float64)
        n_mean = su["n_mean"].to_numpy().astype(np.float64)
        comp = g["completeness"].to_numpy().astype(np.float64)
        cr = self.carriage_rows.filter(pl.col("genome").is_in(rows.tolist()))
        gi = np.searchsorted(rows, cr["genome"].to_numpy())
        ui = np.searchsorted(units, cr["unit"].to_numpy())
        ok = (ui < len(units)) & (units[np.minimum(ui, len(units) - 1)] == cr["unit"].to_numpy())
        carried = np.zeros((len(rows), len(units)), dtype=bool)
        n = np.broadcast_to(n_mean, carried.shape).copy()
        f = np.broadcast_to(f_mean, carried.shape).copy()
        carried[gi[ok], ui[ok]] = True
        n[gi[ok], ui[ok]] = cr["n"].to_numpy()[ok]
        f[gi[ok], ui[ok]] = cr["f"].to_numpy()[ok]
        xt = np.where(carried, 1.0, corrected_carriage(q[None, :], comp[:, None]))
        nb = self.neighbours.filter(pl.col("genome").is_in(rows.tolist()))
        if nb.height == 0:
            edges = np.zeros((1, 2), dtype=np.int64)
            prior = np.ones(1)
        else:
            edges = np.stack(
                [np.searchsorted(rows, nb["genome"].to_numpy()),
                 np.searchsorted(rows, nb["neighbour"].to_numpy())], 1
            )  # fmt: skip
            degree = np.bincount(edges[:, 0], minlength=len(rows))
            prior = 1 / (len(rows) * degree[edges[:, 0]])
        return SpeciesPanel(
            species=species, units=units, q=q, f_mean=f_mean, n_mean=n_mean, m_g=m_g,
            genomes=rows, names=g["name"].to_list(), completeness=comp, carried=carried,
            xt=xt, n=n, f=f, edges=edges, edge_prior=prior,
        )  # fmt: skip
