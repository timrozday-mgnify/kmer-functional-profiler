"""Species model (phase 11, steps 7-9; plan: Genome mode, Species model).

A species index links species to the units their genomes carry: per (species, unit) the
prevalence prior *q* (how often the species' genomes carry the unit, shrunk towards genus,
family and all genomes) and *e*, the expected raw tier-2 hits per genome copy at depth 1
(the mean content of the carriers). Built from either source by one aggregation over
(genome, carried unit, content) rows:

- :func:`species_index_from_catalogue`: an MGnify genome catalogue. Each species'
  pangenome families (Panaroo ``pan-genome.fna``, translated) are annotated as proteins and
  their carriage per conspecific genome read from ``gene_presence_absence.Rtab``; a
  single-genome species' proteins are its representative's ``.faa``. Only the pangenome is
  hashed, not every genome.
- :func:`species_index_from_genomes`: a genome index of ``annotate-genomes`` whose
  ``genomes.tsv`` has a GTDB-style ``taxonomy`` (species from ``s__``; a genome without one
  is its own species).

A genome carries a unit when one of its proteins (families) has it as best unit
(:func:`~kmer_functional_profiler.genomes.best_units`); its content on the unit counts every
protein's hits there. Prevalence is corrected for MAG incompleteness: each genome sees a
carried unit with probability its completeness, so clade sizes are Σ completeness.

Files written to the index directory:

- ``species.tsv``: ``species`` (row number), ``id``, ``name``, ``taxonomy``, ``genomes``,
  ``completeness`` (Σ), ``units`` (pairs), ``content`` (Σ q e);
- ``species_offsets.npy``, ``s_unit.npy``, ``s_q.npy``, ``s_e.npy``: pairs species-major;
  ``unit_offsets.npy``, ``u_species.npy``, ``u_q.npy``, ``u_e.npy``: the same unit-major;
- ``genomes.parquet``: ``genome``, ``name``, ``species`` (id), ``completeness``;
- ``carriage.parquet``, ``clades.parquet``, ``lineage.parquet``: carriage counts per clade
  in kfp-prior's format (it reads them as its carriage table);
- ``held_out.parquet``: units carried by each genome left out with ``exclude`` (benchmark
  truth);
- ``unit_pfam.parquet`` (if the index has it) and ``meta.json``.
"""

import json
from collections.abc import Iterator
from pathlib import Path
from typing import Final

import numpy as np
import polars as pl
from scipy import sparse
from scipy.special import digamma, expit, gammaln
from scipy.stats import gamma

from kmer_functional_profiler import _core
from kmer_functional_profiler.genomes import (
    GTDB_PREFIX,
    MIN_CONTAINMENT,
    MIN_UNITS,
    PROTEINS_PER_BATCH,
    best_units,
    fasta,
    protein_hits,
    taxon_labels,
    unit_functions,
)
from kmer_functional_profiler.index import Index
from kmer_functional_profiler.mask import _sha256
from kmer_functional_profiler.query import _div, _ranges

RANKS: Final = ("family", "genus", "species")  # top down
ALPHAS: Final = np.geomspace(0.01, 100, 21)
MAX_HOLDOUT: Final = 200  # held-out genomes per rank when fitting α
Q_FLOOR: Final = 0.05  # prevalences below this give no pair
CLIP: Final = 1e-6
CATALOGUE_METADATA: Final = "genomes-all_metadata.tsv"


# --- Prevalence and its shrinkage (moved from kfp-prior's carriage table) ---


def lineages(genomes: pl.DataFrame) -> pl.DataFrame:
    """``genome`` and its ``family`` and ``genus`` from a GTDB-style ``taxonomy`` (by
    prefix; null where missing) and its ``species``: the ``species`` column if given, else
    ``s__`` of the taxonomy, else ``genome:<name>`` (its own species)."""
    taxa = genomes["taxonomy"] if "taxonomy" in genomes.columns else [None] * genomes.height
    given = genomes["species"] if "species" in genomes.columns else [None] * genomes.height
    rows = []
    for g, name, taxonomy, species in zip(genomes["genome"], genomes["name"], taxa, given,
                                          strict=True):  # fmt: skip
        parts = [x.strip() for x in (taxonomy or "").split(";")]
        by_rank = {GTDB_PREFIX[x[0]]: x for x in parts if x[1:3] == "__" and len(x) > 3
                   and x[0] in GTDB_PREFIX}  # fmt: skip
        rows.append((g, by_rank.get("family"), by_rank.get("genus"),
                     species or by_rank.get("species") or f"genome:{name}"))  # fmt: skip
    schema = {"genome": pl.UInt32, **dict.fromkeys(RANKS, pl.String)}
    return pl.DataFrame(rows, schema=schema, orient="row")


def _counts(
    carried: pl.DataFrame, lineage: pl.DataFrame, completeness: pl.DataFrame | None = None
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Carriage counts ``rank``, ``clade``, ``unit``, ``n`` and clade sizes ``rank``,
    ``clade``, ``N`` (Σ completeness, 1 per genome without it), per rank and for all
    genomes (rank ``root``, clade "")."""
    long = pl.concat(
        [lineage.select("genome", rank=pl.lit("root"), clade=pl.lit(""))]
        + [lineage.select("genome", rank=pl.lit(r), clade=pl.col(r)) for r in RANKS]
    ).drop_nulls("clade")
    n = carried.join(long, on="genome").group_by("rank", "clade", "unit").agg(n=pl.len())
    weight = (
        long.join(completeness, on="genome", how="left").with_columns(
            pl.col("completeness").fill_null(1.0)
        )
        if completeness is not None
        else long.with_columns(completeness=pl.lit(1.0))
    )
    sizes = weight.group_by("rank", "clade").agg(N=pl.col("completeness").sum())
    return n.cast({"n": pl.Int64}), sizes.cast({"N": pl.Float64})


def _parent(lineage: pl.DataFrame, rank: str) -> pl.DataFrame:
    """Per genome, the nearest clade above ``rank`` (``p_rank``, ``p``), else root."""
    above = RANKS[: RANKS.index(rank)]
    p_rank, p = pl.lit("root"), pl.lit("")
    for r in above:  # top down: the last non-null wins
        p_rank = pl.when(pl.col(r).is_not_null()).then(pl.lit(r)).otherwise(p_rank)
        p = pl.when(pl.col(r).is_not_null()).then(pl.col(r)).otherwise(p)
    return lineage.select("genome", p_rank=p_rank, p=p)


def fit_alpha(
    carried: pl.DataFrame, lineage: pl.DataFrame, rank: str, seed: int = 0
) -> tuple[float, dict[str, float]]:
    """α at ``rank`` by leave-one-genome-out log loss: each held-out genome's carriage of
    the units its parent clade carries (its own clade's, when the parent is root), predicted
    from the rest of its clade shrunk towards the rest of the parent's frequency. Returns
    α and the losses of the best α, of no shrinkage (α → 0) and of the parent alone
    (α → ∞), per prediction. Counts are raw (no completeness weighting)."""
    n, sizes = _counts(carried, lineage)
    eligible = (
        lineage.select("genome", clade=pl.col(rank))
        .join(sizes.filter(pl.col("rank") == rank), on="clade")
        .filter(pl.col("N") >= 2)
    )
    if eligible.height == 0:
        return 1.0, {}
    held = eligible.sample(min(MAX_HOLDOUT, eligible.height), seed=seed).join(
        _parent(lineage, rank), on="genome"
    )
    by_parent = held.filter(pl.col("p_rank") != "root").join(
        n, left_on=["p_rank", "p"], right_on=["rank", "clade"]
    )
    by_clade = held.filter(pl.col("p_rank") == "root").join(
        n.filter(pl.col("rank") == rank), on="clade"
    )
    rows = (
        pl.concat(
            [
                by_parent.select("genome", "clade", "p_rank", "p", "unit"),
                by_clade.select("genome", "clade", "p_rank", "p", "unit"),
            ]
        )  # fmt: skip
        .join(
            n.filter(pl.col("rank") == rank).select("clade", "unit", n_c="n"),
            on=["clade", "unit"],
            how="left",
        )  # fmt: skip
        .join(sizes.filter(pl.col("rank") == rank).select("clade", N_c="N"), on="clade")
        .join(
            n.select(p_rank="rank", p="clade", unit="unit", n_p="n"),
            on=["p_rank", "p", "unit"],
            how="left",
        )  # fmt: skip
        .join(sizes.select(p_rank="rank", p="clade", N_p="N"), on=["p_rank", "p"])
        .join(carried.with_columns(x=pl.lit(1)), on=["genome", "unit"], how="left")
        .with_columns(pl.col("n_c", "n_p", "x").fill_null(0))
    )
    x = rows["x"].to_numpy().astype(np.float64)
    rest = rows["n_c"].to_numpy() - x
    parent = (rows["n_p"].to_numpy() - x) / np.maximum(rows["N_p"].to_numpy() - 1, 1)
    size = rows["N_c"].to_numpy() - 1.0

    def loss(alpha: float) -> float:
        q = np.clip((rest + alpha * parent) / (size + alpha), CLIP, 1 - CLIP)
        return float(-(x * np.log(q) + (1 - x) * np.log1p(-q)).mean())

    losses = [loss(a) for a in ALPHAS]
    best = float(ALPHAS[int(np.argmin(losses))])
    return best, {"best": min(losses), "no_shrinkage": loss(1e-9), "parent_only": loss(1e9)}


def prevalence(
    carried: pl.DataFrame,
    lineage: pl.DataFrame,
    n: pl.DataFrame,
    sizes: pl.DataFrame,
    alpha: dict[str, float],
) -> pl.DataFrame:
    """Per species and candidate unit (any unit its genus carries; its own where it has no
    genus): ``q``, the species' completeness-corrected carriage frequency shrunk down the
    ranks from all genomes, (n + α q_parent) / (N + α) at family, genus and species, capped
    at 1; and ``e``, the mean content of the species' carriers (else its genus', family's,
    all genomes'). Pairs with q < ``Q_FLOOR`` are dropped."""
    species = lineage.unique("species", keep="first").select(RANKS)
    top = pl.coalesce("genus", "species")
    top_rank = (
        pl.when(pl.col("genus").is_not_null()).then(pl.lit("genus")).otherwise(pl.lit("species"))
    )
    pairs = species.with_columns(t_rank=top_rank, t=top).join(
        n.select(t_rank="rank", t="clade", unit="unit"), on=["t_rank", "t"]
    )
    root_n = sizes.filter(pl.col("rank") == "root")["N"][0]
    q = pl.col("n_root") / root_n
    pairs = pairs.join(
        n.filter(pl.col("rank") == "root").select("unit", n_root="n"), on="unit", how="left"
    )
    for r in RANKS:
        counts = n.filter(pl.col("rank") == r).select(pl.col("clade").alias(r), "unit", n_r="n")
        size = sizes.filter(pl.col("rank") == r).select(pl.col("clade").alias(r), N_r="N")
        pairs = pairs.join(counts, on=[r, "unit"], how="left").join(size, on=r, how="left")
        a = alpha[r]
        shrunk = (pl.col("n_r").fill_null(0) + a * q) / (pl.col("N_r") + a)
        q = pl.when(pl.col(r).is_null()).then(q).otherwise(shrunk)
        pairs = pairs.with_columns(q=q).drop("n_r", "N_r")
        q = pl.col("q")
    # mean content of carriers, species first, then up the ranks
    with_lineage = carried.join(lineage, on="genome")
    e = None
    for r in ("species", "genus", "family"):
        mean = with_lineage.drop_nulls(r).group_by(r, "unit").agg(pl.col("c").mean().alias(r[0]))
        pairs = pairs.join(mean, on=[r, "unit"], how="left")
        e = pl.col(r[0]) if e is None else pl.coalesce(e, pl.col(r[0]))
    root = carried.group_by("unit").agg(e_root=pl.col("c").mean())
    pairs = pairs.join(root, on="unit", how="left")
    assert e is not None
    return (
        pairs.select("species", "unit", q=pl.min_horizontal("q", 1.0), e=pl.coalesce(e, "e_root"))
        .filter(pl.col("q") >= Q_FLOOR)
        .sort("species", "unit")
    )


# --- Building the index ---


def _csr(
    pairs: pl.DataFrame, key: str, other: str, n: int, out: Path, prefix: str, offsets: str
) -> None:
    rows = pairs.sort(key, other)
    per = np.bincount(rows[key].to_numpy(), minlength=n)
    np.save(out / offsets, np.r_[0, np.cumsum(per)].astype(np.uint64))
    np.save(out / f"{prefix}{other}.npy", rows[other].to_numpy().astype(np.uint32))
    np.save(out / f"{prefix}q.npy", rows["q"].to_numpy().astype(np.float32))
    np.save(out / f"{prefix}e.npy", rows["e"].to_numpy().astype(np.float32))


def write_species_index(
    out: Path,
    index_dir: Path,
    n_units: int,
    carried: pl.DataFrame,
    genomes: pl.DataFrame,
    exclude: frozenset[str],
    meta: dict[str, object],
) -> dict[str, object]:
    """Aggregate (``genome``, ``unit``, ``c``) carriage rows of ``genomes`` (``genome``,
    ``name``, ``species`` or ``taxonomy``, optional ``taxonomy`` and ``completeness`` in
    [0, 1]) into a species index at ``out`` (module docstring); genomes named in
    ``exclude`` are left out of every count, and their carried units written to
    ``held_out.parquet``."""
    out.mkdir(parents=True, exist_ok=True)
    lineage_all = lineages(genomes)
    held = genomes.filter(pl.col("name").is_in(list(exclude)))
    carried.join(held.select("genome", "name"), on="genome").join(
        lineage_all.select("genome", "species"), on="genome"
    ).select("name", "species", "unit", "c").sort("name", "unit").write_parquet(
        out / "held_out.parquet"
    )
    kept = genomes.filter(~pl.col("name").is_in(list(exclude)))
    lineage = lineage_all.join(kept.select("genome"), on="genome", how="semi")
    carried = carried.join(kept.select("genome"), on="genome", how="semi")
    completeness = kept.select("genome", "completeness") if "completeness" in kept.columns else None
    n, sizes = _counts(carried.select("genome", "unit"), lineage, completeness)
    fits = {r: fit_alpha(carried.select("genome", "unit"), lineage, r) for r in RANKS}
    alpha = {r: a for r, (a, _) in fits.items()}
    pairs = prevalence(carried, lineage, n, sizes, alpha)
    first = lineage.unique("species", keep="first", maintain_order=True)
    taxonomy = (
        kept.select("genome", "taxonomy") if "taxonomy" in kept.columns
        else kept.select("genome", taxonomy=pl.lit(None, pl.String))
    )  # fmt: skip
    species = (
        first.join(taxonomy, on="genome", how="left")
        .join(
            sizes.filter(pl.col("rank") == "species").select(species="clade", completeness="N"),
            on="species",
        )  # fmt: skip
        .join(lineage.group_by("species").agg(genomes=pl.len()), on="species")
        .sort("species")
        .with_row_index("row")
    )
    # A species is named by s__ if it has one, else by its id (catalogue: the representative;
    # genome sets: genome:<name> -> <name>)
    species = species.with_columns(
        name=pl.when(pl.col("species").str.starts_with("genome:"))
        .then(pl.col("species").str.slice(7))
        .otherwise(pl.col("taxonomy").str.extract(r"s__([^;]+)$", 1).fill_null(pl.col("species")))
    )
    ids = species.select("species", row="row")
    pairs = pairs.join(ids, on="species").select(
        species=pl.col("row").cast(pl.UInt32), unit=pl.col("unit").cast(pl.UInt32), q="q", e="e"
    )
    sums = pairs.group_by("species").agg(units=pl.len(), content=(pl.col("q") * pl.col("e")).sum())
    table = (
        species.select(pl.col("row").cast(pl.UInt32).alias("species"), id="species", name="name",
                       taxonomy="taxonomy", genomes="genomes", completeness="completeness")
        .join(sums, on="species", how="left")
        .with_columns(pl.col("units").fill_null(0), pl.col("content").fill_null(0.0))
    )  # fmt: skip
    table.write_csv(out / "species.tsv", separator="\t")
    n_species = table.height
    _csr(pairs, "species", "unit", n_species, out, "s_", "species_offsets.npy")
    _csr(pairs, "unit", "species", n_units, out, "u_", "unit_offsets.npy")
    kept.join(lineage.select("genome", "species"), on="genome").select(
        "genome", "name", "species", *(["completeness"] if completeness is not None else [])
    ).write_parquet(out / "genomes.parquet")
    n.sort("rank", "clade", "unit").write_parquet(out / "carriage.parquet")
    sizes.sort("rank", "clade").write_parquet(out / "clades.parquet")
    lineage.write_parquet(out / "lineage.parquet")
    if (index_dir / "unit_pfam.parquet").exists():
        pl.read_parquet(index_dir / "unit_pfam.parquet").write_parquet(out / "unit_pfam.parquet")
    meta = meta | {
        "index_meta_sha256": _sha256(index_dir / "meta.json"),
        "units": n_units,
        "species": n_species,
        "genomes": kept.height,
        "excluded": held.height,
        "pairs": pairs.height,
        "alpha": alpha,
        "held_out_log_loss": {r: losses for r, (_, losses) in fits.items()},
    }
    (out / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    return meta


def _load_index(index_dir: Path) -> Index:
    index = Index.load(index_dir)
    if index.meta.get("hash") == "sourmash":
        raise ValueError("a species index needs the index's own hash (not a sourmash import)")
    return index


def species_index_from_genomes(
    index_dir: str | Path,
    genome_index: str | Path,
    out: str | Path,
    exclude: frozenset[str] = frozenset(),
) -> dict[str, object]:
    """A species index from an ``annotate-genomes`` genome index (its ``genome_best`` for
    carriage, its content for *e*); ``genomes.tsv`` may have ``completeness`` (percent)."""
    index_dir, genome_index = Path(index_dir), Path(genome_index)
    gmeta = json.loads((genome_index / "meta.json").read_text())
    if gmeta["index_meta_sha256"] != _sha256(index_dir / "meta.json"):
        raise ValueError(f"{genome_index} was annotated with another index than {index_dir}")
    genomes = pl.read_csv(genome_index / "genomes.tsv", separator="\t").cast({"genome": pl.UInt32})
    if "completeness" in genomes.columns:
        genomes = genomes.with_columns(pl.col("completeness").cast(pl.Float64) / 100)
    best = pl.read_parquet(genome_index / "genome_best.parquet").select("genome", "unit").unique()
    offsets = np.load(genome_index / "unit_offsets.npy", mmap_mode="r")
    units = np.unique(best["unit"].to_numpy()).astype(np.int64)
    starts = offsets[units].astype(np.int64)
    lengths = offsets[units + 1].astype(np.int64) - starts
    rows = _ranges(starts, lengths)
    content = pl.DataFrame(
        {
            "unit": np.repeat(units, lengths).astype(np.uint32),
            "genome": np.load(genome_index / "genome.npy", mmap_mode="r")[rows],
            "c": np.load(genome_index / "hits.npy", mmap_mode="r")[rows].astype(np.float64),
        }
    )
    carried = best.join(content, on=["genome", "unit"]).select("genome", "unit", "c")
    meta: dict[str, object] = {
        "source": "genomes",
        "genome_index_meta_sha256": _sha256(genome_index / "meta.json"),
    }
    return write_species_index(
        Path(out), index_dir, gmeta["units"], carried, genomes, exclude, meta
    )


def translate(cds: bytes) -> bytes:
    """A CDS's protein (table 11, frame 0), the terminal stop dropped; internal stops stay
    (k-mers spanning them are never hashed)."""
    protein: bytes = _core.translate_frames(cds.upper())[0]
    return protein[:-1] if protein.endswith(b"*") else protein


def _fasta_named(path: Path) -> Iterator[tuple[str, bytes]]:
    name, seq = None, bytearray()
    with open(path, "rb") as f:
        for line in f:
            if line.startswith(b">"):
                if name is not None:
                    yield name, bytes(seq)
                name, seq = line[1:].split()[0].decode(), bytearray()
            else:
                seq += line.strip()
    if name is not None:
        yield name, bytes(seq)


def _species_dir(catalogue: Path, rep: str) -> Path:
    return catalogue / "species_catalogue" / rep[:-2] / rep


def _pangenome(
    catalogue: Path, rep: str, genome_ids: dict[str, int]
) -> tuple[list[bytes], sparse.csr_matrix, np.ndarray]:
    """A species' family proteins, the family x genome presence matrix and the genome ids
    of its columns: Panaroo's pangenome if the species has one, else its representative's
    proteins, all carried by it."""
    root = _species_dir(catalogue, rep)
    rtab = root / "pan-genome" / "gene_presence_absence.Rtab"
    if not rtab.exists():
        proteins = list(fasta(root / "genome" / f"{rep}.faa"))
        presence = sparse.csr_matrix(np.ones((len(proteins), 1), dtype=np.float64))
        return proteins, presence, np.array([genome_ids[rep]])
    table = pl.read_csv(rtab, separator="\t")
    names = [c for c in table.columns[1:] if c in genome_ids]
    families = dict(_fasta_named(root / "pan-genome" / "pan-genome.fna"))
    table = table.filter(pl.col(table.columns[0]).is_in(list(families)))
    proteins = [translate(families[f]) for f in table[table.columns[0]]]
    presence = sparse.csr_matrix(table.select(names).to_numpy().astype(np.float64))
    return proteins, presence, np.array([genome_ids[n] for n in names])


def species_index_from_catalogue(
    index_dir: str | Path,
    catalogue: str | Path,
    out: str | Path,
    exclude: frozenset[str] = frozenset(),
    species: set[str] | None = None,
) -> dict[str, object]:
    """A species index from an MGnify genome catalogue directory (``CATALOGUE_METADATA``
    and ``species_catalogue/`` as on the FTP site), optionally only the species whose
    representatives are in ``species``."""
    index_dir, catalogue = Path(index_dir), Path(catalogue)
    index = _load_index(index_dir)
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
    parts = []
    # ponytail: every (genome, carried unit) row is held in memory (~12 B each; ~10 GB for
    # human-gut's 289k genomes); aggregate per species and fit α on a sample if a node runs short.
    for rep in genomes["species"].unique(maintain_order=True):
        proteins, presence, columns = _pangenome(catalogue, rep, genome_ids)
        batches = (proteins[i : i + PROTEINS_PER_BATCH]
                   for i in range(0, len(proteins), PROTEINS_PER_BATCH))  # fmt: skip
        hits, lengths = protein_hits(index, batches)
        links = hits.group_by("protein", "unit").agg(c=pl.len())
        best = best_units(index, hits, lengths)
        if best.height == 0:
            continue
        units, unit_col = np.unique(links["unit"].to_numpy(), return_inverse=True)
        f = len(proteins)
        content = sparse.csr_matrix(
            (links["c"].to_numpy().astype(np.float64),
             (links["protein"].to_numpy().astype(np.int64), unit_col)),
            shape=(f, len(units)),
        )  # fmt: skip
        is_best = sparse.csr_matrix(
            (np.ones(best.height),
             (best["protein"].to_numpy().astype(np.int64),
              np.searchsorted(units, best["unit"].to_numpy()))),
            shape=(f, len(units)),
        )  # fmt: skip
        carries = (is_best.T @ presence).tocoo()  # unit x genome: families with it as best
        c = np.asarray((content.T @ presence)[carries.row, carries.col]).ravel()
        parts.append(
            pl.DataFrame(
                {
                    "genome": columns[carries.col].astype(np.uint32),
                    "unit": units[carries.row].astype(np.uint32),
                    "c": c,
                }
            )
        )
    schema = {"genome": pl.UInt32, "unit": pl.UInt32, "c": pl.Float64}
    carried = pl.concat([pl.DataFrame(schema=schema), *parts])
    meta: dict[str, object] = {"source": "catalogue", "catalogue": str(catalogue)}
    return write_species_index(
        Path(out), index_dir, index.units.height, carried, genomes, exclude, meta
    )


class SpeciesIndex:
    """A species index written by :func:`write_species_index`, its pairs memory-mapped."""

    def __init__(self, path: str | Path) -> None:
        path = Path(path)
        self.path = path
        self.meta = json.loads((path / "meta.json").read_text())
        self.species = pl.read_csv(
            path / "species.tsv", separator="\t", schema_overrides={"taxonomy": pl.String}
        )

        def load(name: str) -> np.ndarray:
            array: np.ndarray = np.load(path / f"{name}.npy", mmap_mode="r")
            return array

        self.species_offsets, self.s_unit, self.s_q, self.s_e = map(
            load, ("species_offsets", "s_unit", "s_q", "s_e")
        )
        self.unit_offsets, self.u_species, self.u_q, self.u_e = map(
            load, ("unit_offsets", "u_species", "u_q", "u_e")
        )
        pfam = path / "unit_pfam.parquet"
        self.unit_pfam = pl.read_parquet(pfam) if pfam.exists() else None

    def check_index(self, index_dir: str | Path) -> None:
        """Raise if ``index_dir`` is not the index the species index was built with."""
        if _sha256(Path(index_dir) / "meta.json") != self.meta["index_meta_sha256"]:
            raise ValueError(f"{self.path} was built with another index than {index_dir}")

    def by_unit(self, units: np.ndarray) -> pl.DataFrame:
        """``unit``, ``species``, ``q``, ``e`` pairs of ``units`` (index unit ids)."""
        units = np.asarray(units, dtype=np.int64)
        if len(units) and units.max() >= self.meta["units"]:
            raise ValueError("unit ids beyond the species index: a profile of another index?")
        starts = self.unit_offsets[units].astype(np.int64)
        lengths = self.unit_offsets[units + 1].astype(np.int64) - starts
        rows = _ranges(starts, lengths)
        return pl.DataFrame(
            {
                "unit": np.repeat(units, lengths).astype(np.uint32),
                "species": self.u_species[rows],
                "q": self.u_q[rows].astype(np.float64),
                "e": self.u_e[rows].astype(np.float64),
            }
        )

    def by_species(self, species: np.ndarray) -> pl.DataFrame:
        """``species``, ``unit``, ``q``, ``e`` pairs of ``species`` (all their units)."""
        species = np.asarray(species, dtype=np.int64)
        starts = self.species_offsets[species].astype(np.int64)
        lengths = self.species_offsets[species + 1].astype(np.int64) - starts
        rows = _ranges(starts, lengths)
        return pl.DataFrame(
            {
                "species": np.repeat(species, lengths).astype(np.uint32),
                "unit": self.s_unit[rows],
                "q": self.s_q[rows].astype(np.float64),
                "e": self.s_e[rows].astype(np.float64),
            }
        )


# --- The fit (phase 11, step 8) ---

# ponytail: priors are guesses, to tune on the strain hold-out benchmark (step 10).
PRIOR_PRESENT: Final = 0.01  # ρ: prior probability that a screened species is present
DEPTH_SHAPE: Final = 1.0  # Gamma prior on a present species' depth (weak)
DEPTH_RATE: Final = 1e-3
BG_SHAPE: Final = 0.1  # Gamma prior on a unit's background rate, in hits (sparse)
BG_RATE: Final = 0.01
BG_PRIOR: Final = 0.05  # prior probability that a unit has background hits at all
DETECTED: Final = 0.5  # present_prob reported as detected
TINY: Final = 1e-12


def _log_pois(h: np.ndarray, m: np.ndarray) -> np.ndarray:
    out: np.ndarray = h * np.log(np.maximum(m, TINY)) - m - gammaln(h + 1)
    return out


def _log_nb(h: np.ndarray) -> np.ndarray:
    """Log marginal of h hits from a background rate under its Gamma prior alone."""
    a, b = BG_SHAPE, BG_RATE
    out: np.ndarray = (
        gammaln(h + a) - gammaln(a) - gammaln(h + 1) + a * np.log(b / (b + 1)) - h * np.log1p(b)
    )
    return out


def _explain(h: np.ndarray, m: np.ndarray, nb: np.ndarray) -> np.ndarray:
    """Log probability of h hits from sources at mean m, the background free to add hits:
    either m explains them, or the background does and the sources gave ~none."""
    out: np.ndarray = np.maximum(_log_pois(h, m), nb - m)
    return out


def fit_species(
    prof: pl.DataFrame,
    si: SpeciesIndex,
    *,
    min_containment: float = MIN_CONTAINMENT,
    min_units: int = MIN_UNITS,
    prior: float = PRIOR_PRESENT,
    tol: float = 1e-6,
    max_iter: int = 2000,
) -> tuple[pl.DataFrame, pl.DataFrame, dict[str, float]]:
    """Species presence and depth, fitted jointly with each species' carriage of its units in
    the sample (plan: Species model), by variational EM.

    h_u ~ Poisson(β_u + Σ_s b_s z_{s,u} λ_s e_{s,u}) over the profile's raw unit ``hits``:
    presence b_s ~ Bernoulli(``prior``), depth λ_s ~ Gamma(``DEPTH_SHAPE``, ``DEPTH_RATE``)
    if present, carriage z_{s,u} ~ Bernoulli(q_{s,u}), background β_u = 0, or with
    probability ``BG_PRIOR`` ~ Gamma(``BG_SHAPE``, ``BG_RATE``). Species are screened as
    genome mode's are, on prevalence-weighted content: Σ_{u hit} q e / Σ_u q e
    ≥ ``min_containment`` on ≥ ``min_units`` hit units.

    Each round, with the other sources at their means (μ₋ₛ):

    - carriage r = P(z = 1 | h)
      = q Pois(h; μ₋ₛ + λe) / [q Pois(h; μ₋ₛ + λe) + (1 - q) Pois(h; μ₋ₛ)];
    - the background is on with its posterior given the species' mean M_u (spike against
      slab, :func:`_explain`) and then takes the excess, max(h - M_u, 0); the rest of the
      hits are allocated ∝ p_s r exp(E log λ_s) e;
    - depth given presence: Gamma(a₀ + Σ_u allocated hits as if present, b₀ + Σ_u r e)
      over every unit of the species, hit or not;
    - presence: logit p = logit ρ + Σ_u log[q ℓ(μ₋ₛ + λe) + (1 - q) ℓ(μ₋ₛ)] - log ℓ(μ₋ₛ),
      plus the depth's Laplace term; ℓ lets the background explain a unit's hits instead
      (:func:`_explain`), so hits a species does not need cost it nothing and hits it
      cannot explain earn it nothing.

    Returns per screened species (``species``, ``present_prob``, ``depth``: the mode of λ's
    Gamma posterior given presence less the prior's shape, i.e. allocated hits / Σ r e,
    ``shape``, ``rate``, ``units_hit``, ``containment``), per (species, unit) pair
    (``q``, ``e``, ``hits``, ``carriage_prob``, ``expected_hits`` = r λ e,
    ``hits_assigned``: the allocation, weighted by presence; the rest of a unit's hits is
    the background's); and a report.
    """
    hit = prof.filter(pl.col("hits") > 0).select(
        pl.col("unit").cast(pl.UInt32), pl.col("hits").cast(pl.Float64)
    )
    content = si.species.select(pl.col("species").cast(pl.UInt32), total="content")
    screened = (
        si.by_unit(np.sort(hit["unit"].to_numpy()))
        .group_by("species")
        .agg(units_hit=pl.len().cast(pl.UInt32), on_hit=(pl.col("q") * pl.col("e")).sum())
        .join(content, on="species")
        .with_columns(containment=pl.col("on_hit") / pl.col("total"))
        .filter(pl.col("containment") >= min_containment, pl.col("units_hit") >= min_units)
        .select("species", "units_hit", "containment")
        .sort("species")
    )
    report: dict[str, float] = {"screened_species": screened.height}
    pairs = (
        si.by_species(screened["species"].to_numpy())
        .join(hit, on="unit", how="left")
        .with_columns(pl.col("hits").fill_null(0.0))
    )
    n_s = screened.height
    s = np.searchsorted(screened["species"].to_numpy(), pairs["species"].to_numpy())
    unit_ids, u = np.unique(pairs["unit"].to_numpy(), return_inverse=True)
    n_u = len(unit_ids)
    q = np.clip(pairs["q"].to_numpy(), CLIP, 1 - CLIP)
    lq, l1q = np.log(q), np.log1p(-q)
    e, h = pairs["e"].to_numpy(), pairs["hits"].to_numpy()
    h_u = np.zeros(n_u)
    h_u[u] = h
    nb = _log_nb(h)
    # Start: depth from prevalence-weighted hits, every screened species present.
    rate = DEPTH_RATE + np.bincount(s, q * e, n_s)
    lam = _div(np.bincount(s, q * h, n_s), rate - DEPTH_RATE)
    shape = DEPTH_SHAPE + lam * (rate - DEPTH_RATE)
    p, r = np.ones(n_s), q.copy()
    y_bg = np.zeros(n_u)
    logit_prior = np.log(prior) - np.log1p(-prior)
    logit_bg = np.log(BG_PRIOR) - np.log1p(-BG_PRIOR)
    nb_u = _log_nb(h_u)
    for iteration in range(1, max_iter + 1):  # noqa: B007  (reported after the loop)
        mean_lam, log_lam = shape / rate, digamma(shape) - np.log(rate)
        mu = mean_lam[s] * e
        m = p[s] * r * mu
        others = np.maximum((np.bincount(u, m, n_u) + y_bg)[u] - m, TINY)
        r = expit(lq - l1q - mu + h * np.log1p(mu / others))
        # Background: on with its posterior given the species' mean, taking the excess.
        m = p[s] * r * mu
        mean_u = np.bincount(u, m, n_u)
        on = expit(logit_bg + _explain(h_u, mean_u, nb_u) - _log_pois(h_u, mean_u))
        y_bg = on * np.maximum(h_u - mean_u, 0.0)
        w = r * np.exp(log_lam[s]) * e
        total = np.bincount(u, p[s] * w, n_u)
        left = (h_u - y_bg)[u]
        if_present = _div(left * w, w + total[u] - p[s] * w)
        new_shape = DEPTH_SHAPE + np.bincount(s, if_present, n_s)
        rate = DEPTH_RATE + np.bincount(s, r * e, n_s)
        lam_hat = new_shape / rate
        others = np.maximum((mean_u + y_bg)[u] - m, TINY)
        alone = _explain(h, others, nb)
        with_s = np.logaddexp(lq + _explain(h, others + lam_hat[s] * e, nb), l1q + alone)
        occam = gamma.logpdf(lam_hat, DEPTH_SHAPE, scale=1 / DEPTH_RATE) + 0.5 * np.log(
            2 * np.pi * new_shape / rate**2
        )
        new_p = expit(logit_prior + np.bincount(s, with_s - alone, n_s) + occam)
        change = max(
            float(np.max(np.abs(lam_hat - mean_lam) / np.maximum(lam_hat, TINY), initial=0.0)),
            float(np.max(np.abs(new_p - p), initial=0.0)),
        )
        shape, p = new_shape, new_p
        if change < tol:
            break
    report |= {"iterations": iteration, "max_change": change}
    mean_lam, log_lam = shape / rate, digamma(shape) - np.log(rate)
    w = r * np.exp(log_lam[s]) * e
    total = np.bincount(u, p[s] * w, n_u)
    # depth reported at the posterior mode (the MLE at DEPTH_SHAPE 1): the mean adds a hit
    # (DEPTH_SHAPE), which biases species with few hits upwards
    depth = np.maximum(shape - DEPTH_SHAPE, 0.0) / rate
    species = screened.with_columns(present_prob=p, depth=depth, shape=shape, rate=rate)
    pairs = pairs.with_columns(
        carriage_prob=r,
        expected_hits=r * depth[s] * e,
        hits_assigned=_div((h_u - y_bg)[u] * p[s] * w, total[u]),
    )
    return species, pairs, report


# --- Outputs (phase 11, step 9) ---

CORE: Final = 0.9  # prevalence from which a unit counts as core in the species table


def _intervals(species: pl.DataFrame, pairs: pl.DataFrame) -> tuple[pl.DataFrame, float]:
    """95% ``depth_lo``/``depth_hi`` from each species' Gamma posterior, its variance
    inflated by φ = max(1, Pearson χ² / (units - species)) over the units of detected
    species (mean-field posteriors are too narrow; unit hits are overdispersed)."""
    detected = species.filter(pl.col("present_prob") >= DETECTED).select("species")
    units = pairs.join(detected, on="species", how="semi")["unit"].unique().implode()
    per_unit = (
        pairs.filter(pl.col("unit").is_in(units))
        .join(species.select("species", "present_prob"), on="species")
        .group_by("unit")
        .agg(pl.col("hits").first(), (pl.col("present_prob") * pl.col("expected_hits")).sum())
        .filter(pl.col("present_prob") > 0)
    )
    dof = max(per_unit.height - detected.height, 1)
    chi2 = ((per_unit["hits"] - per_unit["present_prob"]) ** 2 / per_unit["present_prob"]).sum()
    phi = max(1.0, float(chi2) / dof)
    a, b = species["shape"].to_numpy() / phi, species["rate"].to_numpy() / phi
    return species.with_columns(
        depth_lo=gamma.ppf(0.025, a, scale=1 / b), depth_hi=gamma.ppf(0.975, a, scale=1 / b)
    ), phi


def unit_presence(
    prof: pl.DataFrame, pairs: pl.DataFrame, species: pl.DataFrame, unit_pfam: pl.DataFrame | None
) -> tuple[pl.DataFrame, pl.DataFrame | None]:
    """Unit and Pfam presence updated by the species fit, in kfp-prior's format.

    Per unit, carriage by some species: P = 1 - Π_s (1 - p_s r_{s,u}); ``prior`` the same
    with q for r, ``expected_hits`` Σ_s p_s r λ e. Hit units: 1 - (1 - ``present_prob``)
    (1 - P), so a unit no screened species carries keeps ``present_prob`` (taken as 1 if
    the profile has none). Zero-hit units with a prior ≥ ``Q_FLOOR``: P, ``imputed`` if
    ≥ 0.5. ``sibling_hit``: one of the unit's Pfams is observed through another unit. Per
    Pfam, ``present_prob_observed`` and ``present_prob_updated`` (units independent) and
    ``imputed``."""
    weighted = pairs.join(species.select("species", "present_prob"), on="species")
    per_unit = weighted.group_by("unit").agg(
        carried=1 - (1 - pl.col("present_prob") * pl.col("carriage_prob")).product(),
        prior=1 - (1 - pl.col("present_prob") * pl.col("q")).product(),
        expected_hits=(pl.col("present_prob") * pl.col("expected_hits")).sum(),
    )
    observed = prof.filter(pl.col("hits") > 0).select(
        pl.col("unit").cast(pl.UInt32),
        pl.col("hits").cast(pl.Float64),
        present_prob=pl.col("present_prob").cast(pl.Float64)
        if "present_prob" in prof.columns
        else pl.lit(1.0),
    )
    hit = observed.join(per_unit, on="unit", how="left").with_columns(
        pl.col("prior", "expected_hits", "carried").fill_null(0.0)
    )
    hit = hit.with_columns(
        present_prob_updated=1 - (1 - pl.col("present_prob")) * (1 - pl.col("carried"))
    ).drop("carried")
    zero = (
        per_unit.join(observed, on="unit", how="anti")
        .filter(pl.col("prior") >= Q_FLOOR)
        .select(
            "unit",
            hits=pl.lit(0.0),
            present_prob=pl.lit(None, pl.Float64),
            prior="prior",
            expected_hits="expected_hits",
            present_prob_updated="carried",
        )  # fmt: skip
    )
    presence = pl.concat([hit, zero], how="diagonal_relaxed").with_columns(
        imputed=(pl.col("hits") == 0) & (pl.col("present_prob_updated") >= 0.5)
    )
    if unit_pfam is None:
        return presence.with_columns(sibling_hit=pl.lit(None, pl.Boolean)).sort("unit"), None
    labels = unit_pfam.select(pl.col("unit").cast(pl.UInt32), "pfam_accession")
    seen = labels.join(
        presence.filter((pl.col("hits") > 0) & (pl.col("present_prob") >= 0.5)), on="unit"
    )["pfam_accession"].unique()
    sibling = labels.filter(pl.col("pfam_accession").is_in(seen.implode()))["unit"].unique()
    presence = presence.with_columns(
        sibling_hit=(pl.col("hits") == 0) & pl.col("unit").is_in(sibling.implode())
    ).sort("unit")
    by_pfam = (
        presence.join(labels, on="unit")
        .group_by("pfam_accession")
        .agg(
            present_prob_observed=1 - (1 - pl.col("present_prob").fill_null(0.0)).product(),
            present_prob_updated=1 - (1 - pl.col("present_prob_updated")).product(),
        )
        .with_columns(
            imputed=(pl.col("present_prob_observed") < 0.5)
            & (pl.col("present_prob_updated") >= 0.5)
        )
        .sort("pfam_accession")
    )
    return presence, by_pfam


def function_species(
    prof: pl.DataFrame, pairs: pl.DataFrame, species: pl.DataFrame, si: SpeciesIndex
) -> pl.DataFrame:
    """The function x species table: each unit's split hits (``hits_em``) attributed to
    the detected species by the fit's allocation (``hits_assigned`` / ``hits``); the rest
    (background, undetected species, units no screened species carries) is
    ``unclassified``. Species roll up to genus and family by taxonomy (a missing rank takes
    the lowest taxon above it); functions as genome mode's table (Pfam, else unit name).
    Long format: ``function``, ``rank`` (``total``, ``family``, ``genus``, ``species``),
    ``taxon``, ``hits_em``."""
    split = prof.filter(pl.col("hits_em") > 0).select(
        pl.col("unit").cast(pl.UInt32), "hits_em", *(["name"] if "name" in prof.columns else [])
    )
    detected = species.filter(pl.col("present_prob") >= DETECTED).select("species")
    shares = (
        pairs.join(detected, on="species", how="semi")
        .filter(pl.col("hits") > 0)
        .select("unit", "species", share=pl.col("hits_assigned") / pl.col("hits"))
    )
    resp = split.join(shares, on="unit").select(
        "unit", genome="species", hits_em=pl.col("hits_em") * pl.col("share")
    )
    names = si.species.select(
        pl.col("species").cast(pl.UInt32).alias("genome"), "name", "taxonomy",
        group=pl.col("species").cast(pl.UInt32),
    )  # fmt: skip
    labels = (
        taxon_labels(names.join(detected.select(genome="species"), on="genome", how="semi"))
        .filter(pl.col("rank") != "species")
        .with_columns(rank=pl.col("rank").replace("genome", "species"))
    )
    classified = resp.join(labels, on="genome").select("unit", "rank", "taxon", "hits_em")
    ranks = pl.DataFrame({"rank": ["family", "genus", "species"]})
    done = classified.group_by("unit", "rank").agg(done=pl.col("hits_em").sum())
    unclassified = (
        split.select("unit", "hits_em")
        .join(ranks, how="cross")
        .join(done, on=["unit", "rank"], how="left")
        .with_columns(rest=pl.col("hits_em") - pl.col("done").fill_null(0.0))
        .filter(pl.col("rest") > 1e-9 * pl.col("hits_em"))  # rounding, not hits
        .select("unit", "rank", taxon=pl.lit("unclassified"), hits_em="rest")
    )
    total = split.select("unit", rank=pl.lit("total"), taxon=pl.lit(""), hits_em="hits_em")
    return (
        pl.concat([total, classified, unclassified])
        .join(unit_functions(si.unit_pfam, split), on="unit")
        .group_by("function", "rank", "taxon")
        .agg(pl.col("hits_em").sum())
        .sort("function", "rank", "taxon")
    )


def species_profile(
    prof: pl.DataFrame,
    si: SpeciesIndex,
    *,
    min_containment: float = MIN_CONTAINMENT,
    min_units: int = MIN_UNITS,
    prior: float = PRIOR_PRESENT,
) -> dict[str, object]:
    """:func:`fit_species` as reports: ``species`` (detected, ``present_prob`` ≥
    ``DETECTED``: ``id``, ``name``, ``taxonomy``, ``present_prob``, ``depth`` with its
    interval, ``relative_abundance``, ``units_hit``, ``containment``, ``core_hit`` (share of
    units with q ≥ ``CORE`` hit), ``accessory_called`` (units with q < ``CORE`` and
    carriage ≥ 0.5)); ``units`` (their pairs); ``presence`` and ``pfam_presence``
    (:func:`unit_presence`); ``function_species`` (None without ``hits_em``); and
    ``summary``: raw ``hits``, ``explained_fraction`` (hits assigned to detected species
    over all), ``genome_equivalents`` (Σ depth of the detected), counts."""
    species, pairs, report = fit_species(
        prof, si, min_containment=min_containment, min_units=min_units, prior=prior
    )
    species, phi = _intervals(species, pairs)
    counts = pairs.group_by("species").agg(
        core_hit=((pl.col("q") >= CORE) & (pl.col("hits") > 0)).sum()
        / (pl.col("q") >= CORE).sum().clip(1),
        accessory_called=((pl.col("q") < CORE) & (pl.col("carriage_prob") >= 0.5)).sum(),
    )
    names = si.species.select(pl.col("species").cast(pl.UInt32), "id", "name", "taxonomy")
    detected = (
        species.filter(pl.col("present_prob") >= DETECTED)
        .join(counts, on="species")
        .join(names, on="species")
        .with_columns(relative_abundance=pl.col("depth") / pl.col("depth").sum())
        .select(
            "species",
            "id",
            "name",
            "taxonomy",
            "present_prob",
            "depth",
            "depth_lo",
            "depth_hi",
            "relative_abundance",
            "units_hit",
            "containment",
            "core_hit",
            "accessory_called",
        )  # fmt: skip
        .sort("depth", descending=True)
    )
    units = pairs.join(detected.select("species"), on="species", how="semi").select(
        "species", "unit", prevalence="q", e="e", hits="hits", carriage_prob="carriage_prob",
        expected_hits="expected_hits", hits_assigned="hits_assigned",
    )  # fmt: skip
    presence, by_pfam = unit_presence(prof, pairs, species, si.unit_pfam)
    hits = float(prof.filter(pl.col("hits") > 0)["hits"].sum())
    summary = report | {
        "hits": hits,
        "explained_fraction": float(units["hits_assigned"].sum()) / hits if hits else 0.0,
        "genome_equivalents": float(detected["depth"].sum()),
        "species_detected": detected.height,
        "phi": phi,
        "imputed_units": int(presence["imputed"].sum()),
    }
    return {
        "species": detected,
        "units": units,
        "presence": presence,
        "pfam_presence": by_pfam,
        "function_species": function_species(prof, pairs, species, si)
        if "hits_em" in prof.columns
        else None,
        "summary": summary,
    }
