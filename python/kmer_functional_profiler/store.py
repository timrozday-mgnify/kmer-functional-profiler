"""Target sequences by unit (plan, phase 13, step 1): an enhanced output of the index.

The store is one or more members tables (``protein_id``, ``cluster_rep``, ``sequence``)
sorted by ``cluster_rep`` with small row groups, so a fetch reads only the row groups whose
statistics bracket the wanted units: the index's own ``targets.parquet`` when it was built
from a FASTA, or the partitioned build's MERGE output (``members*.parquet``) for MGnify.

:func:`fetch` picks each unit's representative plus up to ``members - 1`` other distinct
member sequences, spread over the cluster's length range; :func:`write_fasta` writes them
with the member id as the first header word, which DIAMOND (``makedb``, with
``--no-parse-seqids`` for ids holding ``|``) and bwa-mem2 (``index``) both take as the
reference name.
"""

from collections.abc import Sequence
from pathlib import Path
from typing import Final

import polars as pl

from kmer_functional_profiler.index import TARGETS

MEMBERS: Final = 3  # sequences per unit: the representative and two others


def store_paths(index_dir: str | Path, store: Sequence[str | Path] = ()) -> list[Path]:
    """The members tables to fetch from: ``store`` if given, else the index's own targets."""
    if store:
        return [Path(p) for p in store]
    own = Path(index_dir) / TARGETS
    if not own.exists():
        raise ValueError(
            f"{index_dir} was not built from a FASTA (no {TARGETS}): pass the members "
            "tables it was built from as the store"
        )
    return [own]


def fetch(
    index_dir: str | Path,
    units: Sequence[int],
    members: int = MEMBERS,
    store: Sequence[str | Path] = (),
) -> pl.DataFrame:
    """Sequences of ``units`` (the index's unit ids): ``unit``, ``cluster_rep``,
    ``protein_id``, ``rep`` and ``sequence``, at most ``members`` distinct sequences per
    unit, the representative first.

    The representative is the member whose id is the ``cluster_rep`` (MGnify), else the
    longest (FASTA targets grouped by name). The others are taken evenly along the rest
    sorted by length.
    """
    # ponytail: spread by length, a cheap proxy for diversity; pick by k-mer distance
    # (farthest point) if reference bias shows up in the pileups
    if members < 1:
        raise ValueError("members must be >= 1")
    wanted = pl.read_parquet(Path(index_dir) / "units.parquet", columns=["unit", "cluster_rep"])
    wanted = wanted.filter(pl.col("unit").is_in(list(units)))
    if wanted.height < len(set(units)):
        raise ValueError("unit ids not in the index")
    reps = wanted["cluster_rep"]
    found = (
        pl.scan_parquet([str(p) for p in store_paths(index_dir, store)])
        .select("protein_id", "cluster_rep", "sequence")
        .filter(pl.col("cluster_rep").is_in(reps.implode()))
        .collect()
    )
    k = members - 1
    ranked = (
        found.with_columns(
            rep=pl.col("protein_id").cast(pl.String) == pl.col("cluster_rep").cast(pl.String),
            length=pl.col("sequence").str.len_bytes(),
        )
        .sort(["cluster_rep", "rep", "length", "protein_id"], descending=[False, True, True, False])
        .unique(["cluster_rep", "sequence"], keep="first", maintain_order=True)
        .with_columns(
            pos=pl.int_range(pl.len()).over("cluster_rep"),
            m=(pl.len() - 1).over("cluster_rep"),
        )
    )
    # Of the m others (longest first), keep min(k, m) evenly spaced, both ends included: the
    # first of each run of floor(i (picks - 1) / (m - 1)) over their index i.
    picks = pl.min_horizontal(pl.lit(k), pl.col("m"))
    i, step, span = pl.col("pos") - 1, picks - 1, (pl.col("m") - 1).clip(1)
    first = (i == 0) | ((i * step) // span != ((i - 1) * step) // span)
    chosen = ranked.filter((pl.col("pos") == 0) | ((i >= 0) & (picks > 0) & first))
    return (
        chosen.with_columns(rep=pl.col("pos") == 0)
        .join(wanted, on="cluster_rep")
        .select("unit", "cluster_rep", "protein_id", "rep", "sequence")
        .sort("unit", "rep", "protein_id", descending=[False, True, False])
    )


def write_fasta(sequences: pl.DataFrame, path: str | Path, ids: str = "protein_id") -> None:
    """``sequences`` as FASTA: ``ids`` as the first header word, then ``unit=`` and
    ``cluster_rep=``, one line per sequence."""
    with open(path, "w") as f:
        for name, unit, rep, seq in sequences.select(
            ids, "unit", "cluster_rep", "sequence"
        ).iter_rows():
            f.write(f">{name} unit={unit} cluster_rep={rep}\n{seq}\n")
