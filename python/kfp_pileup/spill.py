"""Target columns and events spilled to Parquet by reference sequence (plan, phase 13,
``CALL``'s memory, option 5).

``call``'s crude pass writes each batch's rows (hit, seq, col, qpos, residue, kind) to the
partition of their sequence; steps 6b and 8-9 then read one partition at a time, so they
hold one partition's rows instead of every target alignment's. All of a sequence's rows sit
in one partition, so results computed per sequence are the same as from the whole table.
``kind`` 0 is an aligned column (``residue`` its read byte, ``qpos`` < 0 a deletion); any
other ``kind`` is an event (insertion or frameshift byte, ``qpos`` its read-position key).
"""

from collections.abc import Iterator
from pathlib import Path
from typing import Final

import numpy as np
import polars as pl
import pyarrow.parquet as pq  # type: ignore[import-untyped]

PARTS: Final = 64  # partitions, by target length: one partition's rows held at a time
SCHEMA: Final = {
    "hit": pl.UInt32,
    "seq": pl.UInt32,
    "col": pl.Int32,
    "qpos": pl.Int64,
    "residue": pl.UInt8,
    "kind": pl.UInt8,
}


def partitions(length: np.ndarray, parts: int = PARTS) -> np.ndarray:
    """Each sequence's partition, contiguous runs of sequences holding similar shares of
    ``length`` (0 for sequences that get no rows)."""
    total = max(int(length.sum()), 1)
    before = np.cumsum(length) - length
    out: np.ndarray = np.minimum(before * parts // total, parts - 1)
    return out


class Spill:
    """Rows written by sequence partition (``part_of``, per sequence), read back one
    partition at a time."""

    def __init__(self, root: Path, part_of: np.ndarray) -> None:
        self.root, self.part_of = root, part_of
        self._writers: dict[int, pq.ParquetWriter] = {}
        self.rows = 0

    def write(self, rows: pl.DataFrame) -> None:
        rows = rows.select(pl.col(c).cast(t) for c, t in SCHEMA.items())
        part = pl.Series("part", self.part_of[rows["seq"].to_numpy()])
        for (p,), df in rows.with_columns(part).partition_by("part", as_dict=True).items():
            table = df.drop("part").to_arrow()
            if p not in self._writers:
                self._writers[p] = pq.ParquetWriter(self._path(p), table.schema)
            self._writers[p].write_table(table)
        self.rows += rows.height

    def close(self) -> None:
        for w in self._writers.values():
            w.close()

    @property
    def files(self) -> list[Path]:
        return [self._path(p) for p in sorted(self._writers)]

    def _path(self, p: int) -> Path:
        return self.root / f"part{p:04d}.parquet"

    def parts(self) -> Iterator[tuple[int, pl.DataFrame, pl.DataFrame]]:
        """(partition, columns: hit, seq, col, qpos, residue; events: hit, seq, col, kind,
        qpos) of each written partition in turn."""
        for p in sorted(self._writers):
            rows = pl.read_parquet(self._path(p))
            yield (
                p,
                rows.filter(pl.col("kind") == 0).drop("kind"),
                rows.filter(pl.col("kind") != 0).select("hit", "seq", "col", "kind", "qpos"),
            )
