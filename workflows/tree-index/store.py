"""Tree index: a leaf-ordered sequence store (docs/tree-index-plan.md, Progress log: a
leaf-ordered store). One row per leaf position, row = position: the leaf's longest protein
(ties: the first; the representative ``verify.py`` and ``tier2.py`` align against), its
``protein_id`` and ``cluster_rep``. Reading a sample's candidates then touches only their
rows, not the whole members store: candidates come in runs along the tree's leaf order, so
few row groups hold them. It is a reordered copy of the MGnify protein store, published beside
the index, not part of it.

Built per partition part (``build.nf``'s BUILD outputs: each part's positions are one
contiguous range of the merged index, from ``merge.py``'s ``parts.json``), so no step holds the
whole release. Each part writes ``<name>.part.json`` (lo, hi, format, row group size) and, by
``--format``:

- ``parquet``: ``<name>.parquet`` (sequence, protein_id, cluster_rep; ``--row-group`` rows per
  group, zstd), read a row group at a time;
- ``flat`` (the default in ``build.nf``): ``<name>.bin`` (the sequences concatenated),
  ``<name>.offsets.npy`` (u64, a byte offset per row and the end) and ``<name>.ids.npy`` (i64
  protein_id, cluster_rep per row, -1 for none), memory-mapped and sliced: a read costs the
  candidates' own bytes wherever they fall. Parquet reads whole row groups, as cheap only when
  candidates cluster (2,000 of 3 M leaves: 0.04 s in runs of 20, 0.32 s scattered, as much as a
  full scan; flat 0.004-0.007 s), and is ~0.55x the size.

:class:`LeafStore` reads any number of parts from one directory.

Subcommands: ``build`` (one part), ``bench`` (reads a set of positions from a store, timed).
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np
import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq
from numpy.typing import NDArray

FORMATS = ("parquet", "flat")


def longest(proteins: pl.LazyFrame, members: list[str]) -> pl.DataFrame:
    """(pos, sequence, protein_id, cluster_rep): each position's longest protein, by pos."""
    return (
        pl.scan_parquet(members)
        .select("protein_id", "sequence")
        .join(proteins, on="protein_id")
        .with_columns(n=pl.col("sequence").str.len_bytes())
        .sort("pos", "n", "protein_id", descending=[False, True, False])
        .unique("pos", keep="first", maintain_order=True)
        .select("pos", "sequence", "protein_id", "cluster_rep")
        .collect()
    )


def build(
    proteins: Path,
    members: list[str],
    offset: int,
    out: Path,
    name: str,
    fmt: str,
    row_group: int,
) -> dict[str, int | str]:
    """One part's store: ``proteins`` (protein_id, cluster_rep, pos: the part's own
    positions, -1 for none) shifted by ``offset`` to the merged index's, every position from
    the part's least to its greatest one row (an empty sequence where none is held)."""
    prot = pl.scan_parquet(proteins).filter(pl.col("pos") >= 0).with_columns(pl.col("pos") + offset)
    got = longest(prot.select("protein_id", "cluster_rep", "pos"), members)
    lo, hi = int(got["pos"].min()), int(got["pos"].max()) + 1
    rows = (
        pl.DataFrame({"pos": np.arange(lo, hi, dtype=np.int64)})
        .join(got, on="pos", how="left")
        .with_columns(pl.col("sequence").fill_null(""))
    )
    out.mkdir(parents=True, exist_ok=True)
    ids = rows.select("protein_id", "cluster_rep")
    if fmt == "parquet":
        pq.write_table(
            rows.select("sequence", "protein_id", "cluster_rep").to_arrow(),
            out / f"{name}.parquet", row_group_size=row_group, compression="zstd",
        )  # fmt: skip
    else:
        blob = "".join(rows["sequence"].to_list()).encode()
        lengths = rows["sequence"].str.len_bytes().to_numpy().astype(np.uint64)
        offsets = np.zeros(len(lengths) + 1, dtype=np.uint64)
        np.cumsum(lengths, out=offsets[1:])
        (out / f"{name}.bin").write_bytes(blob)
        np.save(out / f"{name}.offsets.npy", offsets)
        # protein_id, cluster_rep per row (-1: none), memory-mapped like the offsets
        np.save(out / f"{name}.ids.npy", ids.fill_null(-1).to_numpy().astype(np.int64))
    part = {"name": name, "lo": lo, "hi": hi, "format": fmt, "row_group": row_group,
            "missing": int(rows["protein_id"].null_count())}  # fmt: skip
    (out / f"{name}.part.json").write_text(json.dumps(part, indent=1) + "\n")
    return part


class LeafStore:
    """The parts in one directory (``*.part.json``), read by position."""

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        self.parts = sorted(
            (json.loads(f.read_text()) for f in self.path.glob("*.part.json")),
            key=lambda p: p["lo"],
        )
        if not self.parts:
            raise ValueError(f"no *.part.json in {self.path}")
        self.lo = np.array([p["lo"] for p in self.parts], dtype=np.int64)
        self._mm: dict[str, tuple[np.memmap, NDArray[np.uint64]]] = {}

    def _groups(
        self, part: dict, rows: NDArray[np.int64], file: Path, columns: list[str]
    ) -> tuple[pa.Table, NDArray[np.int64]]:
        """Columns of the row groups holding ``rows``, and each row's place in them."""
        rg = part["row_group"]
        groups = np.unique(rows // rg)
        table = pq.ParquetFile(file).read_row_groups(groups.tolist(), columns=columns)
        # rows of group g sit after the groups read before it (all full but the file's last)
        start = np.searchsorted(groups, rows // rg) * rg
        return table, start + rows % rg

    def get(self, positions: NDArray[np.int64], ids: bool = False) -> pl.DataFrame:
        """(pos, sequence[, protein_id, cluster_rep]) for each distinct position held, by pos."""
        positions = np.unique(np.asarray(positions, dtype=np.int64))
        part_of = np.searchsorted(self.lo, positions, side="right") - 1
        out = []
        for i in np.unique(part_of[part_of >= 0]).tolist():
            p = self.parts[i]
            sel = positions[(part_of == i) & (positions < p["hi"])]
            rows = sel - p["lo"]
            if p["format"] == "parquet":
                cols = ["sequence", "protein_id", "cluster_rep"] if ids else ["sequence"]
                table, at = self._groups(p, rows, self.path / f"{p['name']}.parquet", cols)
                got = pl.DataFrame(pl.from_arrow(table.take(pa.array(at))))
            else:
                if p["name"] not in self._mm:
                    self._mm[p["name"]] = (
                        np.memmap(self.path / f"{p['name']}.bin", dtype=np.uint8, mode="r"),
                        np.load(self.path / f"{p['name']}.offsets.npy", mmap_mode="r"),
                    )
                blob, off = self._mm[p["name"]]
                a, b = off[rows], off[rows + 1]
                got = pl.DataFrame(
                    {"sequence": [blob[x:y].tobytes().decode() for x, y in zip(a, b, strict=True)]}
                )
                if ids:
                    both = np.load(self.path / f"{p['name']}.ids.npy", mmap_mode="r")[rows]
                    got = got.with_columns(protein_id=both[:, 0], cluster_rep=both[:, 1])
            out.append(got.insert_column(0, pl.Series("pos", sel)))
        if not out:
            return pl.DataFrame(schema={"pos": pl.Int64, "sequence": pl.String})
        return pl.concat(out).filter(pl.col("sequence") != "")


def main() -> None:
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("build")
    p.add_argument("--proteins", required=True, help="a part's (or the index's) proteins.parquet")
    p.add_argument("--members", nargs="+", required=True)
    p.add_argument("--offset", type=int, default=0, help="the part's first merged position")
    p.add_argument("--parts", help="merge.py's parts.json: the offset of --name's part")
    p.add_argument("--name", default="part0000")
    p.add_argument("--format", choices=FORMATS, default="parquet")
    p.add_argument("--row-group", type=int, default=1024)
    p.add_argument("--out", required=True)
    p = sub.add_parser("bench")
    p.add_argument("--store", required=True)
    p.add_argument("--positions", required=True, help=".npy of positions to read")
    a = parser.parse_args()
    if a.command == "build":
        offset = a.offset
        if a.parts:
            offset = {x["name"]: x["offset"] for x in json.loads(Path(a.parts).read_text())}[a.name]
        t0 = time.perf_counter()
        got = build(Path(a.proteins), a.members, offset, Path(a.out), a.name, a.format, a.row_group)
        print(json.dumps(got | {"build_s": time.perf_counter() - t0}))
    else:
        t0 = time.perf_counter()
        got = LeafStore(a.store).get(np.load(a.positions))
        took = time.perf_counter() - t0
        residues = int(got["sequence"].str.len_bytes().sum())
        print(json.dumps({"positions": got.height, "residues": residues, "read_s": took}))


if __name__ == "__main__":
    main()
