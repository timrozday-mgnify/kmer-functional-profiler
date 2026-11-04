"""Tree index: merge part indexes (``build.py`` on ``partition.py``'s parts) into one.

Parts are taken in name order. Positions, node ranges, spanning-tree edges, pieces and protein
positions are offset by the positions before each part; internal nodes are numbered part by part;
postings are re-sorted by (hash, lo). The result has ``build.py``'s format, so ``query.py`` reads
it as is; ``parts.json`` gives each part's first position (``store.py`` builds per part).

Out of core, so a release-sized index merges in bounded memory: tables are copied part by part
into memory-mapped outputs, and postings are merged one hash range (``--ranges``) at a time.
Each part's postings are sorted by hash (``build.py``), and parts are in position order, so a
range's postings sorted stably by (hash, lo) are the whole merge's for that range. The merged
index leaves out the build's files (``kind.npy``, ``pieces.parquet``).
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np
import polars as pl
import query
from numpy.lib.format import open_memmap

ARRAYS = ("hash", "lo", "hi", "sign", "bits", "kcode")  # kind and pieces stay with the parts
TABLES = ("node_level", "pos_size", "node_size", "pos_len", "edge_w")  # concatenated as they are
OFFSET = ("node_lo", "node_hi", "edge_a", "edge_b", "layer_group")  # plus each part's offset


def load(part: Path, name: str) -> np.ndarray:
    return np.load(part / f"{name}.npy", mmap_mode="r")


def tables(parts: list[Path], offs: list[int], out: Path) -> None:
    """Per-position, per-node and per-edge tables, part by part."""
    names = [*TABLES, *OFFSET[:4]] + (
        ["layer_group"] if (parts[0] / "layer_group.npy").exists() else []
    )
    for x in names:
        n = sum(len(load(p, x)) for p in parts)
        dtype = np.uint32 if x in OFFSET else load(parts[0], x).dtype
        dst = open_memmap(out / f"{x}.npy", mode="w+", dtype=dtype, shape=(n,))
        i = 0
        for p, off in zip(parts, offs, strict=True):
            src = load(p, x)
            dst[i : i + len(src)] = src.astype(np.uint64) + off if x in OFFSET else src
            i += len(src)
        dst.flush()
        del dst


def cuts(sorted_hashes: np.ndarray, ranges: int) -> np.ndarray:
    """Each hash range's first row, and the end."""
    inner = np.array([(r << 64) // ranges for r in range(1, ranges)], dtype=np.uint64)
    return np.r_[0, np.searchsorted(sorted_hashes, inner), len(sorted_hashes)]


def postings(parts: list[Path], offs: list[int], out: Path, ranges: int) -> dict:
    """Postings merged by hash range. Each part's postings carry their key's bits and k code;
    a key held in several parts gets their combined bits on every posting. ``kind`` is
    counted, not written (the build's)."""
    pcut = [cuts(load(p, "hash"), ranges) for p in parts]
    n = sum(int(c[-1]) for c in pcut)
    dtypes = {"hash": np.uint64, "lo": np.uint32, "hi": np.uint32, "sign": np.int8,
              "bits": np.uint8, "kcode": load(parts[0], "kcode").dtype}  # fmt: skip
    dst = {x: open_memmap(out / f"{x}.npy", mode="w+", dtype=dtypes[x], shape=(n,)) for x in ARRAYS}
    kinds, distinct, i = np.zeros(3, dtype=np.int64), 0, 0
    for r in range(ranges):
        got: dict[str, list[np.ndarray]] = {x: [] for x in (*ARRAYS, "kind")}
        keys = []  # per part, each hash's bits once
        for p, off, pc in zip(parts, offs, pcut, strict=True):
            a, b = pc[r], pc[r + 1]
            if b == a:
                continue
            for x in (*ARRAYS, "kind"):
                v = load(p, x)[a:b]
                got[x].append(v.astype(np.uint64) + off if x in ("lo", "hi") else np.asarray(v))
            h = got["hash"][-1]
            first = np.flatnonzero(np.r_[True, h[1:] != h[:-1]])
            keys.append(pl.DataFrame({"hash": h[first], "q": got["bits"][-1][first]}))
        if not got["hash"]:
            continue
        merged = {x: np.concatenate(v) for x, v in got.items()}
        order = np.lexsort((merged["lo"], merged["hash"]))  # stable: parts' order on ties
        # 2^-bits = q0 g, and q0 is a k-mer's own: gains across parts add
        k = (
            pl.concat(keys)
            .group_by("hash")
            .agg(bits=-(2.0 ** -(pl.col("q").cast(pl.Float64) * query.BITS_STEP)).sum().log(2))
            .sort("hash")
        )
        h = merged["hash"][order]
        merged["bits"] = query.encode_bits(
            k["bits"].to_numpy()[np.searchsorted(k["hash"].to_numpy(), h)]
        )
        m = len(order)
        for x in ARRAYS:
            dst[x][i : i + m] = merged[x] if x == "bits" else merged[x][order]
        i += m
        kinds += np.bincount(merged["kind"], minlength=3)[:3]
        distinct += k.height
    for x in ARRAYS:
        dst[x].flush()
    dst.clear()
    return {"postings": n, "kinds": kinds, "keys": distinct}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--parts", nargs="+", required=True, help="build.py output directories")
    parser.add_argument("--out", required=True)
    parser.add_argument("--ranges", type=int, default=256, help="hash ranges merged one at a time")
    a = parser.parse_args()
    t0 = time.perf_counter()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    parts = sorted(Path(p) for p in a.parts)
    metas = [json.loads((p / "meta.json").read_text()) for p in parts]
    params = metas[0]["params"]
    if any(m["params"] != params for m in metas):
        raise ValueError("parts were built with different parameters")
    sizes = [m["stats"]["positions"] for m in metas]
    offs = np.r_[0, np.cumsum(sizes)[:-1]].astype(int).tolist()  # each part's first position
    off = int(sum(sizes))
    if off >= 2**32:
        raise ValueError(f"{off} positions do not fit the format's u32 ranges")
    tables(parts, offs, out)
    got = postings(parts, offs, out, a.ranges)
    pos = lambda o: pl.when(pl.col("pos") >= 0).then(pl.col("pos") + o).otherwise(-1)  # noqa: E731
    pl.concat(
        [
            pl.scan_parquet(p / "proteins.parquet").with_columns(pos=pos(o))
            for p, o in zip(parts, offs, strict=True)
        ]
    ).sink_parquet(out / "proteins.parquet")
    prot = pl.scan_parquet(out / "proteins.parquet")
    # a cluster can straddle parts: distinct over all of them, by hash (one u64 per unit)
    units = len(
        np.unique(prot.select(pl.col("cluster_rep").hash()).collect().to_series().to_numpy())
    )
    files = {f.name: f.stat().st_size for f in sorted(out.iterdir()) if f.name != "meta.json"}
    summed = ("tree_leaves", "internal_nodes", "components", "pieces", "floored_pieces")
    kinds = got["kinds"]
    stats = {
        "proteins": prot.select(pl.len()).collect().item(),
        "units": units,
        "positions": off,
        **{k: sum(m["stats"][k] for m in metas) for k in summed},
        "postings": got["postings"],
        "base_postings": int(kinds[0]),
        "floor_postings": int(kinds[1]),
        "topup_postings": int(kinds[2]),
        "keys": got["keys"],
        "bytes": sum(files.values()),
        "bytes_per_posting": sum(files.values()) / max(got["postings"], 1),
        "parts": len(parts),
        "part_build_s": sum(m["stats"]["build_s"] for m in metas),
        "part_peak_rss_bytes": max(m["stats"]["peak_rss_bytes"] for m in metas),
        "merge_s": time.perf_counter() - t0,
    }
    meta = {"format": metas[0]["format"], "params": params, "stats": stats, "files": files}
    (out / "meta.json").write_text(json.dumps(meta, indent=1) + "\n")
    offsets = [
        {"name": p.name, "offset": o, "positions": s}
        for p, o, s in zip(parts, offs, sizes, strict=True)
    ]
    (out / "parts.json").write_text(json.dumps(offsets, indent=1) + "\n")
    print(json.dumps(stats))


if __name__ == "__main__":
    main()
