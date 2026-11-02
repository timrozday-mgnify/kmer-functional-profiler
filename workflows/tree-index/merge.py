"""Tree index: merge part indexes (``build.py`` on ``partition.py``'s parts) into one.

Parts are taken in name order. Positions, node ranges, spanning-tree edges, pieces and protein
positions are offset by the positions before each part; internal nodes are numbered part by part;
postings are re-sorted by (hash, lo). The result has ``build.py``'s format, so ``query.py`` reads
it as is; ``parts.json`` gives each part's first position (``store.py`` builds per part).

Out of core, so a release-sized index merges in bounded memory: tables are copied part by part
into memory-mapped outputs, and postings and keys are merged one hash range (``--ranges``) at a
time. Each part's postings and keys are sorted by hash (``build.py``), and parts are in position
order, so a range's postings sorted stably by (hash, lo) are the whole merge's for that range.
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np
import polars as pl
from numpy.lib.format import open_memmap

ARRAYS = ("hash", "lo", "hi", "sign", "kind")
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
    """Postings and keys, merged by hash range."""
    pcut = [cuts(load(p, "hash"), ranges) for p in parts]
    kcut = [cuts(load(p, "keys"), ranges) for p in parts]
    n = sum(int(c[-1]) for c in pcut)
    dtypes = {
        "hash": np.uint64,
        "lo": np.uint32,
        "hi": np.uint32,
        "sign": np.int8,
        "kind": np.uint8,
    }
    dst = {x: open_memmap(out / f"{x}.npy", mode="w+", dtype=dtypes[x], shape=(n,)) for x in ARRAYS}
    # keys: their number is known only after the merge, so raw files first, .npy at the end
    raw = {x: open(out / f"{x}.raw", "wb") for x in ("keys", "key_bits", "key_k")}  # noqa: SIM115
    kinds, distinct, i = np.zeros(3, dtype=np.int64), 0, 0
    for r in range(ranges):
        got: dict[str, list[np.ndarray]] = {x: [] for x in ARRAYS}
        bits = []
        for p, off, pc, kc in zip(parts, offs, pcut, kcut, strict=True):
            a, b = pc[r], pc[r + 1]
            if b > a:
                for x in ARRAYS:
                    v = load(p, x)[a:b]
                    got[x].append(v.astype(np.uint64) + off if x in ("lo", "hi") else np.asarray(v))
            a, b = kc[r], kc[r + 1]
            if b > a:
                names = (("hash", "keys"), ("bits", "key_bits"), ("k", "key_k"))
                bits.append(pl.DataFrame({x: np.asarray(load(p, y)[a:b]) for x, y in names}))
        if got["hash"]:
            merged = {x: np.concatenate(v) for x, v in got.items()}
            order = np.lexsort((merged["lo"], merged["hash"]))  # stable: parts' order on ties
            m = len(order)
            for x in ARRAYS:
                dst[x][i : i + m] = merged[x][order]
            i += m
            kinds += np.bincount(merged["kind"], minlength=3)[:3]
            distinct += len(np.unique(merged["hash"]))
        if bits:
            # 2^-bits = q0 g, and q0 is a k-mer's own: gains across parts add
            keys = (
                pl.concat(bits)
                .group_by("hash")
                .agg(
                    pl.col("k").first(),
                    bits=-(2.0 ** -pl.col("bits").cast(pl.Float64)).sum().log(2),
                )
                .sort("hash")
            )
            keys["hash"].to_numpy().tofile(raw["keys"])
            keys["bits"].to_numpy().astype(np.float32).tofile(raw["key_bits"])
            keys["k"].to_numpy().tofile(raw["key_k"])
    for x in ARRAYS:
        dst[x].flush()
    dst.clear()
    k_dtype = load(parts[0], "key_k").dtype
    for x, dtype in (("keys", np.uint64), ("key_bits", np.float32), ("key_k", k_dtype)):
        raw[x].close()
        if (out / f"{x}.raw").stat().st_size == 0:  # np.memmap refuses an empty file
            np.save(out / f"{x}.npy", np.empty(0, dtype=dtype))
            (out / f"{x}.raw").unlink()
            continue
        src = np.memmap(out / f"{x}.raw", dtype=dtype, mode="r")
        npy = open_memmap(out / f"{x}.npy", mode="w+", dtype=dtype, shape=src.shape)
        npy[:] = src
        npy.flush()
        del npy, src
        (out / f"{x}.raw").unlink()
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
    shift = lambda p, o, cols: pl.scan_parquet(p).with_columns(cols + o)  # noqa: E731
    pl.concat(
        [
            shift(p / "pieces.parquet", o, pl.col("lo", "hi"))
            for p, o in zip(parts, offs, strict=True)
        ]
    ).sink_parquet(out / "pieces.parquet")
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
