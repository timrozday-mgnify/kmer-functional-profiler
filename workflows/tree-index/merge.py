"""Tree index: merge part indexes (``build.py`` on ``partition.py``'s parts) into one.

Parts are taken in name order. Positions, node ranges, spanning-tree edges, pieces and protein
positions are offset by the positions before each part; internal nodes are numbered part by part;
postings are re-sorted by (hash, lo). The result has ``build.py``'s format, so ``query.py`` reads
it as is.
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np
import polars as pl

ARRAYS = ("hash", "lo", "hi", "sign", "kind")
TABLES = ("node_level", "pos_size", "node_size", "pos_len", "edge_w")  # concatenated as they are


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--parts", nargs="+", required=True, help="build.py output directories")
    parser.add_argument("--out", required=True)
    a = parser.parse_args()
    t0 = time.perf_counter()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    parts = sorted(Path(p) for p in a.parts)
    metas = [json.loads((p / "meta.json").read_text()) for p in parts]
    params = metas[0]["params"]
    if any(m["params"] != params for m in metas):
        raise ValueError("parts were built with different parameters")
    cols: dict[str, list[np.ndarray]] = {
        x: [] for x in (*ARRAYS, *TABLES, "node_lo", "node_hi", "edge_a", "edge_b")
    }
    pieces, proteins, bits = [], [], []
    off = 0
    for p, m in zip(parts, metas, strict=True):
        load = lambda name: np.load(p / f"{name}.npy")  # noqa: B023, E731
        for x in ("hash", "sign", "kind", *TABLES):
            cols[x].append(load(x))
        for x in ("lo", "hi", "node_lo", "node_hi", "edge_a", "edge_b"):
            cols[x].append(load(x).astype(np.uint64) + off)
        names = (("hash", "keys"), ("bits", "key_bits"), ("k", "key_k"))
        bits.append(pl.DataFrame({x: np.load(p / f"{y}.npy") for x, y in names}))
        pieces.append(pl.read_parquet(p / "pieces.parquet").with_columns(pl.col("lo", "hi") + off))
        proteins.append(
            pl.read_parquet(p / "proteins.parquet").with_columns(
                pos=pl.when(pl.col("pos") >= 0).then(pl.col("pos") + off).otherwise(-1)
            )
        )
        off += m["stats"]["positions"]
    if off >= 2**32:
        raise ValueError(f"{off} positions do not fit the format's u32 ranges")
    merged = {x: np.concatenate(v) for x, v in cols.items()}
    order = np.lexsort((merged["lo"], merged["hash"]))
    for x in ARRAYS:
        dtype = np.uint32 if x in ("lo", "hi") else merged[x].dtype
        np.save(out / f"{x}.npy", merged[x][order].astype(dtype))
    for x in ("node_lo", "node_hi", "edge_a", "edge_b", *TABLES):
        dtype = np.uint32 if x in ("node_lo", "node_hi", "edge_a", "edge_b") else merged[x].dtype
        np.save(out / f"{x}.npy", merged[x].astype(dtype))
    # 2^-bits = q0 g, and q0 is a k-mer's own: gains across parts add
    keys = (
        pl.concat(bits)
        .group_by("hash")
        .agg(pl.col("k").first(), bits=-(2.0 ** -pl.col("bits").cast(pl.Float64)).sum().log(2))
        .sort("hash")
    )
    np.save(out / "keys.npy", keys["hash"].to_numpy())
    np.save(out / "key_bits.npy", keys["bits"].to_numpy().astype(np.float32))
    np.save(out / "key_k.npy", keys["k"].to_numpy())
    pl.concat(pieces).write_parquet(out / "pieces.parquet")
    prot = pl.concat(proteins)
    prot.write_parquet(out / "proteins.parquet")
    kinds = np.bincount(merged["kind"], minlength=3)
    files = {f.name: f.stat().st_size for f in sorted(out.iterdir()) if f.name != "meta.json"}
    summed = ("tree_leaves", "internal_nodes", "components", "pieces", "floored_pieces")
    stats = {
        "proteins": prot.height,
        "units": prot["cluster_rep"].n_unique(),
        "positions": off,
        **{k: sum(m["stats"][k] for m in metas) for k in summed},
        "postings": len(order),
        "base_postings": int(kinds[0]),
        "floor_postings": int(kinds[1]),
        "topup_postings": int(kinds[2]),
        "keys": len(np.unique(merged["hash"])),
        "bytes": sum(files.values()),
        "bytes_per_posting": sum(files.values()) / max(len(order), 1),
        "parts": len(parts),
        "part_build_s": sum(m["stats"]["build_s"] for m in metas),
        "part_peak_rss_bytes": max(m["stats"]["peak_rss_bytes"] for m in metas),
        "merge_s": time.perf_counter() - t0,
    }
    meta = {"format": metas[0]["format"], "params": params, "stats": stats, "files": files}
    (out / "meta.json").write_text(json.dumps(meta, indent=1) + "\n")
    print(json.dumps(stats))


if __name__ == "__main__":
    main()
