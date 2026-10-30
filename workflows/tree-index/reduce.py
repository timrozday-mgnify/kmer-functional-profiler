"""Tree index, phase 0: one hash range of every bucket's streaming counts -> per kind, the
histogram of group sizes and the largest groups (``main.nf``'s REDUCE).

Out of core: each bucket's file is read once and split by key into ``--passes``
sub-partitions on local disk; each sub-partition is then summed on its own. Memory is one
bucket's file or one sub-partition, whichever is larger, not the whole range (~2 × 10⁹ rows
of the whole release at 16 ranges, which ran out of memory in one ``group_by``).

A separate script, not a ``phase0.py`` subcommand, so that changing it leaves ``phase0.py``,
and with it the cached SAMPLE and STREAM tasks, untouched.
"""

import argparse
import tempfile
from pathlib import Path

import polars as pl


def reduce(parts: list[str], prefix: str, top: int, passes: int, tmp: str | None) -> None:
    # the files are already split by key % ranges, so sub-partition on other bits
    sub = ((pl.col("key") // 2**32) % passes).alias("sub")
    with tempfile.TemporaryDirectory(dir=tmp) as scratch:
        for i, path in enumerate(parts):
            table = pl.read_parquet(path).with_columns(sub)
            for (s,), part in table.partition_by("sub", as_dict=True).items():
                part.drop("sub").write_parquet(Path(scratch) / f"s{s}.{i}.parquet")
        hists, tops = [], []
        for s in range(passes):
            files = sorted(Path(scratch).glob(f"s{s}.*.parquet"))
            if not files:
                continue
            summed = (
                pl.scan_parquet(files)
                .group_by("kind", "key")
                .agg(pl.sum("count"))
                .collect(engine="streaming")
            )
            hists.append(summed.group_by("kind", size="count").len("n"))
            tops.append(summed.sort("count", descending=True).group_by("kind").head(top))
            for f in files:
                f.unlink()
    schema = {"kind": pl.String, "size": pl.UInt32, "n": pl.UInt32}
    hist = (
        pl.concat(hists).group_by("kind", "size").agg(pl.sum("n"))
        if hists
        else pl.DataFrame(schema=schema)
    )
    best = (
        pl.concat(tops).sort("count", descending=True).group_by("kind").head(top)
        if tops
        else pl.DataFrame(schema={"kind": pl.String, "key": pl.UInt64, "count": pl.UInt32})
    )
    hist.write_parquet(f"{prefix}.hist.parquet")
    best.write_parquet(f"{prefix}.top.parquet")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--parts", nargs="+", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument("--passes", type=int, default=64, help="sub-partitions summed in turn")
    parser.add_argument("--tmp", default=".", help="scratch directory (default: the task's)")
    args = parser.parse_args()
    reduce(args.parts, args.prefix, args.top, args.passes, args.tmp)


if __name__ == "__main__":
    main()
