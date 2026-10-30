"""workflows/tree-index/build.py: the index holds what it should, against brute force."""

import json
import sys
from pathlib import Path

import numpy as np
import polars as pl
from test_tree_phase0 import SCRIPT, SURVEY, run

from kmer_functional_profiler import _core

BUILD = SCRIPT.with_name("build.py")
sys.path.insert(0, str(SCRIPT.parent))
import build  # noqa: E402
import phase0 as p0  # noqa: E402


def test_build_matches_brute_force_and_survey(tmp_path: Path) -> None:
    run("fixture", "--families", "15", "--members", "8", "--buckets", "1", cwd=tmp_path)
    cfg = ["--t-base", "0.02", "--n-min", "16", "--tau", "0.3", "--n-1", "2"]
    run("--members", "members.bucket0.parquet", "--out", "idx", *cfg, cwd=tmp_path, script=BUILD)
    run(
        "--members", "members.bucket0.parquet", "--out", "s", "--t-bases", "0.02",
        "--n-mins", "16", "--taus", "0.3", "--n-singles", "0", "--n-1s", "2",
        cwd=tmp_path, script=SURVEY,
    )  # fmt: skip
    idx = tmp_path / "idx"
    stats = json.loads((idx / "meta.json").read_text())["stats"]
    row = pl.read_csv(tmp_path / "s" / "survey.tsv", separator="\t").row(0, named=True)
    # the index is the survey's configuration written out
    for kind in ("base_postings", "floor_postings", "topup_postings"):
        assert stats[kind] == row[kind], kind

    held = build.holdings(idx)  # raises unless every cover sums to 0 or 1
    members = p0.load_members(tmp_path / "members.bucket0.parquet")
    prot = pl.read_parquet(idx / "proteins.parquet")
    seqs = p0.seq_bytes(members)
    lengths = np.array([len(s) for s in seqs])
    pos = prot["pos"].to_numpy()
    # each tree position's base band is its representative's (longest protein, then first)
    n = stats["tree_leaves"]
    rep = (
        pl.DataFrame({"seq": np.arange(len(seqs)), "pos": pos, "len": lengths})
        .filter((pl.col("pos") >= 0) & (pl.col("pos") < n))
        .sort("pos", "len", "seq", descending=[False, True, False])
        .group_by("pos", maintain_order=True)
        .first()
    )
    base = _core.max_hash(0.02)
    seq, h = p0.hashed([seqs[i] for i in rep["seq"]], 11, base)
    want = pl.DataFrame({"hash": h, "pos": rep["pos"].to_numpy()[seq]}).sort("pos", "hash")
    got = held.filter(pl.col("hash") <= base).filter(pl.col("pos") < n)
    assert got.equals(want.with_columns(pl.col("pos").cast(got["pos"].dtype)))
    # every floored piece's floor is held across it
    t = pl.DataFrame({x: np.load(idx / f"{x}.npy") for x in ("hash", "lo", "hi", "kind")})
    floors = t.filter(pl.col("kind") == build.FLOOR)
    assert floors.height > 0
    spread = floors.select("hash", pos=pl.int_ranges("lo", "hi")).explode("pos")
    assert spread.join(held, on=["hash", "pos"], how="anti").height == 0
    per_pos = held.group_by("pos").len()
    # every protein with a valid k-mer is covered by at least min(n_1, its k-mers)
    covered = dict(per_pos.rows())
    n_valid = _core.distinct_kmers(seqs, np.arange(len(seqs), dtype=np.uint32), 11)
    valid = np.zeros(len(seqs), dtype=np.int64)
    valid[n_valid["group"]] = n_valid["n_kmers"]
    for i in range(len(seqs)):
        if valid[i]:
            assert pos[i] >= 0
            assert covered.get(int(pos[i]), 0) >= min(2, valid[i])
        else:
            assert pos[i] == -1
