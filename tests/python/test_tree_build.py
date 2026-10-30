"""workflows/tree-index/build.py: the index holds what it should, against brute force."""

import json
import sys
from pathlib import Path

import numpy as np
import polars as pl
from test_tree_phase0 import SCRIPT, SURVEY, run

from kmer_functional_profiler import _core

BUILD = SCRIPT.with_name("build.py")
QUERY = SCRIPT.with_name("query.py")
sys.path.insert(0, str(SCRIPT.parent))
import build  # noqa: E402
import phase0 as p0  # noqa: E402
import query  # noqa: E402


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
    survey = pl.read_csv(tmp_path / "s" / "survey.tsv", separator="\t")
    row = survey.filter(pl.col("leaves") == "proteins").row(0, named=True)
    # the index is the survey's configuration written out
    for kind in ("base_postings", "floor_postings", "topup_postings"):
        assert stats[kind] == row[kind], kind

    held = build.holdings(idx)  # raises unless every cover sums to 0 or 1
    members = p0.load_members(tmp_path / "members.bucket0.parquet")
    prot = pl.read_parquet(idx / "proteins.parquet")
    seqs = p0.seq_bytes(members)
    lengths = np.array([len(s) for s in seqs])
    pos = prot["pos"].to_numpy()
    # proteins share a position only when their dense 8-mer sets are equal
    d8 = p0.digests(p0.dense_sets(seqs, 8))
    by_pos = pl.DataFrame({"pos": pos, "d0": d8[:, 0], "d1": d8[:, 1]}).filter(pl.col("pos") >= 0)
    assert by_pos.group_by("pos").agg(pl.struct("d0", "d1").n_unique() == 1)["d0"].all()
    assert by_pos.group_by("d0", "d1").agg(pl.col("pos").n_unique() == 1)["pos"].all()
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


def test_query_counts_match_brute_force(tmp_path: Path) -> None:
    run("fixture", "--families", "15", "--members", "8", "--buckets", "1", cwd=tmp_path)
    cfg = ["--t-base", "0.02", "--n-min", "16", "--tau", "0.3", "--n-1", "2"]
    run("--members", "members.bucket0.parquet", "--out", "idx", *cfg, cwd=tmp_path, script=BUILD)
    idx = tmp_path / "idx"
    members = p0.load_members(tmp_path / "members.bucket0.parquet")
    # a sample: every third protein, some twice, and one unrelated sequence
    picked = [*members["sequence"][::3], *members["sequence"][::9], "MKVLAAGIVGLLLAQ" * 4]
    fasta = tmp_path / "sample.fasta"
    fasta.write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(picked)))
    run("--index", str(idx), "--proteins", str(fasta), "--out", "q", cwd=tmp_path, script=QUERY)
    leaves = pl.read_parquet(tmp_path / "q" / "leaves.parquet")
    nodes = pl.read_parquet(tmp_path / "q" / "nodes.parquet")
    sample = query.sample_counts(query.protein_hashes(query.Index(idx), str(fasta)))
    assert leaves.height > 0 and nodes.height > 0

    # leaves: brute-force holdings joined to the sample
    want = (
        build.holdings(idx)
        .join(sample, on="hash")
        .group_by("pos")
        .agg(hit_kmers=pl.len(), hits=pl.col("n").sum())
    )
    got = leaves.select("pos", "hit_kmers", "hits")
    assert got.sort("pos").equals(want.sort("pos").cast(got.schema))  # type: ignore[arg-type]
    size = np.load(idx / "pos_size.npy")
    held = build.holdings(idx).group_by("pos").len()
    assert (size[held["pos"].to_numpy()] == held["len"].to_numpy()).all()

    # internal nodes: a k-mer is present at a node iff the signs of the postings whose range
    # contains the node's sum to 1 (checked one node at a time)
    t = pl.DataFrame({x: np.load(idx / f"{x}.npy") for x in ("hash", "lo", "hi", "sign")})
    n_lo, n_hi = np.load(idx / "node_lo.npy"), np.load(idx / "node_hi.npy")
    node_size = np.load(idx / "node_size.npy")
    by_node = {r[0]: r[1:] for r in nodes.select("node", "hit_kmers", "hits").rows()}
    for v in range(len(n_lo)):
        state = (
            t.filter((pl.col("lo") <= n_lo[v]) & (pl.col("hi") >= n_hi[v]))
            .group_by("hash")
            .agg(pl.col("sign").cast(pl.Int32).sum())
            .filter(pl.col("sign") == 1)
        )
        assert node_size[v] == state.height
        hit = state.join(sample, on="hash")
        expect = (hit.height, int(hit["n"].sum())) if hit.height else None
        assert by_node.get(v) == expect, v

    # with a leaf per protein set, every sampled protein with a position hits its own leaf
    prot = pl.read_parquet(idx / "proteins.parquet").with_row_index("i")
    own = prot.filter(pl.col("i").is_in(list(range(0, len(members), 3))), pl.col("pos") >= 0)
    assert own.height > 0
    assert own.join(leaves, on="pos", how="anti").height == 0
