"""workflows/tree-index/build.py (and partition.py + merge.py): the index holds what it should,
against brute force; query.py counts it, and report.py reports it."""

import json
import sys
from pathlib import Path

import numpy as np
import polars as pl
from test_tree_phase0 import SCRIPT, SURVEY, run

from kmer_functional_profiler import _core

BUILD = SCRIPT.with_name("build.py")
QUERY = SCRIPT.with_name("query.py")
REPORT = SCRIPT.with_name("report.py")
PARTITION = SCRIPT.with_name("partition.py")
MERGE = SCRIPT.with_name("merge.py")
sys.path.insert(0, str(SCRIPT.parent))
import build  # noqa: E402
import phase0 as p0  # noqa: E402
import query  # noqa: E402
import report  # noqa: E402
from scipy.sparse import coo_array  # noqa: E402
from scipy.sparse.csgraph import connected_components  # noqa: E402

# phase 1's medium tier: base band, floors and top-ups only, which check_index checks in full
MEDIUM = ["--t-base", "0.02", "--n-min", "16", "--tau", "0.3", "--n-1", "2", "--window", "0",
          "--ladder", ""]  # fmt: skip


def fixture(tmp_path: Path) -> Path:
    run("fixture", "--families", "15", "--members", "8", "--buckets", "1", cwd=tmp_path)
    return tmp_path / "members.bucket0.parquet"


def check_index(idx: Path, members_path: Path, n_1: int = 2) -> None:
    """Positions, base band, floors and coverage against brute force."""
    stats = json.loads((idx / "meta.json").read_text())["stats"]
    held = build.holdings(idx)  # raises unless every cover sums to 0 or 1
    members = p0.load_members(members_path)
    seqs = p0.seq_bytes(members)
    prot = members.select("protein_id").join(
        pl.read_parquet(idx / "proteins.parquet"),
        on="protein_id",
        how="left",
        maintain_order="left",
    )
    pos = prot["pos"].to_numpy()
    assert prot.height == stats["proteins"] and not prot["pos"].is_null().any()
    # proteins share a position only when their dense 8-mer sets are equal
    d8 = p0.digests(p0.dense_sets(seqs, 8))
    by_pos = pl.DataFrame({"pos": pos, "d0": d8[:, 0], "d1": d8[:, 1]}).filter(pl.col("pos") >= 0)
    assert by_pos.group_by("pos").agg(pl.struct("d0", "d1").n_unique() == 1)["d0"].all()
    assert by_pos.group_by("d0", "d1").agg(pl.col("pos").n_unique() == 1)["pos"].all()
    # each position's base band is its representative's (longest protein, then first)
    lengths = np.array([len(s) for s in seqs])
    rep = (
        pl.DataFrame({"seq": np.arange(len(seqs)), "pos": pos, "len": lengths})
        .filter(pl.col("pos") >= 0)
        .sort("pos", "len", "seq", descending=[False, True, False])
        .group_by("pos", maintain_order=True)
        .first()
    )
    base = _core.max_hash(0.02)
    seq, h = p0.hashed([seqs[i] for i in rep["seq"]], 11, base)
    want = pl.DataFrame({"hash": h, "pos": rep["pos"].to_numpy()[seq]}).sort("pos", "hash")
    got = held.filter(pl.col("hash") <= base)
    assert got.equals(want.with_columns(pl.col("pos").cast(got["pos"].dtype)))
    # every floored piece's floor is held across it
    t = pl.DataFrame({x: np.load(idx / f"{x}.npy") for x in ("hash", "lo", "hi", "kind")})
    floors = t.filter(pl.col("kind") == build.FLOOR)
    assert floors.height > 0
    spread = floors.select("hash", pos=pl.int_ranges("lo", "hi")).explode("pos")
    assert spread.join(held, on=["hash", "pos"], how="anti").height == 0
    # every protein with a valid k-mer is covered by at least min(n_1, its k-mers)
    covered = dict(held.group_by("pos").len().rows())
    n_valid = _core.distinct_kmers(seqs, np.arange(len(seqs), dtype=np.uint32), 11)
    valid = np.zeros(len(seqs), dtype=np.int64)
    valid[n_valid["group"]] = n_valid["n_kmers"]
    for i in range(len(seqs)):
        if valid[i]:
            assert pos[i] >= 0
            assert covered.get(int(pos[i]), 0) >= min(n_1, valid[i])
        else:
            assert pos[i] == -1


def check_query(tmp_path: Path, idx: Path, members_path: Path) -> None:
    """query.py's leaf and node counts against brute force; every sampled protein hits its
    own leaf."""
    members = p0.load_members(members_path)
    # a sample: every third protein, some twice, and one unrelated sequence
    picked = [*members["sequence"][::3], *members["sequence"][::9], "MKVLAAGIVGLLLAQ" * 4]
    fasta = tmp_path / "sample.fasta"
    fasta.write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(picked)))
    out = tmp_path / f"q-{idx.name}"
    run(
        "--index", str(idx), "--proteins", str(fasta), "--out", str(out), cwd=tmp_path, script=QUERY
    )
    leaves = pl.read_parquet(out / "leaves.parquet")
    nodes = pl.read_parquet(out / "nodes.parquet")
    sample = query.sample_counts(query.protein_hashes(query.Index(idx), str(fasta)))
    assert leaves.height > 0 and nodes.height > 0

    # leaves: brute-force holdings joined to the sample
    held = build.holdings(idx)
    want = (
        held.join(sample, on="hash").group_by("pos").agg(hit_kmers=pl.len(), hits=pl.col("n").sum())
    )
    got = leaves.select("pos", "hit_kmers", "hits")
    assert got.sort("pos").equals(want.sort("pos").cast(got.schema))  # type: ignore[arg-type]
    size = np.load(idx / "pos_size.npy")
    per_pos = held.group_by("pos").len()
    assert (size[per_pos["pos"].to_numpy()] == per_pos["len"].to_numpy()).all()

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

    # evidence: per leaf, the best read's summed key bits less a charge per extra hit
    index = query.Index(idx)
    hs = sample["hash"].to_numpy()
    hs = hs[np.isin(hs, index.keys)]
    rng = np.random.default_rng(0)
    pairs = pl.DataFrame({"read": rng.integers(0, 6, len(hs)).astype(np.uint64), "hash": hs})
    got = query.evidence(index, pairs, 3.0, floor=-1e9)
    bits = pl.DataFrame({"hash": index.keys, "bits": index.key_bits.astype(np.float64)})
    want = (
        held.join(pairs, on="hash")
        .join(bits, on="hash")
        .group_by("read", "pos")
        .agg(score=pl.col("bits").sum() - (pl.len() - 1) * 3.0)
        .group_by("pos")
        .agg(pl.col("score").max())
        .sort("pos")
    )
    assert got["pos"].to_list() == want["pos"].to_list()
    assert np.allclose(got["score"].to_numpy(), want["score"].to_numpy(), atol=1e-3)
    assert (index.key_bits > 0).all() and np.array_equal(index.keys, np.unique(index.hash))

    # with a leaf per protein set, every sampled protein with a position hits its own leaf
    prot = (
        members.select("protein_id")
        .with_row_index("i")
        .join(pl.read_parquet(idx / "proteins.parquet"), on="protein_id")
    )
    own = prot.filter(pl.col("i").is_in(list(range(0, len(members), 3))), pl.col("pos") >= 0)
    assert own.height > 0
    assert own.join(leaves, on="pos", how="anti").height == 0
    check_report(tmp_path, idx, out, held, own)


def check_report(tmp_path: Path, idx: Path, q: Path, held: pl.DataFrame, own: pl.DataFrame) -> None:
    """report.py: AAI bounds, candidates against brute-force single linkage, gather's
    assignment."""
    out = tmp_path / f"r-{idx.name}"
    run("--index", str(idx), "--query", str(q), "--out", str(out), "--max-neighbours", "4",
        cwd=tmp_path, script=REPORT)  # fmt: skip
    leaves = pl.read_parquet(out / "leaves.parquet")
    assert (leaves["aai_lo"] <= leaves["aai"]).all() and (leaves["aai"] <= leaves["aai_hi"]).all()
    full = leaves.filter(pl.col("hit_kmers") == pl.col("size"))
    assert full.height > 0 and (full["aai"] == 1).all()

    # candidates: the seeds and each seed's walk, as ranges; members are their proteins
    cand = pl.read_parquet(out / "candidates.parquet")
    pos = np.concatenate([np.arange(lo, hi) for lo, hi in cand.rows()])
    assert (np.diff(pos) > 0).all() and np.isin(leaves["pos"], pos).all()
    members = pl.read_parquet(out / "members.parquet")
    assert np.isin(members["pos"], pos).all()
    assert members.filter(pl.col("seed"))["pos"].unique().len() == leaves.height

    # each walk, uncapped, is the seed's single-linkage cluster at its radius (brute force:
    # components of the edges <= radius); capped, the nearest of it in minimax distance
    n_pos = len(np.load(idx / "pos_len.npy"))
    tree = report.Spanning(idx, n_pos)
    a, b, w = (np.load(idx / f"edge_{x}.npy").astype(np.int64) for x in "abw")
    grown = 0
    for radius in (0, 5, 20, 80, 10**9):
        keep = w <= radius
        graph = coo_array((np.ones(keep.sum()), (a[keep], b[keep])), shape=(n_pos, n_pos))
        _, comp = connected_components(graph, directed=False)
        for v in leaves["pos"].to_list():
            got, cut = tree.cluster(v, radius, n_pos)
            assert not cut and sorted(got) == np.flatnonzero(comp == comp[v]).tolist()
            near, cut = tree.cluster(v, radius, 4)
            assert near == got[:5] and cut == (len(got) > 5)
            grown += len(got) > 1
    assert grown > 0

    # gather: groups take disjoint k-mers, together every hit k-mer a seed holds; a group's
    # range is its common subtree, which holds its leaf. Every sampled protein is a candidate
    g = pl.read_parquet(out / "gather.parquet")
    hits = pl.read_parquet(q / "hits.parquet")
    seed_kmers = held.join(hits, on="hash").filter(pl.col("pos").is_in(leaves["pos"].implode()))
    assert int(g["kmers"].sum()) == seed_kmers["hash"].n_unique()
    assert g.filter((pl.col("pos") < pl.col("lo")) | (pl.col("pos") >= pl.col("hi"))).height == 0
    assert ((g["tied"] > 1) == (g["node"] >= 0)).all()
    assert np.isin(own["pos"], members["pos"]).all()


def test_build_matches_brute_force_and_survey(tmp_path: Path) -> None:
    members = fixture(tmp_path)
    run("--members", str(members), "--out", "idx", *MEDIUM, cwd=tmp_path, script=BUILD)
    run(
        "--members", str(members), "--out", "s", "--t-bases", "0.02",
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
    check_index(idx, members)
    check_query(tmp_path, idx, members)


def test_partitioned_build_is_sound(tmp_path: Path) -> None:
    members = fixture(tmp_path)
    run("--members", str(members), "--out", "parts", "--part-proteins", "30",
        cwd=tmp_path, script=PARTITION)  # fmt: skip
    info = json.loads((tmp_path / "parts" / "partition.json").read_text())
    parts = sorted((tmp_path / "parts").glob("part*.parquet"))
    assert info["parts"] == len(parts) >= 3
    assert sum(info["part_proteins"]) == pl.read_parquet(members).height
    built = []
    for p in parts:
        out = tmp_path / f"b-{p.stem}"
        run("--members", str(p), "--out", str(out), *MEDIUM, cwd=tmp_path, script=BUILD)
        built.append(str(out))
    run("--parts", *built, "--out", "merged", cwd=tmp_path, script=MERGE)
    idx = tmp_path / "merged"
    stats = json.loads((idx / "meta.json").read_text())["stats"]
    assert stats["parts"] == len(parts)
    check_index(idx, members)
    check_query(tmp_path, idx, members)


def test_window_kmers_match_brute_force() -> None:
    """build.window_kmers: each sequence's least hash per window (shrunk by min_count)."""
    rng = np.random.default_rng(2)
    seqs = [bytes(rng.choice(list(b"ACDEFGHIKLMNPQRSTVWY"), n)) for n in (5, 15, 22, 40, 120, 61)]
    seqs[4] = seqs[4][:50] + b"X" + seqs[4][51:]  # k-mers over X have no hash
    for w, min_count in ((4, 0), (20, 0), (20, 4), (10, 8)):
        seq, h = build.window_kmers(seqs, 11, w, min_count)
        for j, s in enumerate(seqs):
            hs = [_core.hash_proteins([s[i : i + 11]], 11)["hash"] for i in range(len(s) - 10)]
            hs = [int(x[0]) if len(x) else None for x in hs]
            ww = max(1, min(w, len(hs) // min_count)) if min_count else w
            wins = [hs[i : i + ww] for i in range(len(hs) - ww + 1)] if len(hs) >= ww else [hs]
            want = {min(v) for v in ([x for x in win if x is not None] for win in wins) if v}
            assert {int(x) for x in h[seq == j]} == want, (w, min_count, j)


def test_default_modes_query_soundly(tmp_path: Path) -> None:
    """The default build (standard: spread minimizers, k 7 ladder) and ``--mode recall`` (a
    Dayhoff layer on every leaf): covers sum to 0 or 1, and the query and report checks hold
    against brute force, the Dayhoff layer's k-mers included."""
    members = fixture(tmp_path)
    for mode in ("standard", "recall"):
        idx = tmp_path / f"idx-{mode}"
        run(
            "--members", str(members), "--out", str(idx), "--mode", mode, cwd=tmp_path, script=BUILD
        )
        meta = json.loads((idx / "meta.json").read_text())
        assert meta["params"]["layers"] == ("dayhoff:7:8:rare" if mode == "recall" else "")
        codes = set(np.load(idx / "key_k.npy").tolist())
        assert (71 in codes) == (mode == "recall")  # Dayhoff 7: 7 + 32 x 2
        check_query(tmp_path, idx, members)
