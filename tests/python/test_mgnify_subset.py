"""The MGnify subset extraction on the mini release fixture."""

import gzip
import json
import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "workflows" / "mgnify-subset" / "mgnify_subset.py"
RELEASE = ROOT / "tests" / "data" / "mini_release"
GUT = "root:Host-associated:Human:Digestive system"
SHARD = ("--release", str(RELEASE), "--membership", "membership")


def run(cwd: Path, *args: str) -> None:
    subprocess.run([sys.executable, SCRIPT, *args], cwd=cwd, check=True)


def test_extracts_whole_gut_clusters(tmp_path: Path) -> None:
    run(tmp_path, "membership", "--release", str(RELEASE), "--biome", GUT, "--shard-width", "300")
    for i, (lo, hi) in enumerate([(0, 400), (400, 1000)]):
        run(tmp_path, "extract", *SHARD, "--lo", str(lo), "--hi", str(hi), "--prefix", f"s{i}")
    run(tmp_path, "merge", "--pfam", "pfam.parquet", "s0", "s1")

    clusters = pl.read_parquet(RELEASE / "mgy_clusters.parquet")
    gut = clusters.filter(pl.col("cluster_members_biomes").str.contains(GUT, literal=True))
    expected = (
        pl.read_parquet(RELEASE / "mgy_cluster_seqs.parquet")
        .join(gut, on="cluster_rep", how="semi")
        .join(
            pl.read_parquet(RELEASE / "mgy_protein_sequences.parquet"),
            left_on="cluster_member",
            right_on="protein_id",
        )
        .select(
            pl.col("cluster_member").alias("protein_id"), "cluster_rep", "full_length", "sequence"
        )
        .sort("cluster_rep", "protein_id")
    )
    members = pl.read_parquet(tmp_path / "members.parquet")
    assert 0 < gut.height < clusters.height
    assert members.equals(expected)
    pfam = pl.read_parquet(tmp_path / "pfam.parquet")
    assert set(pfam["protein_id"]) <= set(members["protein_id"])
    assert pfam.height > 0


def test_membership_fails_above_max_protein_id(tmp_path: Path) -> None:
    result = subprocess.run(
        [
            sys.executable,
            SCRIPT,
            "membership",
            "--release",
            str(RELEASE),
            "--max-protein-id",
            "500",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0 and "max-protein-id" in result.stderr


def test_buckets_split_whole_clusters_and_subsets_nest(tmp_path: Path) -> None:
    run(tmp_path, "membership", "--release", str(RELEASE))
    for i, (lo, hi) in enumerate([(0, 400), (400, 1000)]):
        shard = ("--lo", str(lo), "--hi", str(hi), "--prefix", f"s{i}")
        run(tmp_path, "extract", *SHARD, *shard, "--buckets", "3", "--no-pfam")
    for b in range(3):
        run(tmp_path, "merge", "--bucket", str(b), "--members", f"b{b}.parquet", "s0", "s1")
    parts = [pl.read_parquet(tmp_path / f"b{b}.parquet") for b in range(3)]
    reps = [set(p["cluster_rep"]) for p in parts]
    assert sum(map(len, reps)) == len(set().union(*reps))  # no cluster in two buckets
    everything = pl.concat(parts).sort("cluster_rep", "protein_id")
    assert everything.height == pl.read_parquet(RELEASE / "mgy_cluster_seqs.parquet").height

    everything.write_parquet(tmp_path / "all.parquet")
    for n in (2, 4):
        run(
            tmp_path,
            "subset",
            "--members",
            "all.parquet",
            "--sample",
            str(n),
            "--out-members",
            f"sub{n}.parquet",
        )
    sub2, sub4 = (pl.read_parquet(tmp_path / f"sub{n}.parquet") for n in (2, 4))
    assert (sub4["cluster_rep"] % 4 == 0).all() and 0 < sub4.height < sub2.height
    assert sub4.join(sub2, on="protein_id", how="anti").is_empty()

    for b in range(3):
        run(
            tmp_path,
            "stats",
            "--members",
            f"b{b}.parquet",
            "--prefix",
            f"b{b}",
            "--k",
            "6",
            "--rate",
            "0.5",
        )
    run(
        tmp_path,
        "combine",
        "b0",
        "b1",
        "b2",
        "--k",
        "6",
        "--sample",
        "1",
        "2",
        "--n-min",
        "4",
        "8",
        "--t-dense",
        "0",
        "0.1",
    )
    clusters = pl.read_parquet(tmp_path / "clusters.parquet")
    assert clusters.height == len(set().union(*reps))
    cost = pl.read_csv(tmp_path / "cost.tsv", separator="\t")
    assert cost.height == 2 * 2 * 2
    assert set(cost.filter(pl.col("sample") == 1)["units"]) == {clusters.height}
    groups = pl.read_csv(tmp_path / "groups.tsv", separator="\t")
    assert set(groups["sample"]) == {1, 2} and (groups["hashes"] > 0).all()


def test_ladder_nests_read_subsets(tmp_path: Path) -> None:
    data = ROOT / "tests" / "data"
    reads = ("--r1", str(data / "reads_1.fastq.gz"), "--r2", str(data / "reads_2.fastq.gz"))
    run(tmp_path, "ladder", *reads, "--pairs", "20", "5", "10")

    def names(n: int, mate: int) -> list[str]:
        with gzip.open(tmp_path / f"reads.{n}_{mate}.fastq.gz", "rt") as f:
            return [line.split("/")[0] for line in f.read().splitlines()[::4]]

    for small, large in [(5, 10), (10, 20)]:
        assert len(names(small, 1)) == small
        assert set(names(small, 1)) < set(names(large, 1))
    assert names(20, 1) == names(20, 2)


def test_em_starts_tells_unidentifiable_splits_from_unconverged_fits(tmp_path: Path) -> None:
    # Units 0 and 1 share a component and swap coverage between starts (their sum holds: an
    # unidentifiable split); unit 2 alone moves 50% (a fit that did not converge); unit 3
    # agrees; unit 4 is not kept by gather (kmers_unique 0) and is ignored.
    base = {"unit": [0, 1, 2, 3, 4], "component": [0, 0, 2, 3, 4], "kmers_unique": [2, 2, 5, 20, 0]}
    starts = {"i.10.0": [3.0, 1.0, 2.0, 4.0, 9.0], "i.10.0.s1": [1.0, 3.0, 3.0, 4.0, 1.0]}
    for name, cov in starts.items():
        pl.DataFrame({**base, "coverage_em": cov}).write_csv(
            tmp_path / f"{name}.profile.tsv", separator="\t"
        )
    run(tmp_path, "em-starts", "i.10.0.profile.tsv", "i.10.0.s1.profile.tsv")
    got = pl.read_csv(tmp_path / "em_starts.tsv", separator="\t").filter(stratum="all")
    unit = got.filter(level="unit").row(0, named=True)
    component = got.filter(level="component").row(0, named=True)
    assert (unit["n"], unit["within_1pct"]) == (4, 0.25)  # only unit 3 agrees
    assert component["n"] == 3 and component["within_1pct"] == pytest.approx(2 / 3, abs=1e-6)
    assert component["max_spread"] == pytest.approx(1 / 2.5)  # unit 2 alone: 2 vs 3, over 2.5
    # half L1 between the normalised profiles: (3,1,2,4)/10 against (1,3,3,4)/11
    want = 0.5 * sum(abs(a / 10 - b / 11) for a, b in zip([3, 1, 2, 4], [1, 3, 3, 4], strict=True))
    assert unit["l1_max"] == pytest.approx(want, abs=1e-6)


def test_query_cost_skips_nested_counts(tmp_path: Path) -> None:
    # The query's stats carry a hits_per_pair histogram (a list), which CSV cannot hold.
    stats = {"counts": {"hit_units": 3, "hits_per_pair": [0, 2, 1]}, "stages": {}}
    (tmp_path / "idx.40.0.json").write_text(json.dumps(stats))
    run(tmp_path, "query-cost", "idx.40.0.json")
    table = pl.read_csv(tmp_path / "query_cost.tsv", separator="\t")
    assert table["hit_units"].to_list() == [3]
    assert "hits_per_pair" not in table.columns


def test_membership_writes_one_file_per_shard(tmp_path: Path) -> None:
    run(tmp_path, "membership", "--release", str(RELEASE), "--shard-width", "300")
    out = tmp_path / "membership"
    assert not list(out.glob("shard=*"))  # the partitioned pieces are compacted away
    files = sorted(out.glob("shard_*.parquet"))
    assert len(files) > 1
    for f in files:  # each holds its own protein_id range, sorted
        ids = pl.read_parquet(f)["protein_id"]
        assert ids.is_sorted() and ids.min() // 300 == ids.max() // 300 == int(f.stem[6:])
    total = pl.read_parquet(RELEASE / "mgy_cluster_seqs.parquet").height
    assert sum(pl.read_parquet(f).height for f in files) == total


def test_extract_refuses_the_old_membership_layout(tmp_path: Path) -> None:
    (tmp_path / "membership" / "shard=0").mkdir(parents=True)
    pl.DataFrame({"cluster_rep": [1], "protein_id": [1]}).write_parquet(
        tmp_path / "membership" / "shard=0" / "data_0.parquet"
    )
    result = subprocess.run(
        [sys.executable, SCRIPT, "extract", *SHARD, "--lo", "0", "--hi", "10", "--prefix", "s"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0 and "old layout" in result.stderr


def test_io_errors_are_retried(monkeypatch: pytest.MonkeyPatch) -> None:
    import importlib.util

    import duckdb

    spec = importlib.util.spec_from_file_location("mgnify_subset", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod.time, "sleep", lambda s: None)
    calls = []

    def flaky(n: int) -> str:
        calls.append(n)
        if len(calls) < 3:
            raise duckdb.IOException("Cannot open file")
        return "ok"

    assert mod.with_io_retries(flaky, 7) == "ok" and calls == [7, 7, 7]

    def broken() -> None:
        calls.append(0)
        raise duckdb.IOException("Cannot open file")

    calls.clear()
    with pytest.raises(duckdb.IOException):
        mod.with_io_retries(broken)
    assert len(calls) == mod.IO_RETRIES
