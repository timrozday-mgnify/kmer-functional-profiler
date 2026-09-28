"""The MGnify subset extraction on the mini release fixture."""

import subprocess
import sys
from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "workflows" / "mgnify-subset" / "mgnify_subset.py"
RELEASE = ROOT / "tests" / "data" / "mini_release"
GUT = "root:Host-associated:Human:Digestive system"
SHARD = ("--release", str(RELEASE), "--membership", "membership.parquet")


def run(cwd: Path, *args: str) -> None:
    subprocess.run([sys.executable, SCRIPT, *args], cwd=cwd, check=True)


def test_extracts_whole_gut_clusters(tmp_path: Path) -> None:
    run(tmp_path, "membership", "--release", str(RELEASE), "--biome", GUT)
    for i, (lo, hi) in enumerate([(0, 400), (400, 1000)]):
        run(tmp_path, "extract", *SHARD, "--lo", str(lo), "--hi", str(hi), "--prefix", f"s{i}")
    run(tmp_path, "merge", "--membership", "membership.parquet", "s0", "s1")

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


def test_merge_fails_on_missing_members(tmp_path: Path) -> None:
    run(tmp_path, "membership", "--release", str(RELEASE))
    run(tmp_path, "extract", *SHARD, "--lo", "0", "--hi", "500", "--prefix", "s0")
    result = subprocess.run(
        [sys.executable, SCRIPT, "merge", "--membership", "membership.parquet", "s0"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0 and "no sequence" in result.stderr
