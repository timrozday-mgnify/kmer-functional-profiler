"""lookup-index's tables, on hmm-annotate's own output for the mini release (mini_annotation)."""

import shutil
import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "workflows" / "lookup-index" / "lookup_index.py"
DATA = ROOT / "tests" / "data"
PFAM = DATA / "mini_annotation" / "pfam"
KO = DATA / "mini_annotation" / "kofam"
HMM = DATA / "mini_fmh" / "Pfam-mini.hmm"
RELEASE = DATA / "mini_reannotate"


def run(cwd: Path, *args: object) -> None:
    subprocess.run([sys.executable, SCRIPT, *map(str, args)], cwd=cwd, check=True)


def groups(cwd: Path, pfam: Path = PFAM, hmm: Path = HMM, n: int = 2) -> pl.DataFrame:
    run(cwd, "groups", "--hmm", hmm, "--pfam", pfam, "--label", f"ko={KO}", "--groups", n,
        "--batch-chunks", 1)  # fmt: skip
    return pl.read_csv(cwd / "groups.tsv", separator="\t")


def domains() -> pl.DataFrame:
    return pl.read_parquet(PFAM / "domains" / "*.parquet")


def test_groups_hold_every_family_once_and_batches_every_chunk(tmp_path: Path) -> None:
    g = groups(tmp_path, n=4)
    expected = domains().filter("passes").group_by("model").len("rows").sort("model")
    assert g.sort("model").select("model", "rows").rows() == expected.rows()
    assert set(g["group"]) == {0, 1, 2, 3}
    b = pl.read_csv(tmp_path / "batches.tsv", separator="\t")
    r = pl.read_csv(PFAM / "ranges.tsv", separator="\t")
    assert b.select("lo", "hi").rows() == r.select("lo", "hi").rows()


def test_groups_fail_on_a_missing_chunk_or_another_hmm_file(tmp_path: Path) -> None:
    shutil.copytree(PFAM, tmp_path / "pfam")
    (tmp_path / "pfam/domains/chunk2.domains.parquet").unlink()
    with pytest.raises(subprocess.CalledProcessError):
        groups(tmp_path, pfam=tmp_path / "pfam")
    (tmp_path / "other.hmm").write_text(HMM.read_text() + "\n")
    with pytest.raises(subprocess.CalledProcessError):
        groups(tmp_path, hmm=tmp_path / "other.hmm")


def two_sets(cwd: Path) -> Path:
    """The annotation as a run of representatives and then members: odd ids, in wider ranges."""
    out = cwd / "pfam2"
    shutil.copytree(PFAM, out)
    r = pl.read_csv(out / "ranges.tsv", separator="\t")
    top = r["hi"].max()
    with open(out / "ranges.members.tsv", "w") as f:
        f.write("chunk\tlo\thi\n")
        for i, (lo, hi) in enumerate([(0, top // 2), (top // 2, top)]):
            f.write(f"m{i}\t{lo}\t{hi}\n")
            d = pl.read_parquet(PFAM / "domains" / "*.parquet")
            d = d.filter(pl.col("protein_id") % 2 == 1, pl.col("protein_id").is_between(lo, hi - 1))
            d.write_parquet(out / "domains" / f"chunkm{i}.domains.parquet")
    for c in r["chunk"]:
        path = out / "domains" / f"chunk{c}.domains.parquet"
        pl.read_parquet(path).filter(pl.col("protein_id") % 2 == 0).write_parquet(path)
    return out


@pytest.mark.parametrize("sets", [1, 2])
def test_split_gives_every_domain_its_cluster_labels_and_group(tmp_path: Path, sets: int) -> None:
    pfam = two_sets(tmp_path) if sets == 2 else PFAM
    g = groups(tmp_path, pfam=pfam)
    clusters = pl.read_parquet(RELEASE / "mgy_cluster_seqs.parquet").rename(
        {"cluster_member": "protein_id"}
    )
    (tmp_path / "membership").mkdir()
    clusters.write_parquet(tmp_path / "membership/shard_000000.parquet")
    got = []
    for batch, lo, hi in pl.read_csv(tmp_path / "batches.tsv", separator="\t").iter_rows():
        run(tmp_path, "split", "--pfam", pfam, "--label", f"ko={KO}", "--groups", "groups.tsv",
            "--membership", "membership", "--lo", lo, "--hi", hi, "--out", f"b{batch}")  # fmt: skip
        for f in sorted((tmp_path / f"b{batch}").glob("group=*/*.parquet")):
            group = int(f.parent.name.removeprefix("group="))
            got.append(
                pl.read_parquet(f, hive_partitioning=False).with_columns(group=pl.lit(group))
            )
    labels = (
        pl.read_parquet(KO / "labels" / "*.parquet")
        .group_by("protein_id")
        .agg(label_ko=pl.col("model").unique().sort().str.join("+"))
    )
    expected = (
        domains()
        .filter("passes")
        .join(clusters, on="protein_id")
        .join(labels, on="protein_id", how="left")
        .join(g, on="model")
        .with_columns(pfam_accession=pl.col("model").str.slice(2).cast(pl.Int32))
    )
    cols = ["protein_id", "pfam_accession", "hmm_from", "cigar", "left", "right", "residues",
            "cluster_rep", "label_ko", "group"]  # fmt: skip
    rows = sorted(pl.concat(got).select(cols).rows(), key=str)
    assert rows == sorted(expected.select(cols).rows(), key=str)
    assert len(rows) == 18 and any("+" in (r[8] or "") for r in rows)


def test_split_fails_on_a_protein_in_no_cluster(tmp_path: Path) -> None:
    groups(tmp_path)
    (tmp_path / "membership").mkdir()
    pl.read_parquet(RELEASE / "mgy_cluster_seqs.parquet").rename(
        {"cluster_member": "protein_id"}
    ).filter(pl.col("protein_id") != 1).write_parquet(tmp_path / "membership/shard_000000.parquet")
    with pytest.raises(subprocess.CalledProcessError):
        run(tmp_path, "split", "--pfam", PFAM, "--groups", "groups.tsv", "--membership",
            "membership", "--lo", 0, "--hi", 30)  # fmt: skip
