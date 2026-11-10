"""The compact route index: exact k-mer codes, and a build whose keys find their units."""

import sys
from pathlib import Path

import numpy as np
import polars as pl
from test_tree_phase0 import run

from kmer_functional_profiler import _core

SCRIPT = Path(__file__).resolve().parents[2] / "workflows" / "tree-index" / "compact.py"
sys.path.insert(0, str(SCRIPT.parent))
import compact  # noqa: E402

AMINO = b"ACDEFGHIKLMNPQRSTVWY"


def test_codes_invert_the_hash() -> None:
    digit = {
        a: int(compact.codes(_core.hash_proteins([bytes([a])], 1)["hash"], 1)[0]) for a in AMINO
    }
    assert sorted(digit.values()) == list(range(20))
    rng = np.random.default_rng(0)
    kmers = [bytes(rng.choice(list(AMINO), 9).astype(np.uint8)) for _ in range(500)]
    got = _core.hash_proteins(kmers, 9)
    want = [sum(digit[kmers[i][j]] * 20 ** (8 - j) for j in range(9)) for i in got["seq"]]
    assert compact.codes(got["hash"], 9).tolist() == want


def test_build_keys_find_their_units(tmp_path: Path) -> None:
    rng = np.random.default_rng(1)
    rows = []
    for u in range(60):  # units of 1-5 members, each member a few substitutions off the rep
        rep = bytes(rng.choice(list(AMINO), 200).astype(np.uint8))
        for j in range(1 + u % 5):
            s = bytearray(rep)
            for i in rng.choice(200, 3 * j, replace=False):
                s[i] = AMINO[rng.integers(20)]
            rows.append((1000 * u + j, 1000 * u, s.decode()))
    members = tmp_path / "members.parquet"
    pl.DataFrame(
        rows, schema=["protein_id", "cluster_rep", "sequence"], orient="row"
    ).write_parquet(members)
    index = tmp_path / "index"
    for s in (0, 1):
        run("shard", "--members", str(members), "--index", str(index), "--shards", "2",
            "--shard", str(s), "--rate", "8", cwd=tmp_path, script=SCRIPT)  # fmt: skip
    run("merge", "--index", str(index), cwd=tmp_path, script=SCRIPT)
    ix = compact.Compact(index)
    reps = ix.reps(np.arange(ix.meta["units"]))
    assert len(reps) == 60 and all(pid % 1000 == 0 for pid, _ in reps.values())
    # every unit's representative k-mers hit that unit, and only keyed ones hit anything
    for u, (_, seq) in reps.items():
        h = _core.hash_proteins([seq], 9)["hash"]
        q, rows_ = ix.lookup(compact.codes(h, 9))
        assert u in set(np.asarray(ix.unit[rows_]).tolist())
    assert ix.meta["postings"] == len(ix.unit)


def test_recruit_calls_reads_to_their_proteins(tmp_path: Path) -> None:
    data = Path(__file__).resolve().parents[1] / "data"
    truth = pl.read_csv(data / "truth.tsv", separator="\t").filter(pl.col("mate") == 0)
    seeded = tmp_path / "seeded.tsv"
    first = truth.head(1)
    pl.DataFrame({"read": [f"pair{first['read'][0]}"], "target": first["protein"]}).write_csv(
        seeded, separator="\t"
    )
    run("--reads", str(data / "reads_1.fastq.gz"), str(data / "reads_2.fastq.gz"),
        "--targets", str(data / "proteins.faa"), "--seeded", str(seeded),
        "--out", str(tmp_path / "r"),
        cwd=tmp_path, script=SCRIPT.with_name("recruit.py"))  # fmt: skip
    got = pl.read_parquet(tmp_path / "r/reads.parquet")
    assert got.height and got["name"].str.starts_with("pair").all()
    top = got.filter(pl.col("weight") == pl.col("weight").max().over("read"))
    want = dict(zip(truth["read"].to_list(), truth["protein"].to_list(), strict=True))
    right = sum(want.get(r) == t for r, t in top.select("read", "target").iter_rows())
    assert right >= 0.8 * top["read"].n_unique()
    assert pl.read_csv(tmp_path / "r/targets.tsv", separator="\t").height == 20
