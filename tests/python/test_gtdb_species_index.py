"""The GTDB species index recipe's genome pick (phase 11, step 11)."""

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("gtdb", ROOT / "workflows/gtdb-species-index/gtdb.py")
assert SPEC is not None and SPEC.loader is not None
gtdb = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gtdb)
FIXTURE = ROOT / "tests/data/mini_gtdb"


def test_ncbi_url() -> None:
    # real assembly (GTDB r226): spaces become _, dashes stay
    assert gtdb.ncbi_url("GCA_903826505.1", "freshwater MAG --- MJ111018B_bin-0298") == (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/903/826/505/"
        "GCA_903826505.1_freshwater_MAG_---_MJ111018B_bin-0298/"
        "GCA_903826505.1_freshwater_MAG_---_MJ111018B_bin-0298_genomic.fna.gz"
    )


def test_pick_keeps_representative_then_best_per_species() -> None:
    base = f"file://{FIXTURE}/ncbi"
    picked = gtdb.pick([str(FIXTURE / "bac120_metadata_r0.tsv")], 2, 50, 10, base=base)
    # contamination 20 drops 8; alpha: rep 1, then 2 (100%) over 3 (80%); beta: rep 4, then
    # 5 over 6 by accession
    assert picked["genome"].to_list() == [f"GCA_00000000{n}.1" for n in (1, 2, 4, 5, 7)]
    assert picked["representative"].to_list() == [True, False, True, False, True]
    assert all(Path(u.removeprefix("file://")).exists() for u in picked["url"])
    one = gtdb.pick([str(FIXTURE / "bac120_metadata_r0.tsv")], 5, 50, 10, {"s__G2 gamma"})
    assert one["genome"].to_list() == ["GCA_000000007.1"]
    assert gtdb.pick([str(FIXTURE / "bac120_metadata_r0.tsv")], 5, 50, 10,
                     max_species=2)["taxonomy"].n_unique() == 2  # fmt: skip
