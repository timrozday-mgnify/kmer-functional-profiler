"""Tree index: the family-complete slice (slice.py) and verify.py's repeat masking."""

import sys
from pathlib import Path

import polars as pl
from test_tree_phase0 import SCRIPT, run

SLICE = SCRIPT.with_name("slice.py")
sys.path.insert(0, str(SCRIPT.parent))
import verify  # noqa: E402


def test_slice_keeps_whole_clusters_of_chosen_families(tmp_path: Path) -> None:
    """Families within the size range are chosen; every member of a cluster holding one of
    their proteins is kept, members of other families' clusters are not."""
    # clusters 1 (proteins 1-3, PF1 on 1), 2 (4-5, PF2 on 4), 3 (6, PF3: too small)
    members = pl.DataFrame(
        {
            "protein_id": [1, 2, 3, 4, 5, 6],
            "cluster_rep": [1, 1, 1, 4, 4, 6],
            "full_length": [True] * 6,
            "sequence": ["MKV"] * 6,
        }
    )
    members.write_parquet(tmp_path / "members.bucket0.parquet")
    pfam = pl.DataFrame(
        {"protein_id": [1, 2, 4, 5, 6], "pfam_accession": ["PF1", "PF1", "PF2", "PF2", "PF3"]}
    )
    pfam.write_parquet(tmp_path / "pfam.parquet")
    fam = tmp_path / "fam"
    run("families", "--pfam", str(tmp_path / "pfam.parquet"), "--families", "1",
        "--min-proteins", "2", "--out", str(fam), cwd=tmp_path, script=SLICE)  # fmt: skip
    chosen = pl.read_csv(fam / "families.tsv", separator="\t")["pfam_accession"].to_list()
    assert len(chosen) == 1 and chosen[0] in ("PF1", "PF2")  # PF3 has one protein
    out = tmp_path / "take.parquet"
    run("take", "--members", str(tmp_path / "members.bucket0.parquet"), "--chosen",
        str(fam / "chosen.parquet"), "--out", str(out), cwd=tmp_path, script=SLICE)  # fmt: skip
    run("concat", "--parts", str(out), "--out", str(tmp_path / "m.parquet"), cwd=tmp_path,
        script=SLICE)  # fmt: skip
    got = pl.read_parquet(tmp_path / "m.parquet")["protein_id"].to_list()
    assert got == ([1, 2, 3] if chosen[0] == "PF1" else [4, 5])


def test_mask_hides_tandem_repeats_only() -> None:
    assert verify.mask(b"MKVLTHTHTHTHAQ") == b"MKVLxxxxxxxxAQ"
    assert verify.mask(b"MKVAAAAGCVCVCVCVQ") == b"MKVxxxxGxxxxxxxxQ"
    assert verify.mask(b"GCVCVCVCQ") == b"GCVCVCVCQ"  # three copies stay
    assert verify.mask(b"MKVLATGHEW") == b"MKVLATGHEW"
