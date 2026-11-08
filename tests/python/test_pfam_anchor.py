"""HMM column mapping and column statistics for the Pfam-anchored unit assignment."""

import importlib.util
import sys
from pathlib import Path

import numpy as np

SCRIPT = Path(__file__).resolve().parents[2] / "workflows" / "pfam-anchor" / "columns.py"
spec = importlib.util.spec_from_file_location("columns", SCRIPT)
assert spec and spec.loader
columns = importlib.util.module_from_spec(spec)
sys.modules["columns"] = columns
spec.loader.exec_module(columns)


def code(text: str) -> list[int]:
    return [columns.AMINO.index(x) + 1 for x in text]


def test_column_residues_skips_insertions_and_leaves_deletions_empty() -> None:
    # Protein residues 2-9 (BCDEFGHI, as ACDEFGHI's letters) on columns 3-10: 2M, 1 inserted, 2M, 2
    # deleted, 3M. Columns: 3 C, 4 D, (E inserted), 5 F, 6 G, 7-8 none, 9 H, 10 I, 11 K.
    codes, positions = columns.column_residues(3, 2, "2M1I2M2D3M", "ACDEFGHIKL", 12)
    assert codes[3:12].tolist() == [*code("CDFG"), 0, 0, *code("HIK")]
    assert positions[3:12].tolist() == [1, 2, 4, 5, -1, -1, 6, 7, 8]
    assert codes[0] == 0
    assert codes[1:3].tolist() == [0, 0]
    assert codes[12] == 0


def test_whole_identity_uses_shared_columns_and_rejects_small_overlap() -> None:
    row = np.array([0, 1, 2, 3, 4], dtype=np.uint8)
    matrix = np.array(
        [
            [0, 1, 2, 3, 4],  # identical
            [0, 1, 2, 5, 6],  # half
            [0, 1, 0, 0, 0],  # shares one column of four
        ],
        dtype=np.uint8,
    )
    identity = columns.whole_identity(matrix, row)
    assert identity.tolist() == [1.0, 0.5, -1.0]


PREPARE = SCRIPT.with_name("prepare.py")
spec = importlib.util.spec_from_file_location("prepare", PREPARE)
assert spec and spec.loader
prepare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare)


def test_orfs_split_at_stops_on_both_strands() -> None:
    # + frame 0: MK*W -> "MK" then "W"; - strand is the reverse complement.
    seq = "ATGAAATAATGG"
    orfs = prepare.orfs_at(seq, 1)
    assert ("MK", "+", 0) in orfs
    assert ("W", "+", 9) in orfs
    rc = seq[::-1].translate(prepare.RC)  # CCATTATTTCAT: P L F H
    assert ("PLFH", "-", 0) in orfs
    assert all(len(p) >= 2 for p, _, _ in prepare.orfs_at(seq, 2))
    assert rc == "CCATTATTTCAT"
