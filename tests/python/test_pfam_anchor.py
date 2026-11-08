"""HMM column mapping and column statistics for the Pfam-anchored unit assignment."""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import polars as pl

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


sys.path.insert(0, str(SCRIPT.parent))
spec = importlib.util.spec_from_file_location("anchor", SCRIPT.with_name("anchor.py"))
assert spec and spec.loader
anchor = importlib.util.module_from_spec(spec)
spec.loader.exec_module(anchor)


def test_diamond_residues_walk_cigar_on_both_strands() -> None:
    # + frame: codons from base 4; 2 matched, 1 query-only, 1 subject-only, 1 matched.
    assert anchor.diamond_residues(4, 1, 10, "2M1I1D1M") == {4: 9, 7: 10, 13: 12}
    # - frame: qstart is the highest base, codons step down.
    assert anchor.diamond_residues(150, -1, 1, "2M") == {150: 0, 147: 1}


def test_family_score_takes_each_units_best_domain() -> None:
    domains = pl.DataFrame(
        {
            "unit": [7, 3, 7],
            "hmm_from": [1, 1, 1],
            "ali_from": [1, 1, 1],
            "cigar": ["4M", "4M", "2M2D"],
            "sequence": ["ACDE", "ACDF", "AC"],
        }
    )
    fam = anchor.Family(4, domains)
    assert fam.units.tolist() == [3, 7]
    codes = columns.column_residues(1, 1, "4M", "ACDE", 4)[0]
    m, s, n = fam.score(codes, "shared")
    assert n == 4
    # unit 3: 3 of 4; unit 7: its full domain (4 of 4) beats the half one (2 of 2 shares too few)
    assert m.tolist() == [3, 4]
    assert s.tolist() == [4, 4]


def test_seed_index_counts_matching_seeds_and_samples_by_first_column() -> None:
    domains = pl.DataFrame(
        {
            "unit": [1, 2, 3],
            "hmm_from": [1, 1, 1],
            "ali_from": [1, 1, 1],
            "cigar": ["6M", "6M", "3M3D"],
            "sequence": ["ACDEFG", "ACDEFW", "ACD"],
        }
    )
    fam = anchor.Family(6, domains)
    seeds = np.arange(1, 5)[:, None] + np.arange(3)  # contiguous k 3 from columns 1-4
    index = anchor.SeedIndex(fam, seeds, anchor.alphabet(None))
    # postings: unit 1 and 2 have 4 seeds each, unit 3 one (ACD only)
    assert len(index.post) == 4 + 4 + 1
    # matches: unit 2's EFW differs from the query's EFG
    query = columns.column_residues(1, 1, "6M", "ACDEFG", 6)[0][None, :]
    assert index.counts(index.post, query, 3).tolist() == [4, 3, 1]
    # s 2 keeps seeds starting at column 2 and 4
    assert index.counts(index.sampled(2), query, 3).tolist() == [2, 1, 0]
    # Murphy-10 groups residues, but G and W stay in different classes: the same counts
    murphy = anchor.SeedIndex(fam, seeds, anchor.alphabet(anchor.MURPHY10))
    assert murphy.counts(murphy.post, query, 3).tolist() == [4, 3, 1]
