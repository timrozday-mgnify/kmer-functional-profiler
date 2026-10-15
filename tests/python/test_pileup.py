"""kfp-pileup (plan, phase 13, steps 2-6, 8-9) on hand-made alignments and the fixture index."""

from pathlib import Path

import polars as pl
import pytest

from kfp_pileup.check import aai_check
from kfp_pileup.pileup import DIAMOND_BITS, blosum62, call
from kfp_pileup.reference import build_reference
from kmer_functional_profiler.index import (
    TARGETS,
    Index,
    IndexParams,
    build_index,
    members_from_fasta,
    read_fasta,
)
from kmer_functional_profiler.query import profile

DATA = Path(__file__).resolve().parents[1] / "data"
READS = (DATA / "reads_1.fastq.gz", DATA / "reads_2.fastq.gz")
A = "MKVLAAGIVLW"  # target, unit 0
B = "MKVLSAGIILW"  # neighbour, unit 1: S at column 4 and I at column 8
C = "MKVLSAGIVLW"  # second target, unit 2: S at column 4 only


def write_refs(path: Path, rows: list[tuple[str, str, int, str]]) -> Path:
    """refs.tsv and ref.fa of (seq_id, role, unit, sequence) rows."""
    path.mkdir()
    pl.DataFrame(
        [(s, r, u, f"u{u}", s, True, 0, len(q)) for s, r, u, q in rows],
        schema=[
            "seq_id",
            "role",
            "unit",
            "cluster_rep",
            "protein_id",
            "rep",
            "component",
            "length",
        ],
        orient="row",
    ).write_csv(path / "refs.tsv", separator="\t")
    (path / "ref.fa").write_text("".join(f">{s}\n{q}\n" for s, _, _, q in rows))
    return path


def diamond_line(read: str, seq: str, bits: float, btop: str) -> str:
    return f"{read}\t{seq}\t{bits}\t1\t33\t1\t11\t1\t{btop}\n"


@pytest.fixture
def protein_refs(tmp_path: Path) -> Path:
    return write_refs(
        tmp_path / "ref",
        [("s0", "target", 0, A), ("s1", "neighbour", 1, B), ("s2", "decoy", 0, A[::-1])],
    )


def test_posteriors_decoys_and_sibling_residues(protein_refs: Path, tmp_path: Path) -> None:
    hits = tmp_path / "r1.tsv"
    hits.write_text(
        diamond_line("r1", "s0", 30, "11")  # A as is: weight 1 / (1 + 2^-3) on A
        + diamond_line("r1", "s1", 27, "4AS3VI2")
        + diamond_line("r2", "s0", 27, "4SA6")  # B as is: best is the neighbour, dropped
        + diamond_line("r2", "s1", 30, "11")
        + diamond_line("r3", "s0", 28, "8LV2")  # a variant (L) no sibling carries
        + diamond_line("r3", "s1", 25, "4AS3LI2")
        + diamond_line("r4", "s0", 29, "4SA6")  # S at 4, as B has: a sibling's residue
        + diamond_line("r4", "s1", 28, "8VI2")
        + diamond_line("r5", "s2", 40, "11")  # best on a reversed decoy: dropped
    )
    call([hits], protein_refs, tmp_path / "out", "diamond", min_depth=0.5, min_freq=0.1)
    assigned = pl.read_parquet(tmp_path / "out" / "assignments.parquet")
    kept = assigned.filter("kept").group_by("read_name").len()
    assert sorted(kept["read_name"]) == ["r1", "r3", "r4"]
    r1 = assigned.filter(pl.col("read_name") == "r1", pl.col("seq_id") == "s0")["posterior"][0]
    assert r1 == pytest.approx(1 / (1 + 2**-3))
    variants = pl.read_csv(tmp_path / "out" / "variants.tsv", separator="\t")
    at = {row["pos"]: row for row in variants.iter_rows(named=True)}
    assert set(at) == {5, 9}  # 1-based: column 4 (S) and column 8 (L)
    assert (at[5]["alt"], at[5]["sibling_share"]) == ("S", 1.0)
    assert (at[9]["alt"], at[9]["sibling_share"]) == ("L", 0.0)
    units = pl.read_csv(tmp_path / "out" / "units.tsv", separator="\t")
    assert units["unit"].to_list() == [0]
    expected = 1 / (1 + 2**-3) + 1 / (1 + 2**-3) + 1 / (1 + 2**-1)
    assert units["reads"][0] == pytest.approx(expected)


def test_leave_one_column_out_weight(tmp_path: Path) -> None:
    refs = write_refs(tmp_path / "ref", [("s0", "target", 0, A), ("s1", "target", 2, C)])
    hits = tmp_path / "r1.tsv"
    # One read, A's residue at column 4, equal scores on both targets.
    hits.write_text(diamond_line("r1", "s0", 30, "11") + diamond_line("r1", "s1", 30, "4AS6"))
    call([hits], refs, tmp_path / "out", "diamond", min_depth=0, min_freq=0)
    pile = pl.read_parquet(tmp_path / "out" / "pileup.parquet").filter(pl.col("col") == 4)
    m = blosum62()
    t_a, t_c = (m[ord("A"), ord(x)] * DIAMOND_BITS for x in "AS")
    expected = 2**-t_a / (2**-t_a + 2**-t_c)  # the column's own evidence removed
    weights = dict(zip(pile["seq_id"], pile["weight"], strict=True))
    assert weights["s0"] == pytest.approx(expected)
    assert weights["s0"] + weights["s1"] == pytest.approx(1.0)
    assert weights["s0"] < 0.5  # without the column, the read leans to neither as much


def test_sam_strands_and_borrowed_sequences(tmp_path: Path) -> None:
    ref = "ACGTTGCAAGGCTTACGATCGATCGGATCCA"
    refs = write_refs(tmp_path / "ref", [("s0", "target", 0, ref), ("s1", "decoy", 0, ref[::-1])])
    read = ref[2:22]
    variant = read[:9] + ("A" if read[9] != "A" else "C") + read[10:]  # alt at column 11
    rc = variant.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    sam = tmp_path / "aln.sam"
    sam.write_text(
        "@SQ\tSN:s0\tLN:31\n"
        f"p1\t65\ts0\t3\t60\t20M\t=\t3\t0\t{variant}\t*\tAS:i:15\n"
        f"p1\t129\ts0\t3\t60\t20M\t=\t3\t0\t{read}\t*\tAS:i:20\n"
        f"p2\t16\ts0\t3\t60\t20M\t*\t0\t0\t{variant}\t*\tAS:i:15\n"  # reverse strand, as stored
        f"p2\t256\ts1\t5\t0\t20M\t*\t0\t0\t*\t*\tAS:i:2\n"  # secondary, its sequence borrowed
    )
    assert rc  # the read p2 came from, before SAM stored it on the reference strand
    call([sam], refs, tmp_path / "out", "bwa", min_depth=1, min_freq=0.1)
    variants = pl.read_csv(tmp_path / "out" / "variants.tsv", separator="\t")
    assert variants["pos"].to_list() == [12]
    assert variants["depth"][0] == pytest.approx(3, abs=1e-3)
    assert variants["alt_freq"][0] == pytest.approx(2 / 3, abs=1e-3)
    assert pl.read_parquet(tmp_path / "out" / "assignments.parquet").height == 4


def test_reference_from_a_fixture_profile(tmp_path: Path) -> None:
    index_dir = tmp_path / "idx"
    index_dir.mkdir()
    members_from_fasta(DATA / "proteins.faa").write_parquet(index_dir / TARGETS)
    build_index(index_dir / TARGETS, index_dir, IndexParams(k=7, t_base=1.0))
    profile(Index.load(index_dir), *READS).write_csv(tmp_path / "p.tsv", separator="\t")
    counts = build_reference(tmp_path / "p.tsv", index_dir, tmp_path / "ref")
    assert counts["target"] == counts["decoy"] > 0
    refs = pl.read_csv(tmp_path / "ref" / "refs.tsv", separator="\t")
    fasta = dict(read_fasta(tmp_path / "ref" / "ref.fa"))
    assert list(fasta) == refs["seq_id"].to_list()
    target = refs.filter(pl.col("role") == "target").row(0, named=True)
    decoy = refs.filter(pl.col("role") == "decoy", pl.col("unit") == target["unit"])
    assert fasta[decoy["seq_id"][0]] == fasta[target["seq_id"]][::-1]
    assert (tmp_path / "ref" / "molecule.txt").read_text() == "protein\n"
    # a withheld unit is in no role
    held = str(target["cluster_rep"])
    build_reference(tmp_path / "p.tsv", index_dir, tmp_path / "held", withhold=[held])
    kept = pl.read_csv(tmp_path / "held" / "refs.tsv", separator="\t", infer_schema=False)
    assert held not in kept["cluster_rep"].to_list()
    assert kept.height == refs.filter(pl.col("unit") != target["unit"]).height


def test_aai_check_bins_and_flags(tmp_path: Path) -> None:
    # 120 units at alignment identity 0.96 whose aai runs 0.03 high: that bin is flagged;
    # 5 at 0.995 with no gap are too few to flag; units with few reads are left out.
    rows = [(0.96, 0.99, 0.98, 1.0, 50.0, 300.0, 3)] * 120
    rows += [(0.995, 0.995, 0.99, 1.0, 50.0, 300.0, 3)] * 5 + [(0.5, 0.99, 0.9, 1.0, 2.0, 300.0, 3)]
    pl.DataFrame(
        rows,
        schema=["identity", "aai", "aai_lo", "aai_hi", "reads", "len_mean", "n_members"],
        orient="row",
    ).with_columns(coverage_zi=pl.col("reads")).write_csv(tmp_path / "u.tsv", separator="\t")
    table = aai_check(tmp_path / "u.tsv")
    by = {(r["breakdown"], r["bin"]): r for r in table.iter_rows(named=True)}
    assert by[("all", "all")]["units"] == 125
    assert by[("identity", "[0.95, 0.97)")]["flagged"]
    assert by[("identity", "[0.95, 0.97)")]["median_gap"] == pytest.approx(0.03)
    assert not by[("identity", "[0.99, inf)")]["flagged"]
    assert by[("identity", "[0.95, 0.97)")]["cover"] == 0.0


def test_em_and_haplotypes_send_a_relatives_reads_to_novel(tmp_path: Path) -> None:
    # Target A (30 aa); a strain of it with one fixed difference (column 2) and a relative
    # missing from the reference that differs at five more columns, linked (6/30).
    ref = "MKVLAAGIVLWQERTYIPASDFGHKLCVNM"
    refs = write_refs(tmp_path / "ref", [("s0", "target", 0, ref), ("s1", "decoy", 0, ref[::-1])])
    relative = {2: "W", 5: "S", 10: "D", 15: "N", 20: "Q", 25: "E"}

    def btop(diffs: dict[int, str]) -> str:
        out, run = "", 0
        for i, s in enumerate(ref):
            if i in diffs:
                out += f"{run or ''}{diffs[i]}{s}"
                run = 0
            else:
                run += 1
        return out + str(run)

    strain = [diamond_line(f"st{i}", "s0", 60, btop({2: "W"})) for i in range(12)]
    rel = [diamond_line(f"rel{i}", "s0", 40, btop(relative)) for i in range(6)]

    def run(lines: list[str], name: str) -> tuple[dict[str, float], Path]:
        hits, out = tmp_path / f"{name}.tsv", tmp_path / name
        hits.write_text("".join(lines))
        call([hits], refs, out, "diamond", min_depth=3, min_freq=0.2)
        a = pl.read_parquet(out / "assignments.parquet").filter(pl.col("role") == "target")
        return dict(zip(a["read_name"], a["em_novel"], strict=True)), out

    # beside the strain, the EM's novel class takes the relative's reads
    novel, out = run(strain + rel, "mixed")
    assert all(novel[f"rel{i}"] > 0.9 for i in range(6))
    assert all(novel[f"st{i}"] < 0.1 for i in range(12))
    # alone, the consensus moves to the relative and its haplotype is reported as novel
    novel, out = run(rel, "alone")
    assert all(novel[f"rel{i}"] > 0.9 for i in range(6))
    haps = pl.read_csv(out / "haplotypes.tsv", separator="\t")
    assert haps.filter("novel")["divergence"].to_list() == [pytest.approx(6 / 30)]
