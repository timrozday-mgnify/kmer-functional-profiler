"""kfp-pileup (plan, phase 13, steps 2-6, 8-9) on hand-made alignments and the fixture index."""

import json
from pathlib import Path

import polars as pl
import pytest

from kfp_pileup.check import aai_check
from kfp_pileup.handoff import handoff
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
        [("s0", "target", 0, A), ("s1", "neighbour", 1, B)],
    )


def test_posteriors_neighbours_and_sibling_residues(protein_refs: Path, tmp_path: Path) -> None:
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
    )
    profile = tmp_path / "profile.tsv"  # a unit outside the reference and an unused column
    profile.write_text("unit\tother\taai\tcoverage_zi\n0\tx\t0.97\t2.0\n9\ty\t0.5\t1.0\n")
    stats = tmp_path / "stats.json"
    call([hits], protein_refs, tmp_path / "out", "diamond", profile, 0.5, 0.1, stats=stats)
    assert set(json.loads(stats.read_text())["stages"]) == {"read", "batches", "refine", "write"}
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
    assert (units["aai"][0], units["coverage_zi"][0], "other" in units.columns) == (
        0.97,
        2.0,
        False,
    )
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
    refs = write_refs(
        tmp_path / "ref", [("s0", "target", 0, ref), ("s1", "neighbour", 1, ref[::-1])]
    )
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
    assert counts["target"] > 0 and set(counts) <= {"target", "neighbour"}  # no decoys
    refs = pl.read_csv(tmp_path / "ref" / "refs.tsv", separator="\t")
    prof = pl.read_csv(tmp_path / "p.tsv", separator="\t")  # targets: present_q <= 0.05
    want = set(prof.filter(pl.col("present_q") <= 0.05)["unit"])
    assert set(refs.filter(pl.col("role") == "target")["unit"]) == want
    # --max-targets: the most hit of them; off (0) by default, as above
    assert len(want) > 2
    build_reference(tmp_path / "p.tsv", index_dir, tmp_path / "cap", max_targets=2)
    capped = pl.read_csv(tmp_path / "cap" / "refs.tsv", separator="\t")
    top2 = prof.filter(pl.col("present_q") <= 0.05).sort(["hits", "unit"], descending=[True, False])
    assert set(capped.filter(pl.col("role") == "target")["unit"]) == set(top2["unit"][:2])
    fasta = dict(read_fasta(tmp_path / "ref" / "ref.fa"))
    assert list(fasta) == refs["seq_id"].to_list()
    target = refs.filter(pl.col("role") == "target").row(0, named=True)
    assert (tmp_path / "ref" / "molecule.txt").read_text() == "protein\n"
    # a withheld unit is in no role
    held = str(target["cluster_rep"])
    build_reference(tmp_path / "p.tsv", index_dir, tmp_path / "held", withhold=[held])
    kept = pl.read_csv(tmp_path / "held" / "refs.tsv", separator="\t", infer_schema=False)
    assert held not in kept["cluster_rep"].to_list()
    assert kept.height == refs.filter(pl.col("unit") != target["unit"]).height
    # neighbours: one target in one component of every unit, at most 2 neighbours, the most hit
    # present_q dropped: the reference computes it from present_prob, as for older profiles
    hit = pl.read_csv(tmp_path / "p.tsv", separator="\t").drop("present_q")
    hit = hit.with_columns(component=pl.lit(0))
    top = hit.sort("present_prob", descending=True)["unit"][0]
    hit = hit.with_columns(present_prob=(pl.col("unit") == top).cast(pl.Float64))
    hit.write_csv(tmp_path / "one.tsv", separator="\t")
    assert hit.height > 3
    assert "neighbour" not in build_reference(
        tmp_path / "one.tsv", index_dir, tmp_path / "n0", neighbours=0
    )
    build_reference(tmp_path / "one.tsv", index_dir, tmp_path / "n2", neighbours=2)
    kept = pl.read_csv(tmp_path / "n2" / "refs.tsv", separator="\t")
    kept = kept.filter(pl.col("role") == "neighbour")["unit"].unique()
    others = hit.filter(pl.col("unit") != top)
    assert kept.len() == 2
    assert others.filter(pl.col("unit").is_in(kept.implode()))["hits"].min() >= (
        others.filter(~pl.col("unit").is_in(kept.implode()))["hits"].max()
    )


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
    refs = write_refs(tmp_path / "ref", [("s0", "target", 0, ref)])
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


def test_spill_partitions_do_not_change_results(tmp_path: Path) -> None:
    # Two targets, each with a strain (a fixed difference: a consensus update) and a novel
    # relative's reads: the same outputs with every column in one partition or split in two.
    ref = "MKVLAAGIVLWQERTYIPASDFGHKLCVNM"
    other = ref[::-1]
    refs = write_refs(tmp_path / "ref", [("s0", "target", 0, ref), ("s1", "target", 1, other)])

    def btop(seq: str, diffs: dict[int, str]) -> str:
        out, run = "", 0
        for i, s in enumerate(seq):
            if i in diffs:
                out += f"{run or ''}{diffs[i]}{s}"
                run = 0
            else:
                run += 1
        return out + str(run)

    rel = {2: "W", 5: "S", 10: "D", 15: "N", 20: "Q", 25: "E"}
    lines = [diamond_line(f"a{i}", "s0", 60, btop(ref, {2: "W"})) for i in range(8)]
    lines += [diamond_line(f"b{i}", "s1", 60, btop(other, {7: "K"})) for i in range(8)]
    lines += [diamond_line(f"r{i}", "s0", 40, btop(ref, rel)) for i in range(4)]
    hits = tmp_path / "r1.tsv"
    hits.write_text("".join(lines))
    outs = {}
    for parts in (1, 2):
        out, stats = tmp_path / f"p{parts}", tmp_path / f"p{parts}.json"
        call(
            [hits], refs, out, "diamond", min_depth=3, min_freq=0.2, stats=stats, spill_parts=parts
        )
        assert json.loads(stats.read_text())["counts"]["spill_parts"] == parts
        assert not [p for p in out.iterdir() if p.name.startswith(".spill")]  # cleaned up
        outs[parts] = out
    for name in ("assignments.parquet", "pileup.parquet"):
        one, two = (pl.read_parquet(outs[p] / name) for p in (1, 2))
        assert one.sort(one.columns).equals(two.sort(two.columns))
    for name in ("sources.tsv", "haplotypes.tsv", "units.tsv", "variants.tsv"):
        assert (outs[1] / name).read_text() == (outs[2] / name).read_text(), name
    haps = pl.read_csv(outs[1] / "haplotypes.tsv", separator="\t")
    assert set(haps["seq_id"]) == {"s0", "s1"}


def test_frameshift_btop_to_nucleotides(tmp_path: Path) -> None:
    # Real DIAMOND 2.2.8 -F 15 alignments (tests/data/frameshift) of 150-nt reads from
    # cds.fna[300:450] with planted indels near read position 75 (cds 375), both strands.
    fs = DATA / "frameshift"
    cds = dict(read_fasta(fs / "cds.fna"))["p"]
    prot = dict(read_fasta(fs / "ref.faa"))["p"]
    refs = write_refs(tmp_path / "ref", [("p", "target", 0, prot)])
    (tmp_path / "r2.tsv").write_text("")
    counts = call(
        [fs / "diamond.tsv", tmp_path / "r2.tsv"], refs, tmp_path / "out", "diamond",
        min_depth=1, min_freq=0.1, reads=[fs / "r.fq", fs / "r.fq"],
    )  # fmt: skip
    assert counts["codon_mismatch"] == 0  # every read codon translates to its residue
    assert counts["frameshift_alignments"] == 7  # all but none and codon_del
    assigned = pl.read_parquet(tmp_path / "out" / "assignments.parquet")
    assert assigned.filter(pl.col("read_name") == "codon_del")["posterior"][0] > 0.99
    pile = pl.read_parquet(tmp_path / "out" / "pileup_nt.parquet")
    indels = {
        (r["col"], r["allele"]): r["weight"]
        for r in pile.filter(~pl.col("allele").is_in(list("ACGT"))).iter_rows(named=True)
    }
    # a base missing from the read: the next column's first base (DIAMOND puts the shift
    # after the codon it garbles); an extra base: inserted before the next codon column
    assert indels[(378, "-")] == pytest.approx(3, abs=0.01)  # del1, del1_rev, off1_del1
    assert indels[(375, "+G")] == pytest.approx(1, abs=0.01)  # ins1 (its planted G)
    assert indels[(375, f"+{cds[375]}")] == pytest.approx(1, abs=0.01)  # ins1_rev
    # codon_del; and del2 (2 nt): a codon out, and the base DIAMOND skipped back in
    assert indels[(375, "-")] == pytest.approx(2, abs=0.01)
    assert indels[(378, f"+{cds[377]}")] == pytest.approx(1, abs=0.01)
    assert indels[(435, "-")] == pytest.approx(1, abs=0.01)  # del1_end: exactly cds 435
    # the consensus is the CDS wherever reads cover it, away from the garbled codons
    cons = dict(read_fasta(tmp_path / "out" / "consensus.fna"))["p"]
    assert len(cons) == 3 * len(prot)
    assert set(cons[:300]) == {"N"}
    covered = [i for i in range(300, 450) if not 372 <= i < 381 and not 435 <= i < 439]
    assert all(cons[i] == cds[i] for i in covered)
    assert "CDS\t1\t900" in (tmp_path / "out" / "consensus.gff").read_text()
    recruited = (tmp_path / "out" / "recruited_1.fq.gz").read_bytes()
    assert recruited  # every read pair has a target alignment


def test_handoff_keeps_assigned_mates(tmp_path: Path) -> None:
    refs = write_refs(
        tmp_path / "ref",
        [("s0", "target", 0, "A" * 30), ("s1", "target", 1, "C" * 30), ("s2", "neighbour", 2, "A")],
    )
    pile = tmp_path / "pileup"
    pile.mkdir()
    pl.DataFrame(
        {
            "read_name": ["p1", "p1", "p1", "p2", "p3"],
            "mate": [0, 1, 0, 0, 0],
            "unit": [0, 0, 1, 0, 2],
            "role": ["target"] * 4 + ["neighbour"],
            "em_posterior": [0.9, 0.9, 0.1, 0.4, 0.9],
        }
    ).write_parquet(pile / "assignments.parquet")
    seq, qual = "ACGTA", "IIIII"
    sam = tmp_path / "in.sam"
    sam.write_text(
        "@SQ\tSN:s0\tLN:30\n@SQ\tSN:s1\tLN:30\n"
        f"p1\t97\ts1\t1\t60\t5M\ts0\t11\t0\t{seq}\t{qual}\tAS:i:5\n"  # unit 1: not assigned
        f"p1\t353\ts0\t2\t0\t2H3M\ts0\t11\t0\t*\t*\tAS:i:4\n"  # secondary on unit 0: kept
        f"p1\t145\ts0\t11\t60\t5M\ts1\t1\t0\t{seq}\t{qual}\tAS:i:5\n"
        f"p2\t65\ts0\t1\t60\t5M\t=\t1\t0\t{seq}\t{qual}\tAS:i:5\n"  # 0.4 < 0.5: dropped
    )
    counts = handoff(pile, refs, sam, tmp_path / "out.sam")
    assert counts == {"pairs": 2, "kept_mates": 2}
    rows = [
        line.split("\t")
        for line in (tmp_path / "out.sam").read_text().splitlines()
        if not line.startswith("@")
    ]
    first, second = rows
    assert (first[2], first[3], first[9], first[10]) == ("s0", "2", seq, qual)  # borrowed
    assert first[5] == "2S3M"  # the hard clip soft, to span the borrowed read
    assert int(first[1]) & 0x900 == 0 and int(first[1]) & 0x2  # primary, proper pair
    assert int(first[1]) & 0x20  # its mate is reverse
    assert first[6:9] == ["=", "11", "14"]
    assert second[4] == first[4] == "10"  # Phred of 1 - 0.9
    assert second[6:9] == ["=", "2", "-14"]
