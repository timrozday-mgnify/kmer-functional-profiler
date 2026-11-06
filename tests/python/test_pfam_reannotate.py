"""The Pfam re-annotation steps on a small release with six mini HMMs."""

import re
import shutil
import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "workflows" / "pfam-reannotate" / "reannotate.py"
RELEASE = ROOT / "tests" / "data" / "mini_reannotate"
HMM = ROOT / "tests" / "data" / "mini_fmh" / "Pfam-mini.hmm"
PERMISSIVE = ["-E", "10", "--domE", "1000000", "-Z", "1000000", "--domZ", "1000000"]


def run(cwd: Path, *args: object) -> None:
    subprocess.run([sys.executable, SCRIPT, *map(str, args)], cwd=cwd, check=True)


def hmmsearch(cwd: Path, faa: str, out: str, *args: str) -> None:
    subprocess.run(
        ["hmmsearch", "--notextw", *args, "--domtblout", out, "-o", f"{out}.txt", HMM, faa],
        cwd=cwd,
        check=True,
    )


def chunks(cwd: Path, *args: object) -> pl.DataFrame:
    run(cwd, "plan", "--release", RELEASE, "--bin-width", 10, *args)
    return pl.read_csv(cwd / "ranges.tsv", separator="\t")


def test_plan_and_fasta_cover_each_target_once(tmp_path: Path) -> None:
    release = pl.read_parquet(RELEASE / "mgy_cluster_seqs.parquet")
    for kind, expected in (
        ("reps", set(release["cluster_rep"])),
        ("all", set(release["cluster_member"])),
    ):
        ranges = chunks(tmp_path, "--targets", kind, "--chunk-proteins", 7)
        assert ranges.height > 1
        got: list[int] = []
        for chunk, lo, hi, n in ranges.iter_rows():
            run(tmp_path, "fasta", "--release", RELEASE, "--targets", kind, "--lo", lo,
                "--hi", hi, "--proteins", n, "--out", f"{chunk}.faa")  # fmt: skip
            got += [
                int(x[1:]) for x in (tmp_path / f"{chunk}.faa").read_text().split() if x[0] == ">"
            ]
        assert sorted(got) == sorted(expected)


def test_fasta_fails_when_the_plan_disagrees(tmp_path: Path) -> None:
    with pytest.raises(subprocess.CalledProcessError):
        run(tmp_path, "fasta", "--release", RELEASE, "--lo", 0, "--hi", 100, "--proteins", 1,
            "--out", "x.faa")  # fmt: skip


@pytest.mark.skipif(shutil.which("hmmsearch") is None, reason="needs HMMER")
def test_above_ga_is_what_cut_ga_reports_and_weak_hits_are_kept(tmp_path: Path) -> None:
    run(tmp_path, "ga", "--hmm", HMM)
    n = chunks(tmp_path, "--targets", "all", "--chunk-proteins", 10_000).item(0, "proteins")
    hi = chunks(tmp_path, "--targets", "all", "--chunk-proteins", 10_000).item(0, "hi")
    run(tmp_path, "fasta", "--release", RELEASE, "--targets", "all", "--lo", 0, "--hi", hi,
        "--proteins", n, "--out", "c.faa")  # fmt: skip
    hmmsearch(tmp_path, "c.faa", "wide.tbl", *PERMISSIVE)
    hmmsearch(tmp_path, "c.faa", "ga.tbl", "--cut_ga")
    ga = pl.read_csv(tmp_path / "ga.tsv", separator="\t")
    assert ga.height == 6
    wide = ("--domtbl", "wide.tbl", "--alignments", "wide.tbl.txt", "--pp")
    run(tmp_path, "parse", *wide, "--ga", "ga.tsv", "--counts", "c.faa.tsv",
        "--out", "hits.parquet", "--stats", "s.tsv")  # fmt: skip
    hits = pl.read_parquet(tmp_path / "hits.parquet")
    expected = sorted(
        (int(c[0]), int(c[4][2:].split(".")[0]), int(c[19]), int(c[20]))
        for c in map(str.split, (tmp_path / "ga.tbl").read_text().splitlines())
        if c[0] != "#" and not c[0].startswith("#")
    )
    above = hits.filter("above_ga").select("protein_id", "pfam_accession", "env_from", "env_to")
    assert sorted(above.rows()) == expected
    assert expected

    # Each alignment rebuilds from the searched sequence, its start and its cigar.
    fasta = (tmp_path / "c.faa").read_text().split()
    sequences = {int(name[1:]): seq for name, seq in zip(fasta[0::2], fasta[1::2], strict=True)}
    aligned = {}
    for block in (tmp_path / "wide.tbl.txt").read_text().split("  == domain ")[1:]:
        lines = block.splitlines()
        i = next(i for i, x in enumerate(lines) if len(x.split()) == 4 and x.split()[3].isdigit())
        t = lines[i + 2].split()
        aligned[(int(t[0]), int(t[1]), int(t[3]))] = t[2]
    assert len(aligned) == hits.height
    for row in hits.iter_rows(named=True):
        position, rebuilt = row["ali_from"] - 1, ""
        for length, op in re.findall(r"(\d+)([MID])", row["cigar"]):
            if op == "D":
                rebuilt += "-" * int(length)
            else:
                rebuilt += sequences[row["protein_id"]][position : position + int(length)]
                position += int(length)
        text = aligned[(row["protein_id"], row["ali_from"], row["ali_to"])]
        assert rebuilt.upper() == text.upper()
        assert len(row["pp"]) == len(text)

    run(tmp_path, "summary", "s.tsv", "s.tsv")
    total = pl.read_csv(tmp_path / "summary.tsv", separator="\t").row(0, named=True)
    assert total["proteins"] == 2 * n
    assert total["domains"] == 2 * hits.height
    assert total["domains_above_ga"] == 2 * len(expected)

    # Hits below GA stay in the table, marked.
    raised = tmp_path / "raised.tsv"
    ga.with_columns(ga_dom=pl.col("ga_dom") + 1000).write_csv(raised, separator="\t")
    run(tmp_path, "parse", *wide, "--ga", raised, "--counts", "c.faa.tsv",
        "--out", "weak.parquet", "--stats", "w.tsv")  # fmt: skip
    weak = pl.read_parquet(tmp_path / "weak.parquet")
    assert weak.height == hits.height
    assert not weak["above_ga"].any()


def test_reverse_flips_sequences_and_compare_bins_both_arms(tmp_path: Path) -> None:
    common = ("--release", RELEASE, "--lo", 0, "--hi", 10, "--targets", "all", "--proteins")
    n = (
        pl.scan_parquet(RELEASE / "mgy_protein_sequences.parquet")
        .filter(pl.col("protein_id") < 10)
        .select(pl.len())
        .collect()
        .item()
    )
    run(tmp_path, "fasta", *common, n, "--out", "f.faa")
    run(tmp_path, "fasta", *common, n, "--reverse", "--out", "r.faa")
    forward = (tmp_path / "f.faa").read_text().split()
    reverse = (tmp_path / "r.faa").read_text().split()
    assert forward[0::2] == reverse[0::2]
    assert [s[::-1] for s in forward[1::2]] == reverse[1::2]

    def hits(name: str, scores: list[float], above: list[bool]) -> None:
        pl.DataFrame(
            {"score": scores, "above_ga": above}, schema_overrides={"score": pl.Float32}
        ).write_parquet(tmp_path / name)

    hits("real.chunk0.parquet", [12.0, 14.9, 31.0], [False, False, True])
    hits("real.chunk1.parquet", [15.0], [False])
    hits("reversed.chunk0.parquet", [11.0, 40.0], [False, False])
    run(
        tmp_path, "compare", "real.chunk0.parquet", "real.chunk1.parquet", "reversed.chunk0.parquet"
    )
    table = pl.read_csv(tmp_path / "compare.tsv", separator="\t")
    assert table.rows() == [
        (10, False, 2, 1),
        (15, False, 1, 0),
        (40, False, 0, 1),
        (30, True, 1, 0),
    ]
    assert table.columns == ["score_from", "above_ga", "real", "reversed"]


def test_parse_fails_on_a_model_with_no_threshold(tmp_path: Path) -> None:
    (tmp_path / "c.tsv").write_text("proteins\tresidues\n1\t100\n")
    (tmp_path / "ga.tsv").write_text("pfam_accession\tga_seq\tga_dom\n2\t25\t25\n")
    (tmp_path / "t.tbl").write_text(
        "1 - 100 Mini1 PF90001.1 8 1e-9 50 0 1 1 1e-9 1e-9 50 0 1 8 1 8 1 8 0.9 -\n"
    )
    (tmp_path / "t.txt").write_text(block("AAAAAAAA", "AAAAAAAA", 1, 8, 1, 8))
    with pytest.raises(subprocess.CalledProcessError):
        run(tmp_path, "parse", "--domtbl", "t.tbl", "--alignments", "t.txt", "--ga", "ga.tsv",
            "--counts", "c.tsv", "--out", "h.parquet", "--stats", "s.tsv")  # fmt: skip


def block(model: str, target: str, m_from: int, m_to: int, t_from: int, t_to: int) -> str:
    return (
        "Accession:   PF90001.1\n"
        "  == domain 1  score: 50 bits;  conditional E-value: 1e-9\n"
        "           xxx CS\n"
        f"   Mini1 {m_from} {model} {m_to}\n"
        "            match line\n"
        f"   1 {t_from} {target} {t_to}\n"
        f"           {'9' * len(model)}\n"
    )


def test_parse_reads_inserts_and_deletions_and_keeps_pp_on_request(tmp_path: Path) -> None:
    (tmp_path / "c.tsv").write_text("proteins\tresidues\n1\t100\n")
    (tmp_path / "ga.tsv").write_text("pfam_accession\tga_seq\tga_dom\n90001\t25\t25\n")
    (tmp_path / "t.tbl").write_text(
        "1 - 100 Mini1 PF90001.1 80 1e-9 50 0 1 1 1e-9 1e-9 50 0 5 11 2 7 2 7 0.9 -\n"
    )
    (tmp_path / "t.txt").write_text(block("AAA.AAAAAA", "AAAaAA--AA", 5, 13, 2, 9))
    parse = ("parse", "--domtbl", "t.tbl", "--alignments", "t.txt", "--ga", "ga.tsv",
             "--counts", "c.tsv", "--stats", "s.tsv")  # fmt: skip
    run(tmp_path, *parse, "--pp", "--out", "with.parquet")
    run(tmp_path, *parse, "--out", "without.parquet")
    row = pl.read_parquet(tmp_path / "with.parquet").row(0, named=True)
    assert row["cigar"] == "3M1I2M2D2M"
    assert row["pp"] == "9" * 10
    assert "pp" not in pl.read_parquet(tmp_path / "without.parquet").columns
    assert pl.read_parquet(tmp_path / "without.parquet")["cigar"].item() == "3M1I2M2D2M"


def test_parse_fails_when_an_alignment_disagrees_with_its_coordinates(tmp_path: Path) -> None:
    (tmp_path / "c.tsv").write_text("proteins\tresidues\n1\t100\n")
    (tmp_path / "ga.tsv").write_text("pfam_accession\tga_seq\tga_dom\n90001\t25\t25\n")
    (tmp_path / "t.tbl").write_text(
        "1 - 100 Mini1 PF90001.1 80 1e-9 50 0 1 1 1e-9 1e-9 50 0 5 11 2 7 2 7 0.9 -\n"
    )
    (tmp_path / "t.txt").write_text(block("AAAA", "AAAA", 5, 9, 2, 5))
    with pytest.raises(subprocess.CalledProcessError):
        run(tmp_path, "parse", "--domtbl", "t.tbl", "--alignments", "t.txt", "--ga", "ga.tsv",
            "--counts", "c.tsv", "--out", "h.parquet", "--stats", "s.tsv")  # fmt: skip
