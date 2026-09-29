"""The fmh-benchmark steps on the mini fixture, with reads from known positions."""

import math
import random
import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest

from kmer_functional_profiler.reference import reverse_complement

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "workflows" / "fmh-benchmark" / "bench.py"
MINI = ROOT / "tests" / "data" / "mini_fmh"
READ, FRAG = 150, 300


def run(cwd: Path, *args: str) -> None:
    subprocess.run([sys.executable, SCRIPT, *args], cwd=cwd, check=True)


def fastq(path: Path, reads: list[str]) -> None:
    path.write_text("".join(f"@r{i}\n{r}\n+\n{'I' * len(r)}\n" for i, r in enumerate(reads)))


def test_sample_truth_and_score(tmp_path: Path) -> None:
    run(tmp_path, "sample", "--genomes-dir", str(MINI / "genomes_extracted_from_kegg"),
        "--n", "2", "--seed", "1")  # fmt: skip
    records = (tmp_path / "sample.fna").read_text().split(">")[1:]
    contigs = {r.split("\n")[0]: r.split("\n")[1] for r in records}
    assert len(contigs) == 2 and all("|CP" in name for name in contigs)

    rng = random.Random(0)
    genes = pl.read_parquet(tmp_path / "genes.parquet")
    kos = pl.read_csv(MINI / "present_genes_and_koids.csv", columns=["gene_id", "ko_id"])
    r1, r2, spans = [], [], []
    for _ in range(40):
        contig = rng.choice(sorted(contigs))
        start = rng.randrange(len(contigs[contig]) - FRAG)
        frag = contigs[contig][start : start + FRAG]
        r1.append(frag[:READ])
        r2.append(reverse_complement(frag[-READ:].encode()).decode())
        spans.append((contig, start, start + FRAG))
    fastq(tmp_path / "r1.fq", r1)
    fastq(tmp_path / "r2.fq", r2)
    hit_genes = {
        g
        for contig, lo, hi in spans
        for g in genes.filter(pl.col("contig") == contig, pl.col("start") < hi, pl.col("end") > lo)[
            "gene_name"
        ]
    }
    expected = set(kos.filter(pl.col("gene_id").is_in(hit_genes))["ko_id"])

    run(tmp_path, "truth", "--fna", "sample.fna", "--genes", "genes.parquet", "--kos",
        str(MINI / "present_genes_and_koids.csv"), "--r1", "r1.fq", "--r2", "r2.fq",
        "--threads", "2")  # fmt: skip
    truth = pl.read_csv(tmp_path / "truth.csv")
    assert set(truth["ko_id"]) == expected
    assert (truth["depth"] > 0).all()

    found = sorted(expected)[:2]
    pl.DataFrame({"name": [*found, "ko:K99999"], "kmers_hit": [3, 1, 5]}).write_csv(
        tmp_path / "p.tsv", separator="\t"
    )
    run(tmp_path, "score", "--truth", "truth.csv", "--profile", "p.tsv", "--sample", "s",
        "--index", "i")  # fmt: skip
    score = pl.read_csv(tmp_path / "score.tsv", separator="\t").row(0, named=True)
    assert (score["tp"], score["n_pred"]) == (2, 3)
    assert score["completeness"] == 2 / len(expected)
    assert score["spearman_tp"] is None  # no abundance column in this profile

    # An exact abundance estimate scores perfectly.
    exact = truth.select(
        name="ko_id",
        kmers_hit=pl.lit(1),
        kmers_unique=pl.lit(1),
        coverage="depth",
        abundance_zi="depth",
        abundance_zi_lo=pl.col("depth") / 2,
        abundance_zi_hi=pl.col("depth") * 2,
    )
    exact.write_csv(tmp_path / "p.tsv", separator="\t")
    run(tmp_path, "score", "--truth", "truth.csv", "--profile", "p.tsv", "--sample", "s",
        "--index", "i")  # fmt: skip
    score = pl.read_csv(tmp_path / "score.tsv", separator="\t").row(0, named=True)
    assert score["spearman_tp"] == pytest.approx(1.0)
    assert score["l1"] == pytest.approx(0.0)
    scores = pl.read_csv(tmp_path / "score.tsv", separator="\t")
    zi = scores.filter(abundance="abundance_zi", min_hits=1).row(0, named=True)
    assert zi["ci_cover"] == 1.0
    assert zi["ci_width"] == pytest.approx(math.log(4))
    # Two false positives, one grouped with a true KO: fp_grouped counts the ungrouped one too.
    first = exact["name"][0]
    grouped = pl.concat(
        [exact, *(exact.head(1).with_columns(name=pl.lit(p)) for p in ("ko:FP1", "ko:FP2"))]
    ).with_columns(
        ambiguity_group=pl.when(pl.col("name").is_in([first, "ko:FP1"])).then(pl.lit(1)),
        group_abundance_zi_lo=pl.lit(0.0),
        group_abundance_zi_hi=pl.lit(1e9),
    )
    grouped.write_csv(tmp_path / "g.tsv", separator="\t")
    run(tmp_path, "score", "--truth", "truth.csv", "--profile", "g.tsv", "--sample", "s",
        "--index", "i", "--out", "g_score.tsv")  # fmt: skip
    g = pl.read_csv(tmp_path / "g_score.tsv", separator="\t")
    assert g.filter(abundance="abundance_zi", min_hits=1)["fp_grouped"].item() == 0.5
    # A score file where a metric is empty (read as String) still stacks with the others.
    empty = pl.read_csv(tmp_path / "score.tsv", separator="\t").with_columns(
        pl.lit(None, dtype=pl.Float64).alias("l1"), sample=pl.lit("t")
    )
    empty.write_csv(tmp_path / "score2.tsv", separator="\t")
    run(tmp_path, "summary", "score.tsv", "score2.tsv")
    assert pl.read_csv(tmp_path / "summary.tsv", separator="\t")["n_samples"].to_list() == [2] * 3
