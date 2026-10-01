"""The fmh-benchmark steps on the mini fixture, with reads from known positions."""

import math
import random
import shutil
import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest

from kmer_functional_profiler.reference import reverse_complement, translate

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
    assert set(truth["label"]) == expected
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
        name="label",
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


def test_tool_profile_and_cost(tmp_path: Path) -> None:
    kos = tmp_path / "kos.csv"
    kos.write_text(",gene_id,ko_id\n0,g1,ko:K00001\n1,g2,ko:K00002\n2,g2,ko:K00003\n")

    def profile(tool: str, files: dict[str, str], *extra: str) -> dict[str, tuple[int, float]]:
        raw = tmp_path / tool
        raw.mkdir()
        for name, text in files.items():
            (raw / name).write_text(text)
        run(tmp_path, "tool-profile", "--tool", tool, "--raw", str(raw), "--kos", str(kos),
            "--out", f"{tool}.tsv", *extra)  # fmt: skip
        out = pl.read_csv(tmp_path / f"{tool}.tsv", separator="\t")
        return {n: (e, a) for n, e, a in out.iter_rows()}

    # DIAMOND: mates of one pair count once; a gene with two KOs counts for both
    hits = "p1/1\tg1|x\t50\t100\t90\np1/2\tg1|x\t50\t100\t90\np2/1\tg2|y\t30\t60\t50\n"
    assert profile("diamond", {"hits.tsv": hits}) == {
        "ko:K00001": (1, 1.0), "ko:K00002": (1, 0.5), "ko:K00003": (1, 0.5)}  # fmt: skip

    # fmh-funprofiler: shared hashes = intersect_bp / scaled; its own abundance
    prefetch = "match_name,intersect_bp,f_match_query\nko:K00001,3000,0.1\nko:K00002,1000,0.3\n"
    ko = "ko_id,abundance\nko:K00001,0.25\nko:K00002,0.75\n"
    assert profile("fmh_funprofiler", {"prefetch.csv": prefetch, "ko.csv": ko}) == {
        "ko:K00001": (3, 0.25), "ko:K00002": (1, 0.75)}  # fmt: skip

    # kMermaid: read pairs, and reads / mean member length
    pl.DataFrame(
        {
            "protein_id": ["a", "b", "c"],
            "cluster_rep": ["ko:K00001", "ko:K00001", "ko:K00002"],
            "sequence": ["M" * 100, "M" * 300, "M" * 50],
        }  # fmt: skip
    ).write_parquet(tmp_path / "members.parquet")
    reads = (
        "seq_name\tcluster_rep\tprot_name\tscore\np1/1\tko:K00001\tx\t5\np1/2\tko:K00001\tx\t4\n"
    )
    assert profile("kmermaid", {"kmermaid.tsv": reads}, "--members", "members.parquet") == {
        "ko:K00001": (1, 2 / 200)}  # fmt: skip

    # HUMAnN: unstratified KO rows only, KEGG's "ko:" prefix added
    table = ("# Gene Family\treads_Abundance-RPKs\nUNMAPPED\t10\nUNGROUPED\t5\nK00001\t2.5\n"
             "K00001|g__Escherichia.s__Escherichia_coli\t2.5\nK00002\t0\n")  # fmt: skip
    assert profile("humann", {"ko.tsv": table}) == {"ko:K00001": (1, 2.5)}
    table4 = ("# Gene Family HUMAnN v4.0.0.alpha.2 Adjusted CPMs\treads\nREADS_UNMAPPED\t10\n"
              "UNGROUPED\t5\nK00002\t7\nK00002|s__Bacteroides_ovatus.t__SGB1871\t7\n")  # fmt: skip
    assert profile("humann4", {"ko.tsv": table4}) == {"ko:K00002": (1, 7.0)}

    trace = ("task_id\tname\tstatus\trealtime\t%cpu\tpeak_rss\n"
             "1\tPROFILE (seed 1 kfp_s100)\tCOMPLETED\t3600000\t100\t2000000000\n"
             "2\tPROFILE (seed 2 kfp_s100)\tCACHED\t7200000\t50\t1000000000\n"
             "3\tDIAMOND (seed 1)\tFAILED\t10\t100\t1\n")  # fmt: skip
    (tmp_path / "trace.tsv").write_text(trace)
    run(tmp_path, "cost", "trace.tsv")
    cost = pl.read_csv(tmp_path / "cost.tsv", separator="\t").row(0, named=True)
    assert (cost["what"], cost["tasks"], cost["hours_mean"], cost["cpu_hours_mean"],
            cost["peak_rss_gb"]) == ("kfp_s100", 2, 1.5, 1.0, 2.0)  # fmt: skip


@pytest.mark.skipif(shutil.which("hmmsearch") is None, reason="needs HMMER")
def test_pfam_domains_truth_and_profile(tmp_path: Path) -> None:
    genomes = str(MINI / "genomes_extracted_from_kegg")
    hmmsearch = ["hmmsearch", "--cut_ga", "--domtblout", "d.tbl", "-o", "/dev/null"]
    run(tmp_path, "pfam-proteins", "--genomes-dir", genomes, "--chunks", "2")
    for i in (0, 1):
        search = [*hmmsearch[:3], f"d{i}.tbl", *hmmsearch[4:], MINI / "Pfam-mini.hmm"]
        subprocess.run([*search, f"proteins.{i}.faa"], cwd=tmp_path, check=True)
    run(tmp_path, "pfam-domains", "--domtbl", "d0.tbl", "d1.tbl", "--proteins",
        "proteins.0.faa", "proteins.1.faa")  # fmt: skip
    domains = pl.read_parquet(tmp_path / "domains.parquet")
    assert domains.height == 6 and domains["gene_name"].n_unique() == 5  # one gene has two
    members = pl.read_parquet(tmp_path / "pfam_members.parquet")
    assert (
        members["sequence"].str.len_chars().to_list()
        == (domains["aa_end"] - domains["aa_start"]).to_list()
    )

    # Domain intervals on the genome translate back to the domain, on either strand.
    run(tmp_path, "sample", "--genomes-dir", genomes, "--n", "3", "--seed", "1")
    records = (tmp_path / "sample.fna").read_text().split(">")[1:]
    contigs = {r.split("\n")[0]: r.split("\n")[1] for r in records}
    genes = pl.read_parquet(tmp_path / "genes.parquet")
    sys.path.insert(0, str(SCRIPT.parent))
    import bench

    features = bench.domain_features(genes, domains)
    assert features.height == 6
    by_name = dict(members.select("cluster_rep", "sequence").iter_rows())
    for contig, start, end, label in features.iter_rows():
        dna = contigs[contig][start:end].encode()
        strand = genes.filter(pl.col("contig") == contig, pl.col("start") <= start,
                              pl.col("end") >= end)["strand"].item()  # fmt: skip
        if strand == "-":
            dna = reverse_complement(dna)
        assert translate(dna).decode() == by_name[label]

    # Truth: reads over a domain make its Pfam present; reads elsewhere on the gene do not.
    contig, start, end, label = features.row(0)
    seq = contigs[contig]
    lo = max(0, start - 100)
    frag = seq[lo : lo + FRAG]
    fastq(tmp_path / "r1.fq", [frag[:READ]])
    fastq(tmp_path / "r2.fq", [reverse_complement(frag[-READ:].encode()).decode()])
    run(tmp_path, "truth", "--fna", "sample.fna", "--genes", "genes.parquet", "--kos",
        str(MINI / "present_genes_and_koids.csv"), "--r1", "r1.fq", "--r2", "r2.fq",
        "--domains", "domains.parquet")  # fmt: skip
    pfam = pl.read_csv(tmp_path / "truth_pfam.csv")
    hit = features.filter(pl.col("contig") == contig, pl.col("start") < lo + FRAG,
                          pl.col("end") > lo)  # fmt: skip
    assert set(pfam["label"]) == set(hit["label"]) and label in set(pfam["label"])

    # A unit profile summed per Pfam: a unit with two Pfams counts for both.
    units = {"unit": [0, 1, 2], "name": ["a", "b", "c"], "kmers_hit": [1, 2, 4],
             "coverage": [0.5, 1.0, 2.0]}  # fmt: skip
    pl.DataFrame(units).write_csv(tmp_path / "u.tsv", separator="\t")
    pl.DataFrame({"unit": [0, 1, 1], "pfam_accession": ["PF1.2", "PF1.2", "PF2.1"],
                  "n_members": [1, 1, 1]}).write_parquet(tmp_path / "up.parquet")  # fmt: skip
    run(tmp_path, "pfam-profile", "--profile", "u.tsv", "--unit-pfam", "up.parquet")
    got = pl.read_csv(tmp_path / "pfam_profile.tsv", separator="\t")
    assert got.select("name", "kmers_hit", "coverage").rows() == [("PF1", 3, 1.5), ("PF2", 2, 1.0)]
    pl.DataFrame({"unit": [0, 1], "pfam_accession": [1007, 42], "n_members": [1, 1]}).write_parquet(
        tmp_path / "up.parquet"
    )
    run(tmp_path, "pfam-profile", "--profile", "u.tsv", "--unit-pfam", "up.parquet")
    got = pl.read_csv(tmp_path / "pfam_profile.tsv", separator="\t")
    assert got["name"].to_list() == ["PF00042", "PF01007"]
