"""The fmh-benchmark steps on the mini fixture, with reads from known positions."""

import gzip
import json
import math
import random
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
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
    per_gene = pl.read_csv(tmp_path / "truth_genes.csv")
    assert set(per_gene["gene_name"]) == hit_genes and (per_gene["depth"] > 0).all()
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
    # A joint query's extra-index (decoy) units are not predictions.
    pl.DataFrame(
        {"name": [*found, "ko:K99999", "7"], "kmers_hit": [3, 1, 5, 9], "source": [0, 0, 0, 1]}
    ).write_csv(tmp_path / "p.tsv", separator="\t")
    run(tmp_path, "score", "--truth", "truth.csv", "--profile", "p.tsv", "--sample", "s",
        "--index", "i")  # fmt: skip
    score = pl.read_csv(tmp_path / "score.tsv", separator="\t").row(0, named=True)
    assert (score["tp"], score["n_pred"]) == (2, 3)

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


def test_mix_host_abundance_and_decoy(tmp_path: Path) -> None:
    def pairs(prefix: str, n: int) -> None:
        for mate in (1, 2):
            with gzip.open(tmp_path / f"{prefix}_R{mate}.fastq.gz", "wt") as f:
                f.writelines(f"@{prefix}{i}/{mate}\nACGT\n+\nIIII\n" for i in range(n))

    pairs("m", 1000)
    pairs("h", 950)
    reads = ["--r1", "m_R1.fastq.gz", "--r2", "m_R2.fastq.gz", "--host-r1", "h_R1.fastq.gz",
             "--host-r2", "h_R2.fastq.gz"]  # fmt: skip
    run(tmp_path, "mix", *reads, "--fraction", "0.9", "--seed", "1")
    texts = [gzip.decompress((tmp_path / f"mixed_R{m}.fastq.gz").read_bytes()) for m in (1, 2)]
    names = [t.decode().splitlines()[::4] for t in texts]
    assert [n[:-2] for n in names[0]] == [n[:-2] for n in names[1]]  # mates stay paired
    host = sum(n.startswith("@h") for n in names[0])
    assert host == 900 and 60 < len(names[0]) - host < 140  # ~100 microbial pairs kept
    with pytest.raises(subprocess.CalledProcessError):  # too few host reads for 99%
        run(tmp_path, "mix", *reads, "--fraction", "0.99")

    with gzip.open(tmp_path / "host.fa.gz", "wt") as f:
        f.write(">chr1 x\n" + "A" * 900 + "\n>chrM\n" + "A" * 10 + "\n>NC_001422.1\nAAAA\n")
    run(tmp_path, "host-abundance", "--fasta", "host.fa.gz", "--mito-copies", "10")
    got = dict(line.split("\t") for line in (tmp_path / "abundance.txt").read_text().splitlines())
    assert float(got["chr1"]) == pytest.approx(0.995 * 0.9)
    assert float(got["chrM"]) == pytest.approx(0.995 * 0.1)
    assert float(got["NC_001422.1"]) == pytest.approx(0.005)

    (tmp_path / "p.faa").write_text(">sp|P1|A\nMKV\nLL\n>sp|P2|B\nMAA\n")
    run(tmp_path, "decoy-members", "--faa", "p.faa")
    decoy = pl.read_parquet(tmp_path / "decoy_members.parquet")
    assert decoy["sequence"].to_list() == ["MKVLL", "MAA"] and decoy["cluster_rep"].n_unique() == 2


def test_split_scores() -> None:
    sys.path.insert(0, str(SCRIPT.parent))
    import bench

    truth = pl.DataFrame({"label": ["a", "b", "c"], "depth": [3.0, 1.0, 5.0]})
    detected = pl.DataFrame(
        {
            "name": ["a", "b", "c", "x", "y"],
            "component": [0, 0, 1, 2, 2],  # c alone; x, y have no truth
            "abundance_zi": [2.0, 2.0, 7.0, 1.0, 1.0],
            "group_size": [2, 2, 1, 1, 1],
            "host_like": [False, False, False, True, False],
        }
    )
    got = bench.split_scores(truth, detected, "abundance_zi")
    # Component 0: shares (0.5, 0.5) against (0.75, 0.25)
    assert got["split_l1"] == pytest.approx(0.5)
    assert got["group_size_mean"] == pytest.approx(1.4)
    assert got["host_like_detected"] == 1


def test_mgnify_genes_and_aai_score(tmp_path: Path) -> None:
    sys.path.insert(0, str(SCRIPT.parent))
    import bench

    # g1: nearest cluster 10 at 96% (and a near hit, 11 at 80%); g2: nearest 20 at 85%;
    # g3: absent from the sample. Two HSPs of one pair: the best counts.
    hits = ("g1\t10\t96.0\t200\t200\t210\t380\ng1\t11\t80.0\t190\t200\t200\t250\n"
            "g1\t11\t70.0\t50\t200\t200\t40\ng2\t20\t85.0\t300\t300\t300\t400\n"
            "g2\t10\t90.0\t300\t300\t210\t300\n"
            "g3\t30\t99.0\t100\t100\t100\t200\n")  # fmt: skip
    (tmp_path / "h.tsv").write_text(hits)
    run(tmp_path, "mgnify-genes", "--hits", "h.tsv")
    units = pl.read_parquet(tmp_path / "gene_units.parquet")
    assert units.filter(gene_name="g1").select("cluster_rep", "rank", "identity").rows() == [
        (10, 1, 0.96), (11, 2, 0.8)]  # fmt: skip
    genes = pl.DataFrame({"gene_name": ["g1", "g2", "g3"], "depth": [3.0, 1.0, 0.0]})
    profile = pl.DataFrame(
        {
            "cluster_rep": [10, 11, 20, 99],
            "kmers_unique": [5, 0, 3, 2],  # 11 explained away; 99 a false detection
            "aai": [0.95, None, 0.80, 0.9],
            "aai_lo": [0.93, None, 0.70, 0.8],
            "aai_hi": [0.99, None, 0.84, 0.95],
            "aai_naive": [0.97, 0.79, 0.86, 0.9],
            "aai_kmers": [6.0, 0.0, 12.0, 2.0],
        }
    )
    got = bench.aai_score(profile, units, genes)
    assert (got["genes_present"], got["genes_in90"]) == (2, 1)
    assert got["completeness_90"] == 1.0  # cluster 10 detected
    assert got["purity_nearest"] == pytest.approx(2 / 3)  # 10 and 20 are nearest; 99 is not
    assert got["recall_0.95"] == 1.0 and got["recall_0.8"] == 1.0  # g2 beyond 90%: 20 found
    assert got["aai_bias_0.95"] == pytest.approx(-0.01) and got["aai_cover_0.95"] == 1.0
    assert got["aai_bias_0.8"] == pytest.approx(-0.05) and got["aai_cover_0.8"] == 0.0
    # union truth: 10 holds g1's (96%) and g2's (90%) k-mers, as one strain at ~97.4%
    union = (1 - (1 - 0.96**11) * (1 - 0.9**11)) ** (1 / 11)
    assert got["aai_bias_union_0.95"] == pytest.approx(0.95 - union)
    assert got["aai_cover_union"] == 0.5  # 10 covered, 20 (85%) not
    assert got["aai_n_kmers5"] == 1 and got["aai_bias_kmers5"] == pytest.approx(0.95 - union)
    assert got["aai_n_kmers10"] == 1 and got["aai_n_kmers0"] == 0
    # single-gene units: 20 (g2 alone, 85%, interval 0.70-0.84 misses); 10 has two genes
    assert got["aai_n_single"] == 1 and got["aai_cover_single"] == 0.0
    assert got["aai_bias_single_0.8"] == pytest.approx(0.80 - 0.85)
    assert got["aai_bias_single_0.95"] is None
    assert got["aai_n_single_kmers10"] == 1 and got["aai_n_single_kmers5"] == 0
    # aai_naive on every hit unit, near hits included: 11 against g1's 80%
    assert got["naive_n"] == 3 and got["naive_bias_0.8"] == pytest.approx(
        (0.86 - 0.85 + 0.79 - 0.8) / 2
    )


def test_nearest_member_truth(tmp_path: Path) -> None:
    sys.path.insert(0, str(SCRIPT.parent))
    import bench

    # Clusters 10 (members 10, 101, 102, 103) and 20 (members 20, 201). g1 hits 10's rep at
    # 85% and member 101 at 93% (102 at 99% but query coverage 0.3, 103 at 97% but subject
    # coverage 0.5: both skipped, as aai-model's pairs); g2 hits 20 at 85% in pass 1 and no
    # member in pass 2 (the rep's identity stays).
    (tmp_path / "h1.tsv").write_text("g1\t10\t85.0\t200\t200\t200\t300\n"
                                     "g2\t20\t85.0\t200\t200\t200\t300\n")  # fmt: skip
    run(tmp_path, "mgnify-genes", "--hits", "h1.tsv", "--out", "reps.parquet")
    members = pl.DataFrame({"protein_id": [10, 101, 102, 103, 20, 201, 30],
                            "cluster_rep": [10, 10, 10, 10, 20, 20, 30],
                            "sequence": ["MK"] * 7})  # fmt: skip
    members.write_parquet(tmp_path / "members.parquet")
    run(tmp_path, "mgnify-members", "--members", "members.parquet", "--gene-units", "reps.parquet")
    clusters = pl.read_parquet(tmp_path / "member_clusters.parquet")
    assert sorted(clusters["protein_id"].to_list()) == [10, 20, 101, 102, 103, 201]  # not 30
    assert (tmp_path / "members.faa").read_text().count(">") == 6
    (tmp_path / "h2.tsv").write_text("g1\t10\t85.5\t200\t200\t200\t300\n"
                                     "g1\t101\t93.0\t200\t200\t230\t350\n"
                                     "g1\t102\t99.0\t60\t200\t60\t100\n"
                                     "g1\t103\t97.0\t200\t200\t400\t350\n")  # fmt: skip
    run(tmp_path, "mgnify-nearest", "--hits", "h2.tsv", "--gene-units", "reps.parquet",
        "--member-clusters", "member_clusters.parquet", "--out", "near.parquet")  # fmt: skip
    near = {
        r["gene_name"]: r for r in pl.read_parquet(tmp_path / "near.parquet").iter_rows(named=True)
    }
    assert near["g1"]["nearest"] == 101 and near["g1"]["identity_nearest"] == pytest.approx(0.93)
    assert near["g1"]["identity"] == pytest.approx(0.85)  # rep identity kept for detection
    assert near["g2"]["nearest"] == 20 and near["g2"]["identity_nearest"] == pytest.approx(0.85)
    assert near["g1"]["member_hit"] and not near["g2"]["member_hit"]
    # aai is scored against the nearest member; coverage also by n_members
    units = pl.read_parquet(tmp_path / "near.parquet")
    genes = pl.DataFrame({"gene_name": ["g1", "g2"], "depth": [2.0, 1.0]})
    profile = pl.DataFrame({"cluster_rep": [10, 20], "kmers_unique": [5, 5], "aai": [0.92, 0.84],
                            "aai_lo": [0.9, 0.8], "aai_hi": [0.95, 0.88], "aai_kmers": [9.0, 9.0],
                            "n_members": [3, 2]})  # fmt: skip
    got = bench.aai_score(profile, units, genes)
    assert got["aai_bias_0.9"] == pytest.approx(0.92 - 0.93) and got["aai_cover_0.9"] == 1.0
    assert got["aai_n_members2"] == 2 and got["aai_cover_members2"] == 1.0
    assert got["aai_n_members1"] == 0
    # rep truth (no member pass): g1's cluster sits at 85%
    rep = bench.aai_score(profile, units.drop("identity_nearest", "nearest"), genes)
    assert rep["aai_bias_0.8"] == pytest.approx(((0.92 - 0.85) + (0.84 - 0.85)) / 2)


def test_aai_score_by_gene_length() -> None:
    sys.path.insert(0, str(SCRIPT.parent))
    import bench

    # four genes at 100/200/400/700 aa, each alone in its own cluster at 98%; the 700-aa
    # one undetected; coverage over depth 2, 2, 1: the short genes read high
    genes = pl.DataFrame({"gene_name": ["a", "b", "c", "d"], "depth": [3.0, 3.0, 3.0, 3.0]})
    genes = genes.with_columns(bases=pl.col("depth") * 3 * pl.Series([100, 200, 400, 700]))
    units = pl.DataFrame(
        {
            "gene_name": ["a", "b", "c", "d"],
            "cluster_rep": [1, 2, 3, 4],
            "identity": 0.98,
            "qcov": 1.0,
            "scov": 1.0,
            "rank": 1,
        }
    )
    profile = pl.DataFrame(
        {"cluster_rep": [1, 2, 3], "kmers_unique": [5, 5, 5], "coverage_em": [6.0, 6.0, 3.0]}
    )
    got = bench.aai_score(profile, units, genes)
    assert [got[f"completeness_90_len{lo}"] for lo in (0, 150, 300, 600)] == [1.0, 1.0, 1.0, 0.0]
    # log2 ratios 1, 1, 0 around their median 1
    assert (got["abund_bias_len0"], got["abund_bias_len150"]) == (0.0, 0.0)
    assert got["abund_bias_len300"] == -1.0 and got["abund_err_len300"] == 1.0
    assert got["abund_bias_len600"] is None


def test_aai_calibration_inverts_a_biased_estimator() -> None:
    sys.path.insert(0, str(SCRIPT.parent))
    import bench

    # aai reads high and compressed (0.6 + 0.4 a) with noise; the inverse map undoes it
    # whatever the identity distribution, and a null aai stays null
    rng = np.random.default_rng(0)
    true = rng.uniform(0.7, 1.0, 20_000)
    aai = 0.6 + 0.4 * true + rng.normal(0, 0.01, true.size)
    pairs = pl.DataFrame({"true": true, "aai": aai, "aai_lo": aai - 0.005, "aai_hi": aai + 0.005})
    cal = bench.fit_aai_calibration(pairs)
    assert np.all(np.diff(cal["aai"]) > 0)
    got = bench.calibrate_aai(pairs, cal)
    err = (got["aai"] - got["true"]).to_numpy()
    assert abs(np.median(err)) < 0.005 and np.median(np.abs(err)) < 0.03
    assert got.select(pl.col("true").is_between("aai_lo", "aai_hi").mean()).item() >= 0.94
    f64 = pl.Float64
    nulls = pl.DataFrame(
        {"aai": [None, 0.9], "aai_lo": [None, 0.89], "aai_hi": [None, 0.91]},
        schema={"aai": f64, "aai_lo": f64, "aai_hi": f64},
    )
    nulls = bench.calibrate_aai(nulls, cal)
    assert nulls["aai"].null_count() == 1 and nulls["aai"][1] == pytest.approx(0.75, abs=0.01)


def test_aai_calibrate_transfers_a_map_between_runs(tmp_path: Path) -> None:
    # Two runs of one biased estimator (0.6 + 0.4 a) with different noise: a map fitted on
    # run A, scored with --map on run B's held-out clusters, removes the bias there; a map
    # from an index built at another k is refused.
    n = 6000
    rng = np.random.default_rng(0)
    true = rng.uniform(0.7, 1.0, n)
    reps = np.arange(n)
    pl.DataFrame(
        {
            "gene_name": [f"g{i}" for i in reps],
            "cluster_rep": reps,
            "identity": true,
            "qcov": 1.0,
            "scov": 1.0,
            "rank": 1,
        }
    ).write_parquet(tmp_path / "gene_units.parquet")
    pl.DataFrame({"gene_name": [f"g{i}" for i in reps], "depth": 1.0}).write_csv(
        tmp_path / "genes.csv"
    )
    params = {"k": 11, "alphabet": "protein", "t_base": 0.02, "n_min": 8, "t_cap": 0.2,
              "t_dense": 0.0}  # fmt: skip
    for run_name, seed, k in (("a", 1, 11), ("b", 2, 11), ("c", 3, 9)):
        aai = 0.6 + 0.4 * true + np.random.default_rng(seed).normal(0, 0.01, n)
        pl.DataFrame(
            {
                "cluster_rep": reps,
                "kmers_unique": 5,
                "aai": aai,
                "aai_lo": aai - 0.005,
                "aai_hi": aai + 0.005,
                "aai_naive": aai,
            }
        ).write_csv(tmp_path / f"{run_name}.tsv", separator="\t")
        (tmp_path / f"idx_{run_name}").mkdir()
        (tmp_path / f"idx_{run_name}" / "meta.json").write_text(
            json.dumps({"params": params | {"k": k}})
        )
    common = ["--genes", "genes.csv", "--gene-units", "gene_units.parquet"]
    run(tmp_path, "aai-calibrate", "--profiles", "a.tsv", *common, "--index", "idx_a",
        "--out", "a.json", "--scores-out", "a_scores.tsv")  # fmt: skip
    run(tmp_path, "aai-calibrate", "--profiles", "b.tsv", *common, "--index", "idx_b",
        "--map", "a.json", "--out", "unused.json", "--scores-out", "a_to_b.tsv")  # fmt: skip
    assert not (tmp_path / "unused.json").exists()
    scores = pl.read_csv(tmp_path / "a_to_b.tsv", separator="\t")
    raw, cal = (scores.filter(method=m).row(0, named=True) for m in ("raw", "calibrated"))
    assert raw["aai_bias_0.7"] > 0.1  # 0.6 + 0.4 * 0.75 reads 0.15 high
    for lo in (0.7, 0.8, 0.9, 0.95):
        assert abs(cal[f"aai_bias_{lo}"]) < 0.01 and cal[f"aai_cover_{lo}"] >= 0.9
    with pytest.raises(subprocess.CalledProcessError):
        run(tmp_path, "aai-calibrate", "--profiles", "c.tsv", *common, "--index", "idx_c",
            "--map", "a.json", "--scores-out", "a_to_c.tsv")  # fmt: skip


def test_divergence_abundance() -> None:
    sys.path.insert(0, str(SCRIPT.parent))
    import bench

    # 12 one-member units at 98% set the scale (coverage_zi = 2 x depth); unit 100 at 85%
    # reads coverage_em at a quarter of 2 x depth, coverage_zi at it, inside its interval;
    # unit 200 has two genes and is left out, but for the strain-mix metrics: against
    # 2 x (3 + 3), coverage_zi reads half and coverage_mix (split) the sum
    n = 12
    hits = pl.DataFrame({
        "gene_name": [f"g{i}" for i in range(n)] + ["d", "s1", "s2"],
        "cluster_rep": list(range(n)) + [100, 200, 200],
        "identity": [0.98] * n + [0.85, 0.9, 0.9],
        "depth": [float(i + 1) for i in range(n)] + [4.0, 3.0, 3.0],
        "rank": [1] * (n + 3),
    })  # fmt: skip
    depth = np.r_[np.arange(1, n + 1), 4.0]
    profile = pl.DataFrame({
        "cluster_rep": list(range(n)) + [100, 200],
        "kmers_unique": [5] * (n + 2),
        "n_members": [1] * (n + 2),
        "coverage_em": list(2 * depth) + [6.0],
        "coverage_zi": list(2 * depth) + [6.0],
        "abundance_zi": list(2 * depth) + [6.0],
        "coverage_zi_lo": list(1.9 * depth) + [5.0],
        "coverage_zi_hi": list(2.1 * depth[:n]) + [9.0, 7.0],
        "coverage_mix": list(2 * depth) + [12.0],
        "mix_rates": [1] * (n + 1) + [2],
    }).with_columns(coverage_em=pl.when(pl.col("cluster_rep") == 100).then(2.0)
                    .otherwise(pl.col("coverage_em")))  # fmt: skip
    got = bench.divergence_abundance(profile, hits)
    assert got["abund_scale"] == pytest.approx(2.0)
    assert got["abund_em_bias_0.8"] == pytest.approx(np.log2(2.0 / 8.0))
    assert got["abund_zi_bias_0.8"] == pytest.approx(0.0) and got["abund_zi_bias_0.95"] == 0
    assert got["abund_zi_cover"] == pytest.approx(13 / 13)  # 8 is inside [5, 9]
    assert got["abund_n_depth2"] == 4 and got["abund_n_depth5"] == 8
    assert got["abund_genes_n_1"] == n + 1 and got["abund_genes_split_1"] == 0
    assert got["abund_genes_zi_bias_2"] == pytest.approx(-1.0)
    assert got["abund_genes_mix_bias_2"] == pytest.approx(0.0) and got["abund_genes_split_2"] == 1
    assert got["abund_genes_n_4"] == 0 and got["abund_genes_mix_bias_4"] is None


def test_study_rebuild_and_ladder(tmp_path: Path) -> None:
    # Base: two MGnify clusters with integer ids and Pfam numbers.
    base = {"protein_id": [1, 2], "cluster_rep": [1, 1], "full_length": [True, True],
            "sequence": ["MAAA", "MAAC"]}  # fmt: skip
    pl.DataFrame(base).write_parquet(tmp_path / "base.parquet")
    pl.DataFrame({"protein_id": [1], "pfam_accession": [1007]}).write_parquet(
        tmp_path / "bp.parquet"
    )
    # Study: g1 is 95% identical to cluster 1, g2 only 80%, g3 has no hit; linclust put all
    # three in g1's cluster.
    pl.DataFrame({"protein_id": ["g1", "g2", "g3"], "cluster_rep": ["g1", "g1", "g1"],
                  "full_length": [True, True, False], "sequence": ["MAAD", "MAAE", "MAAF"]}
                 ).write_parquet(tmp_path / "sm.parquet")  # fmt: skip
    pl.DataFrame({"protein_id": ["g3"], "pfam_accession": ["PF00042"]}).write_parquet(
        tmp_path / "sp.parquet"
    )
    pl.DataFrame({"gene_name": ["g1", "g2", "g1"], "cluster_rep": [1, 1, 2],
                  "identity": [0.95, 0.8, 0.95], "qcov": [1.0, 1.0, 1.0], "scov": [1.0, 1.0, 1.0],
                  "rank": [1, 1, 2]}).write_parquet(tmp_path / "gu.parquet")  # fmt: skip
    run(tmp_path, "study-rebuild", "--members", "base.parquet", "--pfam", "bp.parquet",
        "--study-members", "sm.parquet", "--study-pfam", "sp.parquet", "--gene-units",
        "gu.parquet")  # fmt: skip
    members = pl.read_parquet(tmp_path / "members.parquet")
    assert members.select("protein_id", "cluster_rep").rows() == [
        ("1", "1"), ("2", "1"), ("g1", "1"), ("g2", "g1"), ("g3", "g1")
    ]  # fmt: skip
    assert pl.read_parquet(tmp_path / "pfam.parquet").rows() == [
        ("1", "PF01007"),
        ("g3", "PF00042"),
    ]

    # The joint query recovers 3/4 of the rebuild's completeness gain.
    keys = {"index": "m", "count": "kmers_unique", "abundance": "coverage_em", "min_hits": 1,
            "label": "pfam", "sample": "seed1"}  # fmt: skip
    for i, (arm, completeness) in enumerate((("", 0.5), ("study", 0.65), ("rebuild", 0.7))):
        pl.DataFrame([keys | {"arm": arm, "completeness": completeness, "purity": 0.9}]).write_csv(
            tmp_path / f"s{i}.tsv", separator="\t"
        )
    run(tmp_path, "study-ladder", "--scores", "s0.tsv", "s1.tsv", "s2.tsv")
    ladder = pl.read_csv(tmp_path / "study_ladder.tsv", separator="\t")
    assert ladder["recovered"].to_list() == pytest.approx([0.75])


def test_uhgg_pick_sample_and_score(tmp_path: Path) -> None:
    """The species benchmark's steps on the mini catalogue (tests/data/mini_uhgg)."""
    catalogue = ROOT / "tests" / "data" / "mini_uhgg"
    run(tmp_path, "uhgg-pick", "--metadata", str(catalogue / "genomes-all_metadata.tsv"),
        "--replicates", "2", "--per-sample", "2", "--min-genomes", "3", "--distractors", "1",
        "--depth-min", "1", "--depth-max", "4")  # fmt: skip
    samples = pl.read_csv(tmp_path / "samples.tsv", separator="\t")
    # held-out genomes: non-representatives of >= 3-genome species, >= 90% complete, unique
    assert set(samples["genome"]) <= {"MGYG000000002", "MGYG000000005", "MGYG000000006"}
    assert samples["genome"].n_unique() == samples.height
    assert samples["depth"].is_between(1, 4).all()
    assert (tmp_path / "exclude.txt").read_text().split() == samples["genome"].to_list()
    species = (tmp_path / "species.txt").read_text().split()
    assert (
        set(samples["species"]) < set(species) and len(species) == samples["species"].n_unique() + 1
    )

    genome = samples["genome"][0]
    (tmp_path / "g").mkdir()
    (tmp_path / "g" / f"{genome}.fna").write_text(">c1 x\nACGT\n>c2\nGG\n")
    one = samples.filter(pl.col("genome") == genome).with_columns(sample=pl.lit(1))
    one.write_csv(tmp_path / "one.tsv", separator="\t")
    run(tmp_path, "uhgg-sample", "--samples", "one.tsv", "--sample", "1", "--genomes-dir", "g")
    assert (tmp_path / "sample.fna").read_text() == f">{genome}|c1\nACGT\n>{genome}|c2\nGG\n"
    depth = one["depth"][0]
    assert (
        tmp_path / "coverage.txt"
    ).read_text() == f"{genome}|c1\t{depth}\n{genome}|c2\t{depth}\n"

    # scoring: a species index whose held-out genome carries units 1-3; we call A at the
    # right depth plus a false B, and units 1, 2 and a false 9 present
    si = tmp_path / "si"
    si.mkdir()
    pl.DataFrame({"species": [0, 1], "id": ["SA", "SB"], "name": ["A", "B"]}).write_csv(
        si / "species.tsv", separator="\t"
    )
    pl.DataFrame({"name": ["G"] * 3, "species": ["SA"] * 3, "unit": [1, 2, 3],
                  "c": [1.0] * 3}).write_parquet(si / "held_out.parquet")  # fmt: skip
    pl.DataFrame({"sample": [1], "genome": ["G"], "species": ["SA"], "depth": [2.0]}).write_csv(
        tmp_path / "truth.tsv", separator="\t"
    )
    pl.DataFrame({"id": ["SA", "SB"], "relative_abundance": [0.75, 0.25]}).write_csv(
        tmp_path / "species.tsv", separator="\t"
    )
    pl.DataFrame({"unit": [1, 2, 9, 3], "hits": [5.0, 0.0, 3.0, 0.0],
                  "present_prob": [0.9, None, 0.8, None],
                  "present_prob_updated": [0.95, 0.7, 0.8, 0.2]}).write_csv(
        tmp_path / "presence.tsv", separator="\t")  # fmt: skip
    pl.DataFrame({"species": [0, 0], "unit": [3, 4], "prevalence": [0.5, 0.5], "hits": [0.0, 0.0],
                  "carriage_prob": [0.8, 0.2]}).write_csv(
        tmp_path / "units.tsv", separator="\t")  # fmt: skip
    run(tmp_path, "uhgg-score", "--truth", "truth.tsv", "--species-index", "si", "--sample", "1",
        "--arm", "species", "--pred", "species.tsv", "--presence", "presence.tsv",
        "--units", "units.tsv")  # fmt: skip
    got = dict(pl.read_csv(tmp_path / "species_score.tsv", separator="\t")
               .select("metric", "value").iter_rows())  # fmt: skip
    assert got["purity"] == 0.5 and got["completeness"] == 1.0 and got["l1"] == pytest.approx(0.5)
    assert got["unit_all_completeness_observed"] == pytest.approx(1 / 3)
    assert got["unit_all_completeness_updated"] == pytest.approx(2 / 3)
    assert got["unit_all_purity_updated"] == pytest.approx(2 / 3)
    assert got["carriage_auc"] == 1.0 and got["accessory_pairs"] == 2
    calibration = pl.read_csv(tmp_path / "species_calibration.tsv", separator="\t")
    assert calibration["n"].sum() == 2  # the zero-hit units 2 and 3


def test_uhgg_two_strain_samples_and_lineage_score(tmp_path: Path) -> None:
    """Two-strain samples (phase 11, step 12) and the lineage arm's scoring."""
    catalogue = ROOT / "tests" / "data" / "mini_uhgg"
    run(tmp_path, "uhgg-pick", "--metadata", str(catalogue / "genomes-all_metadata.tsv"),
        "--replicates", "0", "--two-strain-replicates", "1", "--per-sample", "2",
        "--min-genomes", "3", "--distractors", "0")  # fmt: skip
    samples = pl.read_csv(tmp_path / "samples.tsv", separator="\t")
    # only beta has two eligible genomes (alpha's 3 is 80% complete)
    assert sorted(samples["genome"]) == ["MGYG000000005", "MGYG000000006"]

    si = tmp_path / "si"
    si.mkdir()
    pl.DataFrame({"species": [0, 1], "id": ["SA", "SB"], "name": ["A", "B"]}).write_csv(
        si / "species.tsv", separator="\t"
    )
    pl.DataFrame({"name": ["G1", "G2"], "species": ["SA"] * 2, "unit": [1, 2],
                  "c": [1.0] * 2}).write_parquet(si / "held_out.parquet")  # fmt: skip
    pl.DataFrame({"sample": [1, 1], "genome": ["G1", "G2"], "species": ["SA", "SA"],
                  "depth": [2.0, 3.0]}).write_csv(
        tmp_path / "truth.tsv", separator="\t")  # fmt: skip
    pl.DataFrame({"id": ["SA", "SA", "SB"], "lineage": [0, 1, 0],
                  "relative_abundance": [0.5, 0.3, 0.2]}).write_csv(
        tmp_path / "lineages.tsv", separator="\t")  # fmt: skip
    run(tmp_path, "uhgg-score", "--truth", "truth.tsv", "--species-index", "si", "--sample", "1",
        "--arm", "lineage", "--kind", "lineage", "--pred", "lineages.tsv")  # fmt: skip
    got = dict(pl.read_csv(tmp_path / "species_score.tsv", separator="\t")
               .select("metric", "value").iter_rows())  # fmt: skip
    assert got["n_truth"] == 1 and got["completeness"] == 1.0 and got["purity"] == 0.5
    assert got["l1"] == pytest.approx(0.4)
    assert got["two_strain_species"] == 1 and got["two_strain_resolved"] == 1.0
