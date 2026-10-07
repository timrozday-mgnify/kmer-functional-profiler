"""Steps of the fmh-funprofiler benchmark (``main.nf``) on its Zenodo 10055954 inputs.

- ``members``: KEGG proteins grouped by KO -> members table for ``kmer-functional-profiler
  index`` (a gene with several KOs is a member of each).
- ``sample``: draw genomes, write them as one FASTA with ``genome|contig`` names and their
  gene coordinates as ``genes.parquet``.
- ``truth``: map simulated read pairs back to the sample with minimap2 (mappy); a KO is
  present if a primary alignment overlaps one of its genes, as in the paper's CAMISIM
  ground truth. Writes ``label`` (the KO), ``n_reads``, ``bases`` (overlapping aligned
  bases) and ``depth`` (bases / gene length, summed over the KO's genes: read depth times
  copies). With ``--domains``, also the Pfam truth (``--out-pfam``): the same over Pfam
  domains' coordinates on the genome, not whole genes (Benchmark labels).
- ``pfam-proteins``: every gene's protein from the genomes' mapping tables, as FASTA.
- ``pfam-domains``: ``hmmsearch --domtblout`` -> ``domains.parquet`` (gene, Pfam, envelope
  in aa) and ``pfam_members.parquet``, the domains as members of one unit per Pfam (the
  Pfam analogue of the KO index).
- ``reps``: the representative of every MGnify90 cluster in members tables, as FASTA.
- ``mgnify-genes``: DIAMOND blastp of the genomes' proteins against those representatives ->
  ``gene_units.parquet``: per gene its best and near-hit clusters, with identity and
  coverage (the truth of unit-level detection and of containment AAI).
- ``aai-score``: a unit profile against that truth: detection at the MGnify90 level (genes
  whose best cluster passes ``--id``/``--cov``), recall of the nearest cluster beyond it, and
  ``aai`` / ``aai_naive`` against the alignment identities (see :func:`aai_score`). With
  ``--model`` (an ``aai_model.json``) or ``--min-aai-kmers``, ``aai`` is first re-estimated
  from the profile's zero-inflated fit under that survival model and mask
  (:func:`reestimate_aai`).
- ``aai-calibrate``: an inverse map from ``aai`` to alignment identity, fitted on half of
  the clusters of several samples' unit profiles and scored on the other half (see
  :func:`fit_aai_calibration`).
- ``pfam-profile``: a profile against an index with Pfam labels (``unit_pfam.parquet``, e.g.
  MGnify90 clusters) summed per Pfam, each unit counting for each of its Pfams; several
  label tables for a joint query (``--extra-index``), in the query's index order.
- ``study-proteins``: the study ladder's study (phase 10): the proteins of a random
  ``--fraction`` of a sample's genomes, as FASTA, with their Pfam domains.
- ``study-rebuild``: the members and Pfam tables of a rebuild that includes the study: base
  members plus the study's proteins, each in its nearest MGnify90 cluster where it passes
  ``--min-id``/``--min-cov``, else in its own linclust cluster.
- ``study-ladder``: base, base + study index and rebuild scores side by side, the share of
  the rebuild's completeness gain the joint query recovers, and whether base units the
  study does not touch changed.
- ``score``: purity and completeness of one profile against a truth table, one row per
  count present (``kmers_hit``; ``kmers_unique`` after gather; ``kmers_wta`` and
  ``kmers_ufirst`` after winner-take-all and uniqueness-first) and ``--min-hits``
  threshold. Each count is paired with its abundance estimates (``RULES``), scored on the detected
  KOs against ``depth``: Spearman over true positives and L1 between relative abundances,
  and the calibration of posterior intervals where the profile has them.
- ``detected``: one row per KO gather keeps (``kmers_unique`` >= 1), true or false, with
  its evidence (``present_prob``, ``own_evidence``, ...) and its tier-2 hit k-mers summarised:
  median holders, most hits on one k-mer, and the share found in the sample genomes' six
  frames (``in_genome``; null for sourmash-hashed indexes). A false positive's k-mers that
  are in the genomes come from real sequence (another gene or KO); the rest from read errors.
- ``tool-profile``: another tool's output (DIAMOND, fmh-funprofiler, kMermaid, HUMAnN 3/4) as
  a profile SCORE reads: ``name``, ``evidence`` (the tool's detection count) and ``abundance``.
- ``summary``: mean and sd of the scores per index, count and threshold.
- ``cost``: wall time, CPU time and peak memory per step and index or tool, from the raw
  Nextflow trace.
- ``genome-set``, ``genome-truth``, ``genome-score``: genome mode (phase 11): every genome
  of the record as protein FASTA for ``annotate-genomes``; per-sample genome and
  (genome, function) truth from the per-gene truth; and a genome fit (or sylph) scored
  against it.
- ``subsample``, ``prior-score``: kfp-prior's depth ladder (phase 11): read pairs kept at
  a fraction; unit and Pfam presence, observed vs updated, by carrier depth, and the
  calibration of zero-hit units.
- ``iss``: ``iss`` with its arguments, the perfect error model patched (see :func:`iss`).
"""

import argparse
import gzip
import itertools
import json
import random
import sys
from collections.abc import Iterable, Iterator
from multiprocessing.pool import ThreadPool
from pathlib import Path
from typing import IO, Any, Final

import mappy
import numpy as np
import polars as pl

from kmer_functional_profiler import _core
from kmer_functional_profiler.query import (
    AAI_CALIBRATION_PARAMS,
    MIN_AAI_KMERS,
    aai_columns,
    calibrate_aai,
    check_aai_calibration,
    unit_windows,
    ztp_lambda,
)
from kmer_functional_profiler.survival import SurvivalModel

GENOME_COLUMNS = {"gene_name": pl.String, "contig_id": pl.String, "start_position": pl.Int64,
                  "end_position": pl.Int64, "strand": pl.String}  # fmt: skip


def _open(path: str | Path) -> IO[str]:
    """A text file, gzip-compressed if its name ends in .gz."""
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)  # noqa: SIM115


def fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """(header, sequence) records; headers keep everything after ``>``."""
    header, chunks = None, []
    with _open(path) as f:
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def read_kos(path: str | Path) -> pl.DataFrame:
    return pl.read_csv(path, columns=["gene_id", "ko_id"]).drop_nulls().unique()


def members(args: argparse.Namespace) -> None:
    proteins = pl.DataFrame(
        [(h.split("|", 1)[0], s) for h, s in fasta(args.faa)],
        schema=["gene_id", "sequence"],
        orient="row",
    )
    (
        proteins.join(read_kos(args.kos), on="gene_id")
        .select(
            protein_id="gene_id", cluster_rep="ko_id", full_length=pl.lit(True), sequence="sequence"
        )
        .write_parquet(args.out)
    )


def sample(args: argparse.Namespace) -> None:
    genomes = sorted(
        p.name for p in Path(args.genomes_dir).iterdir() if (p / f"{p.name}.fasta").exists()
    )
    chosen = sorted(random.Random(args.seed).sample(genomes, min(args.n, len(genomes))))
    with open("sample.fna", "w") as out:
        for g in chosen:
            for header, seq in fasta(Path(args.genomes_dir) / g / f"{g}.fasta"):
                out.write(f">{g}|{header.split()[0]}\n{seq}\n")
    pl.concat(
        pl.read_csv(
            Path(args.genomes_dir) / g / f"{g}_mapping.csv",
            columns=list(GENOME_COLUMNS),
            schema_overrides=GENOME_COLUMNS,
        ).with_columns(genome=pl.lit(g))  # fmt: skip
        for g in chosen
    ).select(
        "gene_name",
        contig=pl.col("genome") + "|" + pl.col("contig_id"),
        start=pl.min_horizontal("start_position", "end_position") - 1,  # 0-based, half-open
        end=pl.max_horizontal("start_position", "end_position"),
        strand="strand",
    ).write_parquet("genes.parquet")
    Path("genomes.txt").write_text("\n".join(chosen) + "\n")


# (detection count, abundance estimate) pairs scored: naive coverage for every hit, plain and
# zero-inflated EM for the units gather keeps, and the hits winner-take-all / uniqueness-first
# assign.
RULES = (
    ("kmers_hit", "coverage"),
    ("kmers_unique", "coverage_em"),
    ("kmers_unique", "coverage_zi"),
    ("kmers_unique", "coverage_zib"),
    ("kmers_unique", "coverage_zip"),
    ("kmers_unique", "abundance_zi"),
    ("kmers_wta", "coverage_wta"),
    ("kmers_ufirst", "coverage_ufirst"),
    ("evidence", "abundance"),  # other tools (tool_profile)
)
_BIN = 1000  # bp; overlap candidates are only compared within a shared bin


def _map(
    aligner: mappy.Aligner, chunk: list[tuple[int, str, str]]
) -> list[tuple[int, str, int, int]]:
    buf = mappy.ThreadBuffer()
    return [
        (read, h.ctg, h.r_st, h.r_en)
        for read, s1, s2 in chunk
        for h in aligner.map(s1, s2, buf=buf)
        if h.is_primary
    ]


def _binned(df: pl.DataFrame, start: str, end: str) -> pl.DataFrame:
    """One row per ``_BIN``-sized window the [start, end) interval touches."""
    return df.with_columns(
        bin=pl.int_ranges(pl.col(start) // _BIN, (pl.col(end) - 1) // _BIN + 1)
    ).explode("bin")


def _pairs(r1: str, r2: str, size: int = 10_000) -> Iterator[list[tuple[int, str, str]]]:
    chunk = []
    for read, ((_, s1, _), (_, s2, _)) in enumerate(
        zip(mappy.fastx_read(r1), mappy.fastx_read(r2), strict=True)
    ):
        chunk.append((read, s1, s2))
        if len(chunk) == size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def _truth(aligned: pl.DataFrame, features: pl.DataFrame) -> pl.DataFrame:
    """Per ``label``: reads, aligned bases and depth over ``features`` (``contig``, ``start``,
    ``end``, ``label``) that primary alignments overlap; depth sums bases / feature length."""
    # join_where on contig + ranges is a contig hash join then a filter (reads x genes per
    # contig, billions of rows on complete genomes); joining on (contig, bin) keeps it local
    overlaps = (
        _binned(aligned, "r_start", "r_end")
        .join(_binned(features.unique(), "start", "end"), on=["contig", "bin"])
        .filter(pl.col("r_start") < pl.col("end"), pl.col("r_end") > pl.col("start"))
        .drop("bin")
        .unique()
    )
    return (
        overlaps.with_columns(
            bases=pl.min_horizontal("r_end", "end") - pl.max_horizontal("r_start", "start")
        )
        .group_by("label")
        .agg(
            n_reads=pl.col("read").n_unique(),
            bases=pl.col("bases").sum(),
            depth=(pl.col("bases") / (pl.col("end") - pl.col("start"))).sum(),
        )
        .sort("label")
    )


def domain_features(genes: pl.DataFrame, domains: pl.DataFrame) -> pl.DataFrame:
    """Pfam domains (aa envelopes) of ``genes`` as genome intervals: on the reverse strand
    residue i of a protein lies at [end - 3(i + 1), end - 3i)."""
    forward = pl.col("strand") != "-"
    return genes.join(domains, on="gene_name").select(
        "contig",
        start=pl.when(forward)
        .then(pl.col("start") + 3 * pl.col("aa_start"))
        .otherwise(pl.col("end") - 3 * pl.col("aa_end")),
        end=pl.when(forward)
        .then(pl.col("start") + 3 * pl.col("aa_end"))
        .otherwise(pl.col("end") - 3 * pl.col("aa_start")),
        label="pfam",
    )


def truth(args: argparse.Namespace) -> None:
    # mappy releases the GIL while aligning, so threads share one index instead of one per process
    aligner = mappy.Aligner(args.fna, preset="sr")
    with ThreadPool(args.threads) as pool:
        parts = pool.imap(lambda chunk: _map(aligner, chunk), _pairs(args.r1, args.r2))
        hits = [h for part in parts for h in part]
    schema = {"read": pl.Int64, "contig": pl.String, "r_start": pl.Int64, "r_end": pl.Int64}
    aligned = pl.DataFrame(hits, schema=schema, orient="row")
    genes = pl.read_parquet(args.genes)
    kos = genes.join(read_kos(args.kos), left_on="gene_name", right_on="gene_id").select(
        "contig", "start", "end", label="ko_id"
    )
    _truth(aligned, kos).write_csv(args.out)
    gene_features = genes.select("contig", "start", "end", label="gene_name")
    _truth(aligned, gene_features).rename({"label": "gene_name"}).write_csv(args.out_genes)
    if args.domains:
        features = domain_features(genes, pl.read_parquet(args.domains))
        _truth(aligned, features).write_csv(args.out_pfam)


def read_profile(path: str | Path) -> pl.DataFrame:
    """A profile TSV; an empty one (nothing detected) gets numeric columns, not strings."""
    # whole-file inference: a run of 100+ empty aai cells would otherwise infer str
    profile = pl.read_csv(path, separator="\t", infer_schema_length=None)
    # ponytail: the benchmark's only extra index is the host decoy, whose units are not
    # predictions; drop this when `index --role decoy` (plan, phase 10) leaves them out itself
    if "source" in profile.columns:
        profile = profile.filter(pl.col("source") == 0)
    if profile.height:
        return profile
    strings = {"name", "cluster_rep"}
    return profile.cast({c: pl.UInt32 if c == "unit" else pl.Float64
                         for c in profile.columns if c not in strings})  # fmt: skip


def mix(args: argparse.Namespace) -> None:
    """A host spike-in sample at constant depth: each microbial pair kept with probability
    1 - ``fraction``, then host pairs, ``fraction`` x the microbial pairs, from the start of
    the host reads. Truth is recomputed on the mixed reads (host reads map to no genome)."""
    rng = random.Random(args.seed)
    n = 0
    with (
        gzip.open(args.r1, "rt") as m1, gzip.open(args.r2, "rt") as m2,
        gzip.open("mixed_R1.fastq.gz", "wt", compresslevel=1) as o1,
        gzip.open("mixed_R2.fastq.gz", "wt", compresslevel=1) as o2,
    ):  # fmt: skip
        for rec1, rec2 in zip(_records(m1), _records(m2), strict=True):
            n += 1
            if rng.random() >= args.fraction:
                o1.write(rec1)
                o2.write(rec2)
        want = round(n * args.fraction)
        with gzip.open(args.host_r1, "rt") as h1, gzip.open(args.host_r2, "rt") as h2:
            got = 0
            for rec1, rec2 in itertools.islice(zip(_records(h1), _records(h2), strict=True), want):
                o1.write(rec1)
                o2.write(rec2)
                got += 1
    if got < want:
        raise ValueError(f"{got} host pairs, {want} wanted: simulate more host reads")


def _records(f: Iterable[str]) -> Iterator[str]:
    """FASTQ records (four lines each) as strings."""
    lines = iter(f)
    for header in lines:
        yield header + next(lines) + next(lines) + next(lines)


def host_abundance(args: argparse.Namespace) -> None:
    """iss abundance file for a host genome: reads in proportion to record length, the
    mitochondrion at ``--mito-copies`` copies and PhiX at ``--phix`` of the reads."""
    lengths: dict[str, int] = {}
    with _open(args.fasta) as f:
        name = None
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]
                lengths[name] = 0
            elif name is not None:
                lengths[name] += len(line.strip())
    weight = {n: n_bp * (args.mito_copies if n in MITO else 1) for n, n_bp in lengths.items()
              if n not in PHIX}  # fmt: skip
    total = sum(weight.values())
    rows = [(n, (1 - args.phix) * w / total) for n, w in weight.items()]
    rows += [(n, args.phix) for n in lengths if n in PHIX]
    Path(args.out).write_text("".join(f"{n}\t{a}\n" for n, a in rows))


MITO = {"chrM", "NC_012920.1"}
PHIX = {"NC_001422.1"}


def decoy_members(args: argparse.Namespace) -> None:
    """A proteome as a members table, one unit per protein (decoy indexes)."""
    rows = [(i, i, True, seq) for i, (_, seq) in enumerate(fasta(args.faa))]
    schema = ["protein_id", "cluster_rep", "full_length", "sequence"]
    pl.DataFrame(rows, schema=schema, orient="row").write_parquet(args.out)


def read_truth(path: str | Path) -> pl.DataFrame:
    """A truth table, with ``label`` (truth files before Pfam called it ``ko_id``)."""
    return pl.read_csv(path).rename({"ko_id": "label"}, strict=False)


def pfam_proteins(args: argparse.Namespace) -> None:
    """Every gene's protein, dealt round-robin into ``--chunks`` files for parallel search."""
    outs = [open(f"proteins.{i}.faa", "w") for i in range(args.chunks)]  # noqa: SIM115
    n = 0
    for path in sorted(Path(args.genomes_dir).glob("*/*_mapping.csv")):
        genes = pl.read_csv(path, columns=["gene_name", "aa_sequence"])
        for name, seq in genes.unique("gene_name", maintain_order=True).drop_nulls().iter_rows():
            outs[n % args.chunks].write(f">{name}\n{seq}\n")
            n += 1
    for out in outs:
        out.close()


def pfam_domains(args: argparse.Namespace) -> None:
    rows = [
        (f[0], f[4].split(".")[0], int(f[19]) - 1, int(f[20]))  # target, accession, envelope
        for path in args.domtbl
        for line in Path(path).read_text().splitlines()
        if not line.startswith("#")
        for f in [line.split()]
    ]
    schema = {"gene_name": pl.String, "pfam": pl.String, "aa_start": pl.Int64, "aa_end": pl.Int64}
    domains = pl.DataFrame(rows, schema=schema, orient="row").unique().sort(list(schema))
    domains.write_parquet(args.out)
    proteins = pl.DataFrame(
        [r for path in args.proteins for r in fasta(path)],
        schema=["gene_name", "seq"],
        orient="row",
    )
    domains.join(proteins, on="gene_name", maintain_order="left").select(
        protein_id=pl.format("{}:{}-{}", "gene_name", "aa_start", "aa_end"),
        cluster_rep="pfam",
        full_length=pl.lit(True),
        sequence=pl.col("seq").str.slice(pl.col("aa_start"), pl.col("aa_end") - pl.col("aa_start")),
    ).write_parquet(args.out_members)


def reps(args: argparse.Namespace) -> None:
    """Cluster representatives (``protein_id == cluster_rep``) of members tables (files, or
    directories of them, e.g. the MERGE buckets) as FASTA named by ``cluster_rep``."""
    paths = [str(Path(p) / "*.parquet") if Path(p).is_dir() else p for p in args.members]
    # streamed: the whole release has ~1.7e9 clusters
    pl.scan_parquet(paths).filter(pl.col("protein_id") == pl.col("cluster_rep")).select(
        pl.format(">{}\n{}", "cluster_rep", "sequence")
    ).sink_csv(args.out, include_header=False, quote_style="never")


MGNIFY_HIT_SCHEMA = {
    "gene_name": pl.String,
    "cluster_rep": pl.Int64,
    "pident": pl.Float64,
    "length": pl.Int64,
    "qlen": pl.Int64,
    "slen": pl.Int64,
    "bitscore": pl.Float64,  # DIAMOND prints e.g. `99.8`; never infer it from the first rows
}


def mgnify_genes(args: argparse.Namespace) -> None:
    """DIAMOND hits (outfmt 6 ``qseqid sseqid pident length qlen slen bitscore``) -> per gene
    and cluster: identity (0-1), query and subject coverage, bitscore and rank (1 = best)."""
    hits = pl.concat(
        pl.read_csv(p, separator="\t", has_header=False, schema=MGNIFY_HIT_SCHEMA)
        for p in args.hits
    )  # fmt: skip
    (
        hits.group_by("gene_name", "cluster_rep")  # one row per pair: its best HSP
        .agg(pl.all().sort_by("bitscore").last())
        .select(
            "gene_name",
            "cluster_rep",
            identity=pl.col("pident") / 100,
            qcov=pl.col("length") / pl.col("qlen"),
            scov=pl.col("length") / pl.col("slen"),
            bitscore="bitscore",
        )
        .with_columns(
            rank=pl.col("bitscore")
            .rank("ordinal", descending=True)
            .over("gene_name")
            .cast(pl.UInt32)
        )
        .sort("gene_name", "rank")
        .write_parquet(args.out)
    )


IDENTITY_BINS = ((0.95, 1.01), (0.9, 0.95), (0.8, 0.9), (0.7, 0.8), (0.5, 0.7))
LENGTH_BINS = ((0, 150), (150, 300), (300, 600), (600, float("inf")))  # gene length, aa
AAI_KMER_BINS = ((0, 3), (3, 5), (5, 10), (10, float("inf")))  # aai_kmers: 1-2, 3-4, 5-9, 10+


def aai_score(
    profile: pl.DataFrame,
    gene_units: pl.DataFrame,
    genes: pl.DataFrame,
    min_id: float = 0.9,
    min_cov: float = 0.8,
    k: int = 11,
) -> dict[str, float | None]:
    """Unit-level detection and containment AAI against alignment truth.

    ``profile`` is a unit profile (``cluster_rep``, ``kmers_unique``, ``aai``, ``aai_lo``,
    ``aai_hi``, ``aai_naive``); ``gene_units`` the hits of :func:`mgnify_genes`; ``genes``
    the per-gene ``depth`` of the sample (genes with reads are present). A present gene's
    best hit is its *nearest cluster*; it is *in* that cluster at the MGnify90 level if the
    hit has identity >= ``min_id`` and both coverages >= ``min_cov``.

    - ``completeness_90``, ``purity_nearest``: share of the clusters that hold a present
      gene (at the 90% level) that are detected (``kmers_unique`` >= 1); share of detected
      units that are some present gene's nearest cluster (``purity_near``: any hit). By the
      nearest cluster's identity (``IDENTITY_BINS``, ``recall_<lo>``): share of present
      genes whose nearest cluster is detected, also below 90% (resolution beyond it).
    - ``completeness_90_len<lo>``, ``abund_bias_len<lo>``, ``abund_err_len<lo>``: by gene
      length (``LENGTH_BINS``, aa; needs ``bases`` in ``genes``), the share of present genes
      in their cluster at the 90% level whose cluster is detected; and, for detected units
      that are the 90%-level cluster of exactly one present gene and nearest to no other,
      the median and median absolute log2 error of ``coverage_em`` / depth after removing
      the sample's median log2 ratio (a scale-free length bias).
    - ``aai_bias_<lo>``, ``aai_cover_<lo>``, ``aai_spearman``, ``aai_n``: detected units
      with an ``aai``, against the depth-weighted identity of the present genes they are
      nearest to (``aai_mixed``: units nearest to several genes).
    - ``aai_bias_union_<lo>``, ``aai_cover_union_<lo>``, ``aai_cover_union``: the same against
      the union truth of every present gene hitting the unit (:func:`union_truth`), binned
      by it; the interval coverage over all units is the primary ``aai`` metric (step 23).
    - ``aai_bias_kmers<lo>``, ``aai_cover_kmers<lo>``, ``aai_n_kmers<lo>``: against the union
      truth, by the hit k-mers the estimate rests on (``aai_kmers``, ``AAI_KMER_BINS``):
      near the detection limit the reported units are the lucky draws.
    - ``naive_bias_<lo>``, ``naive_within05_<lo>``, ``naive_spearman``, ``naive_n``: every
      profiled unit some present gene hits, ``aai_naive`` against the best identity of a
      present gene to it (near hits included).
    """
    present = genes.filter(pl.col("depth") > 0).select("gene_name", "depth")
    hits = gene_units.join(present, on="gene_name").filter(pl.col("qcov") >= 0.5)
    nearest = hits.filter(pl.col("rank") == 1)
    in90 = nearest.filter(
        pl.col("identity") >= min_id, pl.col("qcov") >= min_cov, pl.col("scov") >= min_cov
    )
    detected = set(profile.filter(pl.col("kmers_unique") >= 1)["cluster_rep"].to_list())
    out: dict[str, float | None] = {
        "genes_present": present.height,
        "genes_with_hit": nearest.height,
        "genes_in90": in90.height,
    }
    clusters90 = set(in90["cluster_rep"].to_list())
    out["completeness_90"] = len(clusters90 & detected) / len(clusters90) if clusters90 else None
    out["purity_nearest"] = (
        len(detected & set(nearest["cluster_rep"].to_list())) / len(detected) if detected else None
    )
    out["purity_near"] = (
        len(detected & set(hits["cluster_rep"].to_list())) / len(detected) if detected else None
    )
    for lo, hi in IDENTITY_BINS:
        part = nearest.filter(pl.col("identity").is_between(lo, hi, closed="left"))
        out[f"recall_{lo}"] = (
            part["cluster_rep"].is_in(list(detected)).mean() if part.height else None  # type: ignore[assignment]
        )
    # by gene length (aa; depth is bases per nt of gene): the stop filter's end loss falls on
    # short genes (plan, Translation and frame detection)
    if "bases" in genes.columns:
        length = genes.filter(pl.col("depth") > 0).select(
            "gene_name", aa=pl.col("bases") / pl.col("depth") / 3
        )
        g90 = in90.join(length, on="gene_name")
        one = (  # units that are one present gene's 90%-level cluster and nearest to no other
            g90.join(nearest.group_by("cluster_rep").len(), on="cluster_rep").filter(
                pl.col("len") == 1
            )
        )
        if "coverage_em" in profile.columns:
            one = one.join(
                profile.filter(pl.col("kmers_unique") >= 1).select("cluster_rep", "coverage_em"),
                on="cluster_rep",
            ).with_columns(ratio=(pl.col("coverage_em") / pl.col("depth")).log(2))
            one = one.with_columns(err=pl.col("ratio") - pl.col("ratio").median())
        for lo, hi in LENGTH_BINS:
            part = g90.filter(pl.col("aa").is_between(lo, hi, closed="left"))
            out[f"completeness_90_len{lo}"] = (
                part["cluster_rep"].is_in(list(detected)).mean() if part.height else None  # type: ignore[assignment]
            )
            if "err" in one.columns:
                e = one.filter(pl.col("aa").is_between(lo, hi, closed="left"))["err"]
                out[f"abund_bias_len{lo}"] = e.median() if e.len() else None  # type: ignore[assignment]
                out[f"abund_err_len{lo}"] = e.abs().median() if e.len() else None  # type: ignore[assignment]
    # aai of detected units, against the identity of the genes they are nearest to
    truth = aai_truth(nearest)
    if "aai" in profile.columns:
        a = (
            profile.filter(pl.col("kmers_unique") >= 1, pl.col("aai").is_not_null())
            .join(truth, on="cluster_rep")
            .join(union_truth(hits, k), on="cluster_rep")
        )
        out["aai_n"] = a.height
        out["aai_mixed"] = int((a["genes"] > 1).sum())
        out["aai_spearman"] = (
            a.select(pl.corr("aai", "true", method="spearman")).item() if a.height > 2 else None
        )

        def bias_cover(part: pl.DataFrame, true: str, key: str) -> None:
            out[f"aai_bias_{key}"] = (part["aai"] - part[true]).median() if part.height else None  # type: ignore[assignment]
            out[f"aai_cover_{key}"] = (
                part.select(pl.col(true).is_between("aai_lo", "aai_hi").mean()).item()
                if part.height and "aai_lo" in part.columns
                else None
            )

        for true, prefix in (("true", ""), ("true_union", "union_")):
            for lo, hi in IDENTITY_BINS:
                bias_cover(
                    a.filter(pl.col(true).is_between(lo, hi, closed="left")), true, f"{prefix}{lo}"
                )
        bias_cover(a, "true_union", "union")
        del out["aai_bias_union"]  # over every identity it says little; its coverage does
        if "aai_kmers" in a.columns:
            for lo, hi in AAI_KMER_BINS:
                part = a.filter(pl.col("aai_kmers").is_between(lo, hi, closed="left"))
                out[f"aai_n_kmers{lo}"] = part.height
                bias_cover(part, "true_union", f"kmers{lo}")
    if "aai_naive" in profile.columns:
        best = hits.group_by("cluster_rep").agg(true=pl.col("identity").max())
        n = profile.join(best, on="cluster_rep")
        out["naive_n"] = n.height
        out["naive_spearman"] = (
            n.select(pl.corr("aai_naive", "true", method="spearman")).item()
            if n.height > 2
            else None
        )
        for lo, hi in IDENTITY_BINS:
            part = n.filter(pl.col("true").is_between(lo, hi, closed="left"))
            err = part["aai_naive"] - part["true"]
            out[f"naive_bias_{lo}"] = err.median() if part.height else None  # type: ignore[assignment]
            out[f"naive_within05_{lo}"] = (err.abs() <= 0.05).mean() if part.height else None
    return out


def union_truth(hits: pl.DataFrame, k: int = 11) -> pl.DataFrame:
    """Per cluster, the identity whose k-mer survival a^k equals that of the union of every
    present gene hitting it (``true_union``): 1 - prod_g (1 - a_g^k), inverted. Reads from
    several genes put their k-mers in one unit, which the nearest-gene truth ignores (plan,
    phase 7, step 23). a^k, not the fitted survival model, so the truth does not depend on
    the model it scores."""
    return hits.group_by("cluster_rep").agg(
        true_union=(1 - (1 - pl.col("identity") ** k).log().sum().exp()) ** (1 / k)
    )


def reestimate_aai(
    profile: pl.DataFrame, model: SurvivalModel, min_kmers: float = MIN_AAI_KMERS
) -> pl.DataFrame:
    """``aai``, ``aai_lo``, ``aai_hi`` and ``aai_kmers`` recomputed from a profile's
    zero-inflated fit as the query would report them with survival model ``model`` and
    ``--min-aai-kmers`` ``min_kmers``, without rerunning it (offline reruns of archived
    profiles, plan, phase 7, step 23). ``pin_sum`` is present x m / ``copies_zi``; the
    clumping is the query's, from tier-2 hits per read; m is ``m_dense`` with a dense tier.
    Units without a fit (``copies_zi`` 0) get a null ``aai`` and ``aai_kmers`` 0."""
    m = profile["m_dense" if "m_dense" in profile.columns else "m_g"].to_numpy().astype(np.float64)
    present = profile["present_zi"].to_numpy()
    copies = profile["copies_zi"].to_numpy()
    fitted = copies > 0
    pin_sum = np.where(fitted, present * m / np.where(fitted, copies, 1.0), 0.0)
    mu = ztp_lambda((profile["hits"] / profile["reads"]).to_numpy())
    mu *= m / profile["m_g"].to_numpy()
    est = aai_columns(
        profile["coverage_zi"].to_numpy(), present, m, pin_sum,
        profile["n_kmers"].to_numpy(), unit_windows(profile, model.k), model, 1 + mu, min_kmers,
    )  # fmt: skip
    est = est.with_columns(
        pl.when(pl.Series(fitted)).then(pl.col(c)).alias(c) for c in ("aai", "aai_lo", "aai_hi")
    )
    return profile.drop("aai", "aai_lo", "aai_hi", "aai_kmers", strict=False).hstack(est)


def aai_truth(nearest: pl.DataFrame) -> pl.DataFrame:
    """Per cluster, the depth-weighted identity (``true``) of the present genes whose nearest
    cluster it is, and how many (``genes``)."""
    return nearest.group_by("cluster_rep").agg(
        true=(pl.col("identity") * pl.col("depth")).sum() / pl.col("depth").sum(),
        genes=pl.len(),
    )


def isotonic(y: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Non-decreasing weighted least-squares fit to ``y`` in its given order (pool adjacent
    violators)."""
    vals: list[float] = []
    wts: list[float] = []
    runs: list[int] = []
    for yi, wi in zip(y.tolist(), w.tolist(), strict=True):
        vals.append(yi)
        wts.append(wi)
        runs.append(1)
        while len(vals) > 1 and vals[-2] > vals[-1]:
            wt = wts[-2] + wts[-1]
            vals[-2] = (vals[-2] * wts[-2] + vals[-1] * wts[-1]) / wt
            wts[-2], runs[-2] = wt, runs[-2] + runs[-1]
            del vals[-1], wts[-1], runs[-1]
    return np.repeat(vals, runs)


def fit_aai_calibration(
    pairs: pl.DataFrame,
    step: float = 0.01,
    min_n: int = 20,
    level: float = 0.95,
    floor: float = 0.7,
) -> dict[str, Any]:
    """Inverse calibration of ``aai`` against alignment identity (Sequence similarity).

    ``pairs`` has ``aai``, ``aai_lo``, ``aai_hi`` and ``true`` per unit. g(a), the median
    ``aai`` of units at true identity a (bins of ``step`` with >= ``min_n`` units, at their
    mean identity), is made
    non-decreasing and inverted: ``identity`` = g^-1(``aai``), linear between knots, a flat
    run of g mapping to its mid identity. Inverting g rather than regressing identity on
    ``aai`` keeps the benchmark's identity distribution out of the map (a regression would
    pull every estimate toward the benchmark's mean identity). ``widen`` is the error the
    sampling interval misses: the ``level`` quantile, over units at identity >= ``floor``,
    of how far ``true`` lies outside the mapped interval. With fewer than two knots there
    is no map: ``widen`` is None and ``calibrate-aai`` refuses it."""
    b = (
        pairs.group_by((pl.col("true") / step).floor().alias("bin"))
        .agg(mid=pl.col("true").mean(), med=pl.col("aai").median(), n=pl.len())
        .filter(pl.col("n") >= min_n)
        .sort("bin")
    )
    g = isotonic(b["med"].to_numpy(), b["n"].to_numpy().astype(np.float64))
    knots, run = np.unique(g, return_inverse=True)
    identity = np.bincount(run, b["mid"].to_numpy()) / np.bincount(run)
    cal: dict[str, Any] = {"aai": knots.tolist(), "identity": identity.tolist()}
    if len(knots) < 2:  # too few units for a map (e.g. a test fixture)
        return cal | {"widen": None}
    mapped = calibrate_aai(pairs.filter(pl.col("true") >= floor), {**cal, "widen": 0.0})
    miss = mapped.select(
        pl.max_horizontal(pl.col("aai_lo") - pl.col("true"), pl.col("true") - pl.col("aai_hi"), 0)
    ).to_series()
    cal["widen"] = float(miss.quantile(level, "higher") or 0.0)
    return cal


def aai_calibrate(args: argparse.Namespace) -> None:
    """Fit :func:`fit_aai_calibration` on half the clusters (by a hash of ``cluster_rep``) of
    every sample's unit profile and score ``aai`` raw and calibrated (:func:`aai_score`) on
    the other half, so no cluster is in both (genes differ by sample, clusters recur). The
    JSON is what ``kmer-functional-profiler calibrate-aai`` attaches to the index; profiles
    of an index already calibrated are refitted on their raw ``aai``. With ``--map`` the
    given map (fitted on another run) is scored on the same held-out half instead of a new
    fit: the transfer test of plan, phase 7, steps 17 and 19."""
    gene_units = pl.read_parquet(args.gene_units)
    parts = []
    for path, genes in zip(args.profiles, args.genes, strict=True):
        truth = pl.read_csv(genes)
        nearest = gene_units.join(
            truth.filter(pl.col("depth") > 0).select("gene_name", "depth"), on="gene_name"
        ).filter(pl.col("qcov") >= 0.5, pl.col("rank") == 1)
        profile = read_profile(path).filter(pl.col("kmers_unique") >= 1)
        if "aai_raw" in profile.columns:  # an index already calibrated: refit on the raw aai
            raw = ("aai_raw", "aai_raw_lo", "aai_raw_hi")
            profile = profile.drop("aai", "aai_lo", "aai_hi").rename(
                dict(zip(raw, ("aai", "aai_lo", "aai_hi"), strict=True))
            )
        parts.append((profile, truth, nearest))
    # a hash, not cluster_rep's parity: nested subsets keep every n-th accession
    held_out = pl.col("cluster_rep").hash(0) % 2 == 1  # stable within a Polars version
    train = pl.concat(
        p.filter(pl.col("aai").is_not_null(), ~held_out)
        .join(aai_truth(n), on="cluster_rep")
        .select("aai", "aai_lo", "aai_hi", "true")
        for p, _, n in parts
    )
    params = {}
    if args.index:
        params = json.loads((Path(args.index) / "meta.json").read_text())["params"]
    if args.map:
        cal = json.loads(Path(args.map).read_text())
        if cal["widen"] is not None:
            check_aai_calibration(cal, params)  # a map from an index built otherwise is refused
    else:
        cal = fit_aai_calibration(train)
        if params:  # the build parameters the map holds for (checked when it is attached)
            cal["params"] = {p: params.get(p) for p in AAI_CALIBRATION_PARAMS}
        cal["fit"] = {"n_train": train.height, "profiles": [Path(p).name for p in args.profiles]}
        Path(args.out).write_text(json.dumps(cal, indent=1))
    rows = []
    for (profile, genes, _), path in zip(parts, args.profiles, strict=True):
        test = profile.filter(held_out)
        fitted = [("calibrated", calibrate_aai(test, cal))] if cal["widen"] is not None else []
        for method, p in [("raw", test), *fitted]:
            got = aai_score(p, gene_units, genes, args.min_id, args.min_cov)
            rows.append({"profile": Path(path).name, "method": method}
                        | {k: v for k, v in got.items() if k.startswith("aai_")})  # fmt: skip
    pl.DataFrame(rows).write_csv(args.scores_out, separator="\t")


def aai_score_step(args: argparse.Namespace) -> None:
    profile = read_profile(args.profile)
    if args.model is not None or args.min_aai_kmers is not None:
        mask = MIN_AAI_KMERS if args.min_aai_kmers is None else args.min_aai_kmers
        fitted = None if args.model is None else json.loads(Path(args.model).read_text())
        profile = reestimate_aai(profile, SurvivalModel.from_json(fitted, args.k), mask)
    row = aai_score(
        profile,
        pl.read_parquet(args.gene_units),
        pl.read_csv(args.genes),
        args.min_id,
        args.min_cov,
        args.k,
    )
    pl.DataFrame([{"sample": args.sample, "index": args.index, "arm": args.arm, **row}]).write_csv(
        args.out, separator="\t"
    )


def pfam_names(labels: pl.DataFrame) -> pl.DataFrame:
    """``unit`` and ``name`` (``PF01007``) of a ``unit_pfam`` table."""
    accession = pl.col("pfam_accession")
    return labels.select(
        "unit",
        # MGnify stores the accession's number (1007 for PF01007); hmmsearch gives PF01007.23
        name=("PF" + accession.cast(pl.String).str.zfill(5))
        if labels.schema["pfam_accession"].is_integer()
        else accession.str.replace(r"\.\d+$", ""),
    )


def pfam_profile(args: argparse.Namespace) -> None:
    """Sum a unit profile per Pfam: counts and point estimates add over the units carrying
    a Pfam; ``present_prob`` is the largest; intervals and per-unit columns are dropped."""
    profile = read_profile(args.profile)
    # several tables: a joint query, whose unit ids are offset by the units of the indexes
    # before (each table's units.parquet beside it)
    offsets = itertools.accumulate(
        (pl.scan_parquet(Path(p).with_name("units.parquet")).select(pl.len()).collect().item()
         for p in args.unit_pfam[:-1]),
        initial=0,
    )  # fmt: skip
    labels = pl.concat(
        pfam_names(pl.read_parquet(p)).with_columns(pl.col("unit").cast(pl.Int64) + offset)
        for p, offset in zip(args.unit_pfam, offsets, strict=True)
    ).cast({"unit": profile.schema["unit"]})
    summed = [c for pair in RULES for c in pair if c in profile.columns] + ["hits"]
    (
        profile.drop("name", strict=False)
        .join(labels, on="unit")
        .group_by("name")
        .agg(
            pl.col(list(dict.fromkeys(c for c in summed if c in profile.columns))).sum(),
            *([pl.col("present_prob").max()] if "present_prob" in profile.columns else []),
        )
        .sort("name")
        .write_csv(args.out, separator="\t")
    )


def study_proteins(args: argparse.Namespace) -> None:
    """A study for the study ladder: ``--fraction`` of the sample's genomes (at least one),
    drawn with ``--seed``, as if recovered as MAGs. Writes their proteins (``study.faa``,
    named by gene), Pfam domains (``study_pfam.parquet``: ``protein_id``,
    ``pfam_accession``) and names (``study_genomes.txt``)."""
    genomes = sorted(
        pl.read_parquet(args.genes)["contig"].str.split("|").list.first().unique().to_list()
    )
    n = max(1, round(args.fraction * len(genomes)))
    chosen = sorted(random.Random(args.seed).sample(genomes, n))
    proteins = pl.concat(
        pl.read_csv(Path(args.genomes_dir) / g / f"{g}_mapping.csv",
                    columns=["gene_name", "aa_sequence"], schema_overrides=GENOME_COLUMNS)
        for g in chosen
    ).drop_nulls().unique("gene_name", keep="first", maintain_order=True)  # fmt: skip
    with open(args.out, "w") as out:
        for name, seq in proteins.iter_rows():
            out.write(f">{name}\n{seq}\n")
    (
        pl.read_parquet(args.domains)
        .filter(pl.col("gene_name").is_in(proteins["gene_name"].implode()))
        .select(protein_id="gene_name", pfam_accession="pfam")
        .unique()
        .sort("protein_id", "pfam_accession")
        .write_parquet(args.out_pfam)
    )
    Path(args.out_genomes).write_text("\n".join(chosen) + "\n")


def _member_paths(paths: list[str]) -> list[str]:
    return [str(Path(p) / "*.parquet") if Path(p).is_dir() else p for p in paths]


def study_rebuild(args: argparse.Namespace) -> None:
    """Members and Pfam tables of the base plus the study (``study-ladder``'s rebuild arm).
    A study protein joins its best MGnify90 cluster (``gene_units`` rank 1) where identity
    and both coverages pass ``--min-id``/``--min-cov``, as a member of that cluster would;
    otherwise it stays in its linclust cluster. Ids become strings, as the study's are."""
    columns = ["protein_id", "cluster_rep", "full_length", "sequence"]
    base = pl.scan_parquet(_member_paths(args.members)).select(columns)
    nearest = (
        pl.read_parquet(args.gene_units)
        .filter(
            pl.col("rank") == 1,
            pl.col("identity") >= args.min_id,
            pl.col("qcov") >= args.min_cov,
            pl.col("scov") >= args.min_cov,
        )
        .select(protein_id="gene_name", nearest=pl.col("cluster_rep").cast(pl.String))
    )
    study = (
        pl.read_parquet(args.study_members)
        .join(nearest, on="protein_id", how="left")
        .with_columns(cluster_rep=pl.coalesce("nearest", "cluster_rep"))
        .select(columns)
    )
    pl.concat(
        [base.with_columns(pl.col("protein_id", "cluster_rep").cast(pl.String)), study.lazy()]
    ).sink_parquet(args.out)
    pfam = pl.read_parquet(args.pfam, columns=["protein_id", "pfam_accession"])
    pl.concat(
        [
            pfam_names(pfam.rename({"protein_id": "unit"})).select(
                protein_id=pl.col("unit").cast(pl.String), pfam_accession="name"
            ),
            pl.read_parquet(args.study_pfam),
        ]
    ).write_parquet(args.out_pfam)


def study_ladder(args: argparse.Namespace) -> None:
    """The study ladder's gate (plan, Additional references: Evaluation). From the scores of
    the plain arm (base), ``study`` (base + study index, joint) and ``rebuild`` per sample,
    index, count, estimate and threshold: each arm's completeness and purity, and
    ``recovered`` = (joint - base) / (rebuild - base) completeness (null without a gain)."""
    keys = ["sample", "index", "count", "abundance", "min_hits"]
    scores = pl.concat(
        [pl.read_csv(p, separator="\t", schema_overrides={"arm": pl.String}) for p in args.scores],
        how="diagonal_relaxed",
    ).with_columns(pl.col("arm").fill_null(""))
    arms = {"": "base", "study": "joint", "rebuild": "rebuild"}
    wide = None
    for arm, tag in arms.items():
        part = scores.filter(pl.col("arm") == arm).select(
            *keys, pl.col("completeness").alias(f"completeness_{tag}"),
            pl.col("purity").alias(f"purity_{tag}"),
        )  # fmt: skip
        wide = part if wide is None else wide.join(part, on=keys, how="inner", nulls_equal=True)
    assert wide is not None
    gain = pl.col("completeness_rebuild") - pl.col("completeness_base")
    wide = wide.with_columns(
        recovered=pl.when(gain > 0).then(
            (pl.col("completeness_joint") - pl.col("completeness_base")) / gain
        )
    )
    wide.sort(keys).write_csv(args.out, separator="\t")


def study_unrelated(args: argparse.Namespace) -> None:
    """From a base and a joint (base + study) unit profile: base units whose components in
    the joint query (tier 2 and dense) hold no study unit (``unrelated``), and how many changed
    ``hits``, ``kmers_unique`` or ``coverage_em`` (``unrelated_changed``; should be 0)."""
    alone, joint = read_profile(args.base), read_profile(args.joint)
    # linked in tier 2 (gather) or, with a dense tier, in the dense hits the EM is fitted on
    links = [c for c in ("component", "component_dense") if c in joint.columns]
    study = joint.filter(pl.col("source") > 0)
    own = joint.filter(
        pl.col("source") == 0,
        *(
            ~pl.col(c).is_in(study[c].drop_nulls().unique().implode()).fill_null(True)
            for c in links
        ),
    )
    compared = ["hits", "kmers_unique", "coverage_em"]
    both = own.select("unit", *compared).join(
        alone.select("unit", *compared), on="unit", how="left", suffix="_alone"
    )
    diff = [
        (pl.col(c) - pl.col(f"{c}_alone")).abs() > 1e-9 * (1 + pl.col(c).abs()) for c in compared
    ]
    changed = int(both.select(pl.any_horizontal(d.fill_null(True) for d in diff).sum()).item())
    pl.DataFrame(
        {"run": [args.run], "unrelated": [own.height], "unrelated_changed": [changed]}
    ).write_csv(args.out, separator="\t")


def abundance_scores(truth: pl.DataFrame, estimate: pl.DataFrame) -> dict[str, float | None]:
    """Spearman over true positives and L1 between relative abundances (0 to 2).

    ``truth`` has ``label`` and ``depth``; ``estimate`` has ``label`` and ``estimate`` for
    the detected KOs. A KO missing from either side has abundance 0 there. If ``estimate``
    also has an interval (``lo``, ``hi``), ``ci_cover`` is the share of true positives whose
    depth, on the estimate's scale (times the median estimate / depth), lies inside it, and
    ``ci_width`` the median log(hi / lo).
    """
    both = truth.select("label", "depth").join(estimate, on="label", how="full", coalesce=True)
    both = both.fill_null(0.0)
    tp = both.filter(pl.col("depth") > 0, pl.col("estimate") > 0)
    rel = [both[c] / both[c].sum() if both[c].sum() > 0 else both[c] for c in ("depth", "estimate")]
    return {
        "spearman_tp": tp.select(pl.corr("depth", "estimate", method="spearman")).item()
        if tp.height > 1
        else None,
        "l1": float((rel[0] - rel[1]).abs().sum()),
        **interval_scores(tp),
    }


def group_scores(truth: pl.DataFrame, detected: pl.DataFrame) -> dict[str, float | None]:
    """How ambiguity groups relate to truth, for profiles that have them.

    ``fp_grouped``: share of false-positive KOs in a group that holds a true KO (the call
    is ambiguous, not wrong). ``group_cover``: share of groups whose ``abundance_zi``
    interval holds the members' true total depth, on the estimate's scale.
    """
    if "ambiguity_group" not in detected.columns:
        return {"fp_grouped": None, "group_cover": None}
    d = detected.join(truth.select("label", "depth"), left_on="name", right_on="label", how="left")
    scale = (
        d.filter(pl.col("depth") > 0, pl.col("abundance_zi") > 0)
        .select((pl.col("abundance_zi") / pl.col("depth")).median())
        .item()
    )
    grouped = d.filter(pl.col("ambiguity_group").is_not_null())
    with_tp = grouped.filter(pl.col("depth").is_not_null())["ambiguity_group"].unique()
    fp = d.filter(pl.col("depth").is_null())
    groups = grouped.group_by("ambiguity_group").agg(
        total=pl.col("depth").fill_null(0).sum() * scale,
        lo=pl.col("group_abundance_zi_lo").first(),
        hi=pl.col("group_abundance_zi_hi").first(),
    )
    return {
        "fp_grouped": fp.select(
            pl.col("ambiguity_group").is_in(with_tp.implode()).fill_null(False).mean()
        ).item()
        if fp.height
        else None,
        "group_cover": groups.select(pl.col("total").is_between("lo", "hi").mean()).item()
        if groups.height
        else None,
    }


def presence_scores(truth: pl.DataFrame, detected: pl.DataFrame) -> dict[str, float | None]:
    """Mean ``present_prob`` of true (``prob_tp``) and false (``prob_fp``) positives, and the
    share of each below 0.5 (``flag_tp``, ``flag_fp``)."""
    if "present_prob" not in detected.columns:
        return dict.fromkeys(("prob_tp", "prob_fp", "flag_tp", "flag_fp"))
    d = detected.select("present_prob", tp=pl.col("name").is_in(truth["label"].implode()))
    out: dict[str, float | None] = {}
    for name, part in (("tp", d.filter("tp")), ("fp", d.filter(~pl.col("tp")))):
        out[f"prob_{name}"] = part["present_prob"].mean() if part.height else None  # type: ignore[assignment]
        out[f"flag_{name}"] = (part["present_prob"] < 0.5).mean() if part.height else None
    return out


def split_scores(
    truth: pl.DataFrame, detected: pl.DataFrame, abundance: str | None
) -> dict[str, float | None]:
    """How the detected units' abundance is split within components (phase 7: k, alphabet).

    ``split_l1``: per component with two or more detected units, the L1 between the
    estimate's and the truth's shares of the component (0 to 2), averaged with the
    component's true depth as weight; components with no true depth are left out.
    ``group_size_mean``: mean ambiguity-group size of the detected units.
    ``host_like_detected``: detected units flagged ``host_like`` (host mask arms).
    """
    out: dict[str, float | None] = dict.fromkeys(
        ("split_l1", "group_size_mean", "host_like_detected")
    )
    if "host_like" in detected.columns:
        out["host_like_detected"] = float(detected["host_like"].cast(pl.Boolean).sum())
    if "group_size" in detected.columns and abundance == "abundance_zi":
        out["group_size_mean"] = detected["group_size"].mean()  # type: ignore[assignment]
    if not abundance or "component" not in detected.columns:
        return out
    units = (
        detected.select("component", "name", estimate=abundance)
        .join(truth.select(name="label", depth="depth"), on="name", how="left")
        .with_columns(pl.col("depth").fill_null(0.0))
        .filter(pl.len().over("component") >= 2)
    )
    share = lambda c: pl.col(c) / pl.col(c).sum().over("component")  # noqa: E731
    per = (
        units.filter((pl.col("depth").sum() > 0).over("component"))
        .filter((pl.col("estimate").sum() > 0).over("component"))
        .group_by("component")
        .agg(l1=(share("estimate") - share("depth")).abs().sum(), weight=pl.col("depth").sum())
    )
    if per.height:
        out["split_l1"] = float((per["l1"] * per["weight"]).sum() / per["weight"].sum())
    return out


def interval_scores(tp: pl.DataFrame) -> dict[str, float | None]:
    if "lo" not in tp.columns or tp.height == 0:
        return {"ci_cover": None, "ci_width": None}
    scaled = pl.col("depth") * (pl.col("estimate") / pl.col("depth")).median()
    return {
        "ci_cover": tp.select(scaled.is_between(pl.col("lo"), pl.col("hi")).mean()).item(),
        "ci_width": tp.select((pl.col("hi") / pl.col("lo")).log().median()).item(),
    }


def score(args: argparse.Namespace) -> None:
    truth = read_truth(args.truth)
    profile = read_profile(args.profile)
    rows = []
    # Profiles from before an estimate existed are scored on detection alone.
    rules = list(dict.fromkeys((c, a if a in profile.columns else None) for c, a in RULES))
    for (count, abundance), min_hits in itertools.product(rules, args.min_hits):
        if count not in profile.columns:
            continue
        detected = profile.filter(pl.col(count) >= min_hits)
        predicted = set(detected["name"])
        present = truth.with_columns(found=pl.col("label").is_in(predicted))
        low = present.filter(pl.col("n_reads") <= pl.col("n_reads").quantile(0.25))
        tp = int(present["found"].sum())
        rows.append(
            {
                "sample": args.sample,
                "label": args.label,
                "index": args.index,
                "arm": args.arm,
                "count": count,
                "abundance": abundance,
                "min_hits": min_hits,
                "n_truth": truth.height,
                "n_pred": len(predicted),
                "tp": tp,
                "purity": tp / max(len(predicted), 1),
                "completeness": tp / max(truth.height, 1),
                "completeness_low25": low["found"].mean(),
                "weighted_completeness": int(present.filter("found")["bases"].sum())
                / max(int(present["bases"].sum()), 1),
                **(
                    abundance_scores(
                        truth,
                        detected.select(
                            label="name",
                            estimate=abundance,
                            **(
                                {"lo": f"{abundance}_lo", "hi": f"{abundance}_hi"}
                                if f"{abundance}_lo" in profile.columns
                                else {}
                            ),
                        ),
                    )
                    if abundance
                    else {"spearman_tp": None, "l1": None, "ci_cover": None, "ci_width": None}
                ),
                **(
                    group_scores(truth, detected)
                    if abundance == "abundance_zi"
                    else {"fp_grouped": None, "group_cover": None}
                ),
                **(
                    presence_scores(truth, detected)
                    if abundance == "abundance_zi"
                    else dict.fromkeys(("prob_tp", "prob_fp", "flag_tp", "flag_fp"))
                ),
                **split_scores(truth, detected, abundance),
            }
        )
    pl.DataFrame(rows).write_csv(args.out, separator="\t")


DETECTED_COLUMNS = ("kmers_unique", "kmers_hit", "hits", "m_g", "n_members", "present_prob",
                    "own_evidence", "ambiguity_group", "abundance_zi")  # fmt: skip


def detected(args: argparse.Namespace) -> None:
    truth = read_truth(args.truth)
    profile = read_profile(args.profile).filter(pl.col("kmers_unique") >= 1)
    kmers = pl.read_parquet(args.kmers).join(profile.select("unit"), on="unit", how="semi")
    meta = json.loads((Path(args.index_dir) / "meta.json").read_text())
    if meta["hash"] == "sourmash" or kmers.height == 0:
        kmers = kmers.with_columns(in_genome=pl.lit(None, pl.Boolean))
    else:
        p = meta["params"]
        wanted = np.unique(kmers["hash"].to_numpy())
        # Contig by contig, keeping only hit hashes: all genome hashes up to the largest hit
        # (t_max 0.2 on floored indexes) took ~10 GB at once.
        found = [np.empty(0, np.uint64)] + [
            h[np.isin(h, wanted)]
            for _, seq in fasta(args.fna)
            for h in [
                _core.hash_dna(
                    [seq.encode()],
                    p["k"],
                    alphabet=p["alphabet"],
                    frames="all",
                    max_hash=int(wanted[-1]),
                )["hash"]
            ]
        ]
        kmers = kmers.with_columns(in_genome=pl.col("hash").is_in(np.unique(np.concatenate(found))))
    per_unit = kmers.group_by("unit").agg(
        holders_median=pl.col("holders").median(),
        hits_max=pl.col("hits").max(),
        in_genome=pl.col("in_genome").mean(),
    )
    (
        profile.join(per_unit, on="unit", how="left")
        .with_columns(
            sample=pl.lit(args.sample),
            index=pl.lit(args.index),
            tp=pl.col("name").is_in(truth["label"].implode()),
        )
        .select(
            "sample",
            "index",
            "name",
            "tp",
            # Every column in every file: detected.tsv stacks them under one header.
            *[c if c in profile.columns else pl.lit(None).alias(c) for c in DETECTED_COLUMNS],
            "holders_median",
            "hits_max",
            "in_genome",
        )
        .write_csv(args.out, separator="\t")
    )


def _read_pair(col: str) -> pl.Expr:
    """Read-pair name: iss mates end in /1 and /2."""
    return pl.col(col).str.replace(r"/[12]$", "")


def tool_profile(args: argparse.Namespace) -> None:
    """One tool's raw output -> ``name``, ``evidence``, ``abundance`` (see each branch).

    Abundances are on each tool's own scale; SCORE uses ranks (Spearman) and relative
    abundances (L1), so only proportionality to depth matters.
    """
    raw = Path(args.raw)
    if args.tool == "diamond":
        # best hit per read (-k 1); a gene with several KOs counts for each. Evidence: read
        # pairs; abundance: aligned / subject length summed over reads, as truth depth is.
        hits = pl.read_csv(raw / "hits.tsv", separator="\t", has_header=False,
                           new_columns=["read", "gene", "length", "slen", "bitscore"])  # fmt: skip
        out = (
            hits.with_columns(
                pair=_read_pair("read"), gene_id=pl.col("gene").str.split("|").list.first()
            )
            .join(read_kos(args.kos), on="gene_id")
            .group_by(name="ko_id")
            .agg(
                evidence=pl.col("pair").n_unique(),
                abundance=(pl.col("length") / pl.col("slen")).sum(),
            )
        )
    elif args.tool == "fmh_funprofiler":
        # evidence: shared hashes (sourmash reports intersect_bp = hashes x scaled);
        # abundance: funcprofiler's own output (normalised f_match_query)
        prefetch = pl.read_csv(raw / "prefetch.csv", infer_schema_length=None)
        out = prefetch.select(
            name="match_name", evidence=(pl.col("intersect_bp") // args.scaled)
        ).join(pl.read_csv(raw / "ko.csv").rename({"ko_id": "name"}), on="name")
    elif args.tool == "kmermaid":
        # one cluster (KO) per read; evidence: read pairs; abundance: reads / mean member
        # length (aa), i.e. proportional to depth
        reads = pl.read_csv(raw / "kmermaid.tsv", separator="\t", quote_char=None)
        lengths = (
            pl.scan_parquet(args.members)
            .group_by(name="cluster_rep")
            .agg(length=pl.col("sequence").str.len_chars().mean())
            .collect()
        )
        out = (
            reads.group_by(name="cluster_rep")
            .agg(evidence=_read_pair("seq_name").n_unique(), n=pl.len())
            .join(lengths, on="name")
            .select("name", "evidence", abundance=pl.col("n") / pl.col("length"))
        )
    elif args.tool in ("humann", "humann4"):
        # unstratified KO rows of the regrouped gene families (RPK in 3.9, adjusted CPM in 4).
        # HUMAnN reports no read counts, so every reported KO has evidence 1 (min_hits > 1
        # scores the same calls).
        table = pl.read_csv(raw / "ko.tsv", separator="\t", quote_char=None)
        # (KO ids as "K00001", KEGG's as "ko:K00001"; UNMAPPED and UNGROUPED do not match)
        out = (
            table.rename(dict(zip(table.columns, ["name", "abundance"], strict=True)))
            .filter(~pl.col("name").str.contains("|", literal=True), pl.col("abundance") > 0)
            .with_columns(
                name="ko:" + pl.col("name").str.extract(r"^(K\d{5})$"), evidence=pl.lit(1)
            )
            .drop_nulls("name")
        )
    else:
        raise ValueError(f"unknown tool {args.tool}")
    out.select("name", "evidence", "abundance").sort("name").write_csv(args.out, separator="\t")


def cost(args: argparse.Namespace) -> None:
    """Per step and index/tool: tasks, mean wall and CPU hours, max peak RSS (GB).

    Needs the raw trace (``trace.raw = true``: ms and bytes). PROFILE tasks are split by
    index (the tag's last word); every other step is its own row.
    """
    trace = pl.read_csv(args.trace, separator="\t", null_values=["-"], infer_schema_length=None)
    step = pl.col("name").str.extract(r"^(\S+)")
    tag = pl.col("name").str.extract(r"\((.*)\)$")
    (
        trace.filter(pl.col("status").is_in(["COMPLETED", "CACHED"]))  # cached: first run's metrics
        .with_columns(
            step=step,
            what=pl.when(step == "PROFILE").then(tag.str.extract(r"(\S+)$")).otherwise(step),
            hours=pl.col("realtime") / 3.6e6,
            cpu_hours=pl.col("realtime") / 3.6e6 * pl.col("%cpu") / 100,
        )
        .group_by("step", "what")
        .agg(
            tasks=pl.len(),
            hours_mean=pl.col("hours").mean(),
            cpu_hours_mean=pl.col("cpu_hours").mean(),
            peak_rss_gb=pl.col("peak_rss").max() / 1e9,
        )
        .sort("step", "what")
        .write_csv(args.out, separator="\t", float_precision=4)
    )


def summary(args: argparse.Namespace) -> None:
    keys = args.keys
    # A metric that is empty in one file (e.g. group_cover when no group formed) reads as
    # String there; every metric is numeric, so cast before stacking.
    scores = pl.concat(
        [
            pl.read_csv(p, separator="\t").with_columns(
                pl.exclude("sample", *keys).cast(pl.Float64, strict=False)
            )
            for p in args.scores
        ],
        how="diagonal_relaxed",
    )
    metrics = [c for c in scores.columns if c not in ("sample", *keys)]
    (
        scores.group_by(keys)
        .agg(
            pl.len().alias("n_samples"),
            pl.col(metrics).mean().name.suffix("_mean"),
            pl.col(metrics).std().name.suffix("_sd"),
        )
        .sort(keys, nulls_last=True)
        .write_csv(args.out, separator="\t")
    )
    scores.sort([*keys, "sample"], nulls_last=True).write_csv(args.scores_out, separator="\t")


def iss(argv: list[str]) -> None:
    """Run ``iss`` with ``argv``, patching two flaws of its perfect model (iss 2.0.1).

    It never sets ``store_mutations``, so ``generate --mode perfect`` fails with an
    AttributeError that iss prints as usage, exiting 0 without output; and its reads are
    125 bp, not the 151 of the NovaSeq model, which would confound the error comparison.
    """
    import iss.app
    from iss.error_models.perfect import PerfectErrorModel

    init = PerfectErrorModel.__init__

    def patched(self: PerfectErrorModel, *args: object, **kwargs: object) -> None:
        init(self, *args, **kwargs)
        self.store_mutations = False
        self.read_length = 151
        for name in ("subst_choices_for", "subst_choices_rev", "ins_for", "ins_rev",
                     "del_for", "del_rev"):  # fmt: skip
            setattr(self, name, getattr(self, name)[:1] * self.read_length)

    PerfectErrorModel.__init__ = patched  # type: ignore[method-assign]
    sys.argv = ["iss", *argv]
    iss.app.main()


def genome_set(args: argparse.Namespace) -> None:
    """Genome mode's reference set (phase 11): every genome of the record, the samples'
    and all others as distractors, as protein FASTA (``--out-dir``/``{genome}.faa``, from
    its mapping table) and ``genomes.tsv`` (``genome``, ``path``) for ``annotate-genomes``."""
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    root = Path(args.genomes_dir)
    names = sorted(p.name for p in root.iterdir() if (p / f"{p.name}.fasta").exists())
    for g in names:
        proteins = (
            pl.read_csv(root / g / f"{g}_mapping.csv", columns=["gene_name", "aa_sequence"],
                        schema_overrides=GENOME_COLUMNS)
            .drop_nulls()
            .unique("gene_name", keep="first", maintain_order=True)
        )  # fmt: skip
        with open(out / f"{g}.faa", "w") as f:
            for name, seq in proteins.iter_rows():
                f.write(f">{name}\n{seq}\n")
    (out / "genomes.tsv").write_text("genome\tpath\n" + "".join(f"{g}\t{g}.faa\n" for g in names))


def genome_truth(args: argparse.Namespace) -> None:
    """Genome mode's truth for one sample, from the per-gene truth (``truth --out-genes``):
    per genome its read ``depth`` (aligned bases over gene bases, all its genes) and
    ``relative_abundance`` (share of Σ depth: cells, not reads), written to ``--out``; and
    per (genome, KO) and (genome, Pfam) the summed depth of its genes carrying it
    (``kind``, ``genome``, ``label``, ``depth``) to ``--out-functions``."""
    genes = pl.read_parquet(args.genes).select(
        "gene_name", genome=pl.col("contig").str.split("|").list.first(),
        length=pl.col("end") - pl.col("start"),
    )  # fmt: skip
    covered = genes.join(pl.read_csv(args.truth_genes), on="gene_name", how="left").with_columns(
        pl.col("bases", "depth").fill_null(0)
    )
    per_genome = covered.group_by("genome").agg(
        depth=pl.col("bases").sum() / pl.col("length").sum()
    )
    per_genome.with_columns(relative_abundance=pl.col("depth") / pl.col("depth").sum()).sort(
        "genome"
    ).write_csv(args.out)
    labels = [read_kos(args.kos).select(gene_name="gene_id", label="ko_id", kind=pl.lit("ko"))]
    if args.domains:
        domains = pl.read_parquet(args.domains).select("gene_name", label="pfam")
        labels.append(domains.unique().with_columns(kind=pl.lit("pfam")))
    (
        covered.join(pl.concat(labels), on="gene_name")
        .group_by("kind", "genome", "label")
        .agg(pl.col("depth").sum())
        .filter(pl.col("depth") > 0)
        .sort("kind", "genome", "label")
        .write_csv(args.out_functions)
    )


def _f1(tp: int, n_pred: int, n_truth: int) -> dict[str, float]:
    purity, completeness = tp / max(n_pred, 1), tp / max(n_truth, 1)
    f1 = 2 * purity * completeness / (purity + completeness) if tp else 0.0
    return {"purity": purity, "completeness": completeness, "f1": f1}


def genome_score(args: argparse.Namespace) -> None:
    """Genome detection and abundance of one fit (``genomes`` output, or a sylph profile
    with ``--tool sylph``) against ``genome-truth``: purity, completeness and F1 over the
    genomes with reads; L1 between relative abundances over their union; Spearman over the
    true positives. With ``--function-taxon``, the function x genome table against the
    (genome, ``--label``) truth: L1 between shares (ours over the total, so unclassified
    costs), F1 over (genome, function) pairs, ``ft_right`` (the share of classified hits
    on true pairs) and ``ft_unclassified``. Genomes of an ambiguity group reported
    together (names joined by commas) match no truth genome."""
    from scipy.stats import spearmanr

    if args.tool == "sylph":
        raw = pl.read_csv(args.genomes, separator="\t")
        pred = raw.select(
            name=pl.col("Genome_file").str.split("/").list.last().str.replace(r"\.fasta$", ""),
            relative_abundance=pl.col("Taxonomic_abundance") / 100,
        )
    else:
        types = {"name": pl.String, "relative_abundance": pl.Float64}  # typed when empty
        pred = pl.read_csv(args.genomes, separator="\t", schema_overrides=types).select(*types)
    truth = pl.read_csv(args.truth).filter(pl.col("depth") > 0)
    both = truth.select("genome", t="relative_abundance").join(
        pred.select(genome="name", p="relative_abundance"), on="genome", how="full",
        coalesce=True,
    ).fill_null(0.0)  # fmt: skip
    tp = both.filter((pl.col("t") > 0) & (pl.col("p") > 0))
    row: dict[str, object] = {
        "sample": args.sample, "label": args.label, "index": args.index,
        "n_truth": truth.height, "n_pred": pred.height, "tp": tp.height,
        **_f1(tp.height, pred.height, truth.height),
        "l1": float((both["t"] - both["p"]).abs().sum()),
        "spearman": float(spearmanr(tp["t"], tp["p"]).statistic) if tp.height > 2 else None,
    }  # fmt: skip
    if args.function_taxon:
        ft = pl.read_csv(args.function_taxon, separator="\t",
                         schema_overrides={"taxon": pl.String, "hits_em": pl.Float64})  # fmt: skip
        total = float(ft.filter(pl.col("rank") == "total")["hits_em"].sum())
        genome = ft.filter(pl.col("rank") == "genome")
        ours = genome.filter(pl.col("taxon") != "unclassified").select(
            genome="taxon", label="function", p=pl.col("hits_em") / max(total, 1e-300)
        )
        want = pl.read_csv(args.truth_functions).filter(pl.col("kind") == args.label)
        want = want.select("genome", "label", t=pl.col("depth") / pl.col("depth").sum())
        pairs = want.join(ours, on=["genome", "label"], how="full", coalesce=True).fill_null(0.0)
        hit = pairs.filter((pl.col("t") > 0) & (pl.col("p") > 0))
        classified = float(ours["p"].sum())
        row |= {f"ft_{k}": v for k, v in _f1(hit.height, ours.height, want.height).items()} | {
            "ft_l1": float((pairs["t"] - pairs["p"]).abs().sum()),
            "ft_right": float(hit["p"].sum()) / classified if classified else None,
            "ft_unclassified": 1 - classified if total else None,
        }
    pl.DataFrame([row]).write_csv(args.out, separator="\t")


def subsample(args: argparse.Namespace) -> None:
    """A rung of the depth ladder (phase 11, kfp-prior): each read pair kept with
    probability ``--fraction``, drawn with ``--seed``, so every genome's depth scales by it."""
    rng = random.Random(args.seed)
    with (
        gzip.open(args.r1, "rt") as i1, gzip.open(args.r2, "rt") as i2,
        gzip.open("sub_R1.fastq.gz", "wt", compresslevel=1) as o1,
        gzip.open("sub_R2.fastq.gz", "wt", compresslevel=1) as o2,
    ):  # fmt: skip
        for rec1, rec2 in zip(_records(i1), _records(i2), strict=True):
            if rng.random() < args.fraction:
                o1.write(rec1)
                o2.write(rec2)


DEPTH_BINS: Final = [0.0, 0.05, 0.1, 0.3, 1.0, 3.0, float("inf")]
PROB_BINS: Final = [i / 10 for i in range(11)]


def _accession(col: str) -> pl.Expr:
    """Pfam accessions as ``PF01007``, from MGnify's numbers or versioned names."""
    return (
        pl.col(col)
        .cast(pl.String)
        .str.replace(r"\.\d+$", "")
        .map_elements(
            lambda a: a if a.startswith("PF") else f"PF{int(a):05d}", return_dtype=pl.String
        )
    )


def _by_depth(truth: pl.DataFrame, predicted: set, level: str) -> list[dict[str, Any]]:
    """Completeness per bin of the depth of the deepest sample genome carrying each true
    item (``item``, ``depth``), observed and updated (``predicted``: {kind: items})."""
    binned = truth.with_columns(
        bin=pl.col("depth").cut(DEPTH_BINS[1:-1], left_closed=True).cast(pl.String)
    )
    rows = []
    for (label,), part in binned.group_by("bin"):
        items = set(part["item"].to_list())
        rows.append({"level": level, "bin": label, "n_truth": len(items)}
                    | {f"completeness_{k}": len(items & v) / len(items)
                       for k, v in predicted.items()})  # fmt: skip
    every = set(truth["item"].to_list())
    rows.append({"level": level, "bin": "all", "n_truth": len(every)}
                | {f"completeness_{k}": len(every & v) / max(len(every), 1)
                   for k, v in predicted.items()}
                | {f"purity_{k}": len(every & v) / max(len(v), 1) for k, v in predicted.items()}
                | {f"n_pred_{k}": len(v) for k, v in predicted.items()})  # fmt: skip
    return rows


def prior_score(args: argparse.Namespace) -> None:
    """kfp-prior on a depth-ladder rung against what the sample's genomes carry: units
    (their proteins' best units in the genome index) and, with ``--pfam-presence``, Pfams
    (their genes' domains). Completeness per bin of carrier depth (the full sample's
    genome depth x ``--fraction``) for the observed (``present_prob`` >= 0.5 with hits) and
    updated (``present_prob_updated`` >= 0.5) calls, and purity over all; and, in
    ``--out-calibration``, the units with zero hits per bin of ``present_prob_updated``:
    their mean prediction and the share truly carried."""
    gi = Path(args.genome_index)
    names = pl.read_csv(gi / "genomes.tsv", separator="\t").select("genome", "name")
    depth = pl.read_csv(args.truth_genomes).select(
        name="genome", depth=pl.col("depth") * args.fraction
    )
    sample = names.join(depth, on="name")
    units = (
        pl.read_parquet(gi / "genome_best.parquet")
        .join(sample.cast({"genome": pl.UInt32}), on="genome")
        .group_by("unit")
        .agg(pl.col("depth").max())
        .select(item=pl.col("unit").cast(pl.Int64), depth="depth")
    )
    presence = pl.read_csv(args.presence, separator="\t")
    called = {
        "observed": set(presence.filter((pl.col("hits") > 0) & (pl.col("present_prob") >= 0.5))
                        ["unit"].to_list()),
        "updated": set(presence.filter(pl.col("present_prob_updated") >= 0.5)["unit"].to_list()),
    }  # fmt: skip
    rows = _by_depth(units, called, "unit")
    if args.pfam_presence:
        genome_of = pl.read_parquet(args.genes).select(
            "gene_name", name=pl.col("contig").str.split("|").list.first()
        )
        pfams = (
            pl.read_parquet(args.domains)
            .join(genome_of, on="gene_name")
            .join(depth, on="name")
            .group_by("pfam")
            .agg(pl.col("depth").max())
            .select(item=_accession("pfam"), depth="depth")
        )
        pp = pl.read_csv(args.pfam_presence, separator="\t").with_columns(
            item=_accession("pfam_accession")
        )
        called = {
            "observed": set(pp.filter(pl.col("present_prob_observed") >= 0.5)["item"].to_list()),
            "updated": set(pp.filter(pl.col("present_prob_updated") >= 0.5)["item"].to_list()),
        }
        rows += _by_depth(pfams, called, "pfam")
    keys = {"sample": args.sample, "fraction": args.fraction, "index": args.index}
    pl.DataFrame([keys | r for r in rows]).write_csv(args.out, separator="\t")
    zero = presence.filter(pl.col("hits") == 0).with_columns(
        carried=pl.col("unit").is_in(units["item"].implode()),
        bin=pl.col("present_prob_updated").cut(PROB_BINS[1:-1], left_closed=True).cast(pl.String),
    )
    (
        zero.group_by("bin")
        .agg(
            n=pl.len(),
            predicted=pl.col("present_prob_updated").mean(),
            carried=pl.col("carried").mean(),
        )  # fmt: skip
        .with_columns(**{k: pl.lit(v) for k, v in keys.items()})
        .select(*keys, "bin", "n", "predicted", "carried")
        .sort("bin")
        .write_csv(args.out_calibration, separator="\t")
    )


def main() -> None:
    if sys.argv[1:2] == ["iss"]:  # iss parses its own arguments
        iss(sys.argv[2:])
        return
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="step", required=True)
    p = sub.add_parser("members")
    p.add_argument("--faa", required=True)
    p.add_argument("--kos", required=True)
    p.add_argument("--out", default="members.parquet")
    p = sub.add_parser("sample")
    p.add_argument("--genomes-dir", required=True)
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--seed", type=int, required=True)
    p = sub.add_parser("truth")
    for name in ("fna", "genes", "kos", "r1", "r2"):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--out", default="truth.csv")
    p.add_argument("--domains", help="domains.parquet (pfam-domains): also write Pfam truth")
    p.add_argument("--out-pfam", default="truth_pfam.csv")
    p.add_argument("--out-genes", default="truth_genes.csv", help="per-gene reads and depth")
    p = sub.add_parser("score")
    for name in ("truth", "profile", "sample", "index"):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--min-hits", type=int, nargs="+", default=[1])
    p.add_argument("--label", default="ko", help="what the truth and profile name: ko, pfam")
    p.add_argument("--arm", default="", help="query arm (frames, read QC, host handling)")
    p.add_argument("--out", default="score.tsv")
    p = sub.add_parser("detected")
    for name in ("truth", "profile", "kmers", "index-dir", "fna", "sample", "index"):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--out", default="detected.tsv")
    p = sub.add_parser("tool-profile")
    tools = ["diamond", "fmh_funprofiler", "kmermaid", "humann", "humann4"]
    p.add_argument("--tool", required=True, choices=tools)
    p.add_argument("--raw", required=True, help="directory with the tool's output files")
    p.add_argument("--kos", help="gene -> KO table (diamond)")
    p.add_argument("--members", help="members.parquet (kmermaid)")
    p.add_argument("--scaled", type=int, default=1000, help="sketch scaled (fmh_funprofiler)")
    p.add_argument("--out", default="profile.tsv")
    p = sub.add_parser("pfam-proteins")
    p.add_argument("--genomes-dir", required=True)
    p.add_argument("--chunks", type=int, default=1, help="writes proteins.{0..chunks-1}.faa")
    p = sub.add_parser("pfam-domains")
    p.add_argument("--domtbl", required=True, nargs="+", help="hmmsearch --domtblout files")
    p.add_argument("--proteins", required=True, nargs="+")
    p.add_argument("--out", default="domains.parquet")
    p.add_argument("--out-members", default="pfam_members.parquet")
    p = sub.add_parser("pfam-profile")
    p.add_argument("--profile", required=True)
    p.add_argument("--unit-pfam", required=True, nargs="+", help="one per index, query order")
    p.add_argument("--out", default="pfam_profile.tsv")
    p = sub.add_parser("study-proteins")
    for name in ("genomes-dir", "genes", "domains"):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--fraction", type=float, default=1.0, help="share of the sample's genomes")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", default="study.faa")
    p.add_argument("--out-pfam", default="study_pfam.parquet")
    p.add_argument("--out-genomes", default="study_genomes.txt")
    p = sub.add_parser("study-rebuild")
    p.add_argument("--members", required=True, nargs="+", help="base members: files or dirs")
    for name in ("pfam", "study-members", "study-pfam", "gene-units"):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--min-id", type=float, default=0.9)
    p.add_argument("--min-cov", type=float, default=0.8)
    p.add_argument("--out", default="members.parquet")
    p.add_argument("--out-pfam", default="pfam.parquet")
    p = sub.add_parser("study-ladder")
    p.add_argument("--scores", required=True, nargs="+")
    p.add_argument("--out", default="study_ladder.tsv")
    p = sub.add_parser("study-unrelated")
    for name in ("base", "joint", "run"):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--out", default="study_unrelated.tsv")
    p = sub.add_parser("reps")
    p.add_argument("--members", required=True, nargs="+", help="parquet files or directories")
    p.add_argument("--out", default="reps.faa")
    p = sub.add_parser("mgnify-genes")
    p.add_argument("--hits", required=True, nargs="+")
    p.add_argument("--out", default="gene_units.parquet")
    p = sub.add_parser("aai-score")
    for name in ("profile", "gene-units", "genes", "sample", "index"):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--arm", default="")
    p.add_argument("--min-id", type=float, default=0.9)
    p.add_argument("--min-cov", type=float, default=0.8)
    p.add_argument("--k", type=int, default=11, help="the index's k")
    p.add_argument("--model", help="re-estimate aai under this aai_model.json")
    p.add_argument("--min-aai-kmers", type=float, help="re-estimate aai with this mask")
    p.add_argument("--out", default="aai_score.tsv")
    p = sub.add_parser("aai-calibrate")
    p.add_argument("--profiles", required=True, nargs="+", help="unit profiles (units.tsv)")
    p.add_argument("--genes", required=True, nargs="+", help="per-gene truth, one per profile")
    p.add_argument("--gene-units", required=True)
    p.add_argument("--index", help="the profiles' index: records its build parameters")
    p.add_argument("--map", help="score this map (another run's --out) instead of fitting one")
    p.add_argument("--min-id", type=float, default=0.9)
    p.add_argument("--min-cov", type=float, default=0.8)
    p.add_argument("--out", default="aai_calibration.json")
    p.add_argument("--scores-out", default="aai_calibration_scores.tsv")
    p = sub.add_parser("mix")
    for name in ("r1", "r2", "host-r1", "host-r2"):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--fraction", type=float, required=True, help="host share of read pairs")
    p.add_argument("--seed", type=int, default=0)
    p = sub.add_parser("host-abundance")
    p.add_argument("--fasta", required=True)
    p.add_argument("--mito-copies", type=float, default=100)
    p.add_argument("--phix", type=float, default=0.005, help="share of host reads")
    p.add_argument("--out", default="abundance.txt")
    p = sub.add_parser("decoy-members")
    p.add_argument("--faa", required=True)
    p.add_argument("--out", default="decoy_members.parquet")
    p = sub.add_parser("cost")
    p.add_argument("trace")
    p.add_argument("--out", default="cost.tsv")
    p = sub.add_parser("summary")
    p.add_argument("scores", nargs="+")
    p.add_argument("--keys", nargs="+", default=["label", "index", "arm", "count", "abundance",
                                                 "min_hits"])  # fmt: skip
    p.add_argument("--out", default="summary.tsv")
    p.add_argument("--scores-out", default="scores.tsv")
    p = sub.add_parser("genome-set")
    p.add_argument("--genomes-dir", required=True)
    p.add_argument("--out-dir", default="genome_set")
    p = sub.add_parser("genome-truth")
    for name in ("genes", "truth-genes", "kos"):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--domains", help="domains.parquet: also (genome, Pfam) truth")
    p.add_argument("--out", default="truth_genomes.csv")
    p.add_argument("--out-functions", default="truth_function_genome.csv")
    p = sub.add_parser("genome-score")
    for name in ("genomes", "truth", "sample", "label", "index"):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--tool", default="kfp", choices=["kfp", "sylph"])
    p.add_argument("--function-taxon")
    p.add_argument("--truth-functions")
    p.add_argument("--out", default="genome_score.tsv")
    p = sub.add_parser("subsample")
    for name in ("r1", "r2"):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--fraction", type=float, required=True)
    p.add_argument("--seed", type=int, default=0)
    p = sub.add_parser("prior-score")
    for name in ("presence", "genome-index", "truth-genomes", "sample", "index"):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--fraction", type=float, default=1.0)
    p.add_argument("--pfam-presence")
    p.add_argument("--genes", help="the sample's genes.parquet (with --pfam-presence)")
    p.add_argument("--domains", help="domains.parquet (with --pfam-presence)")
    p.add_argument("--out", default="prior_score.tsv")
    p.add_argument("--out-calibration", default="prior_calibration.tsv")
    args = parser.parse_args()
    steps = {"members": members, "sample": sample, "truth": truth, "score": score,
             "detected": detected, "summary": summary, "tool-profile": tool_profile,
             "cost": cost, "pfam-proteins": pfam_proteins, "pfam-domains": pfam_domains,
             "pfam-profile": pfam_profile, "mix": mix, "host-abundance": host_abundance,
             "decoy-members": decoy_members, "reps": reps, "mgnify-genes": mgnify_genes,
             "aai-score": aai_score_step, "aai-calibrate": aai_calibrate,
             "study-proteins": study_proteins, "study-rebuild": study_rebuild,
             "study-ladder": study_ladder, "study-unrelated": study_unrelated,
             "genome-set": genome_set, "genome-truth": genome_truth,
             "genome-score": genome_score, "subsample": subsample,
             "prior-score": prior_score}  # fmt: skip
    steps[args.step](args)


if __name__ == "__main__":
    main()
