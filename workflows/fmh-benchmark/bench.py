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
- ``pfam-profile``: a profile against an index with Pfam labels (``unit_pfam.parquet``, e.g.
  MGnify90 clusters) summed per Pfam, each unit counting for each of its Pfams.
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
- ``iss``: ``iss`` with its arguments, the perfect error model patched (see :func:`iss`).
"""

import argparse
import itertools
import json
import random
import sys
from collections.abc import Iterator
from multiprocessing.pool import ThreadPool
from pathlib import Path

import mappy
import numpy as np
import polars as pl

from kmer_functional_profiler import _core

GENOME_COLUMNS = {"gene_name": pl.String, "contig_id": pl.String, "start_position": pl.Int64,
                  "end_position": pl.Int64, "strand": pl.String}  # fmt: skip


def fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """(header, sequence) records; headers keep everything after ``>``."""
    header, chunks = None, []
    with open(path) as f:
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
    if args.domains:
        features = domain_features(genes, pl.read_parquet(args.domains))
        _truth(aligned, features).write_csv(args.out_pfam)


def read_profile(path: str | Path) -> pl.DataFrame:
    """A profile TSV; an empty one (nothing detected) gets numeric columns, not strings."""
    profile = pl.read_csv(path, separator="\t")
    if profile.height:
        return profile
    strings = {"name", "cluster_rep"}
    return profile.cast({c: pl.UInt32 if c == "unit" else pl.Float64
                         for c in profile.columns if c not in strings})  # fmt: skip


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


def pfam_profile(args: argparse.Namespace) -> None:
    """Sum a unit profile per Pfam: counts and point estimates add over the units carrying
    a Pfam; ``present_prob`` is the largest; intervals and per-unit columns are dropped."""
    profile = read_profile(args.profile)
    labels = pl.read_parquet(args.unit_pfam)
    accession = pl.col("pfam_accession")
    labels = labels.select(
        "unit",
        # MGnify stores the accession's number (1007 for PF01007); hmmsearch gives PF01007.23
        name=("PF" + accession.cast(pl.String).str.zfill(5))
        if labels.schema["pfam_accession"].is_integer()
        else accession.str.replace(r"\.\d+$", ""),
    )
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
    keys = ["label", "index", "arm", "count", "abundance", "min_hits"]
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
    scores.sort([*keys, "sample"], nulls_last=True).write_csv("scores.tsv", separator="\t")


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
    p.add_argument("--unit-pfam", required=True)
    p.add_argument("--out", default="pfam_profile.tsv")
    p = sub.add_parser("cost")
    p.add_argument("trace")
    p.add_argument("--out", default="cost.tsv")
    p = sub.add_parser("summary")
    p.add_argument("scores", nargs="+")
    p.add_argument("--out", default="summary.tsv")
    args = parser.parse_args()
    steps = {"members": members, "sample": sample, "truth": truth, "score": score,
             "detected": detected, "summary": summary, "tool-profile": tool_profile,
             "cost": cost, "pfam-proteins": pfam_proteins, "pfam-domains": pfam_domains,
             "pfam-profile": pfam_profile}  # fmt: skip
    steps[args.step](args)


if __name__ == "__main__":
    main()
