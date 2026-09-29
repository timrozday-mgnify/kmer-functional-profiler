"""Steps of the fmh-funprofiler benchmark (``main.nf``) on its Zenodo 10055954 inputs.

- ``members``: KEGG proteins grouped by KO -> members table for ``kmer-functional-profiler
  index`` (a gene with several KOs is a member of each).
- ``sample``: draw genomes, write them as one FASTA with ``genome|contig`` names and their
  gene coordinates as ``genes.parquet``.
- ``truth``: map simulated read pairs back to the sample with minimap2 (mappy); a KO is
  present if a primary alignment overlaps one of its genes, as in the paper's CAMISIM
  ground truth. Writes ``ko_id``, ``n_reads``, ``bases`` (overlapping aligned bases) and
  ``depth`` (bases / gene length, summed over the KO's genes: read depth times copies).
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
- ``summary``: mean and sd of the scores per index, count and threshold.
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
                  "end_position": pl.Int64}  # fmt: skip


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


def truth(args: argparse.Namespace) -> None:
    # mappy releases the GIL while aligning, so threads share one index instead of one per process
    aligner = mappy.Aligner(args.fna, preset="sr")
    with ThreadPool(args.threads) as pool:
        parts = pool.imap(lambda chunk: _map(aligner, chunk), _pairs(args.r1, args.r2))
        hits = [h for part in parts for h in part]
    schema = {"read": pl.Int64, "contig": pl.String, "r_start": pl.Int64, "r_end": pl.Int64}
    aligned = pl.DataFrame(hits, schema=schema, orient="row")
    # join_where on contig + ranges is a contig hash join then a filter (reads x genes per
    # contig, billions of rows on complete genomes); joining on (contig, bin) keeps it local
    overlaps = (
        _binned(aligned, "r_start", "r_end")
        .join(_binned(pl.read_parquet(args.genes), "start", "end"), on=["contig", "bin"])
        .filter(pl.col("r_start") < pl.col("end"), pl.col("r_end") > pl.col("start"))
        .drop("bin")
        .unique()
    )
    (
        overlaps.with_columns(
            bases=pl.min_horizontal("r_end", "end") - pl.max_horizontal("r_start", "start")
        )
        .join(read_kos(args.kos), left_on="gene_name", right_on="gene_id")
        .group_by("ko_id")
        .agg(
            n_reads=pl.col("read").n_unique(),
            bases=pl.col("bases").sum(),
            depth=(pl.col("bases") / (pl.col("end") - pl.col("start"))).sum(),
        )
        .sort("ko_id")
        .write_csv(args.out)
    )


def abundance_scores(truth: pl.DataFrame, estimate: pl.DataFrame) -> dict[str, float | None]:
    """Spearman over true positives and L1 between relative abundances (0 to 2).

    ``truth`` has ``ko_id`` and ``depth``; ``estimate`` has ``ko_id`` and ``estimate`` for
    the detected KOs. A KO missing from either side has abundance 0 there. If ``estimate``
    also has an interval (``lo``, ``hi``), ``ci_cover`` is the share of true positives whose
    depth, on the estimate's scale (times the median estimate / depth), lies inside it, and
    ``ci_width`` the median log(hi / lo).
    """
    both = truth.select("ko_id", "depth").join(estimate, on="ko_id", how="full", coalesce=True)
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
    d = detected.join(truth.select("ko_id", "depth"), left_on="name", right_on="ko_id", how="left")
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
    d = detected.select("present_prob", tp=pl.col("name").is_in(truth["ko_id"].implode()))
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
    truth = pl.read_csv(args.truth)
    profile = pl.read_csv(args.profile, separator="\t")
    rows = []
    # Profiles from before an estimate existed are scored on detection alone.
    rules = list(dict.fromkeys((c, a if a in profile.columns else None) for c, a in RULES))
    for (count, abundance), min_hits in itertools.product(rules, args.min_hits):
        if count not in profile.columns:
            continue
        detected = profile.filter(pl.col(count) >= min_hits)
        predicted = set(detected["name"])
        present = truth.with_columns(found=pl.col("ko_id").is_in(predicted))
        low = present.filter(pl.col("n_reads") <= pl.col("n_reads").quantile(0.25))
        tp = int(present["found"].sum())
        rows.append(
            {
                "sample": args.sample,
                "index": args.index,
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
                            ko_id="name",
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
    truth = pl.read_csv(args.truth)
    profile = pl.read_csv(args.profile, separator="\t").filter(pl.col("kmers_unique") >= 1)
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
            tp=pl.col("name").is_in(truth["ko_id"].implode()),
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


def summary(args: argparse.Namespace) -> None:
    keys = ["index", "count", "abundance", "min_hits"]
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
    p = sub.add_parser("score")
    for name in ("truth", "profile", "sample", "index"):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--min-hits", type=int, nargs="+", default=[1])
    p.add_argument("--out", default="score.tsv")
    p = sub.add_parser("detected")
    for name in ("truth", "profile", "kmers", "index-dir", "fna", "sample", "index"):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--out", default="detected.tsv")
    p = sub.add_parser("summary")
    p.add_argument("scores", nargs="+")
    p.add_argument("--out", default="summary.tsv")
    args = parser.parse_args()
    steps = {"members": members, "sample": sample, "truth": truth, "score": score,
             "detected": detected, "summary": summary}  # fmt: skip
    steps[args.step](args)


if __name__ == "__main__":
    main()
