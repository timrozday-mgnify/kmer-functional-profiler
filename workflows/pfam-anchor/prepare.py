"""The local test case for the Pfam-anchored unit assignment (docs/pfam-anchor-plan.md, phase 0).

Steps run in order; ``all`` runs every one. Everything goes under ``--out``.

- ``annotate``: the ``--families`` HMMs cut from ``--pfam`` -> ``hmms.hmm`` and ``ga.tsv``;
  every ``--members`` protein searched with ``pfam-reannotate``'s thresholds -> ``hits.parquet``
  (``reannotate.py parse``'s schema; held-out proteins included).
- ``reads``: real relatives (``recall.py holdout``) -> ``holdout/``, whose ``members.parquet``
  is the reference; ``domains.parquet``, the reference's hits. Two read sets, each from
  ``spike.py simulate``: ``strains/`` (``--strains`` reference proteins with a GA domain,
  mutated to 100-70% identity) and ``relatives/`` (the held-out proteins unchanged). Each gets
  ``truth.parquet``: strain, source protein, truth unit (the source's unit for strains, its
  nearest remaining relative's for relatives), identity bin.
- ``search``: per read set, six-frame ORFs of both mates (split at stops, ``--min-orf`` aa) ->
  ``orfs.parquet`` (orf_id, read pair, mate, strand, start base on that strand, length), and
  their hits -> ``hits.parquet``; the same for the reads reversed (``null_*``), the null.
- ``diamond``: both mates against the reference, DIAMOND ``blastx`` from the biocontainer,
  best hit per mate -> ``diamond.parquet`` (read, mate, protein_id, pident, length, bitscore).
- ``report``: ``hit`` recall (a read pair has a hit on a family of its source protein, at GA and
  at any score) by read set, identity bin and where the pair lies: ``domain`` (a mate overlaps
  one of the source's GA domains by at least ``--min-orf`` residues), ``outside`` (coding, but
  not that) or ``flank``. Null rows give the share of reversed pairs with any hit. With
  DIAMOND's output, also the share with a DIAMOND hit and its median identity.
  -> ``report.tsv``.
"""

import argparse
import gzip
import itertools
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parents[2]
PY = sys.executable
REANNOTATE = ROOT / "workflows/pfam-reannotate/reannotate.py"
RECALL = ROOT / "workflows/tree-index/recall.py"
SPIKE = ROOT / "workflows/spike-benchmark/spike.py"
# pfam-reannotate's thresholds, and tree-index build.nf's bench_args
HMMSEARCH = ["-E", "10", "--domE", "1000000", "-Z", "1000000", "--domZ", "1000000"]
BENCH = [
    *("--present", "1", "--identities", "1", "0.98", "0.95", "0.9", "0.85", "0.8", "0.75"),
    *("0.7", "--substitution", "blosum62", "--depth-mu", "0.3", "--depth-sigma", "0.6"),
]
DIAMOND = "quay.io/biocontainers/diamond:2.2.8--he361c42_0"
SETS = ("strains", "relatives")

# Six-frame ORFs as kfpfam's measure/real_orfs.py (code 11).
CODON = dict(
    zip(
        ("".join(c) for c in itertools.product("TCAG", repeat=3)),
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        strict=True,
    )
)
RC = str.maketrans("ACGTN", "TGCAN")


MAX_GB = 6.0  # resident memory of a step's processes before the watchdog kills them


def run(*cmd: object, shell: bool = False) -> None:
    """Run a step, killing its process group if its resident memory passes ``MAX_GB``
    (a run with no cap crashed the laptop). ``shell``: one bash command, with pipefail."""
    text = " ".join(map(str, cmd))
    print("+", text, flush=True)
    args = ["bash", "-o", "pipefail", "-c", text] if shell else [str(c) for c in cmd]
    child = subprocess.Popen(args, start_new_session=True)
    while child.poll() is None:
        ps = subprocess.run(["ps", "-A", "-o", "pgid=,rss="], capture_output=True, text=True)
        kb = sum(
            int(r) for g, r in (x.split() for x in ps.stdout.splitlines()) if int(g) == child.pid
        )
        if kb > MAX_GB * 2**20:
            os.killpg(child.pid, signal.SIGKILL)
            sys.exit(f"killed at {kb / 2**20:.1f} GB (MAX_GB {MAX_GB}): {text}")
        time.sleep(0.5)
    if child.returncode:
        sys.exit(f"exit {child.returncode}: {text}")


def orfs_at(seq: str, min_len: int) -> list[tuple[str, str, int]]:
    """(protein, strand, first base on its own strand) per ORF of at least ``min_len`` aa."""
    out = []
    for sign, strand in (("+", seq), ("-", seq[::-1].translate(RC))):
        for frame in range(3):
            prot = "".join(
                CODON.get(strand[i : i + 3], "X") for i in range(frame, len(strand) - 2, 3)
            )
            at = 0
            for p in prot.split("*"):
                if len(p) >= min_len:
                    out.append((p, sign, frame + 3 * at))
                at += len(p) + 1
    return out


def mates(d: Path, reverse: bool = False):
    """(read pair, mate, sequence) in file order; ``reverse``: each read reversed."""
    with gzip.open(d / "reads_1.fq.gz", "rt") as f1, gzip.open(d / "reads_2.fq.gz", "rt") as f2:
        for i, (a, b) in enumerate(
            zip(itertools.islice(f1, 1, None, 4), itertools.islice(f2, 1, None, 4), strict=True)
        ):
            for mate, s in ((1, a.strip()), (2, b.strip())):
                yield i, mate, s[::-1] if reverse else s


def search(a: argparse.Namespace, faa: Path, out: Path, arm: str) -> None:
    """hmmsearch ``faa`` (integer names) against ``hmms.hmm`` -> ``out`` hits Parquet."""
    tmp = out.with_suffix("")
    residues = n = 0
    for line in faa.read_text().splitlines():
        n += line.startswith(">")
        residues += 0 if line.startswith(">") else len(line)
    Path(f"{tmp}.counts.tsv").write_text(f"proteins\tresidues\n{n}\t{residues}\n")
    cmd = (
        f"hmmsearch --notextw --cpu {a.threads} {' '.join(HMMSEARCH)} --domtblout "
        f"'{tmp}.domtbl' '{a.out / 'hmms.hmm'}' '{faa}' | gzip -1 > '{tmp}.txt.gz'"
    )
    run(cmd, shell=True)
    run(
        PY,
        REANNOTATE,
        "parse",
        "--arm",
        arm,
        "--domtbl",
        f"{tmp}.domtbl",
        "--alignments",
        f"{tmp}.txt.gz",
        "--ga",
        a.out / "ga.tsv",
        "--counts",
        f"{tmp}.counts.tsv",
        "--out",
        out,
        "--stats",
        f"{tmp}.stats.tsv",
    )
    for suffix in (".txt.gz", ".domtbl", ".counts.tsv"):
        Path(f"{tmp}{suffix}").unlink()


def annotate(a: argparse.Namespace) -> None:
    a.out.mkdir(parents=True, exist_ok=True)
    wanted = set(pl.read_csv(a.families, separator="\t")["pfam_accession"].to_list())
    with gzip.open(a.pfam, "rt") as f, open(a.out / "hmms.hmm", "w") as out:
        block: list[str] = []
        for line in f:
            block.append(line)
            if line.startswith("//"):
                acc = next(x for x in block if x.startswith("ACC"))
                if int(acc.split()[1][2:].split(".")[0]) in wanted:
                    out.writelines(block)
                block = []
    run(PY, REANNOTATE, "ga", "--hmm", a.out / "hmms.hmm", "--out", a.out / "ga.tsv")
    n = pl.read_csv(a.out / "ga.tsv", separator="\t").height
    if n != len(wanted):
        sys.exit(f"{n} of {len(wanted)} families found in {a.pfam}")
    faa = a.out / "members.faa"
    members = pl.read_parquet(a.members, columns=["protein_id", "sequence"])
    faa.write_text("".join(f">{i}\n{s}\n" for i, s in members.iter_rows()))
    search(a, faa, a.out / "hits.parquet", "members")
    faa.unlink()


def reads(a: argparse.Namespace) -> None:
    held = a.out / "holdout"
    run(
        PY,
        RECALL,
        "holdout",
        "--members",
        a.members,
        "--per-bin",
        a.per_bin,
        "--seed",
        a.seed,
        "--out",
        held,
    )
    ref = pl.read_parquet(
        held / "members.parquet", columns=["protein_id", "cluster_rep", "sequence"]
    )
    hits = pl.read_parquet(a.out / "hits.parquet")
    hits.join(ref.select("protein_id"), on="protein_id").write_parquet(a.out / "domains.parquet")

    pool = ref.join(hits.filter("above_ga").select("protein_id").unique(), on="protein_id")
    chosen = pool.sample(a.strains, seed=a.seed).sort("protein_id")
    (a.out / "strains").mkdir(exist_ok=True)
    targets = a.out / "strains/targets.fa"
    targets.write_text("".join(f">{i}\n{s}\n" for i, _, s in chosen.iter_rows()))
    run(
        PY,
        SPIKE,
        "simulate",
        "--targets",
        targets,
        *BENCH,
        "--seed",
        a.seed,
        "--out",
        a.out / "strains",
    )
    run(
        PY,
        SPIKE,
        "simulate",
        "--targets",
        held / "targets.fa",
        "--present",
        1,
        "--identities",
        1,
        "--depth-mu",
        0.3,
        "--depth-sigma",
        0.6,
        "--seed",
        a.seed,
        "--out",
        a.out / "relatives",
    )

    unit = pl.read_parquet(a.members, columns=["protein_id", "cluster_rep"])
    for name in SETS:
        truth = pl.read_csv(a.out / name / "truth.tsv", separator="\t").select(
            "strain", protein_id=pl.col("member").cast(pl.Int64), identity_bin="identity_bin"
        )
        if name == "relatives":
            rel = pl.read_csv(held / "relatives.tsv", separator="\t")
            truth = truth.drop("identity_bin").join(rel, on="protein_id")
            key = "relative"
        else:
            key = "protein_id"
        truth = truth.join(unit.rename({"protein_id": key, "cluster_rep": "unit"}), on=key)
        truth.select("strain", "protein_id", "unit", "identity_bin").write_parquet(
            a.out / name / "truth.parquet"
        )


def orfs(a: argparse.Namespace) -> None:
    for name, null in itertools.product(SETS, (False, True)):
        d = a.out / name
        prefix = "null_" if null else ""
        rows = []
        with open(d / f"{prefix}orfs.faa", "w") as faa:
            for read, mate, seq in mates(d, null):
                for prot, strand, start in orfs_at(seq, a.min_orf):
                    faa.write(f">{len(rows)}\n{prot}\n")
                    rows.append((len(rows), read, mate, strand, start, len(prot)))
        pl.DataFrame(
            rows, schema=["orf_id", "read", "mate", "strand", "start", "length"], orient="row"
        ).write_parquet(d / f"{prefix}orfs.parquet")
        search(
            a,
            d / f"{prefix}orfs.faa",
            d / f"{prefix}hits.parquet",
            f"{name}{'_null' if null else ''}",
        )
        (d / f"{prefix}orfs.faa").unlink()


def diamond(a: argparse.Namespace) -> None:
    ref = pl.read_parquet(a.out / "holdout/members.parquet", columns=["protein_id", "sequence"])
    (a.out / "reference.faa").write_text("".join(f">{i}\n{s}\n" for i, s in ref.iter_rows()))
    # The container is capped too: Colima's VM may take 20 of the laptop's 24 GB.
    cap = f"{a.max_gb:g}g"
    docker = [
        "docker",
        "run",
        "--rm",
        f"--memory={cap}",
        f"--memory-swap={cap}",
        f"--cpus={a.threads}",
        "-v",
        f"{a.out.resolve()}:/w",
        "-w",
        "/w",
        DIAMOND,
        "diamond",
    ]
    run(*docker, "makedb", "--in", "reference.faa", "--db", "reference", "--threads", a.threads)
    for name in SETS:
        d = a.out / name
        with open(d / "mates.fa", "w") as f:
            f.writelines(f">{r}_{m}\n{s}\n" for r, m, s in mates(d))
        run(
            *docker,
            "blastx",
            "--db",
            "reference",
            "--query",
            f"{name}/mates.fa",
            "--out",
            f"{name}/diamond.tsv",
            "--sensitive",
            "-k",
            1,
            "-b",
            0.5,
            "--threads",
            a.threads,
            "--outfmt",
            6,
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "bitscore",
        )
        hits = pl.read_csv(
            d / "diamond.tsv",
            separator="\t",
            has_header=False,
            new_columns=["q", "protein_id", "pident", "length", "bitscore"],
        )
        hits.unique("q", keep="first", maintain_order=True).with_columns(
            pl.col("q").str.split_exact("_", 1).struct.rename_fields(["read", "mate"])
        ).unnest("q").with_columns(pl.col("read", "mate").cast(pl.Int64)).write_parquet(
            d / "diamond.parquet"
        )
        (d / "mates.fa").unlink()
        (d / "diamond.tsv").unlink()
    (a.out / "reference.faa").unlink()


def pair_hits(d: Path, prefix: str) -> pl.DataFrame:
    """(read, pfam_accession, above_ga): every family a read pair's ORFs hit."""
    o = pl.read_parquet(d / f"{prefix}orfs.parquet", columns=["orf_id", "read"])
    h = pl.read_parquet(
        d / f"{prefix}hits.parquet", columns=["protein_id", "pfam_accession", "above_ga"]
    )
    return (
        h.join(o.rename({"orf_id": "protein_id"}), on="protein_id")
        .group_by("read", "pfam_accession")
        .agg(pl.col("above_ga").any())
    )


def report(a: argparse.Namespace) -> None:
    ga = pl.read_parquet(a.out / "hits.parquet").filter("above_ga")
    rows = []
    for name in SETS:
        d = a.out / name
        truth = pl.read_parquet(d / "truth.parquet")
        pairs = pl.read_parquet(d / "reads.parquet").join(truth, on="strain")
        # Both mates' residue ranges on the source protein (0-based, end exclusive).
        L = pl.col("pair_bases") // 2
        mate_aa = [
            ((pl.col("cds_start") + off) / 3, (pl.col("cds_start") + off + L) / 3)
            for off in (pl.lit(0), pl.col("fragment") - L)
        ]
        dom = pairs.select("read", "protein_id", "cds_start", "fragment", "pair_bases").join(
            ga.select("protein_id", "env_from", "env_to"), on="protein_id"
        )
        overlap = pl.max_horizontal(
            [
                pl.min_horizontal(hi, pl.col("env_to"))
                - pl.max_horizontal(lo, pl.col("env_from") - 1)
                for lo, hi in mate_aa
            ]
        )
        in_domain = dom.filter(overlap >= a.min_orf).select("read").unique()
        source = ga.select("protein_id", "pfam_accession").unique()
        hit = (
            pair_hits(d, "")
            .join(pairs.select("read", "protein_id"), on="read")
            .join(source, on=["protein_id", "pfam_accession"])
            .group_by("read")
            .agg(hit_any=pl.lit(True), hit_ga=pl.col("above_ga").any())
        )
        df = (
            pairs.join(hit, on="read", how="left")
            .join(in_domain.with_columns(dom=pl.lit(True)), on="read", how="left")
            .with_columns(
                pl.col("hit_any", "hit_ga").fill_null(False),
                where=pl.when(pl.col("dom"))
                .then(pl.lit("domain"))
                .when(pl.col("coding") > 0)
                .then(pl.lit("outside"))
                .otherwise(pl.lit("flank")),
            )
        )
        aggs = [pl.len().alias("pairs"), pl.col("hit_ga").mean(), pl.col("hit_any").mean()]
        if (d / "diamond.parquet").exists():
            best = (
                pl.read_parquet(d / "diamond.parquet")
                .sort("bitscore", descending=True)
                .unique("read", keep="first")
                .select("read", "pident")
            )
            df = df.join(best, on="read", how="left")
            aggs += [
                pl.col("pident").is_not_null().mean().alias("diamond_hit"),
                pl.col("pident").median().alias("diamond_pident"),
            ]
        for keys in (["identity_bin", "where"], ["identity_bin"]):
            g = df.group_by(keys).agg(aggs).with_columns(set=pl.lit(name))
            rows.append(g if "where" in keys else g.with_columns(where=pl.lit("all")))
        n = pairs.height
        null = pair_hits(d, "null_").group_by("read").agg(pl.col("above_ga").any())
        rows.append(
            pl.DataFrame(
                {
                    "set": [f"{name}_null"],
                    "where": ["all"],
                    "pairs": [n],
                    "hit_ga": [null["above_ga"].sum() / n],
                    "hit_any": [null.height / n],
                }
            )
        )
    out = (
        pl.concat(rows, how="diagonal_relaxed")
        .select("set", "identity_bin", "where", pl.exclude("set", "identity_bin", "where"))
        .sort("set", "where", "identity_bin", descending=[False, False, True], nulls_last=True)
    )
    out.write_csv(a.out / "report.tsv", separator="\t", float_precision=3)
    with pl.Config(tbl_rows=-1, tbl_cols=-1, float_precision=3):
        print(out)


STEPS = {"annotate": annotate, "reads": reads, "search": orfs, "diamond": diamond, "report": report}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("step", choices=[*STEPS, "all"])
    p.add_argument("--out", type=Path, default=ROOT / "data/pfam-anchor")
    p.add_argument("--members", type=Path, default=ROOT / "data/mini-slice/members.parquet")
    p.add_argument(
        "--families", type=Path, default=ROOT / "results-hpc/tree-index-slice/results/families.tsv"
    )
    p.add_argument("--pfam", type=Path, default=Path.home() / "Downloads/Pfam-A.hmm.gz")
    p.add_argument("--strains", type=int, default=2000)
    p.add_argument("--per-bin", type=int, default=250, help="recall.py holdout's")
    p.add_argument("--min-orf", type=int, default=25)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--max-gb", type=float, default=MAX_GB, help="watchdog's cap per step")
    a = p.parse_args()
    globals()["MAX_GB"] = a.max_gb
    for name, step in STEPS.items():
        if a.step in (name, "all"):
            step(a)


if __name__ == "__main__":
    main()
