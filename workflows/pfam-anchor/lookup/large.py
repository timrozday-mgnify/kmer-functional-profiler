"""The lookup on one very large family at full release size
(hpc/kfp-ablations/pfam-anchor-large/run.sh).

  large.py ids      --release DIR|URL --family 5 --out D      proteins the release annotates with it
      -> ids.u64
  large.py hmm      --pfam Pfam-A.hmm[.gz] --family 5 --out D -> hmms.hmm, ga.tsv (models.txt: its
      version)
  large.py align    --chunk I --out D                         chunkI.faa -> chunkI.parquet
      (hmmsearch)
  large.py queries  --parquets GLOB --family 5 --per-library N --out D
  large.py pack     --out D                                   -> D/bin
  large.py report   --out D                                   -> report.tsv, null.tsv, times.tsv

Reads are re-aligned to the model the targets are aligned to, whatever release the Parquets come
from. Each also gives a
null: the same residues reversed, on the same columns, which shows the score chance reaches at each
family size.
"""

import argparse
import glob
import gzip
import sys
from pathlib import Path

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import prepare  # noqa: E402
from prep import CODE, hmm_lengths_of, pack_queries, rows_of  # noqa: E402


def ids(a: argparse.Namespace) -> None:
    import duckdb

    con = duckdb.connect()
    con.execute(f"SET threads = {a.threads}; SET memory_limit = '{a.memory * 0.75:.0f}GB'")
    con.execute(f"SET temp_directory = '{a.out / 'duckdb_tmp'}'")
    x = con.execute(
        f"SELECT DISTINCT protein_id FROM read_parquet('{a.release}/mgy_proteins_pfam.parquet') "
        f"WHERE pfam_accession = {a.family}"
    ).fetchnumpy()["protein_id"]
    if not len(x):
        sys.exit(f"no protein with family {a.family} in {a.release}")
    (a.out / "ids.u64").write_bytes(np.asarray(x, "<u8").tobytes())
    print(f"{len(x)} proteins with PF{a.family:05d}")


def hmm(a: argparse.Namespace) -> None:
    opener = gzip.open if a.pfam.suffix == ".gz" else open
    block: list[str] = []
    with opener(a.pfam, "rt") as f:
        for line in f:
            block.append(line)
            if line.startswith("//"):
                acc = next((x.split()[1] for x in block if x.startswith("ACC")), "")
                if acc.split(".")[0] == f"PF{a.family:05d}":
                    break
                block = []
        else:
            sys.exit(f"PF{a.family:05d} is not in {a.pfam}")
    (a.out / "hmms.hmm").write_text("".join(block))
    (a.out / "models.txt").write_text(f"{acc}\t{a.pfam}\n")
    prepare.run(
        prepare.PY, prepare.REANNOTATE, "ga", "--hmm", a.out / "hmms.hmm", "--out", a.out / "ga.tsv"
    )
    print(acc)


def search(a: argparse.Namespace, faa: Path, out: Path) -> pl.DataFrame:
    prepare.RUNTIME, prepare.THREADS, prepare.MAX_GB = a.runtime, a.threads, float(a.memory)
    if not out.exists():
        prepare.search(argparse.Namespace(out=a.out, threads=a.threads), faa, out, "members")
    return pl.read_parquet(out)


def align(a: argparse.Namespace) -> None:
    d = search(a, a.out / f"chunk{a.chunk:03d}.faa", a.out / f"chunk{a.chunk:03d}.parquet")
    print(f"chunk {a.chunk}: {d.height} domains, {d['above_ga'].sum()} at GA")


def queries(a: argparse.Namespace) -> None:
    files = sorted(glob.glob(a.parquets))
    if not files:
        sys.exit(f"nothing matches {a.parquets}")
    parts = []
    for p in files:
        run = Path(p).name.split("_")[0]
        d = pl.read_parquet(p, columns=["read", "frame", "domain", "pfam", "residues"]).filter(
            pl.col("pfam").str.starts_with(f"PF{a.family:05d}")
        )
        parts.append(
            d.sample(min(a.per_library, d.height), seed=1).select(
                run=pl.lit(run),
                query=pl.col("read")
                + ":"
                + pl.col("frame").cast(str)
                + ":"
                + pl.col("domain").cast(str),
                sequence="residues",
            )
        )
    q = (
        pl.concat(parts)
        .with_row_index("protein_id")
        .with_columns(pl.col("protein_id").cast(pl.Int64))
    )
    if not q.height:
        sys.exit(f"no read of PF{a.family:05d} in {len(files)} files")
    faa = a.out / "reads.faa"
    faa.write_text(
        "".join(f">{i}\n{s}\n" for i, s in q.select("protein_id", "sequence").iter_rows())
    )
    d = search(a, faa, a.out / "reads.parquet").filter("above_ga")
    best = d.sort("score", descending=True).unique("protein_id", keep="first")
    q = q.join(best, on="protein_id").with_columns(set=pl.lit("real"))
    null = q.with_columns(sequence=pl.col("sequence").str.reverse(), set=pl.lit("null"))
    q = pl.concat([q, null]).with_row_index("qid")
    q.write_parquet(a.out / "queries.parquet")
    print(
        f"{q.height // 2} reads at GA of {faa.read_text().count('>')}, from {q['run'].n_unique()} "
        "libraries, and as many reversed"
    )


def fasta(p: Path):
    name, seq = None, []
    with open(p) as f:
        for line in f:
            if line[0] == ">":
                if name is not None:
                    yield name, "".join(seq)
                name, seq = int(line[1:]), []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def pack(a: argparse.Namespace) -> None:
    """D/bin as prep.py's, for one family: chunk by chunk, so only a chunk's sequences are ever
    held."""
    b = a.out / "bin"
    b.mkdir(exist_ok=True)
    lengths = hmm_lengths_of(a.out / "hmms.hmm")
    ((acc, _),) = lengths.items()
    n_unit = n_row = at = 0
    stored = residues = 0
    with (
        open(b / "reps.seq", "wb") as fs,
        open(b / "reps.rows", "wb") as f1,
        open(b / "reps.pos", "wb") as f2,
        open(b / "reps.unit", "wb") as fu,
        open(b / "reps.ga", "wb") as fg,
        open(b / "reps.off", "wb") as f_off,
    ):
        f_off.write(np.zeros(1, "<u8").tobytes())
        for faa in sorted(a.out.glob("chunk*.faa")):
            t = pl.DataFrame(
                list(fasta(faa)),
                schema=["protein_id", "sequence"],
                orient="row",
                schema_overrides={"protein_id": pl.Int64},
            )
            if not t.height:
                continue
            hits = faa.with_suffix(".parquet")
            if not hits.exists():
                sys.exit(f"{hits} is missing: its align task has not finished")
            off = at + np.cumsum(t["sequence"].str.len_bytes().to_numpy().astype("<u8"))
            f_off.write(off.astype("<u8").tobytes())
            at = int(off[-1])
            for s in t["sequence"]:
                fs.write(CODE[np.frombuffer(s.encode(), np.uint8)].tobytes())
            d = (
                t.with_row_index("u")
                .join(pl.read_parquet(hits), on="protein_id")
                .sort("u", "domain")
            )
            rows_of(d, lengths, f1, f2)
            fu.write((d["u"].to_numpy().astype("<u4") + n_unit).astype("<u4").tobytes())
            fg.write(d["above_ga"].cast(pl.UInt8).to_numpy().tobytes())
            g = d.filter("above_ga").with_columns(n=pl.col("sequence").str.len_chars())
            span = g["ali_to"] - g["ali_from"] + 1
            residues += int(span.sum())
            stored += int(
                sum(
                    min(e + 10, n) - max(s - 10, 1) + 1
                    for s, e, n in g.select("ali_from", "ali_to", "n").iter_rows()
                )
            )
            n_unit += t.height
            n_row += d.height
            print(faa.name, t.height, d.height, flush=True)
    (b / "families.tsv").write_text(f"{acc}\t{lengths[acc]}\t{n_row}\t0\n")
    pack_queries(pl.read_parquet(a.out / "queries.parquet").sort("qid"), b, [acc], lengths)
    ga = int(np.fromfile(b / "reps.ga", "u1").sum())
    size = {
        "clusters": n_unit,
        "rows": n_row,
        "rows_ga": ga,
        "protein_residues": at,
        "domain_residues_ga": residues,
        "stored_residues_ga_flank10": stored,
        "GB_5bit": (stored * 5 / 8 + 10 * ga) / 1e9,
        "GB_1byte": (stored + 10 * ga) / 1e9,
    }
    pl.DataFrame([size]).write_csv(a.out / "size.tsv", separator="\t")
    print(size)


def report(a: argparse.Namespace) -> None:
    q = pl.read_parquet(a.out / "queries.parquet").select(
        pl.col("qid").cast(pl.Int64), "set", "run", length=pl.col("sequence").str.len_chars()
    )
    t = pl.read_csv(a.out / "bin/truthsub.tsv", separator="\t").join(q, on="qid")
    s = pl.read_csv(a.out / "scan.tsv", separator="\t", schema_overrides={"units": pl.String}).join(
        t, on=["sub", "qid"]
    )
    bad = s.filter(pl.col("score") > pl.col("best")).height
    if bad:
        sys.exit(f"{bad} calls score above the exhaustive best: truth and scan disagree")
    s = s.with_columns(
        right=pl.col("score") == pl.col("best"), near=pl.col("score") >= 0.95 * pl.col("best")
    )
    by = [
        pl.len().alias("reads"),
        pl.col("right").mean(),
        pl.col("near").mean(),
        (~pl.col("right")).sum().alias("missed"),
    ]
    r = (
        s.filter(pl.col("set") == "real")
        .group_by("sub", "k", "second")
        .agg(by)
        .sort("sub", "k", "second")
    )
    lib = (
        s.filter((pl.col("set") == "real") & (pl.col("sub") == 0))
        .group_by("k", "second", "run")
        .agg(by)
        .sort("k", "second", "run")
    )
    r.write_csv(a.out / "report.tsv", separator="\t")
    lib.write_csv(a.out / "report_by_library.tsv", separator="\t")
    n = (
        t.group_by("sub", "set")
        .agg(
            pl.len(),
            *(pl.col("best").quantile(x).alias(f"best_{x}") for x in (0.5, 0.95, 0.99, 1.0)),
            length=pl.col("length").mean(),
            ties=pl.col("ties").median(),
        )
        .sort("set", "sub")
    )
    n.write_csv(a.out / "null.tsv", separator="\t")
    pl.Config.set_tbl_rows(60)
    pl.Config.set_tbl_cols(12)
    print(pl.read_csv(a.out / "size.tsv", separator="\t"), r, n)
    for p in sorted(a.out.glob("scan*.times.tsv")):
        x = pl.read_csv(p, separator="\t")
        print(
            p.name,
            x.with_columns(
                rows=pl.col("rows_read") / pl.col("reads"),
                scan_ns=pl.col("scan_s") / pl.col("rows_read") * 1e9,
                top_ns=pl.col("top_s") / pl.col("rows_read") * 1e9,
                ms_a_read=(
                    pl.col("scan_s") + pl.col("top_s") + pl.col("columns_s") + pl.col("gapped_s")
                )
                / pl.col("reads")
                * 1e3,
            ).select("sub", "k", "second", "threads", "rows", "scan_ns", "top_ns", "ms_a_read"),
        )


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("step", choices=["ids", "hmm", "align", "queries", "pack", "report"])
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--family", type=int)
    p.add_argument("--release")
    p.add_argument("--pfam", type=Path)
    p.add_argument("--chunk", type=int)
    p.add_argument("--parquets")
    p.add_argument("--per-library", type=int, default=200)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--memory", type=int, default=16, help="GB")
    p.add_argument("--runtime", choices=["docker", "singularity", "none"], default="none")
    a = p.parse_args()
    a.out = a.out.resolve()
    a.out.mkdir(parents=True, exist_ok=True)
    globals()[a.step](a)


if __name__ == "__main__":
    main()
