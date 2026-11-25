"""The bundle end to end on genomes it has not seen. 16 gut genomes with KOfam and dbCAN labels make
the index; 12,000 windows of
50 residues from the other 5 are the sample, aligned to Pfam at GA as the pipeline does and given in
its schema.

  heldout.py align    proteins and windows against all of Pfam (hmmsearch; about 15 min on 6
      threads)
  heldout.py run      bundles at several collapse levels, profiles, and the labels against the
      genes' own -> report/annot_heldout.tsv
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[1]))
import prepare  # noqa: E402
from prep import OUT  # noqa: E402

D = OUT / "annot"
H = OUT / "annot_heldout"
HMM = Path.home() / "Documents/kfpfam/data/pfam/Pfam-A.hmm"
REP = OUT / "report"
CHUNK = 12_000


def search(name: str, rows: list[tuple[int, str]]) -> pl.DataFrame:
    prepare.RUNTIME, prepare.THREADS, prepare.MAX_GB = "none", 6, 8.0
    H.mkdir(exist_ok=True)
    if not (H / "hmms.hmm").exists():
        (H / "hmms.hmm").symlink_to(HMM)
        prepare.run(prepare.PY, prepare.REANNOTATE, "ga", "--hmm", HMM, "--out", H / "ga.tsv")
    parts = []
    for k in range(0, len(rows), CHUNK):
        out = H / f"{name}{k // CHUNK:02d}.parquet"
        if not out.exists():
            faa = H / f"{name}.faa"
            faa.write_text("".join(f">{i}\n{s}\n" for i, s in rows[k : k + CHUNK]))
            prepare.search(argparse.Namespace(out=H, threads=6), faa, out, "members")
            faa.unlink()
        parts.append(pl.read_parquet(out))
        print(name, k, parts[-1].height, flush=True)
    return pl.concat(parts)


def windows(d: pl.DataFrame, held: list[str]) -> pl.DataFrame:
    rng = np.random.default_rng(5)
    x = d.filter(pl.col("genome").is_in(held) & (pl.col("n") >= 50))
    n = x["n"].to_numpy().astype(int)
    pick = rng.choice(x.height, 12_000, p=(n - 49) / (n - 49).sum())
    start = [int(rng.integers(0, n[i] - 49)) for i in pick]
    s = x["sequence"].to_list()
    pid = x["protein_id"].to_list()
    return pl.DataFrame(
        {
            "w": range(len(pick)),
            "protein_id": [pid[i] for i in pick],
            "sequence": [s[i][a : a + 50] for i, a in zip(pick, start, strict=False)],
        }
    )


def data():
    d = (
        pl.read_parquet(D / "proteins.parquet")
        .with_row_index("protein_id")
        .with_columns(pl.col("protein_id").cast(pl.Int64))
    )
    genomes = d["genome"].unique(maintain_order=True).to_list()
    return d, genomes[::4][:5]


def align() -> None:
    d, held = data()
    print("held out", held)
    w = list(windows(d, held).select("w", "sequence").iter_rows())
    search("windows", w)
    (H / "orfs.faa").write_text(
        "".join(f">{i}\n{s}\n" for i, s in w)
    )  # for kfpfam, in reach.py's run
    search(
        "proteins",
        list(d.filter(~pl.col("genome").is_in(held)).select("protein_id", "sequence").iter_rows()),
    )


def run() -> None:
    d, held = data()
    w = windows(d, held)
    hits = search("windows", list(w.select("w", "sequence").iter_rows())).filter("above_ga")
    s = hits.join(w.select(pl.col("w").alias("protein_id"), "sequence"), on="protein_id")
    s.select(
        read=pl.col("protein_id").cast(pl.String),
        pfam=pl.format("PF{}", pl.col("pfam_accession").cast(pl.String).str.zfill(5)),
        hmm_from="hmm_from",
        cigar="cigar",
        residues=pl.col("sequence").str.slice(
            pl.col("ali_from") - 1, pl.col("ali_to") - pl.col("ali_from") + 1
        ),
    ).write_parquet(H / "sample.parquet")
    ref = d.filter(~pl.col("genome").is_in(held))
    ref.select("protein_id", "sequence").write_parquet(H / "seqs.parquet")
    search("proteins", list(ref.select("protein_id", "sequence").iter_rows())).write_parquet(
        H / "hits.parquet"
    )
    for k in ("ko", "cazy"):
        ref.select("protein_id", pl.col(k).alias("label")).filter(
            pl.col("label").list.len() > 0
        ).explode("label", empty_as_null=True).write_csv(
            H / f"{k}.tsv", separator="\t", include_header=False
        )
    truth = {str(i): p for i, p in w.select("w", "protein_id").iter_rows()}
    lab = {
        k: dict(zip(d["protein_id"].to_list(), d[k].to_list(), strict=False))
        for k in ("ko", "cazy")
    }
    aligned = set(s["protein_id"].cast(pl.String))
    out = []
    for t in (0.0, 0.9, 0.7, 0.5):
        b = H / f"b{t:g}"
        subprocess.run(["rm", "-rf", b], check=True)
        py = [sys.executable, HERE.parent / "bundle.py"]
        kf = Path.home() / "Documents/kfpfam/data"
        extra = (
            ["--kfpfam-index", kf / "index/k4T7.5.tidx", "--kfpfam-prep", kf / "prep"]
            if t == 0.7
            else []
        )
        subprocess.run(
            [
                *py,
                "build",
                "--hits",
                H / "hits.parquet",
                "--seqs",
                H / "seqs.parquet",
                "--hmm",
                HMM,
                "--out",
                b,
                "--collapse",
                str(t),
                "--pfam",
                "38.2",
                "--mgnify",
                "21 gut genomes",
                "--label",
                f"ko={H / 'ko.tsv'}",
                "--label",
                f"cazy={H / 'cazy.tsv'}",
                *extra,
            ],
            check=True,
        )
        subprocess.run(
            [
                *py,
                "profile",
                "--bundle",
                b,
                "--parquet",
                H / "sample.parquet",
                "--out",
                H / f"p{t:g}",
                "--hmm",
                HMM,
            ],
            check=True,
        )
        m = json.loads((b / "manifest.json").read_text())
        r = pl.read_parquet(H / f"p{t:g}.reads.parquet")
        for floor in (0.0, 0.4, 0.6):
            best = (
                r.filter(pl.col("similarity") >= floor)
                .sort("score", "qid", descending=[True, False])
                .unique("read", keep="first", maintain_order=True)
            )
            for k in ("ko", "cazy"):
                has = called = right = has_all = 0
                got = dict(best.select("read", k).iter_rows())
                want, have = {}, {}
                for read, pid in truth.items():
                    own = set(lab[k][pid])
                    has_all += bool(own)
                    if read not in aligned:
                        continue
                    call = set((got.get(read) or "").split("+")) - {""}
                    has += bool(own)
                    called += bool(call)
                    right += bool(own & call)
                    for x in own:
                        want[x] = want.get(x, 0) + 1
                    for x in call:
                        have[x] = have.get(x, 0) + 1
                names = sorted(set(want) | set(have))
                a, c = (
                    np.array([want.get(x, 0) for x in names], float),
                    np.array([have.get(x, 0) for x in names], float),
                )
                out.append(
                    (
                        t,
                        floor,
                        k,
                        m["rows"],
                        m["files"]["rows.bin"]["bytes"],
                        has_all,
                        has,
                        called,
                        right / has,
                        right / called,
                        right / has_all,
                        float(np.corrcoef(a, c)[0, 1]),
                        float(np.abs(a - c).sum() / (a + c).sum()),
                    )
                )
    o = pl.DataFrame(
        out,
        orient="row",
        schema=[
            "collapse",
            "floor",
            "label",
            "rows",
            "bytes",
            "labelled_windows",
            "labelled_aligned",
            "called",
            "recall_aligned",
            "precision",
            "recall_all_windows",
            "pearson",
            "bray_curtis",
        ],
    )
    o.write_csv(REP / "annot_heldout.tsv", separator="\t")
    pl.Config.set_tbl_rows(60)
    pl.Config.set_tbl_cols(16)
    pl.Config.set_tbl_width_chars(200)
    print(f"{len(aligned)} of {w.height} windows aligned at GA")
    print(o)


if __name__ == "__main__":
    {"align": align, "run": run}[sys.argv[1]]()
