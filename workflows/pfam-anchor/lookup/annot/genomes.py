"""Annotation labels of 21 gut genomes (mimicc 20HM+: Pyrodigal, KOfam, dbCAN, Pfam 37.4 at GA), and
what they say about an
annotated MGnify index: how much of a genome carries each label, and how much of that lies in Pfam
domains.

  genomes.py stats   -> report/annot_genomes.json, data for `lookup best` in OUT/annot/
  genomes.py score   after `lookup best OUT/annot`: is the label of the best hit in another genome
      the read's own?
"""

import gzip
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[1]))
from analyze import REP, blosum  # noqa: E402
from prep import CODE, OUT  # noqa: E402

SRC = Path.home() / "Documents/mimicc/mimicc-core-metabolism-report_rendered/20HM+/pipeline-results"
D = OUT / "annot"
WINDOW, N = 50, 6000


def load() -> pl.DataFrame:
    rows = []
    for g in sorted(p for p in SRC.iterdir() if (p / "pyrodigal").is_dir()):
        seqs, name = {}, None
        with gzip.open(next((g / "pyrodigal").glob("*.faa.gz")), "rt") as f:
            for line in f:
                if line[0] == ">":
                    name = line[1:].split()[0]
                    seqs[name] = []
                else:
                    seqs[name].append(line.strip())
        seqs = {k: "".join(v).rstrip("*").replace("*", "X") for k, v in seqs.items()}
        ko, caz, pf = defaultdict(set), {}, defaultdict(list)
        with open(next((g / "kofam/concatenated").glob("*_kofam.txt"))) as f:
            for line in f:
                if line[0] == "*":  # above the KO's own threshold
                    c = line[1:].split()
                    ko[c[0]].add(c[1])
        with open(next((g / "dbcan/concatenated").glob("*_overview.tsv"))) as f:
            overview = f.readlines()[1:]
        for line in overview:
            c = line.rstrip("\n").split("\t")
            if int(c[5]) >= 2:  # dbCAN's own rule: two of three tools
                caz[c[0]] = {
                    x.split("_")[0].split("(")[0] for x in c[6].replace("|", "+").split("+")
                } - {"-"}
        with gzip.open(next((g / "pfam").rglob("*.domtbl.gz")), "rt") as f:
            for line in f:
                if line[0] != "#":
                    c = line.split(maxsplit=22)
                    pf[c[0]].append((c[4].split(".")[0], int(c[17]), int(c[18])))
        for k, s in seqs.items():
            rows.append(
                (
                    g.name,
                    k,
                    s,
                    sorted(ko.get(k, ())),
                    sorted(caz.get(k, ())),
                    [f"{a}:{b}:{e}" for a, b, e in pf.get(k, ())],
                )
            )
    return pl.DataFrame(
        rows,
        schema=["genome", "protein", "sequence", "ko", "cazy", "pfam"],
        orient="row",
        schema_overrides={
            "ko": pl.List(pl.String),
            "cazy": pl.List(pl.String),
            "pfam": pl.List(pl.String),
        },
    )


def mask(length: int, pfam: list[str]) -> np.ndarray:
    m = np.zeros(length, bool)
    for x in pfam:
        _, a, b = x.split(":")
        m[int(a) - 1 : int(b)] = True
    return m


def stats() -> None:
    d = load().with_columns(
        n=pl.col("sequence").str.len_chars(),
        has_ko=pl.col("ko").list.len() > 0,
        has_cazy=pl.col("cazy").list.len() > 0,
        has_pfam=pl.col("pfam").list.len() > 0,
    )
    covered = [int(mask(n, p).sum()) for n, p in d.select("n", "pfam").iter_rows()]
    d = d.with_columns(covered=pl.Series(covered))
    res = d["n"].sum()
    out = {
        "genomes": d["genome"].n_unique(),
        "proteins": d.height,
        "residues": res,
        "mean_length": d["n"].mean(),
    }
    for k in ("pfam", "ko", "cazy"):
        x = d.filter(pl.col(f"has_{k}"))
        out[k] = {
            "proteins": x.height / d.height,
            "residues": x["n"].sum() / res,
            "mean_length": x["n"].mean(),
            "with_pfam": float(x["has_pfam"].mean()),
            "residues_in_pfam": x["covered"].sum() / x["n"].sum(),
        }
    out["p_ko_given_pfam"] = float(d.filter("has_pfam")["has_ko"].mean())
    out["p_ko_given_none"] = float(d.filter(~pl.col("has_pfam"))["has_ko"].mean())
    out["p_cazy_given_pfam"] = float(d.filter("has_pfam")["has_cazy"].mean())
    out["p_cazy_given_none"] = float(d.filter(~pl.col("has_pfam"))["has_cazy"].mean())
    out["kos_per_ko_protein"] = float(d.filter("has_ko")["ko"].list.len().mean())
    out["distinct_ko"] = d["ko"].explode(empty_as_null=True).drop_nulls().n_unique()
    out["distinct_cazy"] = d["cazy"].explode(empty_as_null=True).drop_nulls().n_unique()
    # the family alone: domain-residue-weighted share of a family's commonest KO set
    by = defaultdict(Counter)
    for ko, pf in d.select("ko", "pfam").iter_rows():
        for x in pf:
            a, s, e = x.split(":")
            by[a][tuple(ko)] += int(e) - int(s) + 1
    tot = sum(sum(c.values()) for c in by.values())
    out["family_majority_label"] = sum(c.most_common(1)[0][1] for c in by.values()) / tot
    named = {a: Counter({k: v for k, v in c.items() if k}) for a, c in by.items()}
    t2 = sum(sum(c.values()) for c in named.values())
    out["family_majority_label_ko_proteins_only"] = (
        sum(c.most_common(1)[0][1] for c in named.values() if c) / t2
    )
    out["kos_per_family_mean_weighted"] = sum(sum(c.values()) * len(c) for c in named.values()) / t2
    (REP / "annot_genomes.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))

    D.mkdir(exist_ok=True)
    d.write_parquet(D / "proteins.parquet")
    seqs = [CODE[np.frombuffer(s.encode(), np.uint8)] for s in d["sequence"]]
    (D / "reps.off").write_bytes(
        np.concatenate([[0], np.cumsum([len(s) for s in seqs])]).astype(np.uint64).tobytes()
    )
    (D / "reps.seq").write_bytes(np.concatenate(seqs).tobytes())
    gid = {g: i for i, g in enumerate(d["genome"].unique(maintain_order=True))}
    group = np.array([gid[g] for g in d["genome"]], np.uint32)
    (D / "reps.group").write_bytes(group.tobytes())
    rng = np.random.default_rng(3)
    w = np.clip(d["n"].to_numpy().astype(int) - WINDOW + 1, 0, None).astype(float)  # n is unsigned
    pick = rng.choice(d.height, N, p=w / w.sum())
    q, meta = [], []
    for i in pick:
        s = int(rng.integers(0, len(seqs[i]) - WINDOW + 1))
        q.append(seqs[i][s : s + WINDOW])
        meta.append(
            (
                int(i),
                s,
                int(mask(d["n"][int(i)], d["pfam"][int(i)].to_list())[s : s + WINDOW].sum()),
            )
        )
    (D / "q.off").write_bytes((np.arange(N + 1) * WINDOW).astype(np.uint64).tobytes())
    (D / "q.seq").write_bytes(np.concatenate(q).tobytes())
    (D / "q.group").write_bytes(group[pick].tobytes())
    pl.DataFrame(meta, schema=["protein", "start", "in_pfam"], orient="row").write_parquet(
        D / "windows.parquet"
    )


def score() -> None:
    d = pl.read_parquet(D / "proteins.parquet")
    w = pl.read_parquet(D / "windows.parquet")
    b = pl.read_csv(D / "best.tsv", separator="\t", schema_overrides={"targets": pl.String})
    self_ = np.diag(blosum())
    seq = np.frombuffer((D / "q.seq").read_bytes(), np.uint8).reshape(-1, WINDOW)
    ko, caz = d["ko"].to_list(), d["cazy"].to_list()
    rows = []
    for (p, _, inp), (qid, sc, _ties, targets) in zip(w.iter_rows(), b.iter_rows(), strict=False):
        t = [int(x) for x in targets.split(",")] if targets else []
        frac = sc / max(1, sum(int(self_[c]) for c in seq[qid] if c < 20))
        for name, lab in (("ko", ko), ("cazy", caz)):
            own = set(lab[p])
            got = Counter(tuple(lab[x]) for x in t)
            call = set(got.most_common(1)[0][0]) if got else set()  # the tied hits' commonest label
            rows.append(
                (qid, name, frac, inp >= 25, bool(own), bool(call), bool(own & call), len(t))
            )
    o = pl.DataFrame(
        rows,
        schema=["qid", "label", "frac", "in_pfam", "has", "called", "right", "ties"],
        orient="row",
    )
    o = o.with_columns(
        bin=pl.col("frac")
        .cut(
            [0.4, 0.6, 0.8, 0.95],
            labels=["<.4", ".4-.6", ".6-.8", ".8-.95", "≥.95"],
            left_closed=True,
        )
        .cast(str)
    )
    o.write_csv(REP / "annot_transfer.tsv", separator="\t")
    pl.Config.set_tbl_rows(60)
    agg = [
        pl.len().alias("reads"),
        pl.col("has").sum().alias("labelled"),
        (pl.col("right").sum() / pl.col("has").sum()).alias("recall"),
        (pl.col("right").sum() / pl.col("called").sum()).alias("precision"),
    ]
    print(o.group_by("label", "bin").agg(agg).sort("label", "bin"))
    print(o.group_by("label", "in_pfam").agg(agg).sort("label", "in_pfam"))


if __name__ == "__main__":
    {"stats": stats, "score": score}[sys.argv[1]]()
