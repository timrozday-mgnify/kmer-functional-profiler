"""Steps of the aai-model workflow (``main.nf``): the parameters of the k-mer survival model
that ``aai`` inverts (``survival.SurvivalModel``: gamma rate ``shape`` and mean ``region``
length of a Markov chain of rates along the sequence; plan, phase 7, step 25), fitted on
MGnify protein pairs.

A pair is a protein P and an MGnify90 cluster C it aligns to. Its survival is P's k-mers
found in the union of C's members' k-mers over ``pin_sum``, C's average member's k-mers:
what the query measures when P's gene is in a sample. Its identity is DIAMOND's, P against
C's representative. Survival against identity fixes how much survival exceeds a^k; it
cannot by itself tell a strong rate contrast over short regions from a weaker one over long
regions. Co-survival does: for lags j in ``LAGS``, the share of P's window pairs j apart
whose k-mers are both in C's union, against the model's B_j.

- ``queries``: ``--n`` random members (seeded) whose sequences hold only the 20 standard
  residues (so a k-mer's emission order is its window position), as FASTA in ``--chunks``
  files, and every cluster's representative as ``reps.faa``.
- ``pairs``: DIAMOND hits (outfmt 6 ``qseqid sseqid pident length qlen slen bitscore``) ->
  one row per (P, C), its best HSP, with both coverages >= ``--min-cov``, P not C's
  representative, at most ``--per-bin`` pairs per identity bin of ``--bin`` from
  ``--min-id``: identities spread evenly, so no band dominates the fit.
- ``survival``: per pair, P's k-mers in C's union and ``pin_sum``, and co-survival by lag
  (``both_<j>``, ``window_survival``: P's windows in C's union), every k-mer (no hash
  threshold: sampling under it is unbiased, all k-mers are more precise). When P is one of
  C's members it is left out of C (a gene new to the index). Clusters are processed in
  batches of about ``--batch-residues`` member residues.
- ``fit``: ``aai_model.json`` (``survival``, ``shape``, ``region``, ``categories``, ``k``,
  ``alphabet``, ``fit``), what ``kmer-functional-profiler aai-model`` attaches to an index,
  and ``model_strata.tsv``, the model refitted per identity band and per cluster size, to
  check that one parameter pair holds.
"""

import argparse
import json
import random
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
from scipy.optimize import minimize

from kmer_functional_profiler import _core
from kmer_functional_profiler.survival import CATEGORIES, SurvivalModel

HIT_SCHEMA = {
    "protein_id": pl.Int64,
    "cluster_rep": pl.Int64,
    "pident": pl.Float64,
    "length": pl.Int64,
    "qlen": pl.Int64,
    "slen": pl.Int64,
    "bitscore": pl.Float64,
}
MEMBER_COLUMNS = ["protein_id", "cluster_rep", "full_length", "sequence"]
IDENTITY_BANDS = ((0.6, 0.7), (0.7, 0.8), (0.8, 0.9), (0.9, 1.01))
SIZE_BANDS = ((2, 4), (4, 11), (11, 101), (101, 2**62))  # members, P left out
LAGS = (1, 3, 6, 11, 20, 40, 80, 150)  # co-survival lags, windows
STANDARD = "^[ACDEFGHIKLMNPQRSTVWY]+$"


def _members(paths: list[str]) -> pl.LazyFrame:
    paths = [str(Path(p) / "*.parquet") if Path(p).is_dir() else p for p in paths]
    return pl.scan_parquet(paths).select(MEMBER_COLUMNS)


def _fasta(df: pl.DataFrame, name: str, path: str | Path) -> None:
    df.select(pl.format(">{}\n{}", name, "sequence")).write_csv(
        path, include_header=False, quote_style="never"
    )


def queries(args: argparse.Namespace) -> None:
    members = _members(args.members)
    reps = members.filter(pl.col("protein_id") == pl.col("cluster_rep"))
    reps.select(pl.format(">{}\n{}", "cluster_rep", "sequence")).sink_csv(
        "reps.faa", include_header=False, quote_style="never"
    )
    standard = members.filter(pl.col("sequence").str.contains(STANDARD))
    ids = standard.select("protein_id").collect()["protein_id"].to_list()
    chosen = sorted(random.Random(args.seed).sample(ids, min(args.n, len(ids))))
    picked = standard.filter(pl.col("protein_id").is_in(chosen)).collect()
    for i, part in enumerate(np.array_split(np.arange(picked.height), args.chunks)):
        _fasta(picked[part], "protein_id", f"queries_{i:03d}.faa")


def pairs(args: argparse.Namespace) -> None:
    hits = pl.concat(
        pl.read_csv(p, separator="\t", has_header=False, schema=HIT_SCHEMA) for p in args.hits
    )
    best = (
        hits.group_by("protein_id", "cluster_rep")
        .agg(pl.all().sort_by("bitscore").last())
        .filter(
            pl.col("protein_id") != pl.col("cluster_rep"),
            pl.col("length") / pl.col("qlen") >= args.min_cov,
            pl.col("length") / pl.col("slen") >= args.min_cov,
            pl.col("pident") / 100 >= args.min_id,
        )
        .select("protein_id", "cluster_rep", identity=pl.col("pident") / 100)
        .with_columns(bin=((pl.col("identity") - args.min_id) / args.bin).floor())
    )
    (
        best.sort("protein_id", "cluster_rep")
        .with_columns(order=pl.struct("protein_id", "cluster_rep").hash(args.seed))  # shuffled
        .filter(pl.col("order").rank("ordinal").over("bin") <= args.per_bin)
        .drop("bin", "order")
        .write_parquet(args.out)
    )


def _hashes(df: pl.DataFrame, k: int, alphabet: str) -> pl.DataFrame:
    """(``protein_id``, ``pos``, ``hash``): every k-mer of each sequence in emission order,
    which is its window position for sequences of standard residues."""
    got = _core.hash_proteins([s.encode() for s in df["sequence"]], k, alphabet=alphabet)
    rows = df.select("protein_id", row=pl.int_range(pl.len(), dtype=pl.UInt64))
    return (
        pl.DataFrame({"row": got["seq"], "hash": got["hash"]})
        .with_columns(pos=pl.int_range(pl.len()).over("row"))
        .join(rows, on="row")
        .drop("row")
    )


def pair_survival(
    pairs: pl.DataFrame, members: pl.DataFrame, k: int, alphabet: str
) -> pl.DataFrame:
    """``pairs`` (``protein_id`` P, ``cluster_rep`` C, ``identity``) with ``shared`` (P's
    distinct k-mers in C's union), ``pin_sum``, ``n_members`` and ``survival``, P left out of
    C when it is a member; and, over P's windows in order, ``window_survival`` (the share in
    C's union) and ``both_<j>`` (the share of window pairs j apart both in it). ``members``
    holds C's members and every P. ``pin_sum`` averages over the members ``p_in`` counts
    (full-length ones, or all when none is), as the index does.
    ponytail: no P_IN floor for k-mers only fragments hold; it moves pin_sum by < 1 k-mer."""
    in_c = members.join(pairs.select("cluster_rep").unique(), on="cluster_rep").with_columns(
        counted=pl.col("full_length") | ~pl.col("full_length").any().over("cluster_rep")
    )
    member_kmers = (
        _hashes(in_c, k, alphabet)
        .select("protein_id", "hash")
        .unique()
        .join(in_c.select("protein_id", "cluster_rep"), on="protein_id")
    )
    holders = member_kmers.group_by("cluster_rep", "hash").agg(n=pl.len())
    sizes = (
        member_kmers.group_by("protein_id")
        .agg(n_kmers=pl.len())
        .join(in_c.select("protein_id", "cluster_rep", "counted"), on="protein_id", how="right")
        .with_columns(pl.col("n_kmers").fill_null(0))
    )
    totals = sizes.group_by("cluster_rep").agg(
        total=(pl.col("n_kmers") * pl.col("counted").cast(pl.UInt32)).sum(),
        counted=pl.col("counted").sum(),
        n_members=pl.len(),
    )
    # P's own row in C, if it is a member: left out of the union, the average and the count
    own = pairs.join(sizes, on=["protein_id", "cluster_rep"], how="left").select(
        "protein_id",
        "cluster_rep",
        member=pl.col("counted").is_not_null(),
        own_kmers=(pl.col("n_kmers") * pl.col("counted").cast(pl.UInt32)).fill_null(0),
        own_counted=pl.col("counted").fill_null(False),
    )
    query = members.join(pairs.select("protein_id").unique(), on="protein_id")
    windows = (
        pairs.select("protein_id", "cluster_rep")
        .join(_hashes(query, k, alphabet), on="protein_id")
        .join(holders, on=["cluster_rep", "hash"], how="left")
        .join(own, on=["protein_id", "cluster_rep"])
        .with_columns(found=pl.col("n").fill_null(0) > pl.col("member").cast(pl.UInt32))
        .sort("protein_id", "cluster_rep", "pos")
    )
    key = ["protein_id", "cluster_rep"]
    shared = (
        windows.filter("found").select(*key, "hash").unique().group_by(key).agg(shared=pl.len())
    )
    lagged = windows.group_by(key, maintain_order=True).agg(
        window_survival=pl.col("found").mean(),
        **{
            f"both_{j}": pl.when(pl.col("found").shift(-j).is_not_null())
            .then(pl.col("found") & pl.col("found").shift(-j))
            .mean()
            for j in LAGS
        },
    )
    return (
        pairs.join(own, on=key)
        .join(totals, on="cluster_rep")
        .join(shared, on=key, how="left")
        .join(lagged, on=key, how="left")
        .select(
            *key,
            "identity",
            "window_survival",
            *(f"both_{j}" for j in LAGS),
            shared=pl.col("shared").fill_null(0),
            pin_sum=(pl.col("total") - pl.col("own_kmers"))
            / (pl.col("counted") - pl.col("own_counted").cast(pl.UInt32)),
            n_members=pl.col("n_members") - pl.col("member").cast(pl.UInt32),
        )
        .filter(pl.col("n_members") > 0, pl.col("pin_sum") > 0)
        .with_columns(survival=pl.col("shared") / pl.col("pin_sum"))
    )


def survival_step(args: argparse.Namespace) -> None:
    pairs_df = pl.read_parquet(args.pairs)
    clusters = pairs_df["cluster_rep"].unique()
    wanted = pl.concat([clusters, pairs_df["protein_id"].unique()]).unique()
    # every member of C, and P wherever it sits
    members = (
        _members(args.members)
        .filter(
            pl.col("cluster_rep").is_in(clusters.implode())
            | pl.col("protein_id").is_in(wanted.implode())
        )
        .collect()
    )
    residues = (
        members.filter(pl.col("cluster_rep").is_in(clusters.implode()))
        .group_by("cluster_rep")
        .agg(r=pl.col("sequence").str.len_bytes().sum())
        .sort("cluster_rep")
        .with_columns(batch=pl.col("r").cum_sum() // args.batch_residues)
    )
    parts = []
    for (batch,), group in residues.group_by("batch", maintain_order=True):
        part = pairs_df.join(group.select("cluster_rep"), on="cluster_rep")
        involved = members.filter(
            pl.col("cluster_rep").is_in(group["cluster_rep"].implode())
            | pl.col("protein_id").is_in(part["protein_id"].implode())
        )
        parts.append(pair_survival(part, involved, args.k, args.alphabet))
        print(f"batch {batch}: {part.height} pairs", flush=True)
    pl.concat(parts).write_parquet(args.out)


def predicted(model: SurvivalModel, identity: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """S and B_j (pairs x ``LAGS``) at each identity under ``model``."""
    tab = model.tables
    rows = np.interp(identity, tab["a"][::-1], np.arange(len(tab["a"]))[::-1].astype(float))
    row = np.rint(rows).astype(int)
    surv = tab["S"][row]
    rho = tab["rho"][row][:, [j - 1 for j in LAGS]]
    return surv, surv[:, None] ** 2 + rho * (surv * (1 - surv))[:, None]


def fit_model(df: pl.DataFrame, k: int, categories: int = CATEGORIES) -> dict[str, float]:
    """Least squares of survival and co-survival (``both_<j>``) against the model, over
    (log shape, log region); with ``rmse`` of each and ``rmse_independent`` (a^k)."""
    a = df["identity"].to_numpy()
    s = df["survival"].to_numpy()
    both = df.select(f"both_{j}" for j in LAGS).to_numpy()
    ok = ~np.isnan(both)

    def residuals(model: SurvivalModel) -> tuple[np.ndarray, np.ndarray]:
        surv, pair = predicted(model, a)
        return s - surv, np.where(ok, both - pair, 0.0)

    def loss(x: np.ndarray) -> float:
        model = SurvivalModel(k, float(np.exp(x[0])), 1 + float(np.exp(x[1])), categories)
        r_s, r_b = residuals(model)
        return float((r_s**2).sum() + (r_b**2).sum() / len(LAGS))

    best: Any = min(
        (minimize(loss, x0, method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-9})
         for x0 in ([0.0, np.log(10)], [np.log(0.3), np.log(50)], [np.log(3), np.log(3)])),
        key=lambda r: r.fun,
    )  # fmt: skip
    model = SurvivalModel(k, float(np.exp(best.x[0])), 1 + float(np.exp(best.x[1])), categories)
    r_s, r_b = residuals(model)
    return {
        "shape": model.shape or 0.0,
        "region": model.region,
        "rmse": float(np.sqrt((r_s**2).mean())),
        "rmse_both": float(np.sqrt((r_b[ok] ** 2).mean())) if ok.any() else float("nan"),
        "rmse_independent": float(np.sqrt(((s - a**k) ** 2).mean())),
    }


def fit(args: argparse.Namespace) -> None:
    pairs_df = pl.concat(pl.read_parquet(p) for p in args.survival)
    strata = [("all", pl.lit(True))]
    strata += [
        (f"identity {lo}-{min(hi, 1.0)}", pl.col("identity").is_between(lo, hi, closed="left"))
        for lo, hi in IDENTITY_BANDS
    ]
    strata += [
        (
            f"members {lo}-{hi - 1}" if hi < 2**62 else f"members {lo}+",
            pl.col("n_members").is_between(lo, hi, closed="left"),
        )
        for lo, hi in SIZE_BANDS
    ]
    rows = []
    for name, where in strata:
        part = pairs_df.filter(where)
        if name == "all" or part.height >= args.min_stratum:
            rows.append({"stratum": name, "n": part.height, **fit_model(part, args.k)})
    pl.DataFrame(rows).write_csv(args.strata_out, separator="\t")
    model = {
        "survival": "markov_gamma",
        "shape": rows[0]["shape"],
        "region": rows[0]["region"],
        "categories": CATEGORIES,
        "k": args.k,
        "alphabet": args.alphabet,
        "fit": {key: rows[0][key] for key in ("n", "rmse", "rmse_both", "rmse_independent")},
    }
    Path(args.out).write_text(json.dumps(model, indent=1) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    sub = parser.add_subparsers(dest="step", required=True)
    p = sub.add_parser("queries")
    p.add_argument("--members", required=True, nargs="+", help="parquet files or directories")
    p.add_argument("--n", type=int, default=200_000)
    p.add_argument("--chunks", type=int, default=20)
    p.add_argument("--seed", type=int, default=1)
    p = sub.add_parser("pairs")
    p.add_argument("--hits", required=True, nargs="+")
    p.add_argument("--min-cov", type=float, default=0.8)
    p.add_argument("--min-id", type=float, default=0.6)
    p.add_argument("--bin", type=float, default=0.02)
    p.add_argument("--per-bin", type=int, default=5000)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--out", default="pairs.parquet")
    p = sub.add_parser("survival")
    p.add_argument("--pairs", required=True)
    p.add_argument("--members", required=True, nargs="+")
    p.add_argument("--k", type=int, default=11)
    p.add_argument("--alphabet", default="protein")
    p.add_argument("--batch-residues", type=int, default=50_000_000)
    p.add_argument("--out", default="survival.parquet")
    p = sub.add_parser("fit")
    p.add_argument("--survival", required=True, nargs="+")
    p.add_argument("--k", type=int, default=11)
    p.add_argument("--alphabet", default="protein")
    p.add_argument("--min-stratum", type=int, default=200)
    p.add_argument("--out", default="aai_model.json")
    p.add_argument("--strata-out", default="model_strata.tsv")
    args = parser.parse_args()
    {"queries": queries, "pairs": pairs, "survival": survival_step, "fit": fit}[args.step](args)


if __name__ == "__main__":
    main()
