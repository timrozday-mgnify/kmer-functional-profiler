"""Steps of the aai-model workflow (``main.nf``): the parameters of the k-mer survival model
that ``aai`` inverts (``survival.SurvivalModel``, Markov-beta: ``concentration`` phi of each
region's identity around the pair's, mean ``region`` length of a Markov chain of regions
along the sequence, and ``ends``, the share of windows inside the aligned region; plan,
phase 7, steps 25 and 33), fitted on MGnify protein pairs.

A pair is a protein P and an MGnify90 cluster C it aligns to. Its survival is the share of
P's windows (k-mers in sequence order) found in the union of C's members' k-mers: what the
query sees when P's gene is in a sample. Its identity is DIAMOND's, P against its nearest
member of C (step 33; against C's representative, ``identity_rep``, in step 25, which left
survival in large clusters far above what the identity predicts: P has closer members).
Survival against identity fixes how much survival exceeds a^k; it cannot by itself tell a
strong contrast over short regions from a weaker one over long regions. Co-survival does:
for lags j in ``LAGS``, the share of P's window pairs j apart both in C's union, against
the model's B_j.

- ``queries``: ``--n`` random members (seeded) whose sequences hold only the 20 standard
  residues (so a k-mer's emission order is its window position), as FASTA in ``--chunks``
  files, and every cluster's representative as ``reps.faa``.
- ``candidates``: DIAMOND hits against the representatives (outfmt 6 ``qseqid sseqid
  pident length qlen slen bitscore``) -> one row per (P, C), its best HSP, with both
  coverages >= ``--min-cov``, P not C's representative, at most ``--per-bin`` per bin of
  ``identity_rep``: candidates for the nearest-member pass, spread over identity.
- ``member-db``: the candidates' clusters' members as ``members.faa`` (and
  ``member_clusters.parquet``, member -> cluster), and their P as ``--chunks`` FASTA files.
- ``pairs``: DIAMOND hits of P against those members -> per candidate (P, C), the member
  of C (not P) with the highest identity among hits with both coverages >= ``--min-cov``:
  ``identity`` and ``nearest``; at most ``--per-bin`` pairs per identity bin of ``--bin``
  from ``--min-id``, so no band dominates the fit.
- ``survival``: per pair, P's k-mers in C's union and ``pin_sum``, and co-survival by lag
  (``both_<j>``, ``window_survival``: P's windows in C's union; ``n_windows``), every k-mer
  (no hash threshold: sampling under it is unbiased, all k-mers are more precise). When P
  is one of C's members it is left out of C (a gene new to the index). Clusters are
  processed in batches of about ``--batch-residues`` member residues.
- ``fit``: ``aai_model.json`` (``survival`` "markov_beta", ``concentration``, ``region``,
  ``ends``, ``categories``, ``k``, ``alphabet``, ``fit``), what ``kmer-functional-profiler
  aai-model`` attaches to an index, and ``model_strata.tsv``, the model refitted per
  identity band and per cluster size, to check that one parameter set holds.
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
    "subject": pl.Int64,
    "pident": pl.Float64,
    "length": pl.Int64,
    "qlen": pl.Int64,
    "slen": pl.Int64,
    "bitscore": pl.Float64,
}
MEMBER_COLUMNS = ["protein_id", "cluster_rep", "full_length", "sequence"]
IDENTITY_BANDS = ((0.6, 0.7), (0.7, 0.8), (0.8, 0.9), (0.9, 1.01))
SIZE_BANDS = ((1, 2), (2, 4), (4, 11), (11, 101), (101, 2**62))  # members, P left out
LAGS = (1, 3, 6, 11, 20, 40, 80, 150)  # co-survival lags, windows
STANDARD = "^[ACDEFGHIKLMNPQRSTVWY]+$"
BIN = 0.01  # identity bins the fit averages pairs over


def _members(paths: list[str]) -> pl.LazyFrame:
    paths = [str(Path(p) / "*.parquet") if Path(p).is_dir() else p for p in paths]
    return pl.scan_parquet(paths).select(MEMBER_COLUMNS)


def _fasta(df: pl.DataFrame, name: str, path: str | Path) -> None:
    df.select(pl.format(">{}\n{}", name, "sequence")).write_csv(
        path, include_header=False, quote_style="never"
    )


def _chunks(df: pl.DataFrame, chunks: int, prefix: str) -> None:
    for i, part in enumerate(np.array_split(np.arange(df.height), chunks)):
        _fasta(df[part], "protein_id", f"{prefix}_{i:03d}.faa")


def _hits(paths: list[str], min_cov: float) -> pl.DataFrame:
    """DIAMOND hits with both coverages >= ``min_cov``, ``identity`` a share."""
    hits = pl.concat(
        pl.read_csv(p, separator="\t", has_header=False, schema=HIT_SCHEMA) for p in paths
    )
    return hits.filter(
        pl.col("length") / pl.col("qlen") >= min_cov,
        pl.col("length") / pl.col("slen") >= min_cov,
    ).with_columns(identity=pl.col("pident") / 100)


def _per_bin(df: pl.DataFrame, column: str, args: argparse.Namespace) -> pl.DataFrame:
    """At most ``args.per_bin`` rows (seeded shuffle) per bin of ``args.bin`` of ``column``
    from ``args.min_id``; rows below ``args.min_id`` dropped."""
    return (
        df.filter(pl.col(column) >= args.min_id)
        .sort("protein_id", "cluster_rep")
        .with_columns(
            bin=((pl.col(column) - args.min_id) / args.bin).floor(),
            order=pl.struct("protein_id", "cluster_rep").hash(args.seed),  # shuffled
        )
        .filter(pl.col("order").rank("ordinal").over("bin") <= args.per_bin)
        .drop("bin", "order")
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
    _chunks(standard.filter(pl.col("protein_id").is_in(chosen)).collect(), args.chunks, "queries")


def candidates(args: argparse.Namespace) -> None:
    best = (
        _hits(args.hits, args.min_cov)
        .rename({"subject": "cluster_rep"})
        .filter(pl.col("protein_id") != pl.col("cluster_rep"))
        .group_by("protein_id", "cluster_rep")
        .agg(pl.col("identity").sort_by("bitscore").last().alias("identity_rep"))
    )
    _per_bin(best, "identity_rep", args).write_parquet(args.out)


def member_db(args: argparse.Namespace) -> None:
    cands = pl.read_parquet(args.candidates)
    members = _members(args.members)
    clusters = cands["cluster_rep"].unique().implode()
    in_c = members.filter(pl.col("cluster_rep").is_in(clusters))
    in_c.select("protein_id", "cluster_rep").sink_parquet("member_clusters.parquet")
    in_c.select(pl.format(">{}\n{}", "protein_id", "sequence")).sink_csv(
        "members.faa", include_header=False, quote_style="never"
    )
    query = members.filter(pl.col("protein_id").is_in(cands["protein_id"].unique().implode()))
    _chunks(query.collect(), args.chunks, "nearest")


def pairs(args: argparse.Namespace) -> None:
    cands = pl.read_parquet(args.candidates)
    clusters = pl.read_parquet(args.member_clusters).rename({"protein_id": "subject"})
    nearest = (
        _hits(args.hits, args.min_cov)
        .filter(pl.col("protein_id") != pl.col("subject"))
        .join(clusters, on="subject")
        .join(cands, on=["protein_id", "cluster_rep"])
        .group_by("protein_id", "cluster_rep")
        .agg(
            pl.all().sort_by("identity", "bitscore").last(),
        )
        .select("protein_id", "cluster_rep", "identity", "identity_rep", nearest="subject")
    )
    print(f"{cands.height} candidates, {nearest.height} with a member hit", flush=True)
    _per_bin(nearest, "identity", args).write_parquet(args.out)


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
    """``pairs`` (``protein_id`` P, ``cluster_rep`` C, ``identity``, any other columns) with
    ``shared`` (P's distinct k-mers in C's union), ``pin_sum``, ``n_members`` and
    ``survival``, P left out of C when it is a member; and, over P's windows in order,
    ``n_windows``, ``window_survival`` (the share in C's union) and ``both_<j>`` (the share
    of window pairs j apart both in it). ``members`` holds C's members and every P.
    ``pin_sum`` averages over the members ``p_in`` counts (full-length ones, or all when
    none is, as the index does), after P is left out.
    ponytail: no P_IN floor for k-mers only fragments hold; it moves pin_sum by < 1 k-mer."""
    in_c = members.join(pairs.select("cluster_rep").unique(), on="cluster_rep")
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
        .join(in_c.select("protein_id", "cluster_rep", "full_length"), on="protein_id", how="right")
        .with_columns(pl.col("n_kmers").fill_null(0))
    )
    full = pl.col("full_length").cast(pl.UInt32)
    totals = sizes.group_by("cluster_rep").agg(
        total=pl.col("n_kmers").sum(),
        total_full=(pl.col("n_kmers") * full).sum(),
        n_full=full.sum(),
        n_members=pl.len(),
    )
    # P's own row in C, if it is a member: left out of the union, the averages and the counts
    own = pairs.join(sizes, on=["protein_id", "cluster_rep"], how="left").select(
        "protein_id",
        "cluster_rep",
        member=pl.col("full_length").is_not_null(),
        own_kmers=pl.col("n_kmers").fill_null(0),
        own_full=pl.col("full_length").fill_null(False).cast(pl.UInt32),
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
        n_windows=pl.len(),
        window_survival=pl.col("found").mean(),
        **{
            f"both_{j}": pl.when(pl.col("found").shift(-j).is_not_null())
            .then(pl.col("found") & pl.col("found").shift(-j))
            .mean()
            for j in LAGS
        },
    )
    left = pl.col("n_members") - pl.col("member").cast(pl.UInt32)
    full_left = pl.col("n_full") - pl.col("own_full")
    return (
        pairs.join(own, on=key)
        .join(totals, on="cluster_rep")
        .join(shared, on=key, how="left")
        .join(lagged, on=key, how="left")
        .select(
            *key,
            *(c for c in pairs.columns if c not in key),
            "n_windows",
            "window_survival",
            *(f"both_{j}" for j in LAGS),
            shared=pl.col("shared").fill_null(0),
            # full-length members without P, or all of them when P was the only one (else
            # 0 / 0: step 25's run lost 295 pairs at identity ~1 to NaN this way)
            pin_sum=pl.when(full_left > 0)
            .then((pl.col("total_full") - pl.col("own_kmers") * pl.col("own_full")) / full_left)
            .otherwise((pl.col("total") - pl.col("own_kmers")) / left),
            n_members=left,
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
    """S and B_j (``identity`` x ``LAGS``) at each identity under ``model``."""
    tab = model.tables
    rows = np.interp(identity, tab["a"][::-1], np.arange(len(tab["a"]))[::-1].astype(float))
    row = np.rint(rows).astype(int)
    surv = tab["S"][row]
    rho = tab["rho"][row][:, [j - 1 for j in LAGS]]
    return surv, surv[:, None] ** 2 + rho * (surv * (1 - surv))[:, None]


def binned(df: pl.DataFrame, min_pairs: int = 20) -> pl.DataFrame:
    """Pairs averaged per ``BIN`` of identity (bins of >= ``min_pairs``): ``n``, ``a``,
    ``S`` (window survival), ``both_<j>``, ``n_<j>`` (pairs with windows j apart) and
    ``inv_<j>``, their mean 1 / (``n_windows`` - j), for end loss at lag j.

    Pair-to-pair scatter (~0.08 in S) swamps the model's error, so fitting single pairs
    cannot tell models apart; the bins' means can (step 33)."""
    both = [pl.col(f"both_{j}") for j in LAGS]
    return (
        df.filter(pl.col("window_survival").is_not_nan())
        .group_by(bin=(pl.col("identity") / BIN).floor())
        .agg(
            *(b.mean() for b in both),
            n=pl.len(),
            a=pl.col("identity").mean(),
            S=pl.col("window_survival").mean(),
            **{f"n_{j}": b.is_not_null().sum() for j, b in zip(LAGS, both, strict=True)},
            **{
                f"inv_{j}": (1 / (pl.col("n_windows") - j)).filter(b.is_not_null()).mean()
                for j, b in zip(LAGS, both, strict=True)
            },
        )
        .filter(pl.col("n") >= min_pairs)
        .sort("bin")
    )


def fit_model(
    df: pl.DataFrame,
    k: int,
    categories: int = CATEGORIES,
    min_pairs: int = 20,
    ends: float | None = None,
) -> dict[str, float]:
    """Markov-beta with end loss, by weighted least squares on the identity bins' means
    (:func:`binned`) of survival, c S(a), and co-survival, c_j B_j(a), over (logit c,
    log phi, log(region - 1)). P's aligned share is one block of c n windows, so two
    windows j apart both lie in it with chance c_j = c - j (1 - c) / (n - j). Returns the
    parameters, ``rmse`` and ``rmse_both`` (bins), ``rmse_independent`` (a^k, bins) and
    ``rmse_pairs`` (single pairs, for the scatter). ``ends`` fixes c: within a narrow band
    of identity c trades off against phi, so the strata take the overall fit's."""
    g = binned(df, min_pairs)
    a, s, w = g["a"].to_numpy(), g["S"].to_numpy(), g["n"].to_numpy() / g["n"].sum()
    both = g.select(f"both_{j}" for j in LAGS).fill_null(0).to_numpy()
    nb = g.select(f"n_{j}" for j in LAGS).to_numpy().astype(float)
    wb = nb / max(nb.sum(), 1.0)
    inv = g.select(f"inv_{j}" for j in LAGS).fill_null(0).to_numpy()
    lags = np.array(LAGS, dtype=float)

    def unpack(x: np.ndarray) -> tuple[float, SurvivalModel]:
        if ends is not None:
            x = np.r_[0.0, x]  # a placeholder for c, which is fixed
        c = float(1 / (1 + np.exp(-x[0]))) if ends is None else ends
        model = SurvivalModel(k, None, 1 + float(np.exp(x[2])), categories, float(np.exp(x[1])))
        return c, model

    def residuals(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        c, model = unpack(x)
        surv, pair = predicted(model, a)
        c_j = np.clip(c - lags * (1 - c) * inv, 0.0, 1.0)
        return s - c * surv, both - c_j * pair

    def loss(x: np.ndarray) -> float:
        r_s, r_b = residuals(x)
        return float((w * r_s**2).sum() + (wb * r_b**2).sum())

    starts = ([3.0, np.log(5), np.log(5)], [3.0, np.log(2), np.log(30)],
              [2.0, np.log(15), np.log(2)])  # fmt: skip
    best: Any = min(
        (minimize(loss, x0, method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-10})
         for x0 in (x if ends is None else x[1:] for x in starts)),
        key=lambda r: r.fun,
    )  # fmt: skip
    c, model = unpack(best.x)
    r_s, r_b = residuals(best.x)
    pairs_a = df.filter(pl.col("window_survival").is_not_nan())
    pa, ps = pairs_a["identity"].to_numpy(), pairs_a["window_survival"].to_numpy()
    return {
        "concentration": model.concentration or 0.0,
        "region": model.region,
        "ends": c,
        "n_bins": g.height,
        "rmse": float(np.sqrt((w * r_s**2).sum())),
        "rmse_both": float(np.sqrt((wb * r_b**2).sum() / max(wb.sum(), 1e-12))),
        "rmse_independent": float(np.sqrt((w * (s - a**k) ** 2).sum())),
        "rmse_pairs": float(np.sqrt(((ps - c * predicted(model, pa)[0]) ** 2).mean())),
    }


def _size_band(n_members: np.ndarray) -> np.ndarray:
    """The ``SIZE_BANDS`` index of each pair's cluster size."""
    return np.searchsorted([lo for lo, _ in SIZE_BANDS], n_members, side="right") - 1


def union_fit(
    df: pl.DataFrame,
    base: SurvivalModel,
    k: int,
    union: tuple[float, float] | None = None,
    min_pairs: int = 20,
) -> dict[str, float]:
    """The union term over a one-member ``base`` model (its c, phi and region): weighted
    least squares of the window survival of multi-member pairs, averaged per (``BIN`` of
    identity, size band), against the same pairs' 1 - (1 - S)^m, m = n^(g0 + g1 (a - 0.8))
    (:meth:`SurvivalModel.members`), over (log g0, log g1): both >= 0, so survival falls
    with identity. ``union`` (g0, g1) fixes them, for the error of the overall fit on a
    stratum. Returns ``union``, ``union_slope``, ``n_bins``, ``rmse``,
    ``rmse_independent`` (a^k) and ``rmse_pairs``."""
    df = df.filter(pl.col("window_survival").is_not_nan())
    a, y = df["identity"].to_numpy(), df["window_survival"].to_numpy()
    n = df["n_members"].to_numpy().astype(np.float64)
    group = np.unique(
        np.stack([np.floor(a / BIN), _size_band(n)], axis=1), axis=0, return_inverse=True
    )[1].ravel()
    count = np.bincount(group)
    ok = count >= min_pairs
    w = count[ok] / count[ok].sum()
    obs = np.bincount(group, y)[ok] / count[ok]
    log_miss = np.log1p(-np.minimum(base.survival(a), 1 - 1e-15))  # log(1 - c S), per pair

    def predict(x: np.ndarray) -> np.ndarray:
        g0, g1 = np.exp(x) if union is None else union
        power = np.maximum(g0 + g1 * (a - 0.8), 0.0)
        return np.asarray(-np.expm1(n**power * log_miss))

    def loss(x: np.ndarray) -> float:
        return float((w * (obs - np.bincount(group, predict(x))[ok] / count[ok]) ** 2).sum())

    x = np.zeros(2)
    if union is None:
        starts = (
            [np.log(0.1), np.log(0.3)],
            [np.log(0.3), np.log(0.01)],
            [np.log(0.05), np.log(2)],
        )
        best: Any = min(
            (minimize(loss, x0, method="Nelder-Mead", options={"xatol": 1e-4, "fatol": 1e-12})
             for x0 in starts),
            key=lambda r: r.fun,
        )  # fmt: skip
        x = best.x
    g0, g1 = np.exp(x) if union is None else union
    ak = np.bincount(group, a**k)[ok] / count[ok]
    return {
        "union": float(g0),
        "union_slope": float(g1),
        "n_bins": int(ok.sum()),
        "rmse": float(np.sqrt(loss(x))),
        "rmse_independent": float(np.sqrt((w * (obs - ak) ** 2).sum())),
        "rmse_pairs": float(np.sqrt(((y - predict(x)) ** 2).mean())),
    }


def _band(lo: int, hi: int) -> str:
    return (
        f"members {lo}+"
        if hi >= 2**62
        else f"members {lo}-{hi - 1}"
        if hi > lo + 1
        else f"members {lo}"
    )


def fit(args: argparse.Namespace) -> None:
    """Two stages (step 34). One member: c, phi and region from single-member pairs (end
    loss shows only there: other members cover P's ends), then the union term (g0, g1) from
    multi-member pairs over that model. ``model_strata.tsv``: the one-member model refitted
    per identity band at its c; per cluster size, the union term refitted (``union``,
    ``union_slope``) and the overall model's error there (``rmse_overall``)."""
    pairs_df = pl.concat(pl.read_parquet(p) for p in args.survival)
    one = pairs_df.filter(pl.col("n_members") == 1)
    many = pairs_df.filter(pl.col("n_members") > 1)
    base_fit = fit_model(one, args.k)
    base = SurvivalModel(
        args.k, None, base_fit["region"], CATEGORIES, base_fit["concentration"], base_fit["ends"]
    )
    rows: list[dict[str, Any]] = [{"stratum": "members 1", "n": one.height, **base_fit}]
    print(rows[-1], flush=True)
    for lo, hi in IDENTITY_BANDS:
        part = one.filter(pl.col("identity").is_between(lo, hi, closed="left"))
        if part.height >= args.min_stratum:
            got = fit_model(part, args.k, ends=base_fit["ends"])
            rows.append(
                {"stratum": f"members 1, identity {lo}-{min(hi, 1.0)}", "n": part.height, **got}
            )
            print(rows[-1], flush=True)
    union = union_fit(many, base, args.k) if many.height else {"union": 0.0, "union_slope": 0.0}
    overall = (union["union"], union["union_slope"])
    rows.append({"stratum": "members 2+", "n": many.height, **union})
    for lo, hi in SIZE_BANDS[1:]:
        part = many.filter(pl.col("n_members").is_between(lo, hi, closed="left"))
        if part.height >= args.min_stratum:
            rows.append({
                "stratum": _band(lo, hi), "n": part.height, **union_fit(part, base, args.k),
                "rmse_overall": union_fit(part, base, args.k, overall)["rmse"],
            })  # fmt: skip
            print(rows[-1], flush=True)
    pl.DataFrame(rows, infer_schema_length=None).write_csv(args.strata_out, separator="\t")
    model = {
        "survival": "markov_beta",
        "concentration": base_fit["concentration"],
        "region": base_fit["region"],
        "ends": base_fit["ends"],
        "union": union["union"],
        "union_slope": union["union_slope"],
        "categories": CATEGORIES,
        "k": args.k,
        "alphabet": args.alphabet,
        "fit": {
            "n_one": one.height,
            "rmse_one": base_fit["rmse"],
            "rmse_both_one": base_fit["rmse_both"],
            "n_union": many.height,
            "rmse_union": union.get("rmse", float("nan")),
            "rmse_independent_union": union.get("rmse_independent", float("nan")),
        },
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
    for name, default in (("candidates", 15000), ("pairs", 5000)):
        p = sub.add_parser(name)
        p.add_argument("--hits", required=True, nargs="+")
        if name == "pairs":
            p.add_argument("--candidates", required=True)
            p.add_argument("--member-clusters", required=True)
        p.add_argument("--min-cov", type=float, default=0.8)
        p.add_argument("--min-id", type=float, default=0.6)
        p.add_argument("--bin", type=float, default=0.02)
        p.add_argument("--per-bin", type=int, default=default)
        p.add_argument("--seed", type=int, default=1)
        p.add_argument("--out", default=f"{name}.parquet")
    p = sub.add_parser("member-db")
    p.add_argument("--candidates", required=True)
    p.add_argument("--members", required=True, nargs="+")
    p.add_argument("--chunks", type=int, default=20)
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
    steps = {"queries": queries, "candidates": candidates, "member-db": member_db,
             "pairs": pairs, "survival": survival_step, "fit": fit}  # fmt: skip
    steps[args.step](args)


if __name__ == "__main__":
    main()
