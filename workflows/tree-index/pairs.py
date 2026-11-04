"""Tree index: a k-mer pair filter after tier 1 (docs/tree-index-plan.md, Progress log: a
k-mer pair filter after tier 1).

Tier 1 (``query.py --decoy``) is a recall-first profile: its candidates (``cand_fam_<share>``)
are leaves seeded by any read pair at the decoy bar, a layer key standing for its group. Most
are near or chance relatives, and tier 2 is dense over every leaf of what it is given, so its
work follows them. This step keeps a candidate group only when some read pair holds a *pair*
of its k-mers in line, a cheap stand-in for an ungapped alignment:

- ``--features pair``: k-mer pairs ``--min-gap`` to ``--gap`` apart in order (``--combo 3``:
  triples), the gap not hashed (indels keep a pair), hashed as one feature
  (``_core.pair_route``), so one shared feature is two seeds close together; sampling the
  pair hashes (``--sample mod``: divisible by m; ``min``: window minimizers in position order)
  keeps a pair at the sampling rate, where sampling k-mers keeps it at the rate squared. With
  ``--soft f`` > 0, a group's features are those held by at least f of its leaves (modulo
  sampling, consistent across leaves; a soft intersection, as the layer's composites), not
  its representative's.
- ``--features kmer`` (default): k-mers sampled by hash, two or more votes on a diagonal band
  (``_core.route_reads``).

- **On the fly:** the reads are read again from the FASTQ in batches and translated
  (``_core.pair_frames``); nothing per read is written. Decoys are the same reads reversed, in
  memory.
- **Sparse, sized to memory:** the sampling rate is set so the index fits ``--memory-gb``
  (8 B per pair entry, 12 B per k-mer entry) whatever tier 1 hands over; a read looks up all
  its features.
- **Calibrated:** each candidate's score is the most votes any one read pair gives it on two
  adjacent diagonal bands; the bar is the least score at which decoys (scaled by
  ``--chance-scale``, the tail extended geometrically) are at most ``share`` of the sample's
  candidates, as ``query.decoy_bar``.

Writes to ``--query``: ``cand_pairs_<share>.npy`` (candidate positions kept, for
``recall.py score`` as strategy ``pairs_<share>``) and ``pairs.json`` (index size, modulo,
bars, candidates in and out, timings, peak RSS).
"""

import argparse
import gzip
import json
import math
import resource
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
import query
import verify

from kmer_functional_profiler import _core

ENTRY_BYTES = {"kmer": 12, "pair": 8}  # (hash, rep, position) or (feature, rep), u32 each


def read_batches(paths: list[str], batch: int, reverse: bool = False):  # noqa: ANN201
    """Read pairs' frames, ``batch`` pairs at a time (decoys: each read reversed)."""
    m1: list[bytes] = []
    m2: list[bytes] = []
    with gzip.open(paths[0], "rb") as f, gzip.open(paths[1], "rb") as g:
        for i, (x, y) in enumerate(zip(f, g, strict=True)):
            if i % 4 == 1:
                x, y = x.strip(), y.strip()
                m1.append(x[::-1] if reverse else x)
                m2.append(y[::-1] if reverse else y)
                if len(m1) >= batch:
                    yield _core.pair_frames(m1, m2)
                    m1, m2 = [], []
    if m1:
        yield _core.pair_frames(m1, m2)


def best_votes(
    paths: list[str], target: dict, n_targets: int, a: argparse.Namespace, reverse: bool
) -> tuple[np.ndarray, int]:
    """Each target's best score from any one read pair (votes, or shared pair features; 0:
    none), and the read pairs seen. ``target``: ``seqs`` (representatives) with ``modulo`` /
    ``minimizer``, or a prebuilt ``feature``, ``rep`` index (a soft intersection)."""
    best = np.zeros(n_targets, dtype=np.int64)
    n = 0
    classes = verify.alphabet_classes(a.alphabet)
    for frames in read_batches(paths, a.batch, reverse):
        if a.features == "kmer":
            _, rep, votes = _core.route_reads(
                frames, target["seqs"], a.k, classes, band=a.band, min_votes=2, top=a.top,
                max_occ=a.max_occ, modulo=target["modulo"],
            )  # fmt: skip
        elif "feature" in target:
            _, rep, votes = _core.pair_route_index(
                frames, target["feature"], target["rep"], a.k, classes,
                gap=a.gap, min_gap=a.min_gap, combo=a.combo,
                modulo=target["modulo"], min_hits=1, top=a.top, max_occ=a.max_occ,
                min_bits=a.min_bits,
            )  # fmt: skip
        else:
            _, rep, votes = _core.pair_route(
                frames, target["seqs"], a.k, classes,
                gap=a.gap, min_gap=a.min_gap, combo=a.combo, minimizer=target["minimizer"],
                modulo=target["modulo"], min_hits=1, top=a.top, max_occ=a.max_occ,
                min_bits=a.min_bits,
            )  # fmt: skip
        np.maximum.at(best, rep.astype(np.int64), votes.astype(np.int64))
        n += len(frames)
    return best, n


SOFT_CHUNK = 4000  # leaves a chunk (whole groups), for --soft


def soft_index(
    index: Path, reps: np.ndarray, group: np.ndarray, a: argparse.Namespace, budget: float
) -> tuple[dict, int]:
    """Each candidate group's pair features held by at least ``--soft`` of its leaves (all of
    a one-leaf group's), modulo-sampled so the leaves' features fit the budget; and the
    features before the cut (sizing)."""
    leaves = np.flatnonzero(np.isin(group, reps))
    seqs = verify.rep_seqs(index, leaves, a.members, a.store, own=True)
    leaves = np.array([x for x in leaves.tolist() if x in seqs], dtype=np.int64)
    ls = [seqs[x] for x in leaves.tolist()]
    total = int(
        _core.pair_feature_counts(
            ls,
            a.k,
            verify.alphabet_classes(a.alphabet),
            gap=a.gap,
            min_gap=a.min_gap,
            combo=a.combo,
        ).sum()
    )
    modulo = a.modulo or max(1, math.ceil(total / budget))
    rep_of = np.searchsorted(reps, group[leaves])  # each leaf's target
    o = np.argsort(rep_of, kind="stable")
    leaves, rep_of, ls = leaves[o], rep_of[o], [ls[x] for x in o.tolist()]
    # a chunk of whole groups at a time: only each group's conserved features are kept
    cuts = np.flatnonzero(np.r_[True, rep_of[1:] != rep_of[:-1]])
    bounds = [
        *cuts[:: max(1, len(cuts) // max(1, len(leaves) // SOFT_CHUNK))].tolist(),
        len(leaves),
    ]
    got = []
    for lo, hi in zip(bounds[:-1], bounds[1:], strict=True):
        i, h = _core.pair_features(
            ls[lo:hi],
            a.k,
            verify.alphabet_classes(a.alphabet),
            gap=a.gap,
            min_gap=a.min_gap,
            combo=a.combo,
            modulo=modulo,
        )
        r = rep_of[lo:hi]
        df = pl.DataFrame({"t": r[i.astype(np.int64)], "h": h})
        size = pl.DataFrame({"t": r}).group_by("t").agg(m=pl.len())
        got.append(
            df.group_by("t", "h").agg(c=pl.len()).join(size, on="t")
            .filter(pl.col("c") >= pl.max_horizontal(1.0, a.soft * pl.col("m")))
        )  # fmt: skip
    kept = pl.concat(got)
    target = {"feature": kept["h"].to_numpy().astype(np.uint32),
              "rep": kept["t"].to_numpy().astype(np.uint32), "modulo": modulo}  # fmt: skip
    return target, total


def key_index(index: Path, reps: np.ndarray, a: argparse.Namespace) -> tuple[dict, int]:
    """Each group's pair features (``--combo``) from the ``--k``-mers inside its stored tier 1
    keys only, at their positions on the leaf they are posted on (the group's representative
    position): no sequence beyond the keys' residues decides the features."""
    idx = query.Index(str(index))
    lo, h, code = np.asarray(idx.lo), np.asarray(idx.hash), np.asarray(idx.kcode)
    on = np.isin(lo, reps)
    lo, h, code = lo[on].astype(np.int64), h[on], code[on]
    kk, _, alpha, _ = query.ladder(idx.meta["params"]["layers"])[0]  # the layer's k-mers
    seqs = verify.rep_seqs(index, reps, a.members, a.store, own=True)
    reps = np.array([r for r in reps.tolist() if r in seqs], dtype=np.int64)
    rs, masks = [], []
    for r in reps.tolist():
        s = seqs[r]
        n = max(len(s) - kk + 1, 0)
        o = _core.hash_proteins([s[i : i + kk] for i in range(n)], kk, alphabet=alpha)
        hs = np.zeros(n, dtype=np.uint64)
        hs[o["seq"]] = o["hash"]
        hit = np.zeros(n, dtype=bool)
        for c in np.unique(code[lo == r]):
            hit |= np.isin(query.salted(hs, int(c)), h[(lo == r) & (code == c)])
        m = np.zeros(len(s), dtype=np.uint8)
        for p in np.flatnonzero(hit):
            m[p : p + kk - a.k + 1] = 1  # the k-mers inside the key
        rs.append(s)
        masks.append(m.tobytes())
    i, f = _core.pair_features_masked(
        rs,
        masks,
        a.k,
        verify.alphabet_classes(a.alphabet),
        gap=a.gap,
        min_gap=a.min_gap,
        combo=a.combo,
    )
    return {"feature": f, "rep": i, "modulo": 1, "reps": reps}, len(f)


def vote_bar(
    sample: np.ndarray, decoy: np.ndarray, share: float, scale: float, least: int = 2
) -> int:
    """The least votes at which decoy-seeded candidates (x ``scale``) are at most ``share`` of
    the sample's. Past the decoys observed, their count falls geometrically, by the mean
    ratio of successive counts over the tail observed (ponytail: one rate for the tail)."""
    top = int(max(sample.max(initial=0), decoy.max(initial=0))) + 2
    d = np.array([(decoy >= v).sum() for v in range(top + 1)], dtype=np.float64)
    s = np.array([(sample >= v).sum() for v in range(top + 1)], dtype=np.float64)
    seen = np.flatnonzero(d > 0)
    if scale != 1.0 and len(seen) >= 2:
        hi = seen[-1]
        tail = d[max(least, hi - 3) : hi + 1]
        ratio = float(np.mean(tail[1:] / tail[:-1])) if len(tail) > 1 else 0.5
        ratio = min(max(ratio, 1e-3), 0.9)
        for v in range(hi + 1, top + 1):
            d[v] = d[hi] * ratio ** (v - hi)
    for v in range(least, top + 1):
        if scale * d[v] <= share * max(s[v], 1):
            return v
    return top


def main() -> None:
    p = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    p.add_argument("--index", required=True)
    p.add_argument("--store", help="the index's leaf-ordered store (store.py)")
    p.add_argument("--members", help="the index's members (sequences), without --store")
    p.add_argument("--reads", nargs=2, required=True)
    p.add_argument("--query", required=True, help="query.py --decoy's --out (cand_fam_*)")
    p.add_argument(
        "--from", dest="source", default="fam", help="tier 1's candidates: cand_<this>_<share>"
    )
    p.add_argument("--features", choices=["pair", "kmer"], default="kmer")
    p.add_argument("--k", type=int, default=5, help="k-mers' k (each half of a pair)")
    p.add_argument("--gap", type=int, default=12, help="pair: k-mers paired up to this apart")
    p.add_argument("--combo", type=int, default=2, choices=[2, 3], help="pair: k-mers a feature")
    p.add_argument(
        "--min-gap", type=int, default=4,
        help="pair: and at least this apart (the gap is not hashed: indels keep a pair)",
    )  # fmt: skip
    p.add_argument("--sample", choices=["mod", "min"], default="mod", help="pair: sampling")
    p.add_argument("--soft", type=float, default=0.0, help="pair: > 0, group soft intersection")
    p.add_argument(
        "--alphabet", default="protein",
        help="pair: residue classes for features (reads and targets; e.g. the keys' alphabet)",
    )  # fmt: skip
    p.add_argument(
        "--key-derived", action="store_true",
        help="pair: features from the k-mers inside each group's stored tier 1 keys only",
    )  # fmt: skip
    p.add_argument("--band", type=int, default=8, help="kmer: diagonal band")
    p.add_argument("--top", type=int, default=64, help="representatives kept per read pair")
    p.add_argument("--max-occ", type=int, default=2000)
    p.add_argument("--memory-gb", type=float, default=4.0, help="the sparse index's budget")
    p.add_argument(
        "--min-bits", type=float, default=None,
        help="pair: score shared features by information, log2(targets / targets holding "
        "it), in quarter-bits; drop features below this (promiscuity trim)",
    )  # fmt: skip
    p.add_argument("--modulo", type=int, default=0, help="0: from the budget (dynamic)")
    p.add_argument("--minimizer", type=int, default=0, help="--sample min: 0: from the budget")
    p.add_argument("--batch", type=int, default=200_000)
    p.add_argument("--shares", default="0.01,0.05,0.1")
    p.add_argument("--chance-scale", type=float, default=1.0, help="query.decoy_bar's scale")
    a = p.parse_args()
    t0 = time.perf_counter()
    index, qdir = Path(a.index), Path(a.query)
    info: dict = {}
    for share in a.shares.split(","):
        cand = np.load(qdir / f"cand_{a.source}_{float(share)}.npy").astype(np.int64)
        gfile = index / "layer_group.npy"
        group = (
            np.load(gfile).astype(np.int64)
            if gfile.exists()
            else np.arange(cand.max(initial=0) + 1)
        )
        reps = np.unique(group[cand])
        budget = a.memory_gb * 1e9 / ENTRY_BYTES[a.features]
        t1 = time.perf_counter()
        if a.features == "pair" and a.key_derived:
            target, total = key_index(index, reps, a)
            reps, entries = target.pop("reps"), total
        elif a.features == "pair" and a.soft > 0:
            target, total = soft_index(index, reps, group, a, budget)
            entries = len(target["feature"])
        else:
            seqs = verify.rep_seqs(index, reps, a.members, a.store)
            reps = np.array([r for r in reps.tolist() if r in seqs], dtype=np.int64)
            rs = [seqs[r] for r in reps.tolist()]
            if a.features == "kmer":
                total = sum(max(len(x) - a.k + 1, 0) for x in rs)
                target = {"seqs": rs, "modulo": a.modulo or max(1, math.ceil(total / budget))}
                entries = total // target["modulo"]
            else:
                classes = verify.alphabet_classes(a.alphabet)
                total = int(
                    _core.pair_feature_counts(
                        rs, a.k, classes, gap=a.gap, min_gap=a.min_gap, combo=a.combo
                    ).sum()
                )
                # density: 1/m for modulo, ~2/(w+1) for window minimizers
                need = total / budget
                target = {"seqs": rs, "modulo": 1, "minimizer": 1}
                if a.sample == "mod":
                    target["modulo"] = a.modulo or max(1, math.ceil(need))
                else:
                    target["minimizer"] = a.minimizer or max(1, math.ceil(2 * need - 1))
                entries = int(_core.pair_feature_counts(
                    rs, a.k, classes,
                gap=a.gap, min_gap=a.min_gap, combo=a.combo, minimizer=target["minimizer"],
                    modulo=target["modulo"]).sum())  # fmt: skip
        t2 = time.perf_counter()
        sample, n = best_votes(a.reads, target, len(reps), a, reverse=False)
        t3 = time.perf_counter()
        decoy, _ = best_votes(a.reads, target, len(reps), a, reverse=True)
        least = 2 if a.features == "kmer" else 1
        bar = vote_bar(sample, decoy, float(share), a.chance_scale, least)
        kept_reps = reps[sample >= bar]
        keep = cand[np.isin(group[cand], kept_reps)]
        np.save(qdir / f"cand_pairs_{float(share)}.npy", np.sort(keep))
        info[share] = {
            "candidates_in": len(cand), "groups_in": len(reps), "candidates_out": len(keep),
            "groups_out": len(kept_reps), "features": a.features, "features_total": total,
            "modulo": target.get("modulo"), "minimizer": target.get("minimizer"),
            "index_entries": entries, "index_gb": entries * ENTRY_BYTES[a.features] / 1e9,
            "bar": bar, "read_pairs": n, "index_s": t2 - t1, "sample_s": t3 - t2,
            "decoy_s": time.perf_counter() - t3,
        }  # fmt: skip
        print(json.dumps({share: info[share]}), flush=True)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    info["peak_rss_gb"] = (rss if sys.platform == "darwin" else rss * 1024) / 2**30
    info["total_s"] = time.perf_counter() - t0
    (qdir / "pairs.json").write_text(json.dumps(info, indent=1) + "\n")


if __name__ == "__main__":
    main()
