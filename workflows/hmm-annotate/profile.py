"""Each protein's cost by hmmsearch and by kfpfam (its filter, then HMMER on its candidates), on
the same proteins, for `--profile`. CPU seconds, from each thread's own clock.

- hmmsearch: the protein against every model (pyhmmer's scan, HMMER's pipeline as `--engine
  hmmsearch` runs it), so its passing models are the truth here too.
- kfpfam's filter: `filter join` with FILTER_PER_ORF beside it. Its diagonal, bound and exact
  gapped stages are timed per protein; the gate's time (the log's join less its survivor stage) is
  shared out by k-mer hits, which cost about the same each at any length, and the survivor sort's
  (its survivor stage less the diagonals timed) by diagonals.
- kfpfam's HMMER: each candidate pair alone, at `--forward-slack` as SEARCH runs it.

A call's own cost (the pipeline set up for each model) is paid once per model in a batch, so each
timed search has an empty one's taken off: a 1-residue protein scanned against every model, or
searched with that pair's model. Both are also run as SEARCH runs them, each model on all the
proteins (or on its candidates) at once, and each protein's time is scaled to sum to that
(`*_raw_s` are after the empty call, before the scale; `scale` in the JSON, near 1 if the
correction holds).

Writes per protein --proteins (Parquet), per pair (every kfpfam candidate and every true pair)
--pairs, and --run (JSON).
"""

import argparse
import csv
import json
import math
import re
import struct
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from annotate import forward_cut, read_thresholds, text
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, SequenceFile, TextSequence
from pyhmmer.plan7 import HMMPressedFile, OptimizedProfileBlock, Pipeline

ALPHABET = Alphabet.amino()
LOCAL = threading.local()
INF = float("inf")


def mine(pressed: str) -> OptimizedProfileBlock:
    """This thread's own profiles: a pipeline reconfigures one for each target, so threads can't
    share them."""
    if not hasattr(LOCAL, "block"):
        with HMMPressedFile(pressed) as f:
            LOCAL.block = OptimizedProfileBlock(ALPHABET, list(f))
    return LOCAL.block


def pipe(z: float, f3: float | None = None) -> Pipeline:
    if not hasattr(LOCAL, "pipeline"):
        LOCAL.pipeline = Pipeline(ALPHABET, Z=z, domZ=z, E=INF, domE=INF)
        LOCAL.f3 = LOCAL.pipeline.F3
    p = LOCAL.pipeline
    p.clear()
    p.F3 = LOCAL.f3 if f3 is None else min(LOCAL.f3, f3)
    return p


def passing(hit, seq_t: float, dom_t: float) -> bool:
    return hit.score >= seq_t and any(d.score >= dom_t for d in hit.domains)


def filter_log(path: Path) -> dict:
    """The join's totals, core-s: the gate is join less the survivor stage."""
    line = [x for x in path.read_text().splitlines() if " core-s; postings" in x][-1]
    num = lambda pattern: float(re.search(pattern, line).group(1))  # noqa: E731
    return {
        "join_s": num(r"join ([\d.]+) core-s"),
        "gapped_s": num(r"gapped ([\d.]+) core-s"),
        "survivor_s": num(r"survivor stage ([\d.]+) core-s"),
        "kernel": line.split(":")[0],
    }


def candidates(path: Path) -> dict[int, list[int]]:
    """FILTER_BINARY's pairs: protein (its index in the FASTA) -> candidate models."""
    out: dict[int, list[int]] = {}
    with open(path, "rb") as f:
        first = 0
        while head := f.read(24):
            n, _, size = struct.unpack("<3Q", head)
            pairs = np.frombuffer(f.read(size), "<u4").reshape(-1, 2)
            for o, m in pairs.tolist():
                out.setdefault(first + o, []).append(m)
            first += n
    return out


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--faa", type=Path, required=True)
    p.add_argument("--pressed", required=True)
    p.add_argument("--thresholds", type=Path, required=True)
    p.add_argument(
        "--pairs-bin", type=Path, required=True, help="filter join's FILTER_BINARY output"
    )
    p.add_argument("--per-orf", type=Path, required=True, help="its FILTER_PER_ORF table")
    p.add_argument("--filter-log", type=Path, required=True)
    p.add_argument("--forward-slack", type=float)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--z", type=float, default=1e6)
    p.add_argument("--proteins", type=Path, required=True)
    p.add_argument("--pairs", type=Path, required=True)
    p.add_argument("--run", type=Path, required=True)
    a = p.parse_args()

    thresholds = read_thresholds(a.thresholds)
    index = {t[1]: i for i, t in enumerate(thresholds)}  # HMMER's hits name the model by its NAME
    with SequenceFile(a.faa, digital=True, alphabet=ALPHABET) as f:
        seqs = list(f.read_block())
    ids = [int(text(s.name)) for s in seqs]
    lengths = [len(s) for s in seqs]
    with HMMPressedFile(a.pressed) as f:
        model_length = [x.M for x in f]
    groups = [range(i, min(i + 16, len(seqs))) for i in range(0, len(seqs), 16)]
    pool = ThreadPoolExecutor(a.threads)

    # an empty search's cost: a call's own, which a batch pays once per model, not once per
    # protein or pair
    one = TextSequence(name=b"0", sequence="A").digitize(ALPHABET)

    def empty_scan(_):
        block = mine(a.pressed)
        pipe(a.z)  # built outside the clock: building one is slow
        t = time.thread_time()
        for _ in range(3):
            pipe(a.z).scan_seq(one, block)
        return (time.thread_time() - t) / 3

    empty_scan_s = float(np.median(list(pool.map(empty_scan, range(a.threads)))))

    # hmmsearch, one protein at a time
    def scan(group):
        block = mine(a.pressed)
        out = []
        for i in group:
            t = time.thread_time()
            hits = pipe(a.z).scan_seq(seqs[i], block)
            dt = max(time.thread_time() - t - empty_scan_s, 0.0)
            true = [
                index[text(h.name)]
                for h in hits
                if passing(h, *thresholds[index[text(h.name)]][2:])
            ]
            out.append((i, dt, true))
        return out

    hmmsearch = {i: (dt, true) for part in pool.map(scan, groups) for i, dt, true in part}

    # kfpfam's HMMER, one pair at a time, less that model's empty call
    cand = candidates(a.pairs_bin)
    used = sorted({m for ms in cand.values() for m in ms})

    def empty_calls(models):
        profiles = mine(a.pressed)
        pipe(a.z)
        out = {}
        for m in models:
            t = time.thread_time()
            pipe(a.z).search_hmm(profiles[m], DigitalSequenceBlock(ALPHABET, [one]))
            out[m] = time.thread_time() - t
        return out

    empty = {
        m: v
        for part in pool.map(empty_calls, [used[k::64] for k in range(64)])
        for m, v in part.items()
    }

    def slack(profile, seq_t):
        """SEARCH's Forward cut (annotate.py search_model)"""
        return (
            forward_cut(profile, seq_t - a.forward_slack)
            if a.forward_slack is not None and math.isfinite(seq_t)
            else None
        )

    def pairs(group):
        profiles = mine(a.pressed)
        out = []
        for i in group:
            for m in cand.get(i, []):
                _, _, seq_t, dom_t = thresholds[m]
                cut = slack(profiles[m], seq_t)
                t = time.thread_time()
                hits = pipe(a.z, cut).search_hmm(
                    profiles[m], DigitalSequenceBlock(ALPHABET, [seqs[i]])
                )
                dt = max(time.thread_time() - t - empty[m], 0.0)
                out.append((i, m, dt, any(passing(h, seq_t, dom_t) for h in hits)))
        return out

    kf_pairs = [x for part in pool.map(pairs, groups) for x in part]

    # both as SEARCH runs them: each model on all the proteins (hmmsearch) or on its candidates
    # (kfpfam)
    # the blocks are built off the clock: hmmsearch builds its one block once, and a block's cost is
    # Python's, which for every model on a read's 300,000 ORFs would be most of it
    everything = DigitalSequenceBlock(ALPHABET, seqs)

    def batch(models, every):
        profiles = mine(a.pressed)
        spent = 0.0
        for m in models:
            if every:
                targets = everything
            elif m in by_model:
                targets = DigitalSequenceBlock(ALPHABET, [seqs[i] for i in by_model[m]])
            else:
                continue
            _, _, seq_t, _ = thresholds[m]
            p = pipe(a.z, None if every else slack(profiles[m], seq_t))
            t = time.thread_time()
            p.search_hmm(profiles[m], targets)
            spent += time.thread_time() - t
        return spent

    by_model: dict[int, list] = {}
    for i, ms in cand.items():
        for m in ms:
            by_model.setdefault(m, []).append(i)
    parts = [range(k, len(thresholds), 64) for k in range(64)]
    batch_hmmsearch = sum(pool.map(lambda r: batch(r, True), parts))
    batch_kfpfam = sum(pool.map(lambda r: batch(r, False), parts))

    # the filter's stages, per protein
    log = filter_log(a.filter_log)
    with open(a.per_orf) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    col = lambda k: np.array([float(r[k]) for r in rows])  # noqa: E731
    hits, diagonals = col("hits"), col("diagonals")
    diagonal_s, bound_s, gapped_s = (
        col("diagonal_ns") / 1e9,
        col("bound_ns") / 1e9,
        col("gapped_ns") / 1e9,
    )
    gate_s = max(log["join_s"] - log["survivor_s"], 0) * hits / max(hits.sum(), 1)
    sort_s = max(log["survivor_s"] - diagonal_s.sum(), 0) * diagonals / max(diagonals.sum(), 1)
    filter_s = gate_s + sort_s + diagonal_s + bound_s + gapped_s
    assert len(hits) == len(seqs), f"{a.per_orf} has {len(hits)} proteins, {a.faa} {len(seqs)}"

    kf_hmmer = np.zeros(len(seqs))
    found: dict[int, set] = {}
    for i, m, dt, ok in kf_pairs:
        kf_hmmer[i] += dt
        if ok:
            found.setdefault(i, set()).add(m)

    # What is left after the empty calls, scaled to the batch's total: a protein's search against
    # every model still costs a little more alone, as the profiles don't stay in cache
    raw = np.array([hmmsearch[i][0] for i in range(len(seqs))])
    scale_hmmsearch = batch_hmmsearch / max(raw.sum(), 1e-9)
    scale_kfpfam = batch_kfpfam / max(kf_hmmer.sum(), 1e-9)
    acc = [t[0] for t in thresholds]
    proteins = {
        "protein_id": ids,
        "length": lengths,
        "hmmsearch_s": raw * scale_hmmsearch,
        "hmmsearch_raw_s": raw,
        "true_models": [[acc[m] for m in hmmsearch[i][1]] for i in range(len(seqs))],
        "true_lengths": [[model_length[m] for m in hmmsearch[i][1]] for i in range(len(seqs))],
        "kf_hits": hits.astype(np.int64),
        "kf_diagonals": diagonals.astype(np.int64),
        "kf_bounded": col("bounded").astype(np.int64),
        "kf_gapped": col("gapped").astype(np.int64),
        "kf_gate_s": gate_s,
        "kf_sort_s": sort_s,
        "kf_diagonal_s": diagonal_s,
        "kf_bound_s": bound_s,
        "kf_gapped_s": gapped_s,
        "kf_filter_s": filter_s,
        "kf_candidates": [len(cand.get(i, [])) for i in range(len(seqs))],
        "kf_hmmer_s": kf_hmmer * scale_kfpfam,
        "kf_hmmer_raw_s": kf_hmmer,
        "kf_s": filter_s + kf_hmmer * scale_kfpfam,
        "kf_missed": [
            [acc[m] for m in hmmsearch[i][1] if m not in found.get(i, set())]
            for i in range(len(seqs))
        ],
        "kf_extra": [
            [acc[m] for m in found.get(i, set()) if m not in hmmsearch[i][1]]
            for i in range(len(seqs))
        ],
    }
    pq.write_table(pa.table(proteins), a.proteins, compression="zstd")

    # every candidate pair, and every true pair that was not one
    rows = [(i, m, dt, True, m in hmmsearch[i][1], ok) for i, m, dt, ok in kf_pairs]
    seen = {(i, m) for i, m, *_ in kf_pairs}
    rows += [
        (i, m, None, False, True, False)
        for i in range(len(seqs))
        for m in hmmsearch[i][1]
        if (i, m) not in seen
    ]
    pairs_table = {
        "protein_id": [ids[r[0]] for r in rows],
        "model": [acc[r[1]] for r in rows],
        "protein_length": [lengths[r[0]] for r in rows],
        "model_length": [model_length[r[1]] for r in rows],
        "hmmer_s": [None if r[2] is None else r[2] * scale_kfpfam for r in rows],
        "candidate": [r[3] for r in rows],
        "true": [r[4] for r in rows],
        "found": [r[5] for r in rows],
    }
    pq.write_table(pa.table(pairs_table), a.pairs, compression="zstd")

    run = {
        "proteins": len(seqs),
        "residues": int(sum(lengths)),
        "models": len(thresholds),
        "threads": a.threads,
        "forward_slack": a.forward_slack,
        "hmmsearch_s": {
            "per_protein_sum": float(raw.sum()),
            "empty_scan": empty_scan_s,
            "batch": batch_hmmsearch,
            "scale": scale_hmmsearch,
        },
        "kfpfam_hmmer_s": {
            "per_pair_sum": float(kf_hmmer.sum()),
            "empty_call_median": float(np.median(list(empty.values() or [0]))),
            "batch": batch_kfpfam,
            "scale": scale_kfpfam,
            "pairs": len(kf_pairs),
        },
        "kfpfam_filter_s": {"per_protein_sum": float(filter_s.sum()), **log},
        "true_pairs": sum(len(v[1]) for v in hmmsearch.values()),
        "kfpfam_missed": sum(len(x) for x in proteins["kf_missed"]),
        "kfpfam_extra": sum(len(x) for x in proteins["kf_extra"]),
    }
    a.run.write_text(json.dumps(run, indent=1) + "\n")
    print(json.dumps(run))


if __name__ == "__main__":
    main()
