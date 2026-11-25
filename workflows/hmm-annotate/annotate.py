"""HMMER on kfpfam's candidate (protein, model) pairs, or on all: domain rows, and labels.

Runs in the kfpfam image (pyhmmer and pyarrow). One subcommand per step of ``main.nf``:

- ``press``: an HMM file -> ``OUT.h3{m,i,f,p}``, checked against kfpfam's ``models.tsv``, and the
  file's SHA-256 (``hmm.sha256``), for lookup-index to check its models against.
- ``search``: pairs on stdin, from kfpfam's ``filter join`` run with ``FILTER_BINARY=1`` beside this
  script on the same FASTA (``protein_id`` as the name) -> ``--domains`` Parquet, ``--labels``
  Parquet (protein_id, model; proteins passing the model's thresholds), one-row ``--stats`` TSV.
  ``--more-pairs``: a file of further pairs in the same blocks, searched as well.
  With ``--all-pairs``, every protein against every model as ``hmmsearch`` does, and no stdin.
  ``--near``: also (protein_id, model) within ``--near-slack`` bits of passing, for ``members.py``.
- ``check``: 1 in N of the same proteins against every model, without kfpfam -> per model, what
  passes and how much of it ``search`` has.

A domain passes when its protein's score is at least the model's ``seq`` threshold and its own score
at least the ``dom`` one (``tables.py thresholds``). Only passing domains are written, unless
``--weak-evalue`` also asks for those with an independent E-value at most that. E-values are at a
fixed Z, so they compare across chunks.

HMMER takes at most 100,000 residues: longer proteins are not searched, and counted (``too_long``).

The pair stream's reader and the thread jobs are from the raw-reads pipeline's ``hmmer_pairs.py``.
"""

import argparse
import hashlib
import re
import resource
import struct
import sys
import threading
from array import array
from bisect import bisect_right
from collections import Counter, defaultdict, deque
from collections.abc import Iterator
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pyhmmer
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, SequenceFile
from pyhmmer.plan7 import HMMFile, HMMPressedFile, OptimizedProfile, Pipeline

ALPHABET = Alphabet.amino()
FLANK = 10  # residues kept either side of a domain: bundle.py's FLANK
HMMER_MAX = 100_000  # residues
BATCH = 2_000_000  # streamed pairs searched together
JOB = 512  # pairs per thread job, so that small models don't cost more in hand-over than in HMMER
ROW_GROUP = 65_536  # domains held before they are written
LOCAL = threading.local()
SCHEMA = pa.schema(
    [
        ("protein_id", pa.int64()),
        ("model", pa.string()),
        ("hmm_from", pa.int32()),
        ("hmm_to", pa.int32()),
        ("ali_from", pa.int32()),
        ("ali_to", pa.int32()),
        ("env_from", pa.int32()),
        ("env_to", pa.int32()),
        ("cigar", pa.string()),
        ("score", pa.float32()),
        ("seq_score", pa.float32()),
        ("i_evalue", pa.float64()),
        ("passes", pa.bool_()),
        ("left", pa.int8()),
        ("right", pa.int8()),
        ("residues", pa.string()),
    ]
)
LABELS = pa.schema([("protein_id", pa.int64()), ("model", pa.string())])


def text(value):
    return value.decode() if isinstance(value, bytes) else value


def pid(seq_or_hit) -> int:
    return int(text(seq_or_hit.name))


def cigar(hmm_line: str, target_line: str) -> str:
    """M, I and D against the HMM, from HMMER's aligned lines."""
    ops = "".join(
        "I" if h == "." else "D" if t == "-" else "M"
        for h, t in zip(hmm_line, target_line, strict=True)
    )
    return "".join(f"{len(m[0])}{m[1]}" for m in re.finditer(r"(.)\1*", ops))


def read_thresholds(path: Path) -> list[tuple[str, str, float, float]]:
    """(model, name, seq, dom) in kfpfam's model order."""
    with open(path) as f:
        rows = [line.rstrip("\n").split("\t") for line in f][1:]
    # HMMER's scores are single precision: a threshold is compared as one too
    return [(m, n, array("f", [float(s)])[0], array("f", [float(d)])[0]) for m, n, s, d in rows]


def check_names(names: list[str], expected: list[str], what: object) -> None:
    if names != expected:
        odd = sorted(set(names) ^ set(expected))
        sys.exit(f"{what}: not the models of kfpfam's prep, in its order ({odd[:5]})")


def press(a: argparse.Namespace) -> None:
    with HMMFile(a.hmm) as f:
        n = pyhmmer.hmmer.hmmpress(f, a.out)
    with open(a.models) as f:
        expected = [line.rstrip("\n").split("\t")[1] for line in f][1:]
    with HMMPressedFile(a.out) as f:
        check_names([text(p.name) for p in f], expected, a.hmm)
    with open(a.hmm, "rb") as f:
        Path("hmm.sha256").write_text(hashlib.file_digest(f, "sha256").hexdigest() + "\n")
    print(f"{n} models pressed")


def profiles(pressed: str, thresholds: list) -> list[OptimizedProfile]:
    with HMMPressedFile(pressed) as f:
        found = list(f)
    check_names([text(p.name) for p in found], [t[1] for t in thresholds], pressed)
    return found


def pipeline(z: float) -> Pipeline:
    """This thread's pipeline, cleared: building one is slow. It reports everything past HMMER's
    own filters, and passing is decided from the scores."""
    if not hasattr(LOCAL, "pipeline"):
        LOCAL.pipeline = Pipeline(ALPHABET, Z=z, domZ=z, E=float("inf"), domE=float("inf"))
    LOCAL.pipeline.clear()
    return LOCAL.pipeline


def search_model(profile, block, seqs, threshold, a: argparse.Namespace) -> tuple[list, set]:
    """One model on its candidates: a row per domain kept, and the proteins near passing.
    `seqs`: the block by name, if made."""
    model, _, seq_t, dom_t = threshold
    rows: list[tuple] = []
    near: set[tuple[int, str]] = set()
    hits = pipeline(a.z).search_hmm(profile, block)
    if not len(hits):
        return rows, near
    seqs = seqs or {text(s.name): s for s in block}
    for hit in hits:
        for d in hit.domains:
            passes = hit.score >= seq_t and d.score >= dom_t
            if a.near and hit.score >= seq_t - a.near_slack and d.score >= dom_t - a.near_slack:
                near.add((pid(hit), model))
            if not passes and not d.i_evalue <= a.weak_evalue:
                continue
            ali = d.alignment
            seq = seqs[text(hit.name)]
            left = min(FLANK, ali.target_from - 1)
            right = min(FLANK, len(seq) - ali.target_to)
            start, stop = ali.target_from - 1 - left, ali.target_to + right
            rows.append(
                (
                    pid(hit), model, ali.hmm_from, ali.hmm_to, ali.target_from, ali.target_to,
                    d.env_from, d.env_to, cigar(ali.hmm_sequence, ali.target_sequence),
                    d.score, hit.score, d.i_evalue, passes, left, right,
                    ALPHABET.decode(seq.sequence[start:stop]),
                )
            )  # fmt: skip
    return rows, near


class Output:
    """The domains (none if its path is None) as they come, and the labels held to the end."""

    def __init__(self, domains: Path | None) -> None:
        self.rows: list[tuple] = []
        self.labels: set[tuple[int, str]] = set()
        self.near: set[tuple[int, str]] = set()
        self.domains = self.passing = 0
        self.writer = domains and pq.ParquetWriter(
            domains,
            SCHEMA,
            compression="zstd",
            compression_level=9,
            use_dictionary=["model", "passes"],
        )

    def add(self, rows: list[tuple], near: set[tuple[int, str]]) -> None:
        self.near |= near
        self.domains += len(rows)
        for r in rows:
            if r[12]:
                self.passing += 1
                self.labels.add(r[:2])
        if self.writer:
            self.rows.extend(rows)
            if len(self.rows) >= ROW_GROUP:
                self.flush()

    def flush(self) -> None:
        self.writer.write_table(table(self.rows, SCHEMA))
        self.rows = []

    def close(self, labels: Path, near: Path | None) -> None:
        if self.writer:
            self.flush()
            self.writer.close()
        for found, path in (self.labels, labels), (self.near, near):
            if path:
                t = table(sorted(found), LABELS)
                pq.write_table(t, path, compression="zstd", use_dictionary=["model"])


def table(rows: list[tuple], schema: pa.Schema) -> pa.Table:
    columns = list(zip(*rows, strict=True)) or [[] for _ in schema.names]
    return pa.Table.from_arrays(
        [pa.array(c, f.type) for c, f in zip(columns, schema, strict=True)], schema=schema
    )


def read_block(stream) -> tuple[tuple[int, int], bytes] | None:
    """((proteins, residues), pairs) of the next block; None after the last."""
    data = stream.read(24)
    if not data:
        return None
    if len(data) == 24:
        n, residues, size = struct.unpack("<3Q", data)
        data = stream.read(size)
        if len(data) == size:
            return (n, residues), data
    sys.exit("the filter's output stops inside a block")


def union(streams: list[bytes]) -> bytes:
    """The pairs of any of them, once each, in model order."""
    pairs = np.concatenate([np.frombuffer(s, "<u4").reshape(-1, 2) for s in streams])
    keys = np.unique(pairs[:, 1].astype("u8") << 32 | pairs[:, 0])
    return np.stack([keys & 0xFFFFFFFF, keys >> 32], axis=1).astype("<u4").tobytes()


def read_blocks(streams: list, proteins: SequenceFile) -> Iterator[tuple]:
    """`filter join`'s blocks: their proteins, read forward from the FASTA it read, and their pairs,
    as one flat array of protein in the block, model, protein, model... in model order. Of several
    streams, the pairs in any of them."""
    while any(found := [read_block(s) for s in streams]):
        if not all(found) or len({head for head, _ in found}) > 1:
            sys.exit("the two sets of pairs are not in the same blocks")
        n, residues = found[0][0]
        pairs = array("I")
        pairs.frombytes(found[0][1] if len(found) == 1 else union([p for _, p in found]))
        block = proteins.read_block(sequences=n)
        if (len(block), sum(map(len, block))) != (n, residues):
            sys.exit("the filter's blocks are not of these proteins")
        yield block, pairs
    if proteins.read() is not None:
        sys.exit("the filter's output stops before the last protein")


def batches(blocks, size: int, counts: Counter) -> Iterator[dict[int, list]]:
    """Candidate proteins by model, of as many blocks as make `size` pairs."""
    targets, n = defaultdict(list), 0
    for block, pairs in blocks:
        counts["proteins"] += len(block)
        lengths = [len(s) for s in block]
        counts["residues"] += sum(lengths)
        long = {i for i, x in enumerate(lengths) if x > HMMER_MAX}
        counts["too_long"] += len(long)
        ids, ms = pairs[::2], pairs[1::2]
        i = 0
        while i < len(ms):
            end = bisect_right(ms, ms[i], i)
            targets[ms[i]].extend([block[k] for k in ids[i:end] if k not in long])
            i = end
        n += len(ms)
        counts["pairs"] += len(ms)
        if n >= size:
            yield targets
            targets, n = defaultdict(list), 0
    if targets:
        yield targets


def every_pair(proteins: SequenceFile, optimized: list, counts: Counter) -> Iterator[tuple]:
    """Jobs of every model on every protein, as one block they share: hmmsearch."""
    block = proteins.read_block()
    counts["proteins"], counts["residues"] = len(block), sum(map(len, block))
    kept = DigitalSequenceBlock(ALPHABET, [s for s in block if len(s) <= HMMER_MAX])
    counts["too_long"] = len(block) - len(kept)
    counts["pairs"] = len(kept) * len(optimized)
    seqs = {text(s.name): s for s in kept}
    if seqs:
        yield from ((m, profile, kept, seqs) for m, profile in enumerate(optimized))


def run(jobs, thresholds, pool, output: Output, a: argparse.Namespace) -> None:
    """Every (model number, profile, block of proteins, the block by name or None) of jobs, written
    in their order, so output does not depend on threads. A model may be in jobs once: a profile is
    reconfigured for each target, so two threads can't share one."""

    def work(group):
        return [search_model(p, b, seqs, thresholds[m], a) for m, p, b, seqs in group]

    def collect(future) -> None:
        for found in future.result():
            output.add(*found)

    pending, group, n = deque(), [], 0
    for job in jobs:
        group.append(job)
        n += len(job[2])
        if n >= JOB:
            pending.append(pool.submit(work, group))
            group, n = [], 0
            if len(pending) >= 4 * a.threads:
                collect(pending.popleft())
    if group:
        pending.append(pool.submit(work, group))
    while pending:
        collect(pending.popleft())


def open_proteins(path: Path) -> SequenceFile:
    return SequenceFile(path, format="fasta", digital=True, alphabet=ALPHABET)


def search(a: argparse.Namespace) -> None:
    thresholds = read_thresholds(a.thresholds)
    optimized = profiles(a.pressed, thresholds)
    counts = Counter()
    output = Output(a.domains)
    with open_proteins(a.faa) as proteins, ThreadPoolExecutor(a.threads) as pool:
        if a.all_pairs:
            run(every_pair(proteins, optimized, counts), thresholds, pool, output, a)
        else:
            streams = [sys.stdin.buffer]
            if a.more_pairs:
                streams.append(open(a.more_pairs, "rb"))  # noqa: SIM115, read to the end
            blocks = read_blocks(streams, proteins)
            for targets in batches(blocks, BATCH, counts):
                # built here, not in the workers: across 6 threads blocks ran about 40x slower
                jobs = (
                    (m, optimized[m], DigitalSequenceBlock(ALPHABET, targets[m]), None)
                    for m in sorted(targets)
                    if targets[m]
                )
                run(jobs, thresholds, pool, output, a)
    output.close(a.labels, a.near)
    used = resource.getrusage(resource.RUSAGE_SELF)
    row = {
        **{k: counts[k] for k in ("proteins", "residues", "too_long", "pairs")},
        "domains": output.domains,
        "domains_passing": output.passing,
        "proteins_passing": len({p for p, _ in output.labels}),
        "hmmer_core_s": round(used.ru_utime + used.ru_stime),
    }
    a.stats.write_text("\t".join(row) + "\n" + "\t".join(map(str, row.values())) + "\n")


def check(a: argparse.Namespace) -> None:
    """Columns: model, passing (protein, model) pairs and domains by hmmsearch, how many of each
    `search` has, and how many it has that hmmsearch doesn't (extra: none expected)."""
    thresholds = {n: (m, s, d) for m, n, s, d in read_thresholds(a.thresholds)}
    with open_proteins(a.faa) as f:
        kept = [s for s in f if pid(s) // a.sample % a.check_sample == 0 and len(s) <= HMMER_MAX]
    truth: set[tuple[int, str, int, int]] = set()
    with HMMFile(a.hmm) as hmms:
        for hits in pyhmmer.hmmsearch(
            hmms, DigitalSequenceBlock(ALPHABET, kept), cpus=a.threads, Z=a.z, domZ=a.z,
            E=float("inf"), domE=float("inf"),
        ):  # fmt: skip
            model, seq_t, dom_t = thresholds[text(hits.query.name)]
            for hit in hits:
                if hit.score >= seq_t:
                    truth.update(
                        (pid(hit), model, d.env_from, d.env_to)
                        for d in hit.domains
                        if d.score >= dom_t
                    )
    ids = {pid(s) for s in kept}
    got = pq.read_table(a.domains, columns=["protein_id", "model", "env_from", "env_to", "passes"])
    found = {r[:4] for r in zip(*got.to_pydict().values(), strict=True) if r[4] and r[0] in ids}
    rows: dict[str, Counter] = defaultdict(Counter)
    for name, domains in (
        ("truth", truth), ("found", truth & found), ("extra", found - truth)
    ):  # fmt: skip
        for x in domains:
            rows[x[1]][f"domains_{name}"] += 1
        pairs = {x[:2] for x in domains}
        if name == "extra":
            pairs -= {x[:2] for x in truth}
        for x in pairs:
            rows[x[1]][f"pairs_{name}"] += 1
    columns = [f"{k}_{n}" for k in ("pairs", "domains") for n in ("truth", "found", "extra")]
    with open(a.out, "w") as out:
        out.write("\t".join(["model", *columns]) + "\n")
        for model in sorted(rows):
            out.write("\t".join([model, *(str(rows[model][c]) for c in columns)]) + "\n")
    print(f"{len(kept)} proteins: {len(truth)} passing domains, {len(truth & found)} found")


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = p.add_subparsers(required=True)
    x = sub.add_parser("press")
    x.set_defaults(run=press)
    x.add_argument("--hmm", required=True)
    x.add_argument("--models", type=Path, required=True, help="models.tsv of kfpfam's prep")
    x.add_argument("--out", required=True, help="the pressed files' prefix")
    x = sub.add_parser("search")
    x.set_defaults(run=search)
    x.add_argument("--pressed", required=True, help="press's --out")
    x.add_argument("--domains", type=Path, help="not written if this is left out")
    x.add_argument("--labels", type=Path, required=True)
    x.add_argument("--stats", type=Path, required=True)
    x.add_argument("--near", type=Path)
    x.add_argument("--near-slack", type=float, default=8)
    x.add_argument("--more-pairs", type=Path, help="as members.py fasta writes them")
    x.add_argument("--all-pairs", action="store_true", help="no prefilter: hmmsearch")
    x.add_argument("--weak-evalue", type=float, default=0, help="0: passing domains only")
    c = sub.add_parser("check")
    c.set_defaults(run=check)
    c.add_argument("--hmm", required=True, help="the HMM file itself, or its .h3m")
    c.add_argument("--domains", type=Path, required=True)
    c.add_argument("--sample", type=int, default=1, help="plan's --sample")
    c.add_argument("--check-sample", type=int, default=1, help="1 in this many of the proteins")
    c.add_argument("--out", type=Path, required=True)
    for s in (x, c):
        s.add_argument("--faa", type=Path, required=True)
        s.add_argument("--thresholds", type=Path, required=True)
        s.add_argument("--threads", type=int, default=1)
        s.add_argument("--z", type=float, default=1e6, help="HMMER's Z and domZ")
    a = p.parse_args()
    a.run(a)


if __name__ == "__main__":
    main()
