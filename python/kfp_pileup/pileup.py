"""Read posteriors and crude pileups (plan, phase 13, steps 4-6).

Input: alignments of reads to the per-sample reference (:mod:`kfp_pileup.reference`),
several per read: DIAMOND blastx tabular output with ``DIAMOND_FIELDS`` (one file per mate),
or bwa-mem2 SAM (``mem -a``, mates in one file). Alignments of one read must be adjacent,
as both tools write them.

Each alignment h of a read has a score S_h in bits (DIAMOND's bit score; bwa's ``AS`` times
``score_nats`` / ln 2) and, per aligned reference column, a term t_h,q in the same units:
the substitution score of that column (BLOSUM62 or bwa's match/mismatch), keyed by the read
position q it uses.

- *Read posterior:* w_h = 2^S_h / sum_g 2^S_g over the read's alignments, neighbours
  included. A read whose best alignment is a neighbour is dropped from the pileups.
- *Leave one column out* (the circularity fix): at a target column using read position q
  the weight is 2^(S_h - t_h,q) / sum_g 2^(S_g - t_g,q), every alignment's term at q
  removed, so a read's residue there does not lower its own weight. Deletions (no read
  position) take w_h.
- *Pileups:* per (reference sequence, column, residue) the summed weights. Per column:
  depth, consensus, the top residue other than the reference's with its frequency, and
  ``sibling_share``: the part of that residue's weight from reads whose alignment to another
  unit (target or neighbour) carries the same residue at the same read position (the
  structured-contamination case).
- *Units* (the input of step 6's ``aai`` check): reads (posterior sum), the reads'
  alignment identity over all of the unit's sequences and over its representative, and the
  representative's consensus identity over its covered columns, beside the profile's
  ``aai`` (``ani``) columns.

The crude pass by design: mates are scored independently, insertions count in scores but
not in pileups, and there is no sequencing-error model. :mod:`kfp_pileup.refine` (steps
8-9) runs after it on the alignments kept here.
"""

import math
import re
import tempfile
from collections.abc import Iterator, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Final

import numpy as np
import polars as pl

from kfp_pileup.nucleotide import nucleotide_pileup
from kfp_pileup.reference import REF_FASTA, REFS
from kfp_pileup.refine import NOVEL, refine
from kfp_pileup.spill import PARTS, Spill, partitions
from kmer_functional_profiler._core import DiamondReader
from kmer_functional_profiler.index import read_fasta
from kmer_functional_profiler.query import Timer

LN2: Final = math.log(2)
DIAMOND_FIELDS: Final = "qseqid sseqid bitscore qstart qend sstart send qframe btop"
DIAMOND_BITS: Final = 0.267 / LN2  # bits per BLOSUM62 point (lambda at gaps 11/1, DIAMOND's)
# Nats per bwa score point: a mismatch costs 5 points against a match, and at a 1% error rate
# its log-likelihood ratio is ln(0.99 / (0.01 / 3)) = 5.7 nats.
SCORE_NATS: Final = 1.1
BWA_MATCH, BWA_MISMATCH, BWA_N, GAP_EXTEND = 1, -4, -1, -1  # bwa-mem2's defaults
MIN_DEPTH: Final = 3.0
MIN_FREQ: Final = 0.2
BATCH_READS: Final = 20_000
BYTES: Final = list(range(256))
CHARS: Final = [chr(b) for b in BYTES]
DIAMOND_CHUNK: Final = 1 << 24  # bytes of DIAMOND output read at a time
DROPPED: Final = ("top", "margin", "unit", "floor_reads", "floor_lines")


@dataclass(frozen=True)
class HitFilter:
    """Which of a read's DIAMOND alignments ``call`` keeps, against the read's best (highest
    bit score), decided over all its lines before any is expanded; a line goes by the first
    rule it fails. ``top``: percent below the best (DIAMOND's ``--top``, re-applied because
    routed runs apply it per block); ``max_bits_below``: bits below the best (a ~3-bit
    mismatch costs the EM ~10^3 in likelihood, so 10 bits keeps every alignment it could
    weigh); ``unit_bits_below``: bits below the read's best alignment to the same unit (its
    other members, which the EM does not use and the pileup weighs at 2^-bits);
    ``min_identity``: the best alignment's identity, below which the read is dropped (a
    distant homologue, which the EM sends to the novel class). ``inf`` / 0 keep all.
    ponytail: defaults chosen, not yet checked against spike-protein's truth."""

    top: float = 10.0
    max_bits_below: float = 10.0
    unit_bits_below: float = 5.0
    min_identity: float = 0.7


KEEP_ALL: Final = HitFilter(math.inf, math.inf, math.inf, 0.0)
DEFAULT_FILTER: Final = HitFilter()
ROLES: Final = ("target", "neighbour")
TARGET, NEIGHBOUR = range(2)
GAP: Final = ord("-")
MATE: Final = ["read", "mate"]  # one read's alignments: mates are scored independently

_BLOSUM62: Final = """
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
"""


def blosum62() -> np.ndarray:
    """BLOSUM62 as a 256 x 256 lookup on ASCII bytes (either case); other bytes score as
    ``X``, and a gap (``-``) against anything as one gap extension."""
    lines = _BLOSUM62.strip().splitlines()
    letters = lines[0].split()
    table = np.full((256, 256), -1.0)
    for line in lines[1:]:
        a, *scores = line.split()
        for b, s in zip(letters, scores, strict=True):
            for x in (ord(a), ord(a.lower())):
                for y in (ord(b), ord(b.lower())):
                    table[x, y] = float(s)
    table[GAP, :] = table[:, GAP] = GAP_EXTEND
    return table


def dna_scores() -> np.ndarray:
    """bwa-mem2's scores as a 256 x 256 lookup: match, mismatch, N, and a gap extension."""
    table = np.full((256, 256), float(BWA_MISMATCH))
    for b in b"ACGTacgt":
        table[b, b] = BWA_MATCH
        table[b, b ^ 0x20] = BWA_MATCH  # the other case
    for n in b"Nn":
        table[n, :] = table[:, n] = BWA_N
    table[GAP, :] = table[:, GAP] = GAP_EXTEND
    return table


class _Batch:
    """Alignments of a batch of reads, as columns: one row per alignment, and one per
    aligned reference column (``hit`` its alignment)."""

    def __init__(self) -> None:
        self.names: list[str] = []
        self.read: list[int] = []
        self.mate: list[int] = []
        self.seq: list[int] = []
        self.score: list[float] = []
        self.start: list[int] = []  # first reference column
        self.end: list[int] = []  # past the last reference column
        self.reverse: list[bool] = []
        self.ins: list[int] = []  # read residues inserted (no reference column)
        self.shifts: list[int] = []  # frameshifts (DIAMOND -F): one difference each
        self.hit: list[int] = []
        self.col: list[int] = []
        self.qpos: list[int] = []
        self.rres = bytearray()
        self.sres = bytearray()
        # events without a reference column, for the nucleotide pileup: (hit, the column
        # they precede, kind, read position key): b"I" an inserted residue (its codon),
        # b"\\" an inserted base, b"/" a deleted base (the column's first base)
        self.ev_hit: list[int] = []
        self.ev_col: list[int] = []
        self.ev_kind = bytearray()
        self.ev_key: list[int] = []

    def start_read(self, name: str) -> int:
        self.names.append(name)
        return len(self.names) - 1

    def start_hit(
        self, read: int, mate: int, seq: int, score: float, start: int, reverse: bool
    ) -> int:
        self.read.append(read)
        self.mate.append(mate)
        self.seq.append(seq)
        self.score.append(score)
        self.start.append(start)
        self.end.append(start)
        self.reverse.append(reverse)
        self.ins.append(0)
        self.shifts.append(0)
        return len(self.read) - 1

    def columns(self, hit: int, cols: range, qpos: Sequence[int], rres: bytes, sres: bytes) -> None:
        self.hit.extend([hit] * len(cols))
        self.col.extend(cols)
        self.qpos.extend(qpos)
        self.rres += rres
        self.sres += sres

    def event(self, hit: int, col: int, kind: bytes, key: int) -> None:
        self.ev_hit.append(hit)
        self.ev_col.append(col)
        self.ev_kind += kind
        self.ev_key.append(key)


class _Columns:
    """A batch from :class:`~kmer_functional_profiler._core.DiamondReader`: the columns
    :class:`_Batch` keeps, as arrays."""

    names: list[str]
    read: np.ndarray
    seq: np.ndarray
    score: np.ndarray
    start: np.ndarray
    end: np.ndarray
    reverse: np.ndarray
    ins: np.ndarray
    shifts: np.ndarray
    hit: np.ndarray
    col: np.ndarray
    qpos: np.ndarray
    rres: np.ndarray
    sres: np.ndarray
    ev_hit: np.ndarray
    ev_col: np.ndarray
    ev_kind: np.ndarray
    ev_key: np.ndarray

    dropped_top: int
    dropped_margin: int
    dropped_unit: int
    dropped_floor_reads: int
    dropped_floor_lines: int

    def __init__(self, cols: dict[str, Any], mate: int) -> None:
        for name, value in cols.items():
            setattr(self, name, value)
        self.mate = np.full(len(self.read), mate, np.uint8)


def read_diamond(
    paths: Sequence[str | Path],
    names: Sequence[str],
    refs: list[bytes],
    units: np.ndarray,
    batch_reads: int,
    hit_filter: HitFilter = KEEP_ALL,
) -> Iterator[_Columns]:
    """Batches of DIAMOND alignments (``paths`` the mates' outputs, in mate order) against
    sequences ``names`` with residues ``refs`` of ``units``, each of whole reads (a read's
    alignments are consecutive) and their lines ``hit_filter`` keeps, parsed in Rust from
    chunks of ``DIAMOND_CHUNK`` bytes."""
    reader = DiamondReader(
        [n.encode() for n in names], refs, units.astype(np.uint32), **asdict(hit_filter)
    )
    for mate, path in enumerate(paths):
        with open(path, "rb") as f:
            buf, eof = b"", False
            while True:
                if not eof:
                    more = f.read(DIAMOND_CHUNK)
                    eof, buf = not more, buf + more
                cols, used = reader.parse(buf, batch_reads, eof)
                buf = buf[used:]
                if cols["names"]:
                    yield _Columns(cols, mate)
                elif eof:
                    break
                # else a read longer than the buffer: read more


_CIGAR: Final = re.compile(r"(\d+)([MIDNSHP=X])")
_COMPLEMENT: Final = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")


def _sam_hit(
    batch: _Batch,
    hit: int,
    pos: int,
    cigar: str,
    seq: bytes,
    reverse: bool,
    borrowed: bool,
    ref: bytes,
) -> None:
    """Columns of one SAM alignment; read positions are counted on the original read."""
    length = len(seq)
    col, i = pos - 1, 0
    for num, op in _CIGAR.findall(cigar):
        n = int(num)
        if op in "M=X":
            qpos = range(length - 1 - i, length - 1 - i - n, -1) if reverse else range(i, i + n)
            batch.columns(hit, range(col, col + n), qpos, seq[i : i + n], ref[col : col + n])
            col, i = col + n, i + n
        elif op in "IS" or (op == "H" and borrowed):
            i += n
            batch.ins[hit] += n if op == "I" else 0
        elif op in "DN":
            batch.columns(hit, range(col, col + n), [-1] * n, b"-" * n, ref[col : col + n])
            col += n
    batch.end[hit] = col


def read_sam(
    path: str | Path,
    seq_index: dict[str, int],
    refs: list[bytes],
    batch_reads: int,
    score_bits: float,
) -> Iterator[_Batch]:
    """Batches of bwa-mem2 SAM alignments (mates in one file). A record without its sequence
    (``*``, as secondary alignments may be) borrows the primary's, turned to its strand."""
    batch, reads = _Batch(), 0
    last, read = None, -1
    primary: dict[int, tuple[bytes, bool]] = {}
    with open(path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            flag = int(fields[1])
            if fields[0] != last:
                if reads >= batch_reads:
                    yield batch
                    batch, reads = _Batch(), 0
                last, read, reads, primary = fields[0], batch.start_read(fields[0]), reads + 1, {}
            if flag & 4 or fields[2] == "*":
                continue
            mate, reverse = int(bool(flag & 128)), bool(flag & 16)
            seq, borrowed = fields[9].encode(), False
            if seq == b"*":
                if mate not in primary:
                    continue  # ponytail: no primary before it; bwa writes the primary first
                seq, strand = primary[mate]
                if strand != reverse:
                    seq = seq.translate(_COMPLEMENT)[::-1]
                borrowed = True
            elif not flag & 0x900:
                primary[mate] = (seq, reverse)
            score = next((int(t[5:]) for t in fields[11:] if t.startswith("AS:i:")), 0)
            ref_seq = seq_index[fields[2]]
            hit = batch.start_hit(
                read, mate, ref_seq, score * score_bits, int(fields[3]) - 1, reverse
            )
            _sam_hit(batch, hit, int(fields[3]), fields[5], seq, reverse, borrowed, refs[ref_seq])
    if reads:
        yield batch


def score_batch(
    batch: _Batch | _Columns, refs: pl.DataFrame, matrix: np.ndarray, term_bits: float
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """(assignments, pileup part, per-sequence read stats, target columns, target events)
    of one batch. Assignments keep each alignment's ``hit`` id, reference span, strand,
    read residues (``qlen``) and differences (``mism``, a frameshift one) for
    :mod:`kfp_pileup.refine`; target columns are (hit, seq, col, qpos, residue byte) of
    alignments to targets, and target events (hit, seq, col, kind byte, qpos) their
    insertions and frameshifts."""
    hits = (
        pl.DataFrame(
            {
                "read": batch.read,
                "mate": batch.mate,
                "seq": batch.seq,
                "score": batch.score,
                "start": batch.start,
                "end": batch.end,
                "reverse": batch.reverse,
                "ins": batch.ins,
                "shifts": batch.shifts,
            },
            schema={
                "read": pl.UInt32,
                "mate": pl.UInt8,
                "seq": pl.UInt32,
                "score": pl.Float64,
                "start": pl.Int32,
                "end": pl.Int32,
                "reverse": pl.Boolean,
                "ins": pl.Int32,
                "shifts": pl.Int32,
            },
        )
        .with_row_index("hit")
        .join(refs.select("seq", "role", "unit"), on="seq")
        .with_columns(e=(2.0 ** (pl.col("score") - pl.col("score").max().over(MATE))))
        .with_columns(
            total=pl.col("e").sum().over(MATE),
            kept=(pl.col("role").sort_by("score", descending=True).first().over(MATE) == TARGET),
        )
        .with_columns(posterior=pl.col("e") / pl.col("total"))
    )
    rres, sres = (
        np.frombuffer(bytes(batch.rres), np.uint8),
        np.frombuffer(bytes(batch.sres), np.uint8),
    )
    rows = pl.DataFrame(
        {
            "hit": np.asarray(batch.hit, np.uint32),
            "col": np.asarray(batch.col, np.int32),
            "qpos": np.asarray(batch.qpos, np.int64),
            "residue": rres,  # bytes; the pileup part's residues become strings
            "ref": sres,
            "term": matrix[rres, sres] * term_bits,
        }
    ).join(
        hits.select(
            "hit", "read", "mate", "seq", "role", "unit", "e", "total", "kept", "posterior"
        ),
        on="hit",
    )
    # Leave one column out: the denominator with every alignment's term at the read position
    # removed, sum_g e_g 2^-t_g,q = total + sum over alignments at q of e_g (2^-t_g,q - 1).
    at_q = (
        rows.filter(pl.col("qpos") >= 0)
        .group_by(*MATE, "qpos")
        .agg(adjust=(pl.col("e") * (2.0 ** -pl.col("term") - 1)).sum())
    )
    rows = rows.join(at_q, on=[*MATE, "qpos"], how="left").with_columns(
        weight=pl.when(pl.col("qpos") >= 0)
        .then(pl.col("e") * 2.0 ** -pl.col("term") / (pl.col("total") + pl.col("adjust")))
        .otherwise(pl.col("posterior"))
    )
    piled = rows.filter(pl.col("kept"), pl.col("role") == TARGET)
    others = rows.select(*MATE, "qpos", residue="ref", other="unit")
    sibling = (
        piled.filter(pl.col("residue") != pl.col("ref"), pl.col("qpos") >= 0)
        .join(others, on=[*MATE, "qpos", "residue"])
        .filter(pl.col("other") != pl.col("unit"))
        .select("hit", "col", sibling=pl.lit(True))
        .unique()
    )
    part = (
        piled.join(sibling, on=["hit", "col"], how="left")
        .group_by("seq", "col", "residue")
        .agg(
            weight=pl.col("weight").sum(),
            sibling=pl.col("weight").filter(pl.col("sibling").fill_null(False)).sum(),
        )
        .with_columns(pl.col("residue").replace_strict(BYTES, CHARS, return_dtype=pl.String))
    )
    identity = piled.group_by("hit").agg(identity=(pl.col("residue") == pl.col("ref")).mean())
    stats = (
        hits.filter(pl.col("kept"), pl.col("role") == TARGET)
        .join(identity, on="hit")
        .group_by("seq")
        .agg(
            reads=pl.col("posterior").sum(),
            matched=(pl.col("posterior") * pl.col("identity")).sum(),
        )
    )
    names = pl.Series("read_name", batch.names)
    # read residues and differences of each alignment (insertions count as both)
    counts = rows.group_by("hit").agg(
        qlen=(pl.col("qpos") >= 0).sum().cast(pl.Int32),
        mism=(pl.col("residue") != pl.col("ref")).sum().cast(pl.Int32),
    )
    assignments = hits.join(counts, on="hit", how="left").select(
        "hit",
        read_name=names.gather(hits["read"]),
        mate="mate",
        seq="seq",
        unit="unit",
        role=pl.col("role").replace_strict(dict(enumerate(ROLES))),
        score="score",
        posterior="posterior",
        kept="kept",
        start="start",
        end="end",
        reverse="reverse",
        qlen=pl.col("qlen").fill_null(0) + pl.col("ins"),
        mism=pl.col("mism").fill_null(0) + pl.col("ins") + pl.col("shifts"),
    )
    target_cols = rows.filter(pl.col("role") == TARGET).select(
        "hit", "seq", "col", "qpos", "residue"
    )
    events = pl.DataFrame(
        {
            "hit": np.asarray(batch.ev_hit, np.uint32),
            "col": np.asarray(batch.ev_col, np.int32),
            "kind": np.frombuffer(bytes(batch.ev_kind), np.uint8),
            "qpos": np.asarray(batch.ev_key, np.int64),
        }
    ).join(hits.filter(pl.col("role") == TARGET).select("hit", "seq"), on="hit")
    return assignments, part, stats, target_cols, events


def call(
    alignments: Sequence[str | Path],
    ref_dir: str | Path,
    out_dir: str | Path,
    aligner: str,
    profile: str | Path | None = None,
    min_depth: float = MIN_DEPTH,
    min_freq: float = MIN_FREQ,
    score_nats: float = SCORE_NATS,
    batch_reads: int = BATCH_READS,
    novel: float = NOVEL,
    reads: Sequence[str | Path] = (),
    stats: str | Path | None = None,
    spill_parts: int = PARTS,
    threads: int = 1,
    hit_filter: HitFilter = DEFAULT_FILTER,
) -> dict[str, int]:
    """Read posteriors, pileups and the unit table from ``alignments`` (``aligner``
    ``diamond``: one tabular file per mate; ``bwa``: one SAM) against the reference in
    ``ref_dir``. Writes ``assignments.parquet``, ``pileup.parquet``, ``variants.tsv`` and
    ``units.tsv`` to ``out_dir`` (the crude pass), then :func:`kfp_pileup.refine.refine`'s
    outputs (steps 8-9); with DIAMOND alignments and the ``reads`` (FASTQ, in mate order),
    :func:`kfp_pileup.nucleotide.nucleotide_pileup`'s (step 6b). With ``stats``, each
    stage's time and memory go there as the query's (``Timer``); ``threads`` worker
    processes fit step 9's haplotypes. DIAMOND alignments are kept by ``hit_filter``
    (:class:`HitFilter`; the lines each rule drops are counted in the stats). Returns counts.

    Target columns and events are spilled to Parquet by reference sequence in a temporary
    directory under ``out_dir`` (:mod:`kfp_pileup.spill`), so steps 6b and 8-9 hold one
    partition of them at a time; the alignments table stays in memory."""
    ref_dir, out = Path(ref_dir), Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    timer = Timer(Path(stats) if stats else None, log=stats is not None)
    with timer("read"):
        refs, sequences, prof = _read_inputs(ref_dir, profile)
    with tempfile.TemporaryDirectory(prefix=".spill-", dir=out) as tmp:
        target_len = np.where(refs["role"] == TARGET, refs["length"], 0)
        spill = Spill(Path(tmp), partitions(target_len, spill_parts))
        with timer("batches"):
            hits, parts, stat_parts = _score(
                alignments, refs, sequences, aligner, score_nats, batch_reads, spill,
                hit_filter, timer.counts,
            )  # fmt: skip
            spill.close()
        timer.counts.update(
            alignments=0 if hits is None else hits.height,
            spill_rows=spill.rows,
            spill_parts=len(spill.files),
            spill_bytes=sum(f.stat().st_size for f in spill.files),
        )
        refined, em_units = None, None
        if hits is not None and hits.height:
            with timer("refine"):
                refined, em_units = refine(
                    hits, spill, refs, sequences, out, aligner == "diamond", prof,
                    min_depth, min_freq, novel, timer, threads,
                )  # fmt: skip
            hits = hits.join(refined, on="hit", how="left")
        nt_counts = {}
        if reads and aligner == "diamond" and hits is not None and hits.height:
            with timer("nucleotide"):
                nt_counts = nucleotide_pileup(hits, spill, refs, reads, out)
    with timer("write"):
        counts = _write(
            hits, parts, stat_parts, refs, sequences, prof, em_units, out, min_depth, min_freq
        )
    counts = {**counts, **nt_counts}
    timer.counts.update(counts)
    timer.write()
    return counts


def _read_inputs(
    ref_dir: Path, profile: str | Path | None
) -> tuple[pl.DataFrame, dict[str, str], pl.DataFrame | None]:
    """The reference table and sequences, and the profile's columns ``call`` uses on the
    reference's units only (a real profile holds every hit unit, ~100 columns)."""
    refs = (
        pl.read_csv(
            ref_dir / REFS,
            separator="\t",
            schema_overrides={"cluster_rep": pl.String, "protein_id": pl.String},
        )
        .with_row_index("seq")
        .with_columns(pl.col("unit").cast(pl.UInt32))
        .with_columns(
            pl.col("role").replace_strict(
                {r: i for i, r in enumerate(ROLES)}, return_dtype=pl.UInt8
            )
        )
    )
    sequences = dict(read_fasta(ref_dir / REF_FASTA))
    prof = None
    if profile is not None:
        with open(profile) as f:
            header = f.readline().rstrip("\n").split("\t")
        wanted = [
            c
            for c in header
            if c in ("present_prob", "coverage_zi", "len_mean", "n_members")
            or c.startswith(("aai", "ani"))
        ]
        prof = (
            pl.read_csv(
                profile, separator="\t", columns=["unit", *wanted], infer_schema_length=None
            )
            .with_columns(pl.col("unit").cast(pl.UInt32))
            .filter(pl.col("unit").is_in(refs["unit"].unique().implode()))
        )
    return refs, sequences, prof


def _score(
    alignments: Sequence[str | Path],
    refs: pl.DataFrame,
    sequences: dict[str, str],
    aligner: str,
    score_nats: float,
    batch_reads: int,
    spill: Spill,
    hit_filter: HitFilter,
    dropped: dict[str, int],
) -> tuple[pl.DataFrame | None, list[pl.DataFrame], list[pl.DataFrame]]:
    """The crude pass, batch by batch: alignments with hit ids, the summed pileup's parts
    and per-sequence read stats; target columns and events go to ``spill``; DIAMOND lines
    ``hit_filter`` drops counted in ``dropped``."""
    seq_index = {s: i for i, s in enumerate(refs["seq_id"])}
    ref_bytes = [sequences[s].encode() for s in refs["seq_id"]]
    batches: Iterator[_Batch] | Iterator[_Columns]
    if aligner == "diamond":
        batches = read_diamond(
            alignments, refs["seq_id"].to_list(), ref_bytes, refs["unit"].to_numpy(),
            batch_reads, hit_filter,
        )  # fmt: skip
        matrix, term_bits = blosum62(), DIAMOND_BITS
    elif aligner == "bwa":
        if len(alignments) != 1:
            raise ValueError("bwa: one SAM file with both mates")
        batches = read_sam(alignments[0], seq_index, ref_bytes, batch_reads, score_nats / LN2)
        matrix, term_bits = dna_scores(), score_nats / LN2
    else:
        raise ValueError(f"aligner must be diamond or bwa, not {aligner}")
    parts, stats, n_assigned = [], [], 0
    assigned_parts = []
    for batch in batches:
        if isinstance(batch, _Columns):
            for k in DROPPED:
                dropped[f"dropped_{k}"] = dropped.get(f"dropped_{k}", 0) + getattr(
                    batch, f"dropped_{k}"
                )
        assignments, part, stat, cols, events = score_batch(batch, refs, matrix, term_bits)
        offset = pl.lit(n_assigned, pl.UInt32)
        assigned_parts.append(assignments.with_columns(pl.col("hit") + offset))
        spill.write(
            pl.concat(
                [
                    cols.with_columns(kind=pl.lit(0, pl.UInt8)),
                    events.with_columns(residue=pl.lit(0, pl.UInt8)),
                ],
                how="diagonal_relaxed",
            ).with_columns(pl.col("hit") + offset)
        )
        n_assigned += assignments.height
        parts.append(part)
        stats.append(stat)
        if len(parts) > 32:  # keep the summed pileup small
            parts = [_sum(parts)]
    hits = pl.concat(assigned_parts) if assigned_parts else None
    return hits, parts, stats


def _write(
    hits: pl.DataFrame | None,
    parts: list[pl.DataFrame],
    stats: list[pl.DataFrame],
    refs: pl.DataFrame,
    sequences: dict[str, str],
    prof: pl.DataFrame | None,
    em_units: pl.DataFrame | None,
    out: Path,
    min_depth: float,
    min_freq: float,
) -> dict[str, int]:
    """``assignments.parquet``, ``pileup.parquet``, ``variants.tsv`` and ``units.tsv``."""
    (
        hits.join(refs.select("seq", "seq_id"), on="seq").drop(
            "seq", "hit", "start", "end", "reverse", "qlen", "mism"
        )
        if hits is not None
        else pl.DataFrame()
    ).write_parquet(out / "assignments.parquet")
    pile = (
        _sum(parts)
        if parts
        else pl.DataFrame(
            schema={
                "seq": pl.UInt32,
                "col": pl.Int32,
                "residue": pl.String,
                "weight": pl.Float64,
                "sibling": pl.Float64,
            }
        )
    )
    pile = pile.join(refs.select("seq", "seq_id"), on="seq")
    pile.drop("seq").sort("seq_id", "col", "residue").write_parquet(out / "pileup.parquet")
    columns = _columns(pile, refs, sequences)
    variants = columns.filter(pl.col("depth") >= min_depth, pl.col("alt_freq") >= min_freq)
    variants.write_csv(out / "variants.tsv", separator="\t")
    units = _units(pl.concat(stats) if stats else None, columns, refs, prof)
    if em_units is not None:
        units = units.join(em_units, on="unit", how="left")
    units.write_csv(out / "units.tsv", separator="\t")
    return {
        "alignments": 0 if hits is None else hits.height,
        "columns": columns.height,
        "variants": variants.height,
        "units": units.height,
    }


def _sum(parts: list[pl.DataFrame]) -> pl.DataFrame:
    return pl.concat(parts).group_by("seq", "col", "residue").agg(pl.col("weight", "sibling").sum())


def _columns(pile: pl.DataFrame, refs: pl.DataFrame, sequences: dict[str, str]) -> pl.DataFrame:
    """Per piled column: depth, consensus, and the top non-reference residue."""
    seqs = pl.DataFrame({"seq_id": list(sequences), "sequence": list(sequences.values())})
    pile = (
        pile.join(refs.select("seq", "unit", "cluster_rep", "protein_id", "rep"), on="seq")
        .join(seqs, on="seq_id")
        .with_columns(ref=pl.col("sequence").str.slice(pl.col("col"), 1))
        .drop("sequence")
    )
    keys = ["seq_id", "unit", "cluster_rep", "protein_id", "rep", "col", "ref"]
    per_column = pile.group_by(keys).agg(
        depth=pl.col("weight").sum(),
        consensus=pl.col("residue").sort_by("weight", "residue", descending=[True, False]).first(),
        consensus_freq=pl.col("weight").max() / pl.col("weight").sum(),
    )
    alt = (
        pile.filter(pl.col("residue") != pl.col("ref"))
        .sort("weight", "residue", descending=[True, False])  # ties: the first residue
        .group_by(["seq_id", "col"])
        .agg(
            alt=pl.col("residue").first(),
            alt_weight=pl.col("weight").first(),
            sibling_share=(pl.col("sibling").first() / pl.col("weight").first()),
        )
    )
    return (
        per_column.join(alt, on=["seq_id", "col"], how="left")
        .with_columns(
            alt_freq=(pl.col("alt_weight") / pl.col("depth")).fill_null(0.0),
            pos=pl.col("col") + 1,
        )
        .drop("alt_weight", "col")
        .sort("seq_id", "pos")
    )


def _units(
    stats: pl.DataFrame | None,
    columns: pl.DataFrame,
    refs: pl.DataFrame,
    prof: pl.DataFrame | None,
) -> pl.DataFrame:
    """Per target unit: reads, alignment identity (all sequences, representative) and the
    representative's consensus identity; with the profile's columns for the aai check."""
    targets = refs.filter(pl.col("role") == TARGET)
    per_seq = (
        targets.select("seq", "unit", "cluster_rep", "rep")
        .join(
            stats
            if stats is not None
            else pl.DataFrame(
                schema={"seq": pl.UInt32, "reads": pl.Float64, "matched": pl.Float64}
            ),
            on="seq",
            how="left",
        )
        .with_columns(pl.col("reads", "matched").fill_null(0.0))
    )
    units = per_seq.group_by("unit", "cluster_rep").agg(
        reads=pl.col("reads").sum(),
        identity=pl.col("matched").sum() / pl.col("reads").sum(),
        reads_rep=pl.col("reads").filter("rep").sum(),
        identity_rep=pl.col("matched").filter("rep").sum() / pl.col("reads").filter("rep").sum(),
    )
    rep_columns = (
        columns.filter(pl.col("rep"), pl.col("depth") > 0)
        .group_by("unit")
        .agg(
            rep_covered=pl.len(),
            consensus_identity_rep=(pl.col("consensus") == pl.col("ref")).mean(),
        )
    )
    units = units.join(rep_columns, on="unit", how="left").with_columns(
        pl.col("rep_covered").fill_null(0)
    )
    if prof is not None:
        units = units.join(prof, on="unit", how="left")
    return units.with_columns(pl.col("identity", "identity_rep").fill_nan(None)).sort("unit")
