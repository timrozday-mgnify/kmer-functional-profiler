"""Read posteriors and crude pileups (plan, phase 13, steps 4-6).

Input: alignments of reads to the per-sample reference (:mod:`kfp_pileup.reference`),
several per read: DIAMOND blastx tabular output with ``DIAMOND_FIELDS`` (one file per mate),
or bwa-mem2 SAM (``mem -a``, mates in one file). Alignments of one read must be adjacent,
as both tools write them.

Each alignment h of a read has a score S_h in bits (DIAMOND's bit score; bwa's ``AS`` times
``score_nats`` / ln 2) and, per aligned reference column, a term t_h,q in the same units:
the substitution score of that column (BLOSUM62 or bwa's match/mismatch), keyed by the read
position q it uses.

- *Read posterior:* w_h = 2^S_h / sum_g 2^S_g over the read's alignments, decoys included.
  A read whose best alignment is not a target (a neighbour or a reversed decoy) is dropped
  from the pileups.
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
not in pileups, and there is no sequencing-error model. :mod:`kfp_pileup.refine` (steps 7-9)
runs after it on the alignments kept here.
"""

import math
import re
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import Final

import numpy as np
import polars as pl

from kfp_pileup.reference import REF_FASTA, REFS
from kfp_pileup.refine import NOVEL, refine
from kmer_functional_profiler.index import read_fasta

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
ROLES: Final = ("target", "neighbour", "decoy")
TARGET, NEIGHBOUR, DECOY = range(3)
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
        self.hit: list[int] = []
        self.col: list[int] = []
        self.qpos: list[int] = []
        self.rres = bytearray()
        self.sres = bytearray()

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
        return len(self.read) - 1

    def columns(self, hit: int, cols: range, qpos: Sequence[int], rres: bytes, sres: bytes) -> None:
        self.hit.extend([hit] * len(cols))
        self.col.extend(cols)
        self.qpos.extend(qpos)
        self.rres += rres
        self.sres += sres


_BTOP: Final = re.compile(r"(\d+)|(.)(.)")


def _diamond_hit(batch: _Batch, hit: int, fields: list[str], ref: bytes) -> None:
    """Columns of one DIAMOND alignment from its BTOP (query letter, then subject letter)."""
    qstart, sstart, frame, btop = int(fields[3]), int(fields[5]), int(fields[7]), fields[8]
    # Read position key of the codon of query residue i: twice its first base on the forward
    # strand, plus 1 on the reverse strand (DIAMOND's qstart is then the higher coordinate).
    step = 6 if frame > 0 else -6
    key0 = 2 * (qstart - 1) + (frame < 0)
    col, qi = sstart - 1, 0
    for num, q, s in _BTOP.findall(btop):
        if num:
            n = int(num)
            run = ref[col : col + n]
            batch.columns(
                hit,
                range(col, col + n),
                range(key0 + step * qi, key0 + step * (qi + n), step),
                run,
                run,
            )
            col, qi = col + n, qi + n
        elif s == "-":  # an insertion in the read: no reference column
            qi += 1
            batch.ins[hit] += 1
        elif q == "-":  # a deletion: a reference column without a read residue
            batch.columns(hit, range(col, col + 1), [-1], b"-", s.encode())
            col += 1
        else:
            batch.columns(hit, range(col, col + 1), [key0 + step * qi], q.encode(), s.encode())
            col, qi = col + 1, qi + 1
    batch.end[hit] = col


def read_diamond(
    paths: Sequence[str | Path], seq_index: dict[str, int], refs: list[bytes], batch_reads: int
) -> Iterator[_Batch]:
    """Batches of DIAMOND alignments; ``paths`` are the mates' outputs, in mate order."""
    batch, reads = _Batch(), 0
    for mate, path in enumerate(paths):
        last, read = None, -1
        with open(path) as f:
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if fields[0] != last:
                    if reads >= batch_reads:
                        yield batch
                        batch, reads = _Batch(), 0
                    last, read, reads = fields[0], batch.start_read(fields[0]), reads + 1
                seq = seq_index[fields[1]]
                hit = batch.start_hit(
                    read, mate, seq, float(fields[2]), int(fields[5]) - 1, int(fields[7]) < 0
                )
                _diamond_hit(batch, hit, fields, refs[seq])
    if reads:
        yield batch


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


def _chars(data: bytearray) -> np.ndarray:
    return np.frombuffer(bytes(data), dtype="S1").astype("U1")


def score_batch(
    batch: _Batch, refs: pl.DataFrame, matrix: np.ndarray, term_bits: float
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """(assignments, pileup part, per-sequence read stats, target columns) of one batch.
    Assignments keep each alignment's ``hit`` id, reference span, strand, read residues
    (``qlen``) and differences (``mism``) for :mod:`kfp_pileup.refine`; target columns are
    (hit, col, qpos, residue byte) of alignments to targets."""
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
            "residue": _chars(batch.rres),
            "rbyte": rres,
            "ref": _chars(batch.sres),
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
    others = rows.filter(pl.col("role") != DECOY).select(*MATE, "qpos", residue="ref", other="unit")
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
        mism=pl.col("mism").fill_null(0) + pl.col("ins"),
    )
    target_cols = rows.filter(pl.col("role") == TARGET).select(
        "hit", "col", "qpos", residue="rbyte"
    )
    return assignments, part, stats, target_cols


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
    reads: Sequence[str | Path] = (),
    novel: float = NOVEL,
) -> dict[str, int]:
    """Read posteriors, pileups and the unit table from ``alignments`` (``aligner``
    ``diamond``: one tabular file per mate; ``bwa``: one SAM) against the reference in
    ``ref_dir``. Writes ``assignments.parquet``, ``pileup.parquet``, ``variants.tsv`` and
    ``units.tsv`` to ``out_dir`` (the crude pass), then :func:`kfp_pileup.refine.refine`'s
    outputs (steps 7-9; ``reads``: the FASTQ mates, for qualities); returns counts.

    ponytail: steps 7-9 hold every alignment and its target columns in memory; spill them
    per component to Parquet if ``pileup-cost``'s 100 M pairs do not fit."""
    ref_dir, out = Path(ref_dir), Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
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
    seq_index = {s: i for i, s in enumerate(refs["seq_id"])}
    ref_bytes = [sequences[s].encode() for s in refs["seq_id"]]
    if aligner == "diamond":
        batches = read_diamond(alignments, seq_index, ref_bytes, batch_reads)
        matrix, term_bits = blosum62(), DIAMOND_BITS
    elif aligner == "bwa":
        if len(alignments) != 1:
            raise ValueError("bwa: one SAM file with both mates")
        batches = read_sam(alignments[0], seq_index, ref_bytes, batch_reads, score_nats / LN2)
        matrix, term_bits = dna_scores(), score_nats / LN2
    else:
        raise ValueError(f"aligner must be diamond or bwa, not {aligner}")
    parts, stats, n_assigned = [], [], 0
    assigned_parts, col_parts = [], []
    for batch in batches:
        assignments, part, stat, cols = score_batch(batch, refs, matrix, term_bits)
        offset = pl.lit(n_assigned, pl.UInt32)
        assigned_parts.append(assignments.with_columns(pl.col("hit") + offset))
        col_parts.append(cols.with_columns(pl.col("hit") + offset))
        n_assigned += assignments.height
        parts.append(part)
        stats.append(stat)
        if len(parts) > 32:  # keep the summed pileup small
            parts = [_sum(parts)]
    hits = pl.concat(assigned_parts) if assigned_parts else None
    refined, em_units = None, None
    if hits is not None and hits.height:
        refined, em_units = refine(
            hits, pl.concat(col_parts), refs, sequences, out, aligner == "diamond", profile,
            reads, min_depth, min_freq, novel,
        )  # fmt: skip
        hits = hits.join(refined, on="hit", how="left")
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
    units = _units(pl.concat(stats) if stats else None, columns, refs, profile)
    if em_units is not None:
        units = units.join(em_units, on="unit", how="left")
    units.write_csv(out / "units.tsv", separator="\t")
    return {
        "alignments": n_assigned,
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
        consensus=pl.col("residue").sort_by("weight").last(),
        consensus_freq=pl.col("weight").max() / pl.col("weight").sum(),
    )
    alt = (
        pile.filter(pl.col("residue") != pl.col("ref"))
        .sort("weight", descending=True)
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
    profile: str | Path | None,
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
    if profile is not None:
        prof = pl.read_csv(profile, separator="\t", infer_schema_length=None)
        wanted = [
            c
            for c in prof.columns
            if c in ("present_prob", "coverage_zi", "len_mean", "n_members")
            or c.startswith(("aai", "ani"))
        ]
        units = units.join(
            prof.select(pl.col("unit").cast(pl.UInt32), *wanted), on="unit", how="left"
        )
    return units.with_columns(pl.col("identity", "identity_rep").fill_nan(None)).sort("unit")
