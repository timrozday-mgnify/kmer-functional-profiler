"""Domains and reads on HMM columns, as the flat files that rs/src/data.rs reads."""

import sys
from pathlib import Path
from typing import BinaryIO, Final

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from columns import OP  # noqa: E402

AA: Final = "ACDEFGHIKLMNPQRSTVWY"  # code 0-19; 20: X and others; 21: no residue in the column
CODE: Final = np.full(256, 20, np.uint8)
for i_, a_ in enumerate(AA):
    CODE[ord(a_)] = i_
NONE: Final = 21


def hmm_lengths_of(path: Path) -> dict[int, int]:
    out, acc, length = {}, None, None
    with open(path) as f:
        for line in f:
            if line.startswith("LENG"):
                length = int(line.split()[1])
            elif line.startswith("ACC"):
                acc = int(line.split()[1][2:].split(".")[0])
            elif line.startswith("//"):
                out[acc] = length
    return out


def rows_of(d: pl.DataFrame, lengths: dict[int, int], f_rows: BinaryIO, f_pos: BinaryIO) -> None:
    """Append each domain's row of residue codes, and of protein positions (65535: none)."""
    for acc, hmm_from, ali_from, cigar, seq in d.select(
        "pfam_accession", "hmm_from", "ali_from", "cigar", "sequence"
    ).iter_rows():
        n = lengths[acc]
        row, pos = np.full(n, NONE, np.uint8), np.full(n, 65535, np.uint16)
        codes = CODE[np.frombuffer(seq.encode(), np.uint8)]
        c, p = hmm_from - 1, ali_from - 1
        for k, op in OP.findall(cigar):
            k = int(k)
            if op == "M":
                row[c : c + k] = codes[p : p + k]
                pos[c : c + k] = np.arange(p, p + k)
                c, p = c + k, p + k
            elif op == "I":
                p += k
            else:
                c += k
        f_rows.write(row.tobytes())
        f_pos.write(pos.tobytes())


def pack_queries(q: pl.DataFrame, b: Path, fams: list[int], lengths: dict[int, int]) -> None:
    fam_index = {f: i for i, f in enumerate(fams)}
    (b / "q.fam").write_bytes(
        np.array([fam_index[x] for x in q["pfam_accession"]], np.uint32).tobytes()
    )
    seqs = [CODE[np.frombuffer(s.encode(), np.uint8)] for s in q["sequence"]]
    (b / "q.off").write_bytes(
        np.concatenate([[0], np.cumsum([len(s) for s in seqs])]).astype(np.uint64).tobytes()
    )
    (b / "q.seq").write_bytes(np.concatenate(seqs).tobytes())
    with open(b / "q.rows", "wb") as f1, open(b / "q.pos", "wb") as f2:
        rows_of(q, lengths, f1, f2)
