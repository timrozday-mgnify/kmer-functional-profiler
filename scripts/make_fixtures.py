"""Generate the unit-test fixtures in tests/data/ (deterministic; rerun to regenerate).

Twenty random proteins are back-translated into CDSs followed by random flanking DNA.
Each protein yields two read pairs (300 bp fragment, 2 x 150 bp); R1 gets one of the
test cases below and R2 is the reverse complement of the fragment's far end, so the
true frames cover all six. ``truth.tsv`` records each read's source protein, true
frame and amino-acid span.
"""

import gzip
import random
from pathlib import Path

from kmer_functional_profiler.reference import BASES, CODE_11, reverse_complement

OUT = Path(__file__).resolve().parent.parent / "tests" / "data"
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
CASES = ["clean", "synonymous", "nonsynonymous", "n", "indel", "stop"]
READ_LEN, FRAG_LEN, FLANK_LEN = 150, 300, 300

CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
SYNONYMS = {aa: [c for c, t in zip(CODONS, CODE_11, strict=True) if t == aa] for aa in CODE_11}


def substitute_codon(read: str, offset: int, rng: random.Random, synonymous: bool) -> str:
    """Replace one in-frame codon of ``read`` with a synonymous or missense codon."""
    starts = [
        i
        for i in range(offset, len(read) - 2, 3)
        if len(SYNONYMS[CODE_11[CODONS.index(read[i : i + 3])]]) > 1
    ]
    i = rng.choice(starts)
    aa = CODE_11[CODONS.index(read[i : i + 3])]
    if synonymous:
        choices = [c for c in SYNONYMS[aa] if c != read[i : i + 3]]
    else:
        choices = [c for c, t in zip(CODONS, CODE_11, strict=True) if t not in (aa, "*")]
    return read[:i] + rng.choice(choices) + read[i + 3 :]


def truth_row(pair: int, mate: int, protein: int, cds_start: int, case: str) -> str:
    """TSV row for a read covering the CDS from ``cds_start`` (reverse strand for mate 1)."""
    offset = -cds_start % 3
    aa_start = (cds_start + offset) // 3
    fields = (pair, mate, f"prot{protein}", 3 * mate + offset, aa_start)
    return "\t".join(map(str, (*fields, aa_start + (READ_LEN - offset) // 3, case)))


def fastq(name: str, seq: str) -> str:
    return f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n"


def main() -> None:
    rng = random.Random(20260928)
    OUT.mkdir(parents=True, exist_ok=True)
    proteins, cdss, r1s, r2s = [], [], [], []
    truth = ["read\tmate\tprotein\tframe\taa_start\taa_end\tcase"]
    for p in range(20):
        protein = "M" + "".join(rng.choice(AMINO_ACIDS) for _ in range(rng.randint(120, 400)))
        cds = "".join(rng.choice(SYNONYMS[aa]) for aa in protein) + rng.choice(SYNONYMS["*"])
        context = cds + "".join(rng.choice("ACGT") for _ in range(FLANK_LEN))
        proteins.append(f">prot{p}\n{protein}\n")
        cdss.append(f">prot{p}\n{cds}\n")
        for j in range(2):
            pair = 2 * p + j
            case = CASES[pair % len(CASES)]
            # A stop-case R1 crosses the stop codon; its R2 lies in the flank.
            start = len(cds) - 90 if case == "stop" else rng.randrange(len(cds) - 3 - FRAG_LEN)
            r1 = context[start : start + READ_LEN]
            r2 = reverse_complement(
                context[start + FRAG_LEN - READ_LEN : start + FRAG_LEN].encode()
            ).decode()
            offset = -start % 3
            if case in ("synonymous", "nonsynonymous"):
                r1 = substitute_codon(r1, offset, rng, synonymous=case == "synonymous")
            elif case == "n":
                i = rng.randrange(READ_LEN)
                r1 = r1[:i] + "N" + r1[i + 1 :]
            elif case == "indel":
                i = rng.randrange(30, READ_LEN - 30)
                r1 = r1[:i] + r1[i + 1 :]
            r2_start = start + FRAG_LEN - READ_LEN
            truth.append(truth_row(pair, 0, p, start, case))
            truth.append(truth_row(pair, 1, p, r2_start, "flank" if case == "stop" else "clean"))
            r1s.append(fastq(f"pair{pair}/1", r1))
            r2s.append(fastq(f"pair{pair}/2", r2))

    (OUT / "proteins.faa").write_text("".join(proteins))
    (OUT / "cds.fna").write_text("".join(cdss))
    (OUT / "truth.tsv").write_text("\n".join(truth) + "\n")
    for name, records in (("reads_1.fastq.gz", r1s), ("reads_2.fastq.gz", r2s)):
        with (
            (OUT / name).open("wb") as f,
            gzip.GzipFile(filename="", mode="wb", fileobj=f, mtime=0) as gz,
        ):
            gz.write("".join(records).encode())


if __name__ == "__main__":
    main()
