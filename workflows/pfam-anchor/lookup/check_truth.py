"""The Rust truth against Biopython: the matrix, and 200 (read, target) scores. Fails loudly."""

import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze import blosum  # noqa: E402
from prep import AA, OUT  # noqa: E402
from zymo import scores  # noqa: E402

src = (Path(__file__).parent / "rs/src/data.rs").read_text()
rows = [list(map(int, x.split())) for x in re.findall(r'"(-?\d[-\d ]*)",', src)]
assert (np.array(rows) == blosum()).all()
from Bio.Align import PairwiseAligner, substitution_matrices  # noqa: E402

al = PairwiseAligner(
    mode="local",
    substitution_matrix=substitution_matrices.load("BLOSUM62"),
    open_gap_score=-12,
    extend_gap_score=-1,
)


def load(name):
    off = np.fromfile(OUT / f"bin/{name}.off", "<u8")
    seq = np.memmap(OUT / f"bin/{name}.seq", "u1", "r")
    return lambda i: "".join((AA + "X")[c] for c in seq[int(off[i]) : int(off[i + 1])])


q, t, fam = load("q"), load("reps"), scores("bin")
best = np.fromfile(OUT / "bin/truth.best", np.dtype([("s", "<i2"), ("n", "<u4")]))
rng = np.random.default_rng(0)
for i in rng.choice(len(best), 50, replace=False):
    u, s = fam(int(i))
    assert s.max() <= best["s"][i]
    for k in [int(s.argmax()), *rng.choice(len(u), 3)]:
        assert al.score(t(int(u[k])), q(int(i))) == s[k], (i, u[k], s[k])
print("matrix and 200 scores agree with Biopython")
