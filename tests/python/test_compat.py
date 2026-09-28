"""sourmash-compatible index and query agree with sourmash's own overlap counts."""

import subprocess
from pathlib import Path

import polars as pl
import pytest
import screed
import sourmash
from sourmash.save_load import SaveSignaturesToLocation

from kmer_functional_profiler.compat import import_signatures
from kmer_functional_profiler.index import Index
from kmer_functional_profiler.query import profile

DATA = Path(__file__).resolve().parents[1] / "data"
FMH = Path(__file__).resolve().parents[2] / "data" / "fmh-funprofiler"
K, SCALED = 7, 10


def test_overlap_matches_sourmash(tmp_path: Path) -> None:
    sigs = []
    for rec in screed.open(str(DATA / "proteins.faa")):
        mh = sourmash.MinHash(n=0, ksize=K, is_protein=True, scaled=SCALED)
        mh.add_protein(rec.sequence)
        sigs.append(sourmash.SourmashSignature(mh, name=rec.name))
    with SaveSignaturesToLocation(str(tmp_path / "refs.sig.zip")) as save:
        for sig in sigs:
            save.add(sig)
    import_signatures(tmp_path / "refs.sig.zip", tmp_path / "idx", K)

    reads = (DATA / "reads_1.fastq.gz", DATA / "reads_2.fastq.gz")
    query = sourmash.MinHash(n=0, ksize=K, is_protein=True, scaled=SCALED)
    for path in reads:
        for rec in screed.open(str(path)):
            query.add_sequence(rec.sequence, force=True)
    expected = {
        s.name: overlap for s in sigs if (overlap := len(set(query.hashes) & set(s.minhash.hashes)))
    }
    result = profile(Index.load(tmp_path / "idx"), *reads)
    assert dict(result.select("name", "kmers_hit").rows()) == expected
    assert len(expected) == len(sigs)


@pytest.mark.skipif(not FMH.exists(), reason="fmh-funprofiler data not downloaded")
def test_fmh_funprofiler_demo(tmp_path: Path) -> None:
    """Same KOs and overlaps as ``sourmash prefetch`` in fmh-funprofiler's pipeline."""
    sketches, fastq = FMH / "KOs_sketched_scaled_1000.sig.zip", FMH / "metagenome_example.fastq"
    import_signatures(sketches, tmp_path / "idx", 11)
    ours = profile(Index.load(tmp_path / "idx"), fastq)
    sm = ["uv", "run", "sourmash", "-q"]
    subprocess.run([*sm, "sketch", "translate", "-p", "scaled=1000,k=11,abund", str(fastq),
                    "-o", str(tmp_path / "mg.sig.zip")], check=True)  # fmt: skip
    subprocess.run([*sm, "prefetch", str(tmp_path / "mg.sig.zip"), str(sketches),
                    "-o", str(tmp_path / "p.csv"), "-k", "11", "--scaled", "1000", "--protein",
                    "--threshold-bp", "1000"], check=True)  # fmt: skip
    prefetch = pl.read_csv(tmp_path / "p.csv")
    expected = dict(zip(prefetch["match_name"], prefetch["intersect_bp"] // 1000, strict=True))
    assert dict(ours.select("name", "kmers_hit").rows()) == expected
