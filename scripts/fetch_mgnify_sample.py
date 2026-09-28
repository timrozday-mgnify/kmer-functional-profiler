"""Fetch a small development sample of MGnify Proteins into data/mgnify/ (not the release).

Downloads the release docs and the first ``--mb`` compressed megabytes of the 90%
cluster representatives FASTA with an HTTP range request, keeping only whole records.
Full-release work runs on HPC via the Nextflow pipeline, not this script.
"""

import argparse
import urllib.request
import zlib
from pathlib import Path

RELEASE = "https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/current_release/"
DOCS = ["README.md", "STATS.md", "LICENSE", "md5sum.txt", "mgy30/mgy30README.md"]
OUT = Path(__file__).resolve().parent.parent / "data" / "mgnify"


def fetch(url: str, max_bytes: int | None = None) -> bytes:
    request = urllib.request.Request(url)
    if max_bytes is not None:
        request.add_header("Range", f"bytes=0-{max_bytes - 1}")
    with urllib.request.urlopen(request, timeout=60) as response:
        return bytes(response.read())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mb", type=int, default=20, help="compressed MB of FASTA to read")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    for doc in DOCS:
        (OUT / Path(doc).name).write_bytes(fetch(RELEASE + doc))

    chunk = fetch(RELEASE + "mgy_clusters.fa.gz", args.mb * 2**20)
    text = zlib.decompressobj(wbits=47).decompress(chunk).decode()
    whole_records = text[: text.rindex("\n>") + 1]  # drop the truncated last record
    (OUT / "mgy_clusters_head.faa").write_text(whole_records)
    print(f"{whole_records.count('>')} representatives -> {OUT / 'mgy_clusters_head.faa'}")


if __name__ == "__main__":
    main()
