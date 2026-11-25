# The lookup image: `lookup` (rs/) and bundle.py, for workflows/lookup-index and for the raw-reads
# pipeline's PFAM_LOOKUP. Published nowhere: build it where Docker is, for amd64, and turn it into a
# Singularity image on the HPC, as kfpfam's image is:
#
#   docker build --platform linux/amd64 -t kfp-lookup workflows/pfam-anchor/lookup
#   docker save kfp-lookup:latest -o kfp-lookup.tar
#   singularity build kfp-lookup.sif docker-archive://kfp-lookup.tar      # on the HPC
#
# x86-64-v3 (AVX2): the scan's row loop vectorises; a node older than Haswell would stop at an
# illegal instruction. The Python packages are the versions the repo's tests ran with.
FROM rust:1.92-slim-bookworm AS build
WORKDIR /src
COPY rs/Cargo.toml rs/Cargo.lock rs/
COPY rs/src rs/src
ENV RUSTFLAGS="-C target-cpu=x86-64-v3"
RUN cargo build --release --offline --manifest-path rs/Cargo.toml

FROM python:3.12-slim-bookworm
RUN pip install --no-cache-dir polars==1.44.2 pyarrow==25.0.1 numpy==2.5.3
COPY --from=build /src/rs/target/release/lookup /usr/local/bin/lookup
COPY bundle.py /opt/lookup/bundle.py
ENV LOOKUP_BIN=/usr/local/bin/lookup PYTHONDONTWRITEBYTECODE=1
# only Python is run here: an amd64 build on Apple silicon is emulated, and may lack AVX2
RUN test -x /usr/local/bin/lookup && python3 /opt/lookup/bundle.py --help > /dev/null
