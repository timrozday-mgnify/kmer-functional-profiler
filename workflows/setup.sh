#!/usr/bin/env bash
# Create the repo's Python environment (builds the Rust extension) for the pipelines.
# Needs uv (https://docs.astral.sh/uv/) and a Rust toolchain (rustup) on PATH.
set -euo pipefail
cd "$(dirname "$0")/.."
uv sync --locked --no-dev
.venv/bin/kmer-functional-profiler version
