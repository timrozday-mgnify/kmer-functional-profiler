#!/usr/bin/env bash
#SBATCH --job-name=kfp-hmm-annotate-benchmark
#SBATCH --output=%x.%j.log
#SBATCH --time=14-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# kfpfam against hmmsearch on the same MGnify proteins, on this cluster's CPUs: their CPU time, and the
# domains of hmmsearch that kfpfam misses. Inputs, as run.sh's (PFAM=39.0 SAMPLE=1000 RELEASE=<mirror> SIF=...), and
#   ENGINES="hmmsearch kfpfam-top kfpfam-fast"   hmmsearch first: it is what the others are compared with
#   SBATCH_CONSTRAINT=<feature>                  Slurm's own: one kind of node for every task. Without it the
#                                                engines may run on different CPUs; benchmark.tsv's column says which
# Run:
#   sbatch benchmark.sh
# Read benchmark-<sample>.tsv. At SAMPLE=1000, 5.7 M proteins, hmmsearch is about 113 core-hours at laptop speed.
# The sequences are read from the release once, for all of the engines. Re-submitting resumes, and takes
# finished engines' times from their first run. Not beside run.sh: they share work/ and the -resume cache.
set -euo pipefail
engines=${ENGINES:-hmmsearch kfpfam-top kfpfam-fast}
export SAMPLE=${SAMPLE:-1000} COLLECTION=pfam PFAM=${PFAM:-39.0} CHECK=0 COMPARE=false MEMBERS=false
dirs=()
for e in $engines; do
    ENGINE=$e bash run.sh
    dirs+=("results-pfam$PFAM-$e-$SAMPLE")
done
python=../../../kmer-functional-profiler/.venv/bin/python
"$python" ../../../kmer-functional-profiler/workflows/hmm-annotate/bench.py "${dirs[@]}" | tee "benchmark-$SAMPLE.tsv"
