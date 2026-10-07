#!/usr/bin/env bash
#SBATCH --job-name=kfp-species-benchmark
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# Phase 11, step 10 (dev-species branch): the species model's strain hold-out benchmark on
# the MGnify human-gut catalogue (~2 GB of pangenomes and genomes fetched from the EBI FTP).
# INDEX=<a built index dir> profiles against it (e.g. a release tier's full build); without
# it the 1-in-100 subset is built here, which keeps ~1% of each species' units.
source ../common.sh
if [[ -n ${INDEX:-} ]]; then
    index=(--index "$INDEX")
else
    index=(--members "$SUBSET/members.parquet" --pfam "$SUBSET/pfam.parquet")
fi
nextflow run "$WF/species-benchmark" -c ../site.config "${index[@]}" \
    --outdir results -resume "$@"
