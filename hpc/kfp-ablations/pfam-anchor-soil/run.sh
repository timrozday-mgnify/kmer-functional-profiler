#!/usr/bin/env bash
#SBATCH --job-name=kfp-pfam-anchor-soil
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# pfam-anchor's spikes in a real metagenome: SRR12659369 (Harvard Forest soil, as em-starts-soil). See ../pfam-anchor/run.sh.
BG=SRR12659369 source ../pfam-anchor/run.sh
