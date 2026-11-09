#!/usr/bin/env bash
#SBATCH --job-name=kfp-pfam-anchor-hadza
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# pfam-anchor's spikes in a real metagenome: ERR7738575 (Hadza gut, as tree-index-reads). See ../pfam-anchor/run.sh.
BG=ERR7738575 source ../pfam-anchor/run.sh
