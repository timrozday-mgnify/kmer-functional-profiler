#!/usr/bin/env bash
#SBATCH --job-name=kfp-ablation-reads-miseq
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
source ../common.sh
fmh reads --iss_model miseq
