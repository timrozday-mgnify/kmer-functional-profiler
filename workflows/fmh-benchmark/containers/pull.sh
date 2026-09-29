#!/usr/bin/env bash
# Pull the benchmark's public images into the Singularity cache Nextflow reads, under the
# names it expects, so compute nodes need no registry access. The kMermaid and HUMAnN
# images come from build.sh.
#   bash containers/pull.sh <singularity cache dir>   (the cacheDir in hpc.config)
set -euo pipefail
cache="${1:?usage: pull.sh <singularity cache dir>}"
mkdir -p "$cache"
sing=$(command -v singularity || command -v apptainer)
for image in \
    quay.io/biocontainers/diamond:2.2.8--he361c42_0 \
    quay.io/biocontainers/fmh-funprofiler:1.1.1--pyh106432d_0; do
    name="$(echo "$image" | tr '/:' '--').img"
    [ -e "$cache/$name" ] || "$sing" pull "$cache/$name" "docker://$image"
done
