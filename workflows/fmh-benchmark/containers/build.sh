#!/usr/bin/env bash
# Build the images the benchmark cannot take from biocontainers: kMermaid (no package),
# HUMAnN 3.9 (biocontainer + a bowtie2 that works) and HUMAnN 4 (GitHub only, built on it).
#
#   bash containers/build.sh                 # docker images (local runs)
#   bash containers/build.sh --sif           # also <name>.sif files, for Singularity on HPC
#   bash containers/build.sh --push REGISTRY # also push REGISTRY/<image>
#
# On HPC, pass each image with --kmermaid_container, --humann3_container and
# --humann4_container: the .sif's absolute path, or docker://REGISTRY/<image> after --push.
# Without singularity here, --sif leaves <name>.tar files; copy them to HPC and run there:
#   singularity build <name>.sif docker-archive://<name>.tar
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
# directory, image; in build order (humann4 builds FROM the humann3 image)
images=(
    "kmermaid kfp-kmermaid:edcb4ed"
    "humann3 kfp-humann3:3.9-bt2.5.5"
    "humann4 kfp-humann4:e07b3a3"
)
cp "$here/../kmermaid_kfp.py" "$here/kmermaid/"
trap 'rm -f "$here/kmermaid/kmermaid_kfp.py"' EXIT
sing=$(command -v singularity || command -v apptainer || true)
for entry in "${images[@]}"; do
    read -r dir image <<<"$entry"
    name="kfp-$dir"
    docker build --platform linux/amd64 -t "$image" "$here/$dir"
    case "${1:-}" in
        --sif)
            docker save "$image" -o "$name.tar"
            if [ -n "$sing" ]; then
                "$sing" build "$name.sif" "docker-archive://$name.tar" && rm "$name.tar"
            fi
            ;;
        --push)
            docker tag "$image" "$2/$image"
            docker push "$2/$image"
            ;;
    esac
done
