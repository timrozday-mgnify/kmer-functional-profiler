#!/usr/bin/env bash
# Build the kMermaid image the benchmark uses (the other tools use public biocontainers).
#
#   bash containers/build.sh                 # docker image kfp-kmermaid:edcb4ed (local runs)
#   bash containers/build.sh --sif           # also kfp-kmermaid.sif, for Singularity on HPC
#   bash containers/build.sh --push REGISTRY # also push REGISTRY/kfp-kmermaid:edcb4ed
#
# On HPC, pass the image with --kmermaid_container: the .sif's absolute path, or
# docker://REGISTRY/kfp-kmermaid:edcb4ed after --push. Without docker on the build machine,
# copy kfp-kmermaid.tar (made by --sif) and run there:
#   singularity build kfp-kmermaid.sif docker-archive://kfp-kmermaid.tar
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
image=kfp-kmermaid:edcb4ed
cp "$here/../kmermaid_kfp.py" "$here/kmermaid/"
trap 'rm -f "$here/kmermaid/kmermaid_kfp.py"' EXIT
docker build --platform linux/amd64 -t "$image" "$here/kmermaid"
case "${1:-}" in
    --sif)
        docker save "$image" -o kfp-kmermaid.tar
        if command -v singularity >/dev/null; then
            singularity build kfp-kmermaid.sif docker-archive://kfp-kmermaid.tar
        elif command -v apptainer >/dev/null; then
            apptainer build kfp-kmermaid.sif docker-archive://kfp-kmermaid.tar
        else
            echo "no singularity/apptainer here: copy kfp-kmermaid.tar to HPC and build there" >&2
        fi
        ;;
    --push)
        docker tag "$image" "$2/$image"
        docker push "$2/$image"
        ;;
esac
