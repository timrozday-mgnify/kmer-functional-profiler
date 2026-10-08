#!/usr/bin/env bash
# Fetch what the species benchmark needs from an MGnify genome catalogue (FTP layout).
#   fetch.sh species URL SPECIES.txt OUT CPUS   each species' genome/ (faa, fna) and, if it has
#                                               one, pan-genome/ (Rtab, pan-genome.fna) into
#                                               OUT/species_catalogue/<prefix>/<rep>/
#   fetch.sh genomes URL SAMPLES.tsv OUT CPUS   each sampled genome's DNA (the FASTA section of
#                                               its all_genomes GFF) as OUT/<genome>.fna
set -euo pipefail
mode=$1 url=$2 list=$3 out=$4 cpus=$5
mkdir -p "$out"

one_species() {  # url out rep
    local d="species_catalogue/${3:0:${#3}-2}/$3"
    mkdir -p "$2/$d/genome"
    curl -fsSL --retry 3 -o "$2/$d/genome/$3.faa" "$1/$d/genome/$3.faa"
    curl -fsSL --retry 3 -o "$2/$d/genome/$3.fna" "$1/$d/genome/$3.fna"
    if curl -fsSL --retry 3 -o "$2/$d/genome/rtab.tmp" "$1/$d/pan-genome/gene_presence_absence.Rtab" 2>/dev/null; then
        mkdir -p "$2/$d/pan-genome"
        mv "$2/$d/genome/rtab.tmp" "$2/$d/pan-genome/gene_presence_absence.Rtab"
        curl -fsSL --retry 3 -o "$2/$d/pan-genome/pan-genome.fna" "$1/$d/pan-genome/pan-genome.fna"
    else
        rm -f "$2/$d/genome/rtab.tmp"  # a one-genome species has no pangenome
    fi
}

one_genome() {  # url out genome species
    curl -fsSL --retry 3 "$1/all_genomes/${4:0:${#4}-2}/$4/genomes1/$3.gff.gz" \
        | gzip -dc | sed -n '/^##FASTA/,$p' | tail -n +2 > "$2/$3.fna"
    [[ -s "$2/$3.fna" ]] || { echo "no sequence for $3" >&2; exit 1; }
}

export -f one_species one_genome
case $mode in
    species) xargs -P "$cpus" -I{} bash -c 'one_species "$0" "$1" {}' "$url" "$out" < "$list" ;;
    genomes) tail -n +2 "$list" | cut -f2,3 | xargs -P "$cpus" -L1 bash -c 'one_genome "$0" "$1" "$2" "$3"' "$url" "$out" ;;
    *) echo "mode: species or genomes" >&2; exit 1 ;;
esac
