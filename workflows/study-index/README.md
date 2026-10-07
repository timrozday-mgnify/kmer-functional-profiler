# study-index

Builds a small index from a study's own proteins (MAG gene calls, assembly gene calls),
to be queried jointly with a base index such as MGnify90 (plan: Additional references,
step 4). Study units then compete with base units for the reads they share, and the
base's rows are unchanged where they share nothing.

```text
STUDY_GENES        contigs or MAGs -> pyrodigal -p meta; protein ids get the file name as a
                   prefix (<file>_<contig>_<n>), so contig names may repeat across MAGs
STUDY_PROTEINS     protein FASTA instead (--proteins): gene calling skipped, ids as given
STUDY_CLUSTER      MMseqs2 linclust (--cluster_args, MGnify90's 90% rule by default)
STUDY_MEMBERS      members.parquet (protein_id, cluster_rep, full_length, sequence);
                   full_length is false for genes Prodigal marks partial
STUDY_PFAM_SEARCH  hmmsearch --cut_ga over --pfam_chunks shares of the proteins
STUDY_PFAM         pfam.parquet (protein_id, pfam_accession)
STUDY_INDEX_BUILD  kmer-functional-profiler index --like BASE --pfam pfam.parquet ->
                   <outdir>/study/<name>/ and <name>.stats.json
```

`--like BASE` copies every build parameter of the base index (k, alphabet, *t\_base*,
floor, dense tier, ...), which a joint query requires. `trace.tsv` records each task's time
and memory: the plan's gate asks that one study of ~10⁶ proteins extends the base in
minutes on one node (the index build; the Pfam search is the long step).

## Run

Needs Java 17+, Nextflow and the project's venv (`bash workflows/setup.sh`); pyrodigal,
MMseqs2 and HMMER run in biocontainers (`-profile docker` or `singularity`).

```bash
nextflow run workflows/study-index -profile singularity,slurm --name my_study --contigs 'mags/*.fa.gz' --base /path/to/mgnify_index --pfam_hmm /path/to/Pfam-A.hmm
```

Set `--pfam_hmm` (or `--pfam_url`) to the Pfam release the base was labelled with, so
both indexes use the same accessions. Then query both together:

```bash
kmer-functional-profiler query /path/to/mgnify_index reads_1.fq.gz reads_2.fq.gz --extra-index results/study/my_study
```

Study units are reported with `source` 1 and ids offset by the base's units.

**Decoys.** The same recipe builds a decoy index from a host or contaminant proteome:
`--proteins proteome.fa.gz --role decoy`. A decoy competes like any extra index, but is
reported as one row per decoy index (`name` `decoy`).

## Test

```bash
.venv/bin/kmer-functional-profiler index tests/data/mini_fmh/mgnify_members.parquet workflows/study-index/results-test/base --k 11 --t-base 0.5 --n-min 0
```

```bash
cd workflows/study-index && nextflow run . -profile test,docker
```

The fixture's three genomes become a 75-protein study index under `results-test/study/`.
