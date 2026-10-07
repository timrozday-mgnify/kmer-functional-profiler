# gtdb-species-index

Builds a species index (plan: Genome mode, Species model; phase 11, step 11) from a GTDB
release, through the genome-set path of `species-index`. Species and their taxonomy are
GTDB's; carriage and content come from each genome's proteins annotated with the index.

```text
METADATA       bac120/ar53 metadata of the release (--gtdb_url, --metadata), once (storeDir)
PICK           gtdb.py: genomes of CheckM2 completeness >= 50, contamination <= 10, at most
               --max_per_species per species (representative first, then completeness - 5
               contamination); --species or --max_species for a subset. genomes.tsv (with
               each genome's NCBI URL) and annotate.tsv
FETCH          each shard's genomic FASTA from NCBI (--ncbi_url; --max_downloads at once)
GENES          pyrodigal -p single (meta for genomes too short to train on), table 11
ANNOTATE       kmer-functional-profiler annotate-genomes INDEX annotate.tsv
SPECIES_INDEX  kmer-functional-profiler species-index INDEX --genomes -> <outdir>/species_index/
```

Every genome, representatives included, is gene-called by pyrodigal from NCBI's DNA, so all
proteins come from one caller. (GTDB ships representatives' Prodigal proteins only.)

## Run

Needs Java 17+, Nextflow and the project's venv (`bash workflows/setup.sh`); pyrodigal runs
in a biocontainer (`-profile docker` or `singularity`). `--index` is the index profiles are
made with: `species` refuses a profile from another.

```bash
nextflow run workflows/gtdb-species-index -profile singularity,slurm --index /path/to/index --max_species 1000
```

Then fit species with it:

```bash
kmer-functional-profiler species profile.tsv results/species_index species.tsv
```

## Licence

GTDB's taxonomy and metadata are CC BY-SA 4.0. The species index carries GTDB names
(`species.tsv`, `lineage.parquet`, `carriage.parquet`), so a shipped one is a derivative:
distribute it under CC BY-SA 4.0, credit GTDB (Parks et al., Nucleic Acids Res. 2022), and
keep it apart from the GPL code. NCBI places no restrictions on its genome data; the index
holds unit counts derived from it, not the sequences.

## Test

```bash
.venv/bin/kmer-functional-profiler index tests/data/mini_fmh/mgnify_members.parquet workflows/gtdb-species-index/results-test-gtdb/index --k 11 --t-base 0.5 --n-min 0 --pfam tests/data/mini_fmh/mgnify_pfam.parquet
```

```bash
cd workflows/gtdb-species-index && nextflow run . -profile test,docker
```

`tests/data/mini_gtdb` is the mini_uhgg genomes as a GTDB release with an NCBI mirror: three
species, five genomes kept (two per species at most; one genome fails contamination).
