# lookup-index

Builds the MGnify90 lookup bundle (`workflows/pfam-anchor/lookup/bundle.py`) at release scale, from
hmm-annotate's Pfam domains of every protein, representatives and members. The raw-reads pipeline's
`PFAM_LOOKUP` profiles samples with it.

```text
MEMBERSHIP  the release's cluster table by protein_id range   (mgnify-subset's; it is sorted by representative)
GROUPS      checks the annotations; families in ~--groups groups of equal rows, the largest alone
SPLIT       per batch of annotation chunks: cluster and labels joined, rows partitioned by group
BUILD       per group: bundle.py build-part, its families whole
ASSEMBLE    bundle.py assemble -> bundle/
CHECK       bundle.py verify, and profile of --check_parquet
```

`bundle.py` is taken from the repo; the image gives it the `lookup` binary and Python packages.

## Run

Test profile (hmm-annotate's output for the mini release, checked in), after
`cargo build --release` in `workflows/pfam-anchor/lookup/rs`:

```bash
nextflow run workflows/lookup-index -profile test
```

On HPC: `hpc/kfp-ablations/lookup-index/run.sh`. The image is built where Docker is:

```bash
docker build --platform linux/amd64 -t kfp-lookup workflows/pfam-anchor/lookup
```

```bash
docker save kfp-lookup:latest -o kfp-lookup.tar
```

```bash
singularity build kfp-lookup.sif docker-archive://kfp-lookup.tar
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--pfam` | | hmm-annotate's output directory for Pfam, `--targets all --sample 1` |
| `--hmm` | | The `Pfam-A.hmm[.gz]` it was made with (checked by SHA-256) |
| `--labels` | | `ko=DIR,ec=DIR`: hmm-annotate outputs whose models are the labels |
| `--groups` | `400` | `BUILD` tasks, about |
| `--batch_chunks` | `50` | Annotation chunks per `SPLIT` task |
| `--pfam_version`, `--mgnify_version` | `39.0`, `2026_07` | Recorded in the bundle's manifest |
| `--check_parquet` | | A sample's `*_pfam_alignments.parquet` to profile with the bundle |
| `--lookup_sif` | | The lookup image |

## Output

- **`bundle/`**: the index; `bundle.py`'s header has the format. `manifest.json` records rows, member
  rows, member rows whose representative has no row (unreached), and rows whose label is past the
  65,535 commonest of its family (`rows_label_dropped`).
- **`groups.tsv`, `batches.tsv`**: rows per family and its group; protein_id range per `SPLIT` task.
- **`check.*`**: `bundle.py verify`'s line, and the sample's profile if one was given.

`BUILD` holds a family whole: 64 GB, 128 GB on the retry. PF00005 sets it (11.5 M representatives).
