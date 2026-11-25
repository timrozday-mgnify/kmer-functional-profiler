# hmm-annotate

Annotates MGnify proteins with any HMM collection (Pfam, KOfam, ...): kfpfam's k-mer prefilter picks
candidate models, and HMMER (pyhmmer) scores only those. It gives the domain rows with alignments
that `lookup-index` builds from, and a label table (protein, model) for any collection.

```text
PREP           HMM file -> kfpfam prep and tiled index, pressed models, thresholds.tsv
PLAN           ids of the targets -> ranges.tsv, protein_id ranges of ~--chunk_proteins   (pfam-reannotate's)
FASTA          per range: the sequences                                                   (pfam-reannotate's)
SEARCH         per range: filter join | annotate.py search -> domains/ and labels/
CHECK          per range: 1 in --check_sample proteins against every model, no prefilter  (--check_sample N)
COMPARE        per range: the passing domains against the release's own Pfam table     (--compare_release true)
SUMMARY        -> summary.tsv, one total row
CHECK_SUMMARY  -> check.tsv and compare.tsv, recall against each, worst model first
```

## Setup

Needs Java 17+, Nextflow, [uv](https://docs.astral.sh/uv/) and kfpfam's image.

```bash
bash workflows/setup.sh   # creates .venv in the repo root
```

kfpfam's image is published nowhere. Build it where Docker is, from kfpfam's `prep-any-hmm` branch
or later, and convert it on the HPC. Older ones need `ACC` and `GA` lines, which KOfam's profiles
lack, and panic on a protein over 65,535 residues.

```bash
docker build -f docker/filter.Dockerfile --platform linux/amd64 -t kfpfam-filter <kfpfam checkout>
```

```bash
docker save kfpfam-filter:latest -o kfpfam-filter.tar
```

```bash
singularity build kfpfam-filter.sif docker-archive://kfpfam-filter.tar
```

## Run

Test profile (six mini HMMs, 84 proteins; needs a kfpfam checkout with the filter built, and a Python
with numpy, pyhmmer and pyarrow):

```bash
nextflow run workflows/hmm-annotate -profile test --kfpfam_filter <kfpfam>/filter/target/release/filter --kfpfam_python <python> --kfpfam_prep_py <kfpfam>/measure/prep.py
```

On HPC, `hpc/kfp-ablations/hmm-annotate/run.sh` runs 1 in 1,000 proteins with the check by default.

| Option | Default | Meaning |
| --- | --- | --- |
| `--hmm` | | The collection: one HMM file, plain or gzipped |
| `--cutoff` | `gathering` | The HMM file's own thresholds: `gathering`, `trusted` or `noise` |
| `--thresholds` | | Or a TSV with a header: model, threshold, `full` or `domain`. KOfam's `ko_list` as it is |
| `--targets` | `all` | `all`: every protein; `reps`: MGnify90 representatives |
| `--sample` | `1000` | 1 in N targets by `protein_id` |
| `--check_sample` | `0` | 1 in N of those also against every model, without the prefilter; `0`: none |
| `--compare_release` | `false` | Compare with the release's `mgy_proteins_pfam.parquet`, made with Pfam 38.2 |
| `--domains` | `true` | Write the domain rows; labels are always written |
| `--weak_evalue` | `0` | Also domains under the thresholds, at this E-value (Z = 10⁶); `0`: none |
| `--chunk_proteins` | `500000` | Proteins per `SEARCH` task |
| `--kfpfam_sif` | | kfpfam's image |
| `--kfpfam_bits`, `_bar`, `_gate`, `_min_gapped` | 6.5, 17, 16.5, 21.25 | kfpfam's Mode A |
| `--kfpfam_block` | `80000` | Proteins the filter holds at a time |

A model in the HMM file with no threshold is an error. A `-` in the table is not: that model never
passes, and shows only with `--weak_evalue`.

## Output

- **`domains/chunkN.domains.parquet`**: one row per passing domain, in the HMM file's model order:
  `protein_id`, `model` (accession without its version, or the name), `hmm_from`, `hmm_to`,
  `ali_from`, `ali_to`, `env_from`, `env_to`, `cigar`, `score`, `seq_score`, `i_evalue`, `passes`,
  and `left`, `right`, `residues`: the aligned residues with up to 10 either side, and how many
  those are. `cigar` is `M`, `I` and `D` against the model.
- **`labels/chunkN.labels.parquet`**: `protein_id`, `model`, for proteins passing the model's
  thresholds.
- **`thresholds.tsv`**: the thresholds used (`inf`: never passes; `-inf`: not used).
- **`hmm.sha256`**: of the HMM file, unzipped. `lookup-index` checks its models against it.
- **`summary.tsv`**: proteins, residues, proteins over 100,000 residues (HMMER's limit, not
  searched), candidate pairs, domains and proteins passing, and pyhmmer's core-seconds.
- **`check.tsv`**: per model, passing (protein, model) pairs and domains by full `hmmsearch`, how many
  `SEARCH` found, and how many it has that `hmmsearch` doesn't (`extra`, none expected).

- **`compare.tsv`**: the same columns with the release's own Pfam table as the truth; a domain is the same if
  its envelope is. Use it with Pfam 38.2 only, and `--cutoff gathering`, as the table was made.

`work/` holds the FASTA of every chunk until deleted: about 1.1 TB for every protein of the release.

## Measured

On the laptop (M4 Pro), 19,655 proteins of the 40-family slice (4.9 M residues) against its 40 models:
13,293 of the 13,318 domains at GA (0.9981), worst family 0.975 (PF14306), no extra domain, 64 bytes a
domain row. This is whole proteins of families that are there; the recall over all of Pfam, the recall on
KOfam (whose models are longer) and the cost per residue are for the 1-in-1,000 run to give.
