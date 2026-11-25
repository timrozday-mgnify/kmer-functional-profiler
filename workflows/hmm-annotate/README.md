# hmm-annotate

Annotates MGnify proteins with any HMM collection (Pfam, KOfam, ...): kfpfam's k-mer prefilter picks
candidate models, and HMMER (pyhmmer) scores only those. Or HMMER scores every model, with
`--engine hmmsearch`. It gives the domain rows with alignments
that `lookup-index` builds from, and a label table (protein, model) for any collection.

```text
PREP           HMM file -> kfpfam prep, tiled index and gapped bars, pressed models, thresholds.tsv
PLAN           ids of the targets -> ranges.tsv, protein_id ranges of ~--chunk_proteins   (pfam-reannotate's)
FASTA          per range: the sequences                                                   (pfam-reannotate's)
SEARCH         per range: filter join | annotate.py search -> domains/ and labels/
               (--engine hmmsearch: annotate.py search --all-pairs, no index)
MEMBER_PAIRS   with --members: the models each representative nearly passes -> its members' pairs, by member id
MEMBER_FASTA   per --member_width ids: those members' sequences and pairs
MEMBER_SEARCH  per --member_width ids: HMMER on the pairs -> domains/ and labels/
CHECK          per range: 1 in --check_sample proteins against every model, no prefilter  (--check_sample N)
COMPARE        per range: the passing domains against the release's own Pfam table     (--compare_release true)
SUMMARY        -> summary.tsv, one total row
CHECK_SUMMARY  -> check.tsv and compare.tsv, recall against each, worst model first
```

## Setup

Needs Java 17+, Nextflow (not 26.04, whose `-resume` reuses tasks run with other parameters), [uv](https://docs.astral.sh/uv/) and kfpfam's image.

```bash
bash workflows/setup.sh   # creates .venv in the repo root
```

kfpfam's image is published nowhere. Build it where Docker is, from kfpfam's `prep-any-hmm` branch
(`e54c21a`) or later, and convert it on the HPC. Older ones need `ACC` and `GA` lines, which KOfam's profiles
lack, panic on a protein over 65,535 residues, and are very slow at these settings.

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
| `--engine` | `kfpfam` | `kfpfam`: its prefilter, then HMMER on the candidates; `hmmsearch`: HMMER on every model |
| `--cutoff` | `gathering` | The HMM file's own thresholds: `gathering`, `trusted` or `noise` |
| `--thresholds` | | Or a TSV with a header: model, threshold, `full` or `domain`. KOfam's `ko_list` as it is |
| `--targets` | `all` | `all`: every protein; `reps`: MGnify90 representatives |
| `--members` | `false` | With `--targets reps`: their members too, so every protein. A member (90% identity) is searched with HMMER against only the models its representative passes to within `--member_slack` (8) bits |
| `--sample` | `1000` | 1 in N targets by `protein_id`; with `--members`, 1 in N clusters |
| `--check_sample` | `0` | 1 in N of those also against every model, without the prefilter; `0`: none |
| `--compare_release` | `false` | Compare with the release's `mgy_proteins_pfam.parquet`, made with Pfam 38.2 |
| `--domains` | `true` | Write the domain rows; labels are always written |
| `--weak_evalue` | `0` | Also domains under the thresholds, at this E-value (Z = 10⁶); `0`: none |
| `--chunk_proteins` | `500000` | Proteins per `SEARCH` task |
| `--kfpfam_sif` | | kfpfam's image |
| `--kfpfam_mode` | `top` | `top`: as near every domain as it goes; `fast`: gapped bars per model from its threshold, 99.8% of domains. Settings: `kfpfam_modes` in `nextflow.config` |
| `--kfpfam_band`, `_diagonals` | 15, `any` | A band around every diagonal at the bar, until one reaches the gapped bar |
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
- **`ranges.tsv`**, and **`ranges.members.tsv`** with `--members` (chunks `mN`): the `protein_id` range of each chunk.
  `lookup-index` reads both.
- **`thresholds.tsv`**: the thresholds used (`inf`: never passes; `-inf`: not used).
- **`hmm.sha256`**: of the HMM file, unzipped. `lookup-index` checks its models against it.
- **`summary.tsv`**: proteins, residues, proteins over 100,000 residues (HMMER's limit, not
  searched), candidate pairs, domains and proteins passing, and pyhmmer's core-seconds.
- **`check.tsv`**: per model, passing (protein, model) pairs and domains by full `hmmsearch`, how many
  `SEARCH` found, and how many it has that `hmmsearch` doesn't (`extra`, none expected). `check.members.tsv`: of the members.

- **`compare.tsv`**: the same columns with the release's own Pfam table as the truth; a domain is the same if
  its envelope is. Use it with Pfam 38.2 only, and `--cutoff gathering`, as the table was made.

`work/` holds the FASTA of every chunk until deleted: about 1.1 TB for every protein of the release.

## Measured

On the laptop (M4 Pro), 19,655 proteins of the 40-family slice (4.9 M residues) against its 40 models:
13,293 of the 13,318 domains at GA (0.9981), worst family 0.975 (PF14306), no extra domain, 64 bytes a
domain row. This is whole proteins of families that are there; the recall over all of Pfam, the recall on
KOfam (whose models are longer) and the cost per residue are for the 1-in-1,000 run to give. That was
with kfpfam's Mode A, before `--kfpfam_mode`.

Every fourth model of Pfam 38.2 (7,534), M4 Pro, 4 threads, CPU-s, native ARM, 2026-11-25. Truth is
`--engine hmmsearch`, which gives `hmmsearch`'s domains exactly; no mode found one it doesn't.

| 80,000 proteins (13.6 M residues) | filter | HMMER | against `hmmsearch` | domains missed, of 13,880 |
| --- | ---: | ---: | ---: | ---: |
| `hmmsearch` | | 404 | 1 | 0 |
| top as for reads (band 16, every diagonal) | 344 | 105 | 0.9× | 0 |
| **`top`** (band 15, `any`) | 150 | 105 | 1.6× | 0 |
| HMMER on the passing pairs alone: the floor | | 7.5 | 54× | 0 |

`fast` on 108,471 other proteins (18.0 M residues): 26 + 30 CPU-s against 530, **9.4×**, 38 of 17,667
domains missed (0.2%). Its bars were chosen on the first 80,000, and were 0.7 bits higher then, with
a query length in them.

`fast` assumes no query length. The same settings, (protein, model) pairs passing:

| | filter | HMMER | against `hmmsearch` | pairs missed |
| --- | ---: | ---: | ---: | ---: |
| 80,000 whole proteins | 21 | 26 | 8.6× of 403 | 41 of 12,867 (0.32%) |
| the same in 290,248 pieces of 50 residues | 9 | 4.5 | 32× of 438 | 13 of 16,070 (0.08%) |

Why this is not the 30× of reads:
- the same residues cut into 50-residue pieces take top's filter 47 to 54 CPU-s, against 344 (band 16)
  or 150 (band 15): a diagonal at the bar is scored over the whole protein, and each hit's window
  is full;
- top's bars are for reads, so it sends HMMER 2.2 M pairs, 170 times the passing ones;
- `hmmsearch` is at its best on this CPU: on those pieces, all of it, top is 6.8× faster. Codon's x86
  (AVX-512 for kfpfam, SSE for HMMER) should favour kfpfam. Not measured on proteins.

`--members`, 49,998 proteins in 8,953 clusters of the 40-family slice, the same models, through the
workflow; domains missed against `hmmsearch` on every protein (357 CPU-s), none extra:

| Representatives by | representatives missed, of 5,115 | members missed, of 16,143 |
| --- | ---: | ---: |
| `hmmsearch` | 0 | 12 (0.07%) |
| `kfpfam`, `top` | 1 | 12 (0.07%) |
| `kfpfam`, `fast` | 6 (0.12%) | 228 (1.4%), 196 of them in one cluster of 574 |

Members cost about 9 CPU-s, 35× less than `hmmsearch` on them, whatever the representatives' engine. A
representative missed loses its cluster's domain, so `fast`'s loss is the same rate in expectation but
comes in clusters; looser bars (`slack` 14 or 18) did not change it. The slice is of chosen families with
large clusters: 1 in 1,000 clusters and `check.members.tsv` give it for the release. Estimated for the
release at 38.2 and laptop speed: 63,600 core-h by `hmmsearch`; with `--members` about 16,000 with
`hmmsearch` on the representatives (23% of the residues) and 2,500 with `fast`. 39.0 has 1.8× the models.
