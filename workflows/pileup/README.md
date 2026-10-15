# pileup

Read recruitment and crude pileups downstream of a kfp profile (plan: Read recruitment and
variant pileups, phase 13; the nucleotide follow-up of phase 14, step 4). kfp profiles the
sample and shrinks the reference; an existing aligner does the alignment; `kfp-pileup` (a
separate tool, plan: Tools and scope) turns the alignments into read posteriors and pileups.

Per sample:

1. `QUERY`: `kmer-functional-profiler query --aai` against the index (built from
   `--targets` by `INDEX` if no `--index` is given).
2. `REFERENCE`: `kfp-pileup reference`: units with `present_prob >= --min_present` as
   targets, the other hit units of their components as neighbour decoys, the targets
   reversed as null decoys; up to `--members` sequences per unit, from kfp's `fetch` (the
   index's `targets.parquet`, or `--store` members tables, e.g. the partitioned build's MERGE
   output).
3. The aligner the index's molecule calls for:
   - amino-acid index: `DIAMOND` (`makedb`, then `blastx` per mate, every hit within
     `--diamond_top` % of the read's best, tabular with BTOP);
   - nucleotide index: `BWAMEM2` (`index`, then `mem -a`, SAM).
4. `CALL`: `kfp-pileup call`: read posteriors (leave one column out), decoy drop, pileups,
   variants and per-unit identity.

Outputs under `<outdir>/<sample>/`:

| File | Contents |
| --- | --- |
| `profile.tsv` | kfp's profile |
| `ref/ref.fa`, `ref/refs.tsv` | the per-sample reference and its roles |
| `pileup/assignments.parquet` | every alignment: read, mate, reference sequence, unit, role, score (bits), posterior, `kept` (best hit a target) |
| `pileup/pileup.parquet` | summed weights per (sequence, column, residue), with the part from sibling residues |
| `pileup/variants.tsv` | columns with depth >= `--min_depth` and an alternative residue at >= `--min_freq`, with `sibling_share` |
| `pileup/units.tsv` | per target unit: reads, alignment identity (all sequences, representative), the representative's consensus identity, and the profile's `aai` (`ani`) columns, for the `aai` check |

`trace.tsv` records each task's time: the aligner's against `QUERY`'s is step 3's rule (add
a dense pre-filter if the aligner takes more than ~5x).

```bash
nextflow run workflows/pileup -profile docker --samples samples.csv --index idx/
nextflow run workflows/pileup -profile test,docker              # fixture proteins, DIAMOND
nextflow run workflows/pileup -profile test_nucleotide,docker   # fixture CDS, bwa-mem2
```

The steps are the subworkflow `PILEUP` (reads, index, store), which `workflows/read-xref`
reuses.

`samples.csv` has columns `sample,r1,r2`; relative paths are read from the CSV's directory.
DIAMOND and bwa-mem2 run in biocontainers (`-profile docker` locally, `singularity` on HPC);
kfp and kfp-pileup run on the host, in the venv `workflows/setup.sh` makes.
