# route-dropped

What read routing drops, searched against every profile unit (plan, phase 13, routing). Reads
are single-end: each mate is routed on its own (`kfp-pileup route`'s k-mer table over the
per-sample reference, k = 8, min_kmers = 3 by default) and dropped when it shares fewer
than `min_kmers` distinct k-mers with every block.

1. `DROPPED` (`dropped.py dropped`): the dropped reads, a seeded subsample of `--reads`
   (0 = all), as one FASTA; `dropped.parquet` (each read's total and shared k-mers) and
   `dropped.json` (reads routed and dropped, dropped reads whose partner is routed, dropped
   reads by total k-mers).
2. `FULLREF` (`dropped.py fullref`): every unit of the profile, with 1 and with 3 members,
   and each sequence's role: target or neighbour (in the per-sample reference), profile
   only, or withheld (`--withhold`, a spike's left-out units).
3. `DIAMOND`: the dropped reads against each full reference, as `workflows/pileup` aligns
   (`--sensitive -F 15 --top 10`).

Then `dropped.py summarise` (not a process: it reads `trace.tsv`, written at the end) gives
per arm and per total-k-mer bin the reads aligned, their identity, the role of the unit hit
(routing's misses: high identity on a target or neighbour), how the unrouted baseline
treated them (`--assignments`) and DIAMOND's CPU-h.

```bash
nextflow run workflows/route-dropped -profile local,docker --ref <sample>/ref \
    --profile <sample>/profile.tsv --index /abs/index --r1 r1.fq.gz --r2 r2.fq.gz --reads 0
.venv/bin/python workflows/route-dropped/dropped.py summarise --dropped results/dropped \
    --arm m1=results/m1 --arm m3=results/m3 --trace results/trace.tsv
```

On HPC: `hpc/kfp-ablations/route-dropped/run.sh` (with `../site.config`).
