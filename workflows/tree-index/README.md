# tree-index

Phase 0 of [the tree-index plan](../../docs/tree-index-plan.md): measure the k-mer tree on a
sample of MGnify before any index is built. It also streams a few counts over the whole
release. Its gate decides whether phases 1–7 go ahead.

```text
SAMPLE   per members bucket: clusters with cluster_rep % --sample == 0, their stored k-mers
CORE     -> sample.parquet, core_hashes.npy
STREAM   per bucket of the whole release: set digests, pigeonhole signatures and k-mer holders
         in hash slices, pre-aggregated and split by hash range; representatives outside the
         sample sharing >= --min-shared stored k-mers with it                (--stream false: skip)
REDUCE   per hash range: group sizes
COMBINE  -> stream.tsv, stream.json
EXTRA    -> extra.parquet (at most --max_extra representatives, sharing most first)
MEASURE  the sample (core) and the sample plus representatives (all) -> measure-core/, measure-all/
TIER1    today's standard tier on the sample's proteins -> tier1/meta.json
REPORT   -> gate.tsv
```

Stored sets follow the plan. Each sequence's set is its k-mers with hash ≤ *t\_v*, where
*t\_v* = max(*t\_base*, *n\_min* / (*L* − *k* + 1)). That rate is capped at *t\_cap* and
rounded up to one of 8 log-spaced levels. The defaults are k 11, *t\_base* 0.02,
*n\_min* 16 and *t\_cap* 0.2, which match today's standard tier.

The Rust kernels (`crates/core/src/tree.rs`) do the heavy work:

- candidate pairs from an inverted index (exact, uncapped);
- Kruskal with the virtual empty root, and the merge forest;
- each diff scheme's postings per k-mer, computed on virtual trees of the k-mer's holders.

The scheme postings are checked against a full-tree brute force on random forests.

## Setup

Needs Java 17+, Nextflow, [uv](https://docs.astral.sh/uv/) and a Rust toolchain:

```bash
bash workflows/setup.sh
```

## Run

Test profile: synthetic families in two buckets, about a minute.

```bash
nextflow run workflows/tree-index -profile test
```

On a laptop, point `--store` at a local members table. The `laptop` profile runs one job at
a time with two threads:

```bash
nextflow run workflows/tree-index -profile laptop --store "$PWD/members.parquet" --sample 1 \
    --measure_args '--dense-seqs 20000 --sensitive-seqs 20000'
```

Full scale runs from the release's members buckets on Slurm. See
`hpc/kfp-ablations/tree-index-phase0/run.sh`.

## Outputs

| Output | Content |
| --- | --- |
| `gate.tsv` | The gate. Postings (Fitch) ≤ 1.5× today's tier-1 postings, and leaves ≤ 2× units, at matched *t*. Also the ratios without floors and with floors on components only, and the diff-scheme decision. |
| `measure-core/measure.json` | Every scalar below, per variant (`floors`, `no_floors`). |
| `collapse` | Empty sets, equal-set collapses, contained sets, leaves, and the containment rate of partial and full-length sequences. |
| `postings` | Fitch, majority and add-only postings on the merge forest; the sequence-node scheme (`sequence_mst`, whole set differences along the spanning tree); flat per-leaf postings; build times. |
| `staged` | Stages 3–4 against the exact tree: candidates and false positives per stage, the largest signature group, capped postings, total weight over the exact tree's, and Fitch on the staged tree. |
| `cuts.tsv` | Postings saved by dropping k-mers with more than 2, 4 or 8 gains, and by today's cut (held by more than 16 units). |
| `changes.tsv`, `gains.tsv` | Histograms of Fitch changes and gains per k-mer. The consistency index is 1 / changes. |
| `identity.tsv` | Edge cost (weight / larger stored set; 1 = as costly as storing fresh) by overlap-alignment identity, for spanning-tree edges and for pairs sharing any stored k-mer. |
| `prefix_postings.tsv` | Stage 4's prefix posting lengths: holders among representatives per prefix k-mer. |
| `chunks_stored` | Duplicated postings, chunk count and the share of forest components over budget, at each `--chunk-sizes` budget (stored k-mers). |
| `density` | Distinct over total k-mers at *t* 0.06, 0.2 and 1, and at the stored rates. |
| `one_topology` | On a subsample (`--dense-seqs`), Fitch at *t* 1 (k 11, and k 8) along the stored-rate tree, against a tree built at those settings. Dense chunking is included. |
| `sensitive` | On one representative per storage component of a subsample (`--sensitive-seqs`), the exact k 6 MST against the seeded build. Seeds are 20 lowest murphy10 10-mers, with centre and sort-neighbour candidates; weights come from a *t* 0.25 sketch. Reported: weights in divergence units, the share of exact MST edges among the seeded candidates, the share of pairs within *R* that the seeded graph does not connect at ≤ *R*, and k 11 storage along the sensitive topology against the storage tree. |
| `measure-all/` | The same on the sample plus the cross-cluster representatives. `mst_edges_across_units` shows cross-cluster sharing. |
| `stream.json`, `stream.tsv` | Whole-release counts. Per kind (`digest`, `signature`, `slice<t>`, `dense<t>`): groups, members, the largest group, and distinct over total. Also `equal_set_collapse_share`, the 20 largest groups per kind, and totals. |

Not covered here:

- the sparse-against-dense comparison on `presence-calibration`'s misses, which needs a
  topology over the full index's units (phase 5);
- how many chunks the spike benchmark's candidate subtrees touch;
- NN-descent refinement (stage 5);
- dense tie-breaks.
