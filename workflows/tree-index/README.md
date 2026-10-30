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
TIER1_MATCHED  today's index at the survey's other rates (--tier1_variants)   (--survey false: skip)
SURVEY   the floor survey (survey.py) on core and all -> survey-core/, survey-all/
INDEX_SIZES    file sizes of the release's standard tier (--release_index)
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

### Floor survey

`survey.py` measures the design that followed phase 0's failed gate (the plan's Progress log:
forest floors, pieces by coherence, floors on every piece). One topology per sample and leaf
mode is kept whole: the storage tree on per-sequence floored sets. Its leaves (`--leaves`) are:

- `collapsed`: phase 0's, distinct sets at the topology rates after the equal-set and
  containment collapses, each stored as its longest member's k-mers;
- `proteins`: a leaf per distinct dense set (every 8-mer, `--topology-k`, at rate 1,
  `--topology-t`), so proteins share a leaf only when those sets are equal. There is no
  containment collapse, and each leaf stores its own k-mers. The topology is the storage tree
  (Kruskal with the virtual empty root) under |A Δ B| of the dense sets, from an exact pair
  join (`--topology-max-holders` caps common k-mers, then weighs pairs exactly).

On each tree, each configuration stores:

- **base band:** each leaf's k-mers at or below *t\_base*, diffed by Fitch over the whole tree;
- **piece floors:** top-down from each component, a node whose floor is at least *τ* coherent
  becomes a piece (*τ* 0 keeps the components). Its floor is tier 1's: up to *n\_min* k-mers
  above *t\_base*, those held by most of its leaves, posted once at its root. Coherence is the
  mean share of the floor a leaf holds;
- **single-leaf pieces:** a floor of their own *n\_single* lowest hashes above *t\_base*;
- **leaf minimum:** each leaf topped up to *n\_1* stored k-mers of its own, and each protein
  without a leaf (no k-mer up to *t\_cap*) given one.

The default grid is both leaf modes, *t\_base* 0.01, 0.02 and 0.05, *n\_min* 8, 16 and 32,
*τ* 0, 0.3, 0.6 and 0.9, *n\_single* 0, 2, 4 and 16, and *n\_1* 0, 1, 2 and 4
(`--survey_args` to change it). On a quarter of the sample, both modes with *τ* 0 and 0.3 and
*n\_single* 0 took 2.7 min and 2.5 GB.

| Output | Content |
| --- | --- |
| `survey-*/survey.tsv` | One row per leaf mode and configuration. Pieces, single-leaf pieces, mean coherence; base, floor and top-up postings, and their total; distinct k-mers; proteins with no stored k-mer (`blind_proteins`); `over_standard` (the standard tier on the same proteins) and `over_matched` (tier 1 at the same *t\_base* and *n\_min*, where built); release-scale projections of postings, keys and pieces (× release proteins / sample proteins). |
| `survey-*/survey.json` | Proteins and units; per leaf mode, leaves, internal nodes, components and, per *t\_base*, how proteins without an own base k-mer reach the hierarchy (component with relatives, sharing a leaf, alone, no leaf, no valid k-mer). Tier 1's own floor coherence per (*t\_base*, *n\_min*). |
| `tier1-*.json` | Tier 1 at each `--tier1_variants` rate, built on the sample's proteins. |
| `index_sizes.tsv` | File sizes of the release's standard tier, every file under its directory, so posting tables and unit columns can be sized separately. |

Not covered here:

- the sparse-against-dense comparison on `presence-calibration`'s misses, which needs a
  topology over the full index's units (phase 5);
- how many chunks the spike benchmark's candidate subtrees touch;
- NN-descent refinement (stage 5);
- dense tie-breaks.

## Building the tiers (`build.nf`)

`build.nf` builds the three tiers chosen from the per-protein survey (the plan's Progress log)
from one members table, in parts:

```
FIXTURE    synthetic members when --members is '' (the test profile)
PARTITION  partition.py: proteins linked by a possible k-11 storage edge (union-find during
           the pair join, _core.overlap_components) or an equal sequence; components packed
           into parts of ~--part_proteins (2 M) -> partition.json
BUILD      build.py per part and tier: a leaf per distinct dense 8-mer set, the dense 8-mer
           topology, base band by Fitch, floors per piece, per-leaf top-ups
MERGE      merge.py per tier: parts concatenated, positions offset -> index-<tier>/
```

```bash
nextflow run workflows/tree-index/build.nf -profile test,laptop     # the fixture, 5 parts
nextflow run workflows/tree-index/build.nf -profile slurm --members /path/members.parquet
```

`--tiers` holds `name:build.py options` pairs, `;`-separated. The defaults are quick
(*t\_base* 0.01, *n\_min* 8, *τ* 0, *n\_1* 1), medium (0.02, 16, 0.3, 2) and large (0.05, 32,
0.3, 2). Each `index-<tier>/` is in `build.py`'s format (`tree-v0`), which `query.py` reads.
Its `meta.json` has postings by kind, bytes, positions, parts and the parts' build time and peak
RSS.

Measured locally on phase 0's core sample (1.07 M proteins, 4 threads): PARTITION 37 s and
2.6 GB (~2.4 KB per protein); BUILD 68 s and 3.4 GB per 262 k proteins (the medium tier).

Parts lose the dense-topology edges between them, so a k-mer shared across parts is posted in
each. Cuts fall between k-11 storage components, which few such edges cross.
