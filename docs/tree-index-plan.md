# A k-mer tree index: every sequence, no clustering threshold

Plan, 2026-10-24 (the user's request): explore a tier-1 index built as a tree of k-mer diffs over every sequence of a set (all of MGnify Proteins, not only 90% cluster representatives or unions), instead of fixed-identity units; its feasibility, precedent, and a build order. A separate design from [kmer-functional-profiler-plan.md](kmer-functional-profiler-plan.md) (the main plan), which it would replace at tier 1 only if the gates below pass. No code yet.

## Verdict

Feasible, and worth a measurement phase before any build. The encoding is sound and has close precedent. Its cost is known exactly: the number of postings is the parsimony length of the k-mer presence matrix on the tree. Near-identical sequences then cost almost nothing at kfp's sampling rate. Four points decide it:

- **Compactness.** At *t\_base* = 0.02 a substitution changes ~22 k-mers and so ~0.44 sampled postings. A member at 99% identity to its tree neighbour (187 aa, the release mean) costs ~1 posting, and often 0, in which case it collapses into its neighbour's node. Leaf resolution therefore costs roughly as much as the per-unit unions the standard tier already stores. It does not cost the ~2 × 10¹⁰ postings a flat per-protein index would need. This is an estimate: phase 0 measures it.
- **No identity threshold, but resolution is still set by something.** It is set by the sampling density (*t*, k) and by the data. Nodes are the distinct sampled k-mer sets, so what the index can tell apart is exactly what its k-mers can tell apart. Reported units are subtrees chosen per sample, not fixed clusters. For MGnify specifically this is the main gain; it is also why the method generalises to other sets of short sequences.
- **The query keeps its cost; the reporting changes.** Lookups are the same: one table probe per sampled read k-mer. Per-leaf hit counts come from an interval sweep over leaf order, not from per-unit posting expansion, and k-mers shared by a whole clade become one posting instead of a large unit set (today such k-mers are dropped above 16 units). The new work is the assignment model over nested subtrees (greedy, recall-first candidates, or a placement-style model), which today's per-unit gather and EM do not cover.
- **Build is the expensive part, but k-mer joins keep it to sorts.** The edge weight is defined on stored k-mer sets, so the tree can be built exactly from them, with no alignment (Building the tree): set digests, pigeonhole signatures and a rarity-ordered prefix join, then Kruskal by integer weight. That is estimated at hours of HPC, against today's reuse of MGnify90 membership, which is free, and it needs no clustering of the input.

Not shown yet, and the reason for phase 0: whether the tree's postings and node table fit in the standard tier's budget (~108–196 GB), and whether 82% partial sequences collapse into their full-length containers as expected.

## The structure

**Tree.** Sequences are the leaves of a binary tree whose internal nodes are virtual ancestors. That tree is the Kruskal merge tree of the spanning-tree build (Building the tree): each internal node is the merge of two single-linkage clusters, and the root is the virtual empty set. The spanning tree's edge weight between sequences *a* and *b* is |S\_a Δ S\_b|, the size of the symmetric difference of their sampled k-mer sets. Leaves are numbered in the merge tree's leaf order, so every internal node, every single-linkage cluster and every chunk is a contiguous leaf range.

*Amended 2026-10-24:* the first version made every node a sequence, with diffs against the parent sequence's set. That scheme is kept as the phase-0 baseline (Diff schemes).

**States and diffs.** Each k-mer *x* has a state (present or absent) at every node:

- Leaves take their sequences' stored sets.
- Internal nodes' states are chosen per k-mer by Fitch parsimony, to minimise changes.
- A posting is a change along an edge: (node, +) where *x* is gained, (node, −) where it is lost. A node's posting covers its leaf range.
- A leaf's set is the k-mers with a net gain along its root path. As an inverted index, each k-mer maps to the nested leaf ranges where it is gained or lost.

**Exact cost.** With internal states chosen by Fitch, a k-mer's postings are its Fitch parsimony score on the tree: the minimum number of changes, counting the gain from the empty root. So total postings = Σ\_x parsimony(*x*), the parsimony length of the k-mer presence/absence matrix. A flat per-sequence index costs Σ\_x |S\_x| instead; today's tier 1 costs Σ\_x |units(*x*)| before set deduplication.

*Corrected 2026-10-24:* the first version claimed this cost for a tree whose nodes are sequences. There every state is fixed, so the cost is Σ over edges of |S\_a Δ S\_b|, at least the parsimony score on the same topology and often more.

**When a diff pays.** With uniform substitutions at identity *p*, two sequences share ~*p*¹¹ of their 11-mers, so an edge costs ~2(1 − *p*¹¹) of the child's k-mers: 0.21 at 99%, 0.57 at 97%, 0.86 at 95%, 1.0 at ~94%, 1.37 at 90%. Below ~94% (at k 11), storing a sequence fresh is cheaper than a diff. Groups then merge only at the virtual empty root, as a forest of components under it, and k-mers shared across components are posted in each, as unit sets post them today. That cut falls out of the storage cost; it is not imposed as a clustering threshold.

**Collapses (generic, no threshold):**

- *Equal sets:* sequences whose sampled sets are equal share one leaf, which lists several sequences. Leaf count = distinct sampled sets, so resolution follows *t*.
- *Contained sequences:* a sequence whose k-mers are (almost) all inside another's, such as a fragment of a full-length protein, joins the container's leaf as a member. Giving it a leaf of its own would post every k-mer outside the fragment as a loss. In MGnify 4.7 × 10⁹ of 5.74 × 10⁹ proteins are partial, so this collapse sets the node table's size. Phase 0 measures its rate.

## Variable thresholds (FracMinHash by length)

Today *t\_g* = max(*t\_base*, *n\_min* / (*L\_g* − *k* + 1)) per unit, and the query keeps every k-mer with hash ≤ *t\_max*. The nesting carries over to nodes unchanged: S\_v = kmers(*v*) ∩ {h ≤ *t\_v*}, and a hit counts for *v* only if h ≤ *t\_v*.

- **Coordinated floors.** A threshold set by length is close to a bottom-*n* sketch. Near-identical neighbours of similar length pick nearly the same low-hash k-mers, so their floors share postings high in the tree. Floors therefore cost far less per node than *n\_min* fresh k-mers. Fragments would break this, but they collapse into their containers.
- **Threshold differences cost nothing.** A leaf whose threshold is below a k-mer's hash cannot hold it, which says nothing about its sequence. In Fitch that leaf is missing data ("?"), assigned whichever state needs fewer changes, so differing thresholds never force losses. Thresholds are quantised to a few levels (e.g. 8, log-spaced from *t\_base* to *t\_cap*). Each posting carries its level, and the query sweeps each level band separately: leaf *v*'s count = Σ over bands *j* ≤ level(*v*). So a leaf is never credited with a k-mer above its threshold. This is exact, and 8 sweeps over the touched ranges are cheap.
- **Where floors apply.** Today the floor skips singletons to fit memory. The tree's equivalent: the floor applies to nodes inside a forest component (an edge cheaper than a fresh set), not to isolated sequences. This is a cost rule, not an identity rule. Phase 0 sizes both choices.

## Building the tree from k-mers

The edge weight is defined on the stored k-mer sets themselves: storage cost, not similarity. So k-mer methods compute the objective exactly, rather than approximating sequence identity. Three facts make the build cheap.

- **Fresh storage is an edge to the empty set.** Add a virtual root whose set is empty: its edge to any *v* costs |S\_v|, which is exactly storing *v* fresh. The minimum spanning tree over every sequence plus that root is the storage-optimal tree under symmetric-difference cost, rooted at the empty set, with the forest cuts built in. (It is within 2× of the Steiner optimum, which may also invent ancestral sets.)
- **Only pairs sharing at least half the smaller set can be tree edges.** On the cycle *a*–*b*–root, edge (*a*, *b*) is excluded whenever |S\_a Δ S\_b| > max(|S\_a|, |S\_b|), i.e. when |S\_a ∩ S\_b| < min(|S\_a|, |S\_b|) / 2. Candidate generation is therefore an overlap set-similarity join with a fixed threshold, a problem with exact, fast algorithms. In particular, pairs sharing no stored k-mer never matter.
- **Weights are small integers.** With the floor, most stored sets hold ~16 k-mers (*t\_g* = 16 / (*L* − 10) binds below ~800 aa), so Δ ∈ {0, …, 32}. Kruskal over integer weights is a bucket queue plus union-find, linear in the candidate edges, and it can run level by level: Δ = 0, then 1, 2, …. Each level only needs the edges of that weight, so later levels see the components that earlier levels already merged, not raw sequences.

**Staged build (each stage exact for its weights, unless noted):**

1. **Hash once.** Every sequence's stored set S\_v (its own level), with today's Rust kernels, partitioned by sequence id. The same pass writes the inverted index (k-mer → sequences) partitioned by hash range, which is today's partitioned build with unit = sequence.
2. **Δ = 0: group by set hash.** Identical stored sets become one node. Linear (one sort or hash of 64-bit set digests).
3. **Δ ≤ d (d ≈ 3) by pigeonhole signatures** (PartEnum, Arasu et al. VLDB 2006). Split hash space into d + 1 parts by a second hash. If |S\_a Δ S\_b| ≤ d, the d differing k-mers fall in at most d parts, so a and b agree exactly on at least one part's restriction. Each node emits d + 1 (part, digest of S ∩ part) signatures, and an equi-join on signatures gives every pair within d, plus false positives removed by an exact Δ (a merge of two ~16-element sorted arrays, ~20 ns). Kruskal levels 1…d follow.
   - *Empty parts* are skipped (they match everything). A pair agreeing only on empty parts is missed here and found in stage 4.
   - *Large signature groups* (many nodes agreeing on one part) are not expanded into all pairs. Members are sorted by full digest and linked to their sort neighbours and to the group's longest member (Linclust's centre). This is approximate, and only for groups above a cap.
   - This stage covers the near-identical bulk: pairs above ~99% identity differ by about one stored k-mer.
4. **Δ > d: prefix-filtered overlap join** (Chaudhuri et al. ICDE 2006; AllPairs, Bayardo et al. WWW 2007; PPJoin, Xiao et al. WWW 2008), over one representative per component left by stage 3. Order k-mers globally by rarity, rarest first. A pair sharing at least half the smaller set must share a k-mer within the first ⌊|S|/2⌋ + 1 of the smaller set's sorted k-mers, so only those prefix postings are joined, and overlaps are counted only for the pairs they generate.
   - Rarity order is what makes this cheap. Unique k-mers (one sequence) fill prefixes and generate no pairs, and the frequent k-mers of large near-identical clades are already gone, because stage 3 merged those clades into single representatives. That is the reason for the staging.
   - Posting lists above a cap fall back to sorted-neighbourhood candidates: sort by the bottom few hashes in a few random orders and link window neighbours, which is MinHash LSH by sorting.
   - Kruskal levels d + 1 … 32 follow; the virtual root's edges close the forest.
5. **Optional refinement: NN-descent** (Dong, Charikar, Li, WWW 2011). Neighbours of neighbours are candidates; one or two rounds over the tree's own edges catch what the stage-3 caps and empty parts missed.
6. **Fitch states per k-mer** (on the merge tree from stages 2–5), by hash range from stage 1's inverted index:
   - A k-mer's holders, in leaf order, with their pairwise lowest common ancestors, form a small virtual tree of O(h) nodes for h holders (LCAs by range-minimum over the merge tree).
   - Subtrees with no holder hang off it as fixed "absent", and leaves whose level is below the k-mer's hash as missing data.
   - Fitch's bottom-up and top-down passes on that virtual tree give the states and the postings.
   - O(Σ h log h) over all k-mers, i.e. about linear in stored occurrences; parallel by hash range. The merge tree's MST edges are kept with their weights (N − 1 rows) for updates.

**Why it is fast:**

- Stages 1–3 are sorts and equi-joins over ~(d + 2) × 5.74 × 10⁹ 64-bit digests (~3 × 10¹⁰ records at d = 3, a few hundred GB), done by hash range as the partitioned build already does: hours of HPC, no alignment, and no pairwise comparison beyond the candidates.
- Stage 4 runs on the representatives stage 3 leaves, over prefix postings that are mostly length 1.

**Free by-products:**

- *A dendrogram.* Kruskal's merge order at levels Δ = 0, 1, 2, … is a single-linkage hierarchy, so a partition into clusters at any stored-k-mer distance is available without another clustering. This is useful for reporting and for benchmarks against fixed-identity units.
- *Containment.* Stage 4's overlap counts flag a set entirely inside another (overlap = |S\_a|) for the containment collapse.
- *Components.* Stage 4's join graph is also the source of cross-component links (relatives below the fresh-cost cut) for grouped routing and round 2.

**Thresholds in weights.** Weights only shape the topology; storage takes threshold differences as missing data (Variable thresholds). So an edge's weight counts the difference over k-mers both ends could hold: |(S\_a Δ S\_b) ∩ {h ≤ min(*t\_a*, *t\_b*)}|. This is symmetric, so the direction problem of the sequence-node scheme (and Edmonds' arborescence) does not arise.

**Two alternatives kept for comparison:**

- *Self-build by insertion.* Insert sequences longest first, in batches. Each batch is queried against the index built so far, with the query's own interval sweep, which gives |S\_new ∩ S\_v| for every touched leaf in one pass. Each sequence is attached beside the leaf minimising |S\_new| + |S\_v| − 2|∩| (as in Updates), or collapsed into it if contained. This is one code path for build and query, and order-dependent, like CD-HIT, but with the right objective. A good fallback if stages 3–4's joins are too large.
- *Dense tie-breaks.* Many nodes tie at the stored density. Ties can be broken by Δ at a denser nested rate (S at *t\_store* ⊂ S at *t\_dense*, computed only for tied candidates), so subtrees also group relatives sensibly for reporting. Storage cost is unchanged.

## Diff schemes

How a k-mer's presence is written on the tree (the user's question, 2026-10-24). Four schemes:

| Scheme | Tree | A k-mer is posted where it… | Compactness on a given tree | Build | Use |
| --- | --- | --- | --- | --- | --- |
| **Sequence MST** (the first version; phase-0 baseline) | nodes are sequences | is gained or lost against the parent sequence's set | Σ edge Δ, at least Fitch's; pays repeatedly on star-like groups: with a real centre, every member deletes the centre's private k-mers (members × privates), where a virtual consensus pays once | simplest: per-edge set differences, independent per edge, no per-k-mer pass | ± sweep; needs a second numbering for single-linkage clusters |
| **Add-only** (ProPhyle) | leaves are sequences, internal nodes virtual | covers a largest subtree whose leaves all hold it; no losses | worst for point mutations: a k-mer missing from one leaf of an n-leaf clade costs ~depth postings (one per sibling subtree along the path) against 2 with losses | per-k-mer pass over the holders' virtual tree; the simplest of the three | + only, a slightly simpler sweep, but more postings |
| **Majority** | the same | changes state, an internal node being present if most of its descendants hold it | close to Fitch, but wrong on near-even splits: weighting by descendant counts lets one large clade outvote several small ones | per-k-mer pass, the same machinery as Fitch | ± sweep |
| **Fitch parsimony** (the design) | the same | changes state, states chosen to minimise changes | the minimum on that tree: never worse than add-only or majority | per-k-mer pass (Building the tree, stage 6), the same cost as majority | ± sweep; exact CI; range postings |

**Why Fitch on the merge tree:**

- *Compactness:* the minimum on its tree, and virtual consensus ancestors fix the sequence MST's star-like groups.
- *Exact cost:* the plan's parsimony cost becomes true, and the consistency-index cut becomes exact.
- *One numbering:* in its leaf order every posting, single-linkage cluster and chunk is a contiguous range, so the neighbour rule and chunking need no second numbering.
- *No extra structure:* the merge tree is the Kruskal order the staged build already produces.

The merge tree's shape can be worse than the spanning tree's, so Fitch on it is not guaranteed below the sequence MST's Σ edge Δ. Phase 0 measures all four schemes on the same sample.

**Decision rule:** Fitch on the merge tree, unless the sequence MST is within ~10% of its postings and clearly cheaper to build. Even then, the single leaf order favours the merge tree.

## Uninformative k-mers

The tree gives each k-mer a natural score: its number of independent gains, compared with the one gain a clade marker needs. That is the phylogenetic consistency index (CI = minimum changes / observed changes). With Fitch states the change count is the true minimum, so CI is exact; schemes with fixed or heuristic internal states overstate changes.

- *Low CI* (many separate gains across unrelated subtrees): convergent or chance k-mers, low-complexity runs, compositional bias, common motifs. These are the k-mers that cost the most postings and carry the least information, so one cut on gains (or CI) removes both cost and noise. It replaces today's "held by more than 16 units".
- *One gain high in the tree:* a clade marker, kept however many sequences hold it (one posting). Today such k-mers in dense families are dropped by the 16-unit cut, so the tree keeps signal the current tier loses.
- *Gained once, then lost many times:* a variable-site k-mer. It weighs like today's low *p\_in*, which can be computed per subtree from the tree (the share of a subtree's sequences holding *x*).

## Querying a metagenome

1. **Lookup** (unchanged): sampled read k-mers (h ≤ *t\_max*, sketch-first and merge-join as today), each probing the k-mer table for its posting list.
2. **Counting:** each hit adds +1 over its gain ranges and −1 over its loss ranges, per level band. A sweep over the sorted range endpoints gives piecewise-constant hit counts over leaf order. The cost is O(P log P) for P hit postings, independent of the 10⁹ leaves; leaf set sizes |S\_v| and levels are read only for touched ranges.
3. **What to report (the user's three options):**
   - *Recall-first candidates (first to build):* every subtree with hits above a background, its members fetched as one contiguous range (the sequence store reordered by leaf order), into a dense sample-specific index or the pileup's reference. This is the main plan's round-2 candidate stage, with subtrees in place of components. It also gives grouped routing its similarity groups for free, and fixes the scattered fetch (members of one MGnify90 cluster spread over the whole id space).
   - *Greedy (gather):* take the node with the largest remaining hits normalised by |S\_v|, remove its k-mers' hits, repeat. Nodes tied within noise are reported as their smallest common subtree, like an LCA but ranked by evidence rather than by taxonomy.
   - *Placement-style model (later):* for each touched subtree, the best node by containment (hits / |S\_v| at its level) gives AAI as containment^(1/k), with an interval as today. The reported unit is the subtree of nodes whose containment is within that interval of the best: resolution adapts to the data, coarse for divergent strains and fine for close ones. Presence uses today's Poisson background, per reported subtree.
4. **Roll-ups:** each node's sequences carry their annotations (MGnify ships Pfam for every protein), so subtree → Pfam or KOfam roll-ups are direct. Mapping nodes to MGnify90 membership gives a like-for-like comparison with today's tier 1.

## Against today's tier 1

| | Tier 1 today (MGnify90 units) | Tree index |
| --- | --- | --- |
| Units | 1.66 × 10⁹ fixed 90% clusters | leaves = distinct sampled sets, after collapses (phase 0: expected ~1–2.5 × 10⁹); reported units = subtrees per sample |
| Postings | Σ\_x \|units(x)\|, set-deduplicated; k-mers in > 16 units dropped | Σ\_x parsimony(x) ≈ per-unit unions + within-unit deletions − cross-unit sharing; low-CI k-mers dropped. Expected within ~0.5–2× (phase 0) |
| Clade-wide k-mers | dropped above 16 units | one posting |
| Lookup | one probe per sampled k-mer | the same |
| Hit expansion | posting set → units | interval sweep per level band |
| Resolution | fixed 90% | follows the data and *t* |
| AAI | per unit, through *p\_in* | to the best node (nearest sequence), plus per-subtree *p\_in* |
| Build | from MGnify's membership (free) | the staged k-mer build (sorts and joins) plus a Fitch pass per k-mer: hours of HPC, no clustering of the input |
| Generic sets | needs a clustering made for the set | any set; parameters k, *t*, a CI cut |
| Fetch for round 2 / pileup | scattered over the store | one contiguous range per subtree |

It may be more efficient at query time: fewer, shorter posting lists for clade k-mers, and no per-unit expansion of large sets. It is less efficient to build. Index size has to be measured.

## Density

- **No dense global tier.** *Measured 2026-10-30 (phase 0): dense distinct over total is 0.14 over the whole release, not the ~0.5 assumed below, so ~1.5 × 10¹¹ distinct dense 11-mers rather than ~5 × 10¹¹; the TB figures below shrink by ~3.5×. A dense global tier is still ~10× the standard tier.* Every distinct k-mer needs at least one gain, so the tree's postings are at least the distinct k-mer count, ~5 × 10¹¹ dense 11-mers if about half of MGnify's ~10¹² occurrences are distinct (the main plan's guess). With deletions that is ~5–8 TB, ~2× smaller than a flat dense index but ~50× the standard tier.
- **Tier 1 can be somewhat denser at the same size.** Today's postings are ~1.7 units per k-mer; the tree posts a clade's k-mer once. Deletions, the node table and the k-mers today's 16-unit cut drops take part of that back, so the room is ~1.2–2× the sampling rate at matched size (phase 0). The query does not get denser: it already probes every hash ≤ *t\_max* = 0.2, so per-node rates under the cap add hits, not lookups.
- **Resolution follows density.** At tier 1's rates (effectively ~0.06, with floors), near-identical members share nodes. Dense, almost every sequence is its own node (MGnify is 100%-deduplicated), and only contained fragments collapse. Per-sequence resolution and presence on a candidate's own k-mers therefore still need a dense second tier.
- **One topology, many densities.** Sampling is uniform, so a tree near-optimal at one rate is near-optimal at another. Dense, near-dense and other-k diffs can be computed along the edges of the tree built once at tier 1's rate (plus dense tie-breaks), as a linear pass over edges with no new spanning tree. Phase 0 checks this against a tree built at the dense rate.

## Uses beyond tier 1 (open question)

The same structure may serve the per-sample dense stages, which today each build a flat table. Each use has to be checked for **build speed** and **use speed** against what it replaces, not only for size.

| Use | Today | With the tree | Build must be | Use must be |
| --- | --- | --- | --- | --- |
| **Tier 2 / round-2 dense presence** (main plan, round-2 candidate stage) | dense flat table over candidates' members | the candidate subtrees' members (contiguous leaf ranges), Fitch states at the dense rate on tier 1's merge tree; per-leaf hits by interval sweep | minutes per sample: fetch + dense hashing + a Fitch pass over the candidates' k-mers; no spanning tree | ≤ today's dense pass |
| **Routing** (`kfp-pileup route`; grouped routing) | flat table of every k 8 k-mer of the reference's members → block | the same merge tree at k 8 (all six frames on the read side); blocks and similarity groups = subtrees (leaf ranges); per-(pair, subtree) counts by sweep | ≤ today's table build (~1 min per group at `hadza10M`) | ≤ today's pass (~1.5 min at `hadza10M`) |
| **Read reclassification** (per-read candidate units for an aligner; main plan amendment 2026-10-21) | routing table with unit values | per-read hits per node, each read's candidates = the nodes or smallest subtree its k-mers support | as routing | per-read sweeps on few postings: to measure |

Two layouts for the dense diffs, to compare:

- **Per sample:** fetch the candidate subtrees' members and compute their Fitch states on the stored merge tree. Nothing extra on disk; the cost is hashing and the Fitch pass per sample.
- **Precomputed, node-major:** dense postings for every node stored on disk in leaf order (~5–8 TB), with candidate ranges read with `pread` and inverted into an in-memory k-mer table, as sylph's stage 2 does. No hashing per sample; it costs disk and one sequential read per range.

Where it should help most: candidates are clades of near-identical members, so most of a flat dense table over them is repeated k-mers. For `hadza10M`'s ~8 M hit units with all members (~10¹⁰ residues), a flat dense table is ~10¹⁰ postings, and the tree should be 3–5× fewer (estimate). Where it may not help: k 8 shares more k-mers across relatives (0.9⁸ ≈ 0.43 against 0.31 at k 11), so routing's diffs are larger, and per-read sweeps add work that a flat lookup does not have.

## Chunks by similarity

An index too large to handle at once, such as a dense tier 2, is split along the tree into chunks of similar proteins, and each chunk is queried in turn (the user's idea, 2026-10-24).

- **Chunks come from the merge tree.** Every merge-tree subtree is a contiguous leaf range, so a chunk is one contiguous block on disk, and its members one contiguous range of the store.
  - Top-down from the root: a subtree with ≤ *B* leaves (by stored size) is a chunk; a larger one is split into its two children. Chunks adjacent in leaf order are then packed up to *B*.
  - The highest merges are the heaviest spanning-tree edges, so cuts fall where relatives diverge most and near-identical competitors stay together. The virtual root's children (the forest components) are free boundaries. No separate tree-partitioning step (Kundu & Misra) is needed, and no chunk-major renumbering.
- **A chunk decodes alone.** A chunk's root state is stored in full: one set per chunk.
- **Cost:** the k-mers shared across each cut edge, posted in both chunks. This is zero at forest cuts, and small when *B* is much larger than a typical component (phase 0 measures duplication against *B*).

**Querying chunk by chunk:**

1. The sample is sketched once at the dense rate: distinct hashes and counts, sorted, on disk. Each chunk pass is a sequential merge-join of that sketch against the chunk's sorted k-mers, with no read pass per chunk. Read-level outputs (reclassification, routing) read the stored translated frames (main plan: `kmer-functional-profiler translate`) instead.
2. Only chunks holding candidate subtrees are loaded. This is the main plan's grouped routing ("built as needed"), with the tree giving the groups.
3. Each pass gives per-node hit counts, complete for its nodes, since each node is in one chunk. Merged counts go to gather, EM or presence, so results are identical to an unchunked run: chunking changes memory, not results. Models that need a read's evidence across chunks keep (read, node, hits) rows per chunk and merge them at the end.

Memory is one chunk's table plus one range of the sample sketch, both set by *B*.

*Precedent:* MiniPhy groups genomes into batches of up to ~10k by their phylogeny, and Phylign searches one batch at a time (phylogenetic compression, below); HIBF lays bins out by similarity.

## Recall-first candidates and neighbours

Tier 1's job before tier 2 is not to miss: every unit the sample holds should reach the dense stage, even when tier 1 hit only a relative of it (the user's request, 2026-10-24).

**The property.** In a minimum spanning tree, the largest edge on the path between two nodes is no larger than their direct distance (the minimax path property). So the single-linkage cluster at level *R* containing *v* (everything reachable from *v* by edges ≤ *R*) holds every node within *R* of *v*. The Kruskal merge tree is that hierarchy, so the cluster at any level is a lookup (the highest ancestor with merge level ≤ *R*). Leaves are numbered in its leaf order (The structure), so each cluster is one contiguous range.

**The rule.** For each hit node *v*, take its single-linkage cluster at *R\_v* = *d̂⁺*(X, *v*) + δ:

- *d̂⁺*(X, *v*) is the upper bound of the sample protein X's divergence from *v*, from *v*'s containment interval, in the tree's distance units.
- δ is the tolerance: how close a sequence must be to X to count as one not to miss.
- **Guarantee:** by the triangle inequality, d(*v*, *w*) ≤ d(*v*, X) + d(X, *w*) ≤ *R\_v*, so every indexed *w* within δ of X is a candidate whenever any relative of it is hit and its bound holds. It holds for any distance that is a metric; Jaccard distance is one.
- *It adapts on its own:* a close hit expands little, a divergent one more, and a hit with few k-mers (a wide interval) more again, the right direction for recall.
- *Budget:* a cluster over a cap is cut nearest-first in minimax distance from *v* (the order in which Prim's algorithm from *v* visits nodes), so chaining drops the farthest neighbours, not arbitrary ones.
- *Built from a candidate graph rather than all pairs,* the guarantee holds wherever the candidate graph has some chain of edges ≤ *R* between the two nodes (not necessarily their direct edge); phase 0 measures how often it does not.

**It needs a topology that reaches past the storage cuts.** The storage tree (k 11, tier 1's rate) splits below ~94% identity, so on it the rule never reaches family-level neighbours. There are two ways to reach further, compared in the next section: a sparse topology from tier 1's own k-mers with the overlap requirement relaxed to any shared k-mer (reach ~85%), or a dense sensitive topology (reach ~60–70%). Fitch states are computed on whichever merge tree is built.

**What it cannot fix.** If no relative of X is hit, no neighbour rule helps. The main plan's deep misses (90–95% strains of units with 1–5 sampled k-mers) need tier-1 density or sensitivity, for which the tree's ~1.2–2× room is the lever.

**The alternative:** candidates as the credible region of X's placement on the tree (pplacer and RAPPAS style). It is better founded statistically, but it is more code and gives no clean guarantee.

## A sensitive topology at MGnify scale

**Separating seeds from weights** is what makes it feasible. Candidate pairs need selective seeds; edge weights need sensitive, dense k-mers; neither needs all pairs.

1. **Nodes.** The storage collapses first (equal and contained sets at k 11, from the staged build), then one representative per storage forest component (sets ≥ ~94% apart). At least ~1.7 × 10⁹, since MGnify90's 1.2 × 10⁹ singletons are separate components by definition; perhaps ~2 × 10⁹ (phase 0). The storage tree resolves everything inside a component; the sensitive topology links components.
2. **Candidate edges** as Linclust (which clusters 1.6 × 10⁹ sequences to 50% in ~10 h on one server): m ≈ 20 lowest-hash seeds per node on a reduced alphabet (k ~ 10–12), grouped by seed. Each node is linked to its group's centre and its sort neighbours: ~20–40 edges per node, ~4–8 × 10¹⁰ edges. Sorting ~10¹¹ records (~1 TB) by hash range is hours of HPC.
3. **Weights** are Jaccard distance on dense short-k sets (k 6: 70% identity still shares 0.7⁶ ≈ 12% of 6-mers; chance overlap between unrelated 180-mers is ~5 × 10⁻⁴ of 6-mers). They are estimated from a per-node sketch: FracMinHash at *t* ≈ 0.25 (~45 hashes) or bottom-128, ~0.2–0.5 KB per node, ~0.5–1 TB for 2 × 10⁹ nodes. Weighting is two streaming joins of the edges (sorted by each endpoint) against the sketch table; compute is negligible (an edge is a merge of two short sorted arrays), I/O is ~2 sequential TB passes.
4. **Kruskal** over the weighted edges: an external sort by weight (~1 TB) and a union–find over 2 × 10⁹ nodes (8 GB). Hours; or a distributed Borůvka.
5. **Storage on it.** k 11 Fitch states are computed on the sensitive topology's merge tree. Within close clades the two trees nearly coincide, so storage should cost little more than on the storage-optimal tree (phase 0).

**Estimate:** a few hundred to ~1,000 core-hours and a few TB of scratch, comparable to MGnify's own cascaded clustering and to today's full index build. HPC only. An exact dense MST (all pairs) is not feasible and not needed. *Not yet known:* the node count after collapses, the largest seed groups, and how much the sketch's noise in weights at low Jaccard (≤ 0.05) costs the guarantee.

**Against a sparse topology from tier 1** (2026-10-24, estimates, scaled against today's full index build, ~1,300 core-h: ~20 h of one node):

| | Sparse (tier 1's k 11 sets, ~16 k-mers per protein; edges: any shared stored k-mer) | Dense (reduced-alphabet seeds, k 6 sketch weights) |
| --- | --- | --- |
| Hashing | free (the tier-1 build's) | a seed pass over ~2 × 10⁹ representatives (~4 × 10¹¹ residues): a few core-h |
| Candidate edges | digests and pigeonhole signatures (~3 × 10¹⁰ records, ~250 GB), and the prefix join relaxed to one shared k-mer on representatives | ~10¹¹ seed records (~1 TB) sorted; group centres and sort neighbours |
| Weights | merges of the stored ~16-element sets: negligible | a sketch table (0.5–1 TB), two streaming joins |
| Kruskal | ~3 × 10¹⁰ edges | ~4–8 × 10¹⁰ edges, ~1 TB external sort |
| **Total** | **~50–300 core-h, ~0.3–0.5 TB scratch** | **~300–1,000 core-h, ~2–4 TB scratch** |
| Against today's build | +5–20% | +25–75% |
| Storage edges (≥ ~94%) | exact (this is the storage tree) | the same |
| Neighbour reach | ~85%: ~16 k-mers share 1.4 at 80% identity on average and 0.3 at 70% (16 × *p*¹¹), so most 70% relatives share none | ~60–70%: at 70% about 12% of 6-mers are shared (~21 of 175; ~5–8 in the sketch) |
| Weight noise | coarse: Jaccard in steps of ~1/16 | finer: steps of ~1/45 to 1/128 |
| Neighbour guarantee | δ up to ~0.1–0.15, noisy there | δ up to ~0.3 |
| Relation to the main plan | in effect its index-level neighbour graph (recall-first pileup, next step 0) | new |

Neither changes query time: the topology is a build cost only. The two are not exclusive. The sparse storage tree is built anyway, so the decision is whether to spend ~250–700 core-h more for neighbours between ~85% and ~60–70%.

**Decision rule (proposed):** build the dense topology only if it reaches substantially more of tier 1's misses than the sparse one, for its extra cost. Specifically, only if a substantial share of `presence-calibration`'s spiked units with no own hit have their nearest hit relative below ~85% (phase 0). The main plan's deep misses were mostly 90–95% strains of units with 1–5 sampled k-mers, inside the sparse reach, so sparse is the default until phase 0 says otherwise.

## Updates

Adding proteins needs no re-clustering, unlike today's tier 1, whose units follow MGnify90 membership. The difficulty is keeping the leaf order contiguous, which a delta index handles (the user's question, 2026-10-24).

- **One protein** (the build's insertion path): query it against the index; the interval sweep gives its overlap with every touched leaf. Then:
  - an equal or contained stored set joins that leaf as a member (nothing else changes);
  - otherwise it becomes a new leaf beside the leaf *u* minimising |S\_new| + |S\_u| − 2|∩|. A new internal node splits *u*'s edge and takes *u*'s state, so the new leaf stores S\_new Δ S\_u as its own gains and losses. If no leaf beats its own size, it becomes a fresh component under the root.
- **Delta index (LSM-style).** The main index stays frozen. New leaves go to a small delta: each has its attachment (a main or delta leaf) and its diffs, in a small k-mer table.
  - Queries stay exact. The main sweep runs as now. A delta leaf's count is its attachment leaf's count (a lookup in the sweep's breakpoints), minus hits on its losses, plus hits on its gains, per level band. Delta leaves are evaluated in topological order, so chains work.
- **Compaction** folds the delta in: renumber leaves, recompute Fitch states for the k-mers the new leaves touch, rewrite postings in one sequential merge, split chunks above *B*.
- **Staying near-optimal.** Attaching leaves one at a time drifts from the MST, since a new protein may be a better neighbour for existing leaves. This costs storage, not correctness, and the batch fix is exact:
  - the MST of the enlarged set uses only the old MST's edges (kept with the merge tree) plus the new proteins' candidate edges, so a batch update is Kruskal over those, linear in their size, and gives the new merge tree;
  - new candidate edges come from the staged joins (digests, pigeonhole signatures, prefix join). Existing leaves' signatures are recomputed from their sets, reconstructed from the postings along each leaf's root path, without re-hashing their sequences.
- **Refreshed at compaction:** per-k-mer gain counts (kept current, with the uninformative cut applied here, since a k-mer can cross it), per-subtree *p\_in*, chunk boundaries, the background calibration. Removed proteins are tombstoned until then.
- **A new MGnify release** is then: hash only the new proteins; join them against existing nodes; Kruskal over old edges plus new ones; one compaction. Small user additions (a lab's own proteins) go to the delta with no rebuild.

## Prior art

The closest precedents, in order of closeness:

- **Mantis-MST** (Almodaresi, Pandey, Ferdman, Johnson, Patro; RECOMB 2019, JCB 2020): colour classes of a coloured de Bruijn graph stored as deltas along a minimum spanning tree, using the graph to find near neighbours (>11× smaller than RRR on 10,000 samples). The same MST-of-deltas encoding, over sample colour vectors rather than sequences. [RECOMB](https://link.springer.com/10.1007/978-3-030-17083-7_1), [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6612864)
- **RowDiff** (Danciu, Karasikov, Mustafa, Kahles, Rätsch; Bioinformatics 2021): annotations stored as diffs between adjacent graph vertices, 30% smaller than Mantis-MST on 10,000 RNA-seq sets. Diffs along a topology, built in linear time. [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8346655/)
- **ProPhyle** (Břinda et al.): k-mers propagated bottom-up on a reference tree, each stored at the highest node whose descendants all share it, then queried through a BWT index; lossless, so accurate on close genomes. This is the tree index with additions only, on a given tree. [PyPI](https://pypi.org/project/prophyle/0.3.1.0/)
- **Phylogenetic compression** (Břinda et al., Nature Methods 2025; MiniPhy): ordering genomes by an estimated tree compresses assemblies, de Bruijn graphs and k-mer indexes 10–100×, and makes search over all sequenced bacteria possible on a desktop. Evidence that tree order is what makes such collections small. [GitHub](https://github.com/karel-brinda/miniphy), [PDF](https://baymlab.hms.harvard.edu/pdfs/BrindaEtAlnmeth.pdf)
- **panCAKES / CAKES / CLAM** (Prior, Howard, Light, Ishaq, Daniels; IEEE Big Data 2024): a divisive hierarchical clustering tree that both compresses (each point as differences from a centre) and searches without full decompression, for any distance whose value tracks encoding cost (edit distance, Jaccard). The same "compression structure = search index" idea, as a generic metric method. [arXiv](https://arxiv.org/abs/2409.12161v2)
- **CaBLASTP** (Daniels, Gallant, Peng, Cowen, Baym, Berger; ISMB 2013): protein databases compressed to a coarse set of unique segments plus links, searched coarse-first, at a cost that scales with unique data. [Bioinformatics](https://academic.oup.com/bioinformatics/article/29/13/i283/189542)
- **Sequence Bloom trees and HowDe-SBT** (Solomon & Kingsford 2016; Harris & Medvedev 2020): a tree built by clustering filters, internal nodes summarising descendants, split into "determined/how" bits; queries descend. The reverse layout (top-down membership tests); at 10⁹ leaves an inverted index with intervals scales better. [Bioinformatics](https://academic.oup.com/bioinformatics/article/36/3/721/5553093), [GitHub](https://github.com/medvedevgroup/HowDeSBT)
- **Fulgor and Mac-dBG** (Fan, Pibiri, Patro et al., WABI 2023 onwards): colour sets factored across the collection ("meta-colours", NP-hard to optimise, solved heuristically) and repetition-aware colour compression (2.5–22× smaller). The best current practice for the posting-list side. [WABI 2023](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.WABI.2023.18), [Mac-dBG](https://par.nsf.gov//servlets/purl/10546249), [WABI 2025](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.WABI.2025.6)
- **HIBF / Raptor** (Mehringer et al., Genome Biology 2023): a similarity-based hierarchical layout of bins, built 211× faster and up to 14× smaller. Hierarchy by similarity for membership queries. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10230713/)
- **Linclust** (Steinegger & Söding, Nature Communications 2018): linear-time grouping by m lowest-hash k-mers, comparing each sequence only to its group's centre (1.6 × 10⁹ sequences in 10 h on one server). The candidate-edge step for the spanning tree. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6026198)
- **Placement and adaptive resolution:** Kraken's LCA (Wood & Salzberg 2014), the lossy single-node version of a posting list; RAPPAS phylo-k-mers and App-SpaM, alignment-free read placement on a tree ([bioRxiv](https://www.biorxiv.org/content/10.1101/328740v5), [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.10.19.344986v2.full.pdf)); TreeCluster (Balaban et al., PLoS ONE 2019), clusters chosen on a tree by criteria rather than by one global threshold ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6705769)).
- **Kaiju** (Menzel, Ng, Krogh, Nature Communications 2016): protein-level classification against a whole reference protein set through an FM index, with LCA for ties. Indexing every sequence works at NR scale. At MGnify's 1.07 × 10¹² residues a full-text index is ~TB, which is what sampled k-mer diffs avoid. [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4833860/)

*What I found no precedent for:* this combination of a k-mer diff tree over every sequence of a protein set, nested FracMinHash thresholds per node, a CI-based k-mer cut, and abundance or presence at adaptive subtree resolution for metagenome profiling. Each part has precedent, but the combination appears new.

## Benchmark scoring without fixed units

Without MGnify90 units, "the right unit" is no longer defined, so accuracy is scored by distance to the truth rather than by unit identity (the user's question, 2026-10-24).

- **Primary: threshold-free nearest-match curves.**
  - *Recall side:* for each true spiked sequence, the identity of the nearest reported sequence, as a distribution by depth and by strain divergence.
  - *Precision side:* for each reported unit (leaf or subtree), the identity of its nearest true sequence in the sample.
  - This mirrors how phylogenetic placement is scored (node distance to the true placement in EPA-ng and RAPPAS benchmarks). It separates near misses ("missed at 99%, within 95%") from failures, which fixed thresholds hide.
- **Threshold read-outs:** the shares at ≥ 99%, 95%, 90% (and 80%) identity on both curves, which give recall and precision per level. They are measured directly against the truth rather than through clusters at those levels, so no pair is split by a cluster boundary.
- **Abundance:** earth mover's distance between true and estimated abundance, with 1 − AAI between sequences as the ground distance (weighted UniFrac on a tree, as CAMI's OPAL does on taxonomies). One number per sample; credit given to a 98% neighbour costs little, to a 60% relative a lot.
- **Distances independent of the method.** The tree is built from the method's own k-mers, so its distances would favour it. Identities are computed by alignment (DIAMOND) on the sequences involved (truth plus reported, small per sample). `presence-calibration`'s and `spike-benchmark`'s `truth.tsv` already hold each strain's coding sequence (`cds`), its source `member` and its `identity`, so no new simulation is needed.
- **Today's tier 1 scored the same way.** Each reported MGnify90 unit is expanded to its member sequences and scored by the same curves, a like-for-like comparison. The existing roll-ups (MGnify90 via membership, Pfam, KOfam) are kept for continuity with past results.
- **Pitfall:** in dense families a reported sequence can be near the truth by chance, so precision inflates at loose levels (90%). The whole curve, not only the 90% read-out, is reported.

## Risks

- **Postings above budget.** If MGnify's members are less redundant than assumed (members at 90–94% to their nearest neighbour post more than a fresh set), the tree degenerates into today's per-unit storage plus a node table. Phase 0 decides.
- **Leaf and node tables.** One row per leaf (level, |S\_v|, length, first member) and one per internal node (leaf range, merge level; N − 1 rows). At ~2 × 10⁹ leaves and ~12 B per row each that is ~48 GB, memory-mapped and read only for touched ranges; ~140 GB if partials do not collapse. Postings that carry their leaf range directly (two leaf ids) would remove the internal-node table at ~4 B more per posting. Phase 0 sizes both.
- **Ties at the stored density.** The stored sets define the storage objective exactly (Building the tree), but with ~16 k-mers per set many candidate parents tie, and the tree's topology among relatives is then arbitrary. That does not matter for storage, but it does for reported subtrees, so ties are broken at a denser nested rate.
- **Single-linkage chaining.** A spanning tree chains: a subtree can span distant ends through intermediates. This matters for what a reported subtree means (a function), not for storage. The reporting model has to bound a subtree's diameter, as TreeCluster does, rather than take whole subtrees.
- **Assignment model.** Gather and EM over nested subtrees are new code. Today's component-wise EM does not carry over directly.
- **Build at scale.** The staged k-mer build (Building the tree) is sorts and joins over ~3 × 10¹⁰ digests plus a prefix join on stage 3's representatives. Its unknowns are the size of stage 3's largest signature groups and stage 4's posting lengths after merging, both measured in phase 0. MGnify90 membership is not needed; it is kept only as a check (does the tree recover the clusters?).

## Implementation plan

Laptop work runs one job at a time with capped threads; full-scale runs go to HPC as Nextflow (SLURM + Singularity) for the user to run.

0. **Measure, no index (gate for everything after).** *Run 2026-10-30: **gate failed** on postings (3.8× tier 1 with per-sequence floors; 0.64× without floors); leaves 1.98× units (passes). Stopped here pending a decision on floors; see Progress log. Survey 2026-10-30: with floors per piece and a per-leaf minimum, **passes** on the full sample (1.11× tier 1, every protein with a valid 11-mer covered).* On a sample of MGnify: a few hundred thousand MGnify90 clusters with all their members, plus the representatives that share k-mers with them across clusters (to see cross-cluster sharing). Laptop on a biome subset, HPC for the full-scale sample. At *t\_base* 0.02, *n\_min* 16, k 11, with today's hashing (`crates/core` `protein_kmers`):
   - spanning tree under |S\_a Δ S\_b| with the virtual empty root: exact (all candidate pairs from the inverted index, Kruskal) at this size, and the staged build (stages 1–4) beside it. Report the staged tree's total weight against the exact one, each stage's candidates, largest signature group, and prefix posting lengths;
   - equal-set and containment collapses;
   - *diff schemes* (Diff schemes): postings and build time for the sequence MST, add-only, majority and Fitch on the same merge tree, on the same sample. **Decision:** Fitch on the merge tree unless the sequence MST is within ~10% of its postings and clearly cheaper to build;
   - nodes against sequences and units; the containment-collapse rate of partials;
   - Σ parsimony (postings) against today's tier-1 postings for the same proteins (built by today's `index`), with and without floors, and with floors on components only;
   - CI distribution, and postings saved by cuts at 2, 4 and 8 gains, against the 16-unit cut;
   - edge-cost histogram by identity (does the ~94% cut appear?);
   - distinct over total k-mers at *t* 0.06 (tier 1's effective rate), 0.2 and 1, which sets every density estimate above;
   - *one topology, many densities:* diffs at *t* 1 and at k 8 along the tree built at tier 1's rate, against a tree built at those settings;
   - *sensitive topology:* on the sample, the exact MST at k 6 Jaccard (all candidate pairs from the dense inverted index) against the seeded build (reduced-alphabet seeds, sketch weights): total weight, the share of exact edges recovered, and the share of node pairs within *R* (0.05–0.3) that lack a chain of edges ≤ *R* in the seeded candidate graph; node counts after collapses; storage cost of k 11 diffs along the sensitive topology against the storage-optimal tree;
   - *sparse against dense topology* (A sensitive topology at MGnify scale): both built on the sample, with build CPU-h and scratch for each. On `presence-calibration`'s spiked units with no own hit (their profiles and `truth.tsv`): the identity of each miss's nearest hit relative, and the share each topology's neighbour rule reaches at δ 0.05, 0.1 and 0.15, with the candidate count. **Decision:** dense only if it reaches substantially more misses than sparse, for its extra cost;
   - *chunking:* duplicated postings against chunk size *B* (at *t* 1 and at tier 1's rate), the share of forest components larger than *B*, and how many chunks the spike benchmark's candidate subtrees touch;
   - **gate:** postings ≤ 1.5× today's, and nodes ≤ 2× units, at matched *t*. If it fails, stop here and record why.
1. **Build prototype** (Python, Rust kernels as the existing tier, on the laptop scale of phase 0). *Started 2026-10-30: `workflows/tree-index/build.py` builds one tier with the survey's design on the exact topology, checked against brute force; see Progress log.*
   - the staged build: set digests (Δ = 0), pigeonhole signatures (Δ ≤ d), the rarity-ordered prefix join on representatives (Δ > d), Kruskal by integer level with the virtual empty root, and dense tie-breaks; collapses;
   - the merge tree and leaf order; Fitch states per k-mer (stage 6, missing data for leaves below a k-mer's level); range postings with level bands; leaf and node tables; a sequence store reordered by leaf order.
   - Checked against a brute-force reconstruction: every sequence's sampled set from its root path equals its own; and against brute-force Fitch on small trees (equal minimum change counts).
   - *Updates:* a held-out batch (e.g. 10% of phase 0's sample) added through the delta index. Query counts must be identical to a rebuild containing the batch, and identical again after compaction. The batched Kruskal update's total weight must equal a from-scratch tree's, under equal candidate edges. Report the drift (weight) of leaf-only insertion without the batched update.
2. **Query: counting.** Lookup as today (sketch, merge-join), interval sweep per band, node counts and |S\_v|. Checked against brute-force per-sequence containment on fixtures. *Done 2026-10-30 (`workflows/tree-index/query.py`); see Progress log.*
3. **Reporting, recall first.** Candidate subtrees above background, members as leaf ranges into the main plan's round-2 dense stage and the pileup's `REFERENCE`. Then greedy gather with LCA ties, and per-node containment AAI. *Done 2026-10-30 (`workflows/tree-index/report.py`), checked locally; see Progress log. Calibration and the δ sweep are phase 4's.*
4. **Benchmark against tier 1**, at matched index size and *t*:
   - `fmh-benchmark` and `presence-calibration` (spiked truth), scored as in Benchmark scoring without fixed units:
     - nearest-match identity curves for recall and precision, with read-outs at 99, 95, 90 and 80%;
     - earth mover's distance for abundance, with alignment identities (DIAMOND) as ground distances;
     - the same for today's tier 1, its units expanded to members;
     - and the roll-ups to MGnify90 (via membership), Pfam and KOfam;
   - AAI calibration against `presence-calibration`'s truth;
   - *neighbour rule* (Recall-first candidates and neighbours): δ at 0.05, 0.1 and 0.15. Report the share of spiked units with no own hit that reach the candidate set, by depth and identity; the candidate count against "every hit node"; and how often the true divergence exceeds *d̂⁺*. Against the main plan's index-level neighbour graph (recall-first pileup, next step 0), by that step's rule: worth it if it reaches a substantial share of unhit spikes for < ~2× the candidates;
   - cost: index size, build CPU-h, query wall and peak RSS per stage;
   - **gate (proposed):** recall at matched precision ≥ tier 1's at every read-out level (99–80%) and every roll-up level, for ≤ 1.2× query cost and ≤ 1.5× index size. Below that, the tree stays a round-2 candidate generator only.
4b. **Tier 2, routing and reclassification** (Uses beyond tier 1): on `spike-benchmark` locally, then `pileup-fdr`'s `hadza10M` candidate set on HPC. Each use is built from the candidate subtrees along tier 1's topology, per sample and node-major precomputed, and compared with the flat table it replaces:
   - build time and peak RSS;
   - use time (dense pass, routing pass, per-read candidates);
   - table size;
   - and the same truth metrics: recall and precision for tier 2; for routing, `route.json`'s routed and dropped pairs and the spike benchmark's `read_precision_em` and `sibling_caught_em`.
   - *chunked:* each use run chunk by chunk (Chunks by similarity) at two budgets *B*. Outputs must be identical to the unchunked run; report peak RSS and wall against *B*;
   - **Gate per use:** build and use each ≤ 1.2× today's at equal truth metrics; a use that fails keeps its flat table.
5. **Full scale** (HPC, the user's run): the staged build and Fitch pass on all of MGnify (MGnify90 membership only as a check), then the phase-4 benchmark at full size.
6. **Generic sets:** a set with no published clustering (e.g. CARD resistance proteins, or viral proteins) built and queried with the same defaults, to show nothing in it depends on MGnify's clustering.
7. **Later:** the placement-style model with adaptive subtree resolution and bounded diameter; posting-list compression in the Fulgor/Mac-dBG manner.

## Progress log

* **Plan (2026-10-24, the user's request).** Written from the main plan's tier-1 design (*t\_g* floors, *p\_in*, unit sets, 16-unit promiscuity cut, standard tier ~108–196 GB) and a literature search. The estimates above are analytic: 11-mer sharing ~*p*¹¹ under uniform substitutions, 22 changed k-mers per substitution, MGnify 2026\_07 counts. No data measured yet; phase 0 is next.

* **Amendment (2026-10-24, the user's point: k-mer sharing is the objective itself, not a proxy for similarity).** Building the tree added: a virtual empty root makes fresh storage an edge, so the MST over sequences plus that root is the storage-optimal tree. Only pairs sharing ≥ half the smaller set can be edges (an overlap set-similarity join). Integer weights allow Kruskal level by level. Stages: set digests (Δ = 0), PartEnum pigeonhole signatures (Δ ≤ d), a rarity-ordered prefix join on the remaining representatives (Δ > d), optional NN-descent; self-build by insertion kept as the alternative. This replaces the Linclust-plus-Borůvka sketch, and MGnify90 membership is no longer needed to seed it.

* **Amendment (2026-10-24, the user's question: is the tree also useful for tier 2, routing or reclassification, and is it fast enough to build and use?).** Density added: no dense global tier (postings ≥ distinct k-mers, ~5–8 TB dense); tier 1 perhaps ~1.2–2× denser at matched size; tier 2 still needed for resolution. Uses beyond tier 1 added: tier 2 / round-2 presence, routing and read reclassification built per sample from the candidate subtrees along tier 1's topology ("one topology, many densities"), per sample or node-major precomputed. Phase 0 now measures distinct-over-total k-mers and the reuse of one topology across densities; phase 4b gates each use on build and use time against its flat table.

* **Amendment (2026-10-24, the user's idea: split an index too large for one go into similarity chunks from the tree and query them in turn).** Chunks by similarity added: whole forest components packed into chunks of size ≤ *B*, larger components split at their heaviest edges (Kundu & Misra), chunk roots stored in full, chunk-major numbering. Querying merge-joins one stored sample sketch against each needed chunk, and merges per-node counts before the model, so results are identical to an unchunked run. Phase 0 measures duplication against *B*; phase 4b runs each use chunked, with identical outputs required.

* **Amendment (2026-10-24, the user's question: can the index be updated with new proteins efficiently?).** Updates added: insertion by the index's own query; an LSM-style delta index (new nodes with parent and diffs, counts from the parent's count plus own diffs, exact); compaction (renumber, one merge, chunk splits, gain counts and the uninformative cut refreshed); an exact batched MST update by Kruskal over the old tree's edges plus new candidate edges; a release update that hashes only new proteins. Phase 1 checks delta queries and compaction against a rebuild, and the batched update against a from-scratch tree.

* **Amendment (2026-10-24, the user's request: err on recall in tier 1's output for tier 2, with a clean neighbour rule; then, is a dense topology feasible?).** Recall-first candidates and neighbours added: each hit node's single-linkage cluster at *d̂⁺*(X, *v*) + δ, with a guarantee from the MST minimax property and the triangle inequality; capped in Prim order. A sensitive topology added, since the k 11 storage tree splits below ~94%: Linclust-style reduced-alphabet seeds for candidate edges, k 6 Jaccard weights from per-node sketches, Kruskal, and storage diffs along it. Estimated at a few hundred to ~1,000 core-h and a few TB of scratch on HPC. Phase 0 compares it with the exact MST on the sample; phase 4 sweeps δ.

* **Amendment (2026-10-24, the user's request: compare a dense topology's cost with tier 1's sparse one).** Comparison added (estimates against today's ~1,300 core-h build):
  - *Sparse* (tier 1's k 11 sets, edges from any shared stored k-mer): ~50–300 core-h and ~0.3–0.5 TB scratch (+5–20%); reach ~85%; coarse weights; guarantee to δ ~0.1–0.15. In effect the main plan's index-level neighbour graph.
  - *Dense* (reduced-alphabet seeds, k 6 sketch weights): ~300–1,000 core-h and ~2–4 TB (+25–75%); reach ~60–70%; guarantee to δ ~0.3.
  - Neither costs query time.
  - *Decision rule:* sparse by default, since the main plan's deep misses (90–95% strains) are inside its reach; dense only if phase 0 shows a substantial share of hit-less spikes whose nearest hit relative is below ~85%, reached by dense and not sparse, for its extra cost. Phase 0 now builds both on the sample and measures this.

* **Amendment (2026-10-24, the user's question: how does the diff work at each node, and should it be majority-based or add-only?).** The design changes from sequence nodes to the Kruskal merge tree, with sequences as leaves, virtual internal nodes, and Fitch states per k-mer:
  - *Fitch over the others:* majority is a heuristic for Fitch at the same build cost, and Fitch is never worse than majority or add-only on a given tree. Add-only pays ~depth postings per k-mer for each point mutation.
  - *Corrected:* the first version's "Exact cost" (postings = parsimony length) was wrong for sequence nodes, whose states are fixed: their cost is Σ edge Δ. It is right for Fitch states on virtual nodes.
  - *Changed with it:* range postings over one leaf order, which also orders single-linkage clusters and chunks; threshold differences as Fitch missing data, with the per-band sweep kept; edge weights over k-mers both ends could hold (symmetric, so no Edmonds); stage 6 (Fitch per k-mer from the inverted index, O(Σ h log h)); chunking top-down on the merge tree (no Kundu & Misra); updates by attaching new leaves beside their nearest leaf, with Fitch recomputed at compaction; leaf and node tables (~48 GB).
  - Diff schemes added, with a decision in phase 0: Fitch on the merge tree unless the sequence MST is within ~10% of its postings and clearly cheaper to build. The sequence MST is kept as the baseline.

* **Amendment (2026-10-25, the user's question: how to benchmark without 90% units, e.g. at several AAI thresholds?).** Benchmark scoring without fixed units added:
  - Threshold-free nearest-match identity curves (recall: truth to nearest reported; precision: reported to nearest truth), with read-outs at 99, 95, 90 and 80% as the threshold figures.
  - Earth mover's distance for abundance, with 1 − AAI as the ground distance.
  - Identities by alignment (DIAMOND) on the sequences involved, independent of the method's k-mers, from `truth.tsv`'s `cds`.
  - Today's tier 1 scored the same way (units expanded to members); the MGnify90, Pfam and KOfam roll-ups kept.
  - Phase 4's gate now applies at every read-out level.

* **Phase 0 set up (2026-10-30, the user's request: start phase 0, with the whole-release streaming counts).** Code on `tree-index-phase0`; not run at full scale.
  - *What runs where.* Everything that needs a tree runs on the sample: 1 in 5,000 MGnify90 clusters (~330 k, ~1.2 M proteins), taken from release-tiers' 256 members buckets. The sample's proteins are also built as today's standard tier, to compare against. One streaming pass over all of MGnify counts what changes with N without building a tree: equal-set groups, pigeonhole signature groups, and k-mer holders in three hash slices (stored and all k-mers). That pass also collects the representatives outside the sample that share ≥ 4 stored k-mers with it (at most 10⁶, sharing most first), for `measure-all`. Workflow `workflows/tree-index`; run `hpc/kfp-ablations/tree-index-phase0`.
  - *Kernels* (`crates/core/src/tree.rs`):
    - candidate pairs from an inverted index, exact and uncapped;
    - Kruskal with the virtual empty root, and the merge forest;
    - Fitch, majority and add-only postings per k-mer, on the virtual tree of its holders, with leaves below the k-mer's level as missing data. Checked against a full-tree brute force on 3,000 random forests;
    - overlap-alignment identity for the edge-cost histogram.
  - *Choices made in setting it up:*
    - **Containment** is exact on the common hash range: every k-mer the smaller set could share is in the larger one. A contained set joins its closest container. The plan's "almost all inside" is not implemented; add a tolerance if the measured rate is low because of single k-mer differences.
    - **The node gate counts leaves.** Internal nodes are reported beside it (`nodes_per_unit`). Counting both would put every tree near 2× its leaves, so the gate would mean nothing.
    - **The sequence-MST baseline** counts each spanning-tree edge's whole set difference, plus root edges. Each node's state must equal its stored set, so a parent's k-mers above the child's threshold are deletions too. The threshold-restricted weight that builds the tree is a lower bound, reported as `spanning_tree_weight`; with one level the two are equal (tested).
    - **The staged build's stage 4** probes each representative's rarity-ordered prefix against every representative's postings, kept when the prober is the smaller set. Postings over 1,000 holders fall back to sorted-neighbourhood candidates (4 orders, window 2). Stage 3 enumerates groups up to 64 and links larger ones to sort neighbours and their largest member. The exact tree is the reference for both.
    - **Sensitive-topology weights** are in divergence units, 1 − (2J / (1 + J))^(1/6), an identity estimate from k 6 Jaccard J. The order is still Jaccard's, so the trees are unchanged, and the radii *R* 0.05–0.3 mean identity distance. Raw Jaccard distance at *R* ≤ 0.3 would only reach pairs above ~95% identity.
    - **Distinct over total at uniform rates** is the same at any *t* in expectation, since sampling is uniform. The sample reports 0.06, 0.2 and 1 directly. The stream estimates it from narrow hash slices (width 1/1024) at 0, 0.06 and 0.2, for all k-mers and for stored ones, with exact holder counts; HyperLogLog is not needed. Pigeonhole signatures are counted in full (`--sig-every 1`), as the plan asked for the largest group.
  - *Not covered* (pending): the sparse-against-dense reach on `presence-calibration`'s misses, which needs a topology over the full index's units (phase 5); how many chunks the spike benchmark's candidate subtrees touch; NN-descent (stage 5); dense tie-breaks.
  - *Fixture only, not evidence:* on synthetic families (60 roots, members at 80–100%, half with fragments), Fitch < add-only < majority < sequence MST with floors (9,686, 10,169, 10,953, 11,831 postings). The postings ratio against tier 1 is 2.9 with floors and 0.93 without, because the fixture's many short fragments raise the floors. The edge-cost curve crosses 1 near 92% identity. The real numbers come from the HPC run.
  - *Amended 2026-10-30 (HPC run: REDUCE killed, exit 137).* One `group_by` held a whole hash range of the release at once, ~2 × 10⁹ rows at 16 ranges. REDUCE is now `workflows/tree-index/reduce.py`. It splits each bucket's file by key into `--reduce_passes` (64) sub-partitions in the task directory, then sums them one at a time. Memory is about one sub-partition (~3 × 10⁷ rows); scratch disk is about the range's input, freed as it goes. Outputs match the old REDUCE on the fixture, apart from the order of tied largest groups. It is a separate script so `phase0.py` is unchanged and `-resume` keeps SAMPLE and STREAM.

* **Phase 0 results (2026-10-30, the user's HPC run of `hpc/kfp-ablations/tree-index-phase0`; REDUCE re-run out of core after exit 137).** Sample: 1 in 5,000 clusters, 308,997 units, 1,067,329 proteins. Tier 1: the standard tier (`tbase02_all16`) built on the same proteins, 1,727,926 postings. MEASURE took 92 s and 8.7 GB on the sample (382 s and 24 GB with the representatives); STREAM about 7–9 min and up to 14 GB per bucket.
  - **Gate: failed on postings.** Fitch with per-sequence floors: 6,600,680 postings, **3.82×** tier 1 (limit 1.5×). Floors on forest components only: 2.23×. Without floors: **0.64×**. Leaves: 612,857, **1.98×** units (limit 2×, passes narrowly); internal nodes 293,064.
  - **Why: the floors, not the tree.** By hash band, k-mers at or below *t\_base* (0.02) cost 1,252,535 postings for 988,458 distinct k-mers. That is already below tier 1's total. K-mers above it, held only because of per-sequence floors, cost 5,348,145 postings for 4,552,777 distinct k-mers (1.7 holders each).
    - Every distinct k-mer needs at least one gain, so no tree can share these away.
    - Each sequence's floor samples its own ~16 low-hash k-mers. Members of a 90% cluster share only ~p¹¹ of their 11-mers, so most floor k-mers are new.
    - Tier 1 floors per unit instead: 16 k-mers per floored non-singleton unit (77,242 of them), chosen by *p\_in* / *n\_groups*, so mostly consensus k-mers, plus the base stratum.
    - The design's "coordinated floors share postings high in the tree" (Variable thresholds) does not hold: near-identical neighbours do share floor k-mers, but most members are not near-identical at 11-mers.
    - The comparison is therefore not at matched *t*: tier 1's floors are per cluster and selected, the tree's per sequence and complete.
  - **Without floors the tree is smaller than tier 1, but loses short proteins.** 1,106,270 postings (0.64×), but 166,704 proteins (15.6%) have no k-mer at or below 0.02 and so no leaf. Tier 1's singletons (226,659 of its units) are also unfloored, so it loses many of the same proteins; that overlap was not measured.
  - **Diff schemes:** Fitch 6.60 M; add-only 1.065×; majority 1.11×; sequence MST 1.16×. Fitch's pass took 0.78 s against 0.06 s for the sequence MST, both negligible. **Decision: Fitch on the merge tree.**
  - **Staged build:** total weight 1.019× the exact tree's, and Fitch on it 0.4% above. Stage 3: 514,307 candidates, 88% of them false positives; largest signature group 215 in the sample, 35,155 over the whole release. Stage 4: 572,613 representatives, 3.04 M prefix candidates, no capped posting in the sample (88 with the representatives). The exact join on the sample took 3 s (3.3 × 10⁸ probes), so phase 1 needs the staged build only at full scale.
  - **Collapses:** 2,073 empty sets. Equal sets collapse 13.5% of proteins in the sample and 15.3% over the whole release (largest group 4,314). Contained sets: 308,057 of 920,914 distinct sets. **Only 42.7% of partial proteins are contained** (29.2% of full-length ones), far below the "partials collapse into their containers" the node estimate assumed. Strict containment on the common hash range may be too strict (a tolerance is the open option); the node gate still passes.
  - **Uninformative k-mers:** dropping k-mers with more than 2, 4 or 8 gains saves 9.8%, 5.1% and 2.4% of postings. The 16-unit cut drops 88 k-mers in a 1-in-5,000 sample, which says nothing about full scale.
  - **Edge cost by identity:** pairs sharing a stored k-mer have median cost 1 (as costly as storing fresh) at 89–90% identity, not the ~94% estimated from *p*¹¹. Spanning-tree edges stay below 1 down to 80%.
  - **Density:** distinct over total is 0.26 on the sample at every uniform rate (0.06, 0.2, 1), as expected for uniform sampling. Over the whole release it is 0.14 for all k-mers and 0.34 for stored ones in the sample, 0.19 and 0.37 in the whole release's 0.06 and 0.2 slices. The largest dense holder count in a slice is 494,986 proteins. Stored k-mers over the whole release: 9.15 × 10¹⁰, ~16 per protein.
  - **One topology, many densities:** on 49,994 proteins at *t* 1, Fitch along the stored-rate tree is 6.1% above a tree built at that rate (k 11), and 7.1% above at k 8. **Passes:** dense and routing diffs can reuse the stored topology.
  - **Chunking:** at the stored rate, duplicated postings are 2.5% at a budget of 10³ stored k-mers, 0.12% at 10⁴, and 0 from 10⁵ (no component is larger). At *t* 1: 64%, 18%, 2.3% and 0 at 10³–10⁶. Budgets of ≥ 10⁵ dense k-mers per chunk are cheap.
  - **Sensitive topology: inconclusive at this sample.** It ran on 11,431 component representatives from random clusters, which rarely have relatives in the sample: only 9–809 pairs fall within *R* 0.05–0.3.
    - The exact k 6 MST is mostly chance edges (205 components, against 9,026 for the seeded graph); the seeded candidates hold only 9% of its edges.
    - Within *R*, candidate generation loses at most 2.3% of pairs (at *R* 0.3) with exact weights, but sketch weights lose 6–33%. The noise of the *t* 0.25 sketch, not the seeds, is the weak point.
    - k 11 storage along the sensitive topology (263,239) is no worse than along the storage tree (264,553).
    - A sample built around families (whole Pfam families, or `presence-calibration`'s misses with their relatives) is needed to decide it.
  - **Cross-cluster sharing** (`measure-all`; the 10⁶ representatives cap was reached): 55% of spanning-tree edges link different MGnify90 clusters, against 1.7% in the sample alone. Fitch 8.74 M (Fitch without floors 1.85 M) over 1,308,997 units.
  - **Interpretation.** The tree's encoding works as designed:
    - Fitch is the cheapest scheme;
    - the staged build is near-exact;
    - one topology serves other densities;
    - chunking is cheap;
    - without floors it is 0.64× tier 1.

    What fails is per-sequence floors. They multiply distinct k-mers by ~4.6 (988 k to 5.5 M), and no topology can share those. By the phase-0 rule this stops the plan as written. Continuing needs a floor design that is not per sequence, measured against tier 1 at matched sensitivity: for example, floors per subtree (forest component or single-linkage cluster) chosen by consensus as tier 1's are, or a lower per-sequence *n\_min*. That is the user's decision (pending).

* **Forest floors explored (2026-10-30, the user's idea: a forest rather than one tree, so the hierarchy that forces per-sequence floors is not kept).** Local exploration, not a workflow yet. Run on a quarter of phase 0's sample: clusters with `cluster_rep % 20000 == 0`, 77,317 units, 262,391 proteins. On this quarter, tier 1 has 432,440 postings and per-sequence floors 1,642,628 (3.80×, as on the whole sample).
  - **Design tried.**
    - The topology is phase 0's exact storage tree, built on floored sets.
    - Its merge tree is cut into pieces at merge level > *c* (stored-k-mer weight); *c* = ∞ keeps the forest components.
    - Leaves keep only base-band k-mers (≤ *t\_base*), diffed by Fitch inside each piece.
    - Each piece gets one floor in tier 1's manner: *t\_p* = max(*t\_base*, min(*t\_cap*, *n\_min* / distinct k-mers of its union)), keeping the *n\_min* candidates above *t\_base* held by most members. These are posted once at the piece's root.
    - Singleton pieces get no floor, as tier 1's singleton units get none.
  - **Results** (postings over tier 1's; blind = proteins whose piece has no k-mer at all):

    | Cut *c* | Pieces | With ≥ 2 members | Base Fitch | Piece floors | Total | / tier 1 | Blind proteins |
    | --- | --- | --- | --- | --- | --- | --- | --- |
    | 0 | 151,163 | 30,365 | 497,145 | 301,946 | 799,091 | 1.85× | 22,463 |
    | 2 | 145,785 | 29,682 | 473,031 | 294,928 | 767,959 | 1.78× | 22,126 |
    | 4 | 131,199 | 26,829 | 418,021 | 267,960 | 685,981 | 1.59× | 21,062 |
    | 6 | 117,315 | 24,383 | 377,654 | 243,841 | 621,495 | 1.44× | 19,885 |
    | 8 | 106,298 | 22,733 | 352,125 | 225,687 | 577,812 | 1.34× | 19,018 |
    | 16 | 84,565 | 20,452 | 315,185 | 196,217 | 511,402 | 1.18× | 17,213 |
    | ∞ (components) | 80,880 | 20,647 | 309,853 | 196,385 | 506,238 | **1.17×** | 16,936 |
    | tier 1 | 77,317 units | 20,648 | | | 432,440 | 1× | 16,075 |

  - **Reading.**
    - With floors on pieces rather than sequences, the tree comes within 1.17× of tier 1 at almost the same blindness (16,936 against 16,075 proteins), inside the phase-0 limit of 1.5×. It also keeps per-sequence resolution in the base band, which tier 1 does not have.
    - **Cutting finer does not help.** Every extra piece needs its own floor and loses base-band sharing across the cut, so cost falls monotonically as cuts coarsen. The useful forest is the storage forest (components), with floors at component roots.
    - The fine cuts (≤ 8) buy finer floor resolution at 1.3–1.9×.
    - Flooring singleton pieces too (each with its own *n\_min* candidates) costs 562,540 more postings at components (2.47× in all) and leaves no protein blind.
  - **Caveats.**
    - Random 1-in-20,000 clusters rarely have relatives in the sample, so components here are close to clusters. At full scale, cross-cluster links (55% of spanning-tree edges in `measure-all`) make fewer, larger components: fewer floors, but coarser ones.
    - Floor candidates are ranked by holders only; tier 1 also divides by *n\_groups*.
    - Measured on a quarter sample; the workflow does not yet compute it.
  - *Next (pending the user's decision):* add piece floors as a variant of `measure` (cut ladder, with and without singleton floors) and re-run phase 0 on HPC, including `measure-all`, to re-check the gate at full sample size.
  - **Both: the full hierarchy kept, floors lower down (the user, same day: "the hierarchy beyond this might still be useful for neighbours but not for the diffs"; "can I have both?").** Postings and topology are independent.
    - The whole merge tree and its single leaf order are kept: the neighbour rule's single-linkage clusters at any level, chunks, and reported subtrees are all contiguous ranges of it.
    - The sparse or sensitive linking of components above the storage forest is kept too.
    - Only *where postings sit* is restricted. Base-band k-mers are diffed by Fitch over the whole tree; above a component no base-band k-mer is shared anyway, so this costs the same as the forest. Floor k-mers are posted at piece roots at cut *c*.
    - The two bands hold disjoint k-mers (≤ *t\_base* and > *t\_base*), so the costs add. A floor k-mer chosen in several pieces is counted once per piece, an upper bound.
    - Merge nodes above the cuts carry no postings, only their N − 1 rows in the node table, which exists anyway.

    | Floor cut *c* | Floor pieces | Base Fitch (whole tree) | Piece floors | Total | / tier 1 | Blind proteins |
    | --- | --- | --- | --- | --- | --- | --- |
    | 0 | 151,163 | 309,853 | 301,946 | 611,799 | 1.41× | 22,463 |
    | 4 | 131,199 | 309,853 | 267,960 | 577,813 | 1.34× | 21,062 |
    | 8 | 106,298 | 309,853 | 225,687 | 535,540 | 1.24× | 19,018 |
    | 16 | 84,565 | 309,853 | 196,217 | 506,070 | 1.17× | 17,213 |
    | ∞ (components) | 80,880 | 309,853 | 196,385 | 506,238 | 1.17× | 16,936 |

    Every cut now passes 1.5×: the cost of finer floors is only the floors themselves. The cut sets floor resolution (and blindness: finer pieces of one member leave short proteins unfloored), not diff cost or neighbour reach.
  - **Pieces by a metric, not a fixed cut (the user, same day: "can the point the diffs start be based on some metric (e.g. some information metric) rather than an arbitrary threshold").**
    - Postings alone always choose the coarsest cut (components), so the metric has to value resolution. The natural one is **floor coherence**: the mean share of a piece's floor k-mers that a member holds (the floor's mean *p\_in*), about the chance that a member's reads hit its floor.
    - Rule: top-down from each component, a node whose floor is at least *τ* coherent is a piece; otherwise its children are tried. Pieces grow where members agree and split where they do not, so the cut level adapts per clade.
    - *τ* can be set from data rather than chosen: today's tier-1 units have a floor coherence of 0.62 (mean; median 0.625) by the same measure, over 19,187 units with ≥ 2 members on this quarter sample. Quarter sample, base band 309,853 postings over the whole tree:

    | *τ* | Pieces | Floored | Floor postings | Total | / tier 1 | Floored pieces' mean coherence | Blind proteins |
    | --- | --- | --- | --- | --- | --- | --- | --- |
    | 0.3 | 83,471 | 20,284 | 211,822 | 521,675 | 1.21× | 0.66 | 17,127 |
    | 0.5 | 90,761 | 23,067 | 245,779 | 555,632 | 1.28× | 0.70 | 17,666 |
    | 0.6 | 97,378 | 24,357 | 263,447 | 573,300 | 1.33× | 0.73 | 18,404 |
    | 0.7 | 104,986 | 24,250 | 263,742 | 573,595 | 1.33× | 0.77 | 19,629 |
    | 0.9 | 117,064 | 21,183 | 231,862 | 541,715 | 1.25× | 0.80 | 21,630 |
    | tier 1 | 77,317 units | 19,187 | | 432,440 | 1× | 0.62 | 16,075 |

    - **Reading.** At *τ* = 0.3, floors are already more coherent than tier 1's (0.66 against 0.62), for 1.21× its postings, with about tier 1's blindness. That beats the fixed cut at level 8 (1.24×, 19,018 blind).
    - Raising *τ* buys coherence (0.80 at *τ* 0.9) but splits more clades into single members, which get no floor, so blindness rises. Above 0.7 floor postings fall for that reason.
    - Matching tier 1's coherence (*τ* ≈ 0.6) costs 1.33×. All are within the 1.5× gate.
    - **Alternatives considered:**
      - Minimum description length (MDL), counting each member that lacks a floor k-mer as a correction posting. With no weight on the corrections it is the postings-only answer (components); with full weight it is per-sequence floors again (3.8×). It needs a weight, which is *τ* in other units.
      - A detection-probability rule: *τ* set so that a member's reads give ≥ *m* floor hits with probability *P* at a given depth. Same rule, with *τ* derived from the query's needs.
    - Floor candidates here are ranked by holders only, and the union's k-mer count is estimated from hashes ≤ *t\_cap*.
  - **Proteins without a k-mer of their own (the user, same day: "that should never be"; would they still reach tier 2 through a higher node?).** On the quarter sample, 40,199 of 262,391 proteins (15.3%) have no 11-mer at or below *t\_base* 0.02. By how the hierarchy reaches them:
    - 18,589 sit in a storage component with other leaves. A hit on a relative reaches them by the neighbour rule (its single-linkage cluster at *R\_v*), so a recall-first tier 1 sends them to the sample's dense tier 2 whenever a relative is hit within reach.
    - 4,174 share a leaf (equal or contained set), so they come with that leaf.
    - 16,936 are alone: their own component and their own leaf. The storage hierarchy cannot reach them. Linking components above the storage forest (the sparse or sensitive topology, built from more k-mers than are stored, at no posting cost) would, if they have any relative; otherwise only an own k-mer can.
    - 500 have no leaf at all: no k-mer at or below even *t\_cap* 0.2. Of those, 38 have no valid 11-mer (too short or ambiguous).
    - Tier 1 today has the same gap: 16,075 proteins on this quarter sample. On the whole sample, 65,004 of 308,997 units (21%) have no posting.
  - **Guaranteeing own k-mers** (a per-protein bottom-*m* floor: each protein's *m* lowest hashes, coordinated across near-identical proteins). Extra stored k-mers on the quarter sample, an upper bound on their postings:

    | Own k-mers ≥ | Extra | / tier 1 (432,440) | Proteins needing it |
    | --- | --- | --- | --- |
    | 1 | 40,161 | +9% | 40,161 |
    | 2 | 128,107 | +30% | 87,966 |
    | 4 | 424,518 | +98% | 165,987 |
    | 8 | 1,294,562 | +299% | 236,950 |

    **≥ 1 own k-mer costs +9%** and leaves only the 38 proteins with no valid 11-mer. Those need a different route (a shorter *k*, or an exact-sequence table). The hierarchy route, together with ≥ 1–2 own k-mers, gives every protein some way into tier 2.
  - **No unfloored pieces (the user, same day: single-member exceptions make sense for clusters, not for a tree).**
    - The exception was carried over from tier 1, where singleton units go unfloored to fit memory. It is dropped: every piece gets a floor, single-leaf pieces included.
    - A single-leaf piece has no one to share its floor with, so each floor k-mer is a fresh posting. The floor size for single leaves (*n\_1*) is therefore the lever.
    - Quarter sample, total postings over tier 1. Blind proteins become 0 in every column (only the 38 with no valid 11-mer remain):

    | *τ* | Multi-leaf pieces only (blind) | *n\_1* = 16 | 8 | 4 | 2 | 1 |
    | --- | --- | --- | --- | --- | --- | --- |
    | 0.3 | 1.21× (17,127) | 2.55× | 2.17× | 1.74× | **1.49×** | **1.35×** |
    | 0.6 | 1.33× (18,404) | 2.93× | 2.46× | 1.95× | 1.65× | **1.49×** |
    | 0.9 | 1.25× (21,630) | 3.46× | 2.80× | 2.09× | 1.68× | 1.47× |

    - Within the 1.5× gate, every protein is covered at *τ* 0.3 with *n\_1* ≤ 2, or at *τ* 0.6 with *n\_1* = 1.
    - These are upper bounds: a single leaf's floor here is up to *n\_1* candidates above *t\_base* even when it has base k-mers already. A top-up rule (floor only to reach *n\_1* stored k-mers in all) is cheaper.
    - The rule becomes uniform in tree terms: a floor of up to *n\_min* k-mers per piece, chosen by coherence, with a minimum of *n\_1* stored k-mers per leaf.
  - **Survey set up for HPC (2026-10-30, the user's request).**
    - `workflows/tree-index/survey.py` implements the design above in the workflow: SURVEY on core and all, TIER1_MATCHED for tier 1 at *t\_base* 0.01 and *n\_min* 8, and INDEX_SIZES for the release's standard-tier file sizes. It reuses the run's cached tasks, since `phase0.py` is unchanged.
    - **Grid:** *t\_base* {0.01, 0.02} × *n\_min* {8, 16} × *τ* {0, 0.3, 0.6, 0.9} × *n\_single* {0, 2, 4, 16} × *n\_1* {0, 1, 2, 4}.
    - **Two differences from the exploration scripts:**
      - Coherence and floors are now per leaf (one representative per leaf) rather than per protein, so pieces of one leaf with several identical proteins count as single-leaf pieces.
      - Proteins without a leaf (no k-mer up to *t\_cap*) get a leaf of their own under the minimum.
    - **On the quarter sample** (17 s, 2.1 GB): at components with no single-leaf floor and *n\_1* = 1, **1.10× the standard tier** with only the 38 proteins that have no valid 11-mer uncovered. *n\_1* = 2: 1.19×; *n\_single* = 2 with *n\_1* = 1: 1.35×.
  - **Survey results (HPC, 2026-10-30).** `tree-index-phase0.tar.gz`: SURVEY on core (1,067,329 proteins, 308,997 units) and all (+10⁶ cross-cluster representatives), TIER1_MATCHED, INDEX_SIZES. SURVEY took 3.3 min / 5.6 GB (core) and 6 min / 21 GB (all).
    - **Core confirms the quarter sample.** At *t\_base* 0.02, *n\_min* 16, no single-leaf floor (*n\_single* 0), totals over the standard tier (1,727,926):

      | *τ* | *n\_1* = 0 (blind) | 1 | 2 | 4 |
      | --- | --- | --- | --- | --- |
      | 0 (components) | 1.06× (91,752) | **1.11×** | 1.20× | 1.49× |
      | 0.3 | 1.08× (90,434) | **1.12×** | 1.21× | 1.50× |
      | 0.6 | 1.11× (97,041) | 1.16× | 1.26× | 1.56× |
      | 0.9 | 0.85× (127,500) | 0.91× | 1.03× | 1.42× |

      With *n\_1* ≥ 1 only 172 proteins stay uncovered (no valid 11-mer). Every *τ* ≤ 0.6 with *n\_1* ≤ 2 passes the 1.5× gate.
    - *τ* 0.9 is cheaper only because it splits clades into single leaves, which get no floor at *n\_single* 0; the top-up covers them with fewer k-mers. It trades floors for top-ups rather than buying resolution.
    - **Single-leaf floors are the expensive part.** *n\_single* 2 / 4 / 16 cost 1.36× / 1.63× / 2.47× (τ 0, *n\_1* 1). Only *n\_single* ≤ 2 passes.
    - **Lower rate.** At *t\_base* 0.01, *n\_min* 8 the tree is 1.23× tier 1 built at the same rate (869,021; *τ* 0, *n\_1* 1), and 0.62× the standard tier. Whether that rate is sensitive enough is phase 4's question.
    - **All (cross-cluster sharing).** The extra 10⁶ representatives add 0.90 M base postings (+72%) and 41 k floor postings for +94% proteins. Its ratios in `survey.tsv` (1.65× at the settings above) are **not comparable**: the denominator is tier 1 on core only. A matched check needs tier 1 built on core + extra.
    - **Release scale.** Projected linearly by protein count (5,378× from core): 1.03 × 10¹⁰ postings (τ 0, *n\_1* 1) against 9.3 × 10⁹ for the standard tier on the same scaling; from all (2,776×) 7.9 × 10⁹. Linear scaling ignores the growth of sharing with depth, so these are upper bounds.
    - **Gap:** `index_sizes.tsv` has no tier-1 file (INDEX_SIZES lists one directory level only), so the release's tier-1 bytes are still unknown.
    - **Reading.** The phase-0 gate passes on the full sample with floors per piece (components or *τ* 0.3), no single-leaf floor, and *n\_1* = 1: ~1.11× tier 1 with every valid protein covered. Pending the user's decision: proceed to phase 1 on these settings, or first re-run with a tier 1 on core + extra for the matched `all` check.
  - **Tree tiers proposed (2026-10-30, the user's request: ~64, ~128 and ~256 GB; quick profile, good profile followed by the dense tier, and a performance boost).** Sizes use the standard tier's all-in compact cost, ~11.6 B per posting (~108 GB for ~9.3 × 10⁹ postings, unit table included). Postings are projected from all (lower bound) and from core (upper bound). Every tier has a floor on each multi-leaf piece and a per-leaf minimum *n\_1*, so no protein with a valid 11-mer is blind; none has single-leaf floors (*n\_single* 0), since those cost the most.

    | Tier | *t\_base* | *n\_min* | *τ* | *n\_1* | Postings (release) | ~Size | Why |
    | --- | --- | --- | --- | --- | --- | --- | --- |
    | quick (~64 GB) | 0.01 | 8 | 0 (components) | 1 | 4.3–5.7 × 10⁹ | 50–67 GB | slim's rates; cheapest floors; one own k-mer per leaf |
    | medium (~128 GB) | 0.02 | 16 | 0.3 | 2 | 8.4–11.3 × 10⁹ | 98–130 GB | the standard tier's rates; floors more coherent than tier 1's for +1%; two own k-mers so a leaf can reach `min_hits` 2 before the dense tier |
    | large (~256 GB) | 0.05 | 16 | 0.3 | 4 | ~1.5–1.9 × 10¹⁰ (estimate) | ~175–220 GB | the large tier's rate, which calibrates `aai`; floors rarely bind at 0.05, so a higher *n\_1* is cheap |

    - The large tier is not measured: the survey stopped at *t\_base* 0.02. Base postings doubled from 0.01 to 0.02 (623 k to 1.25 M on core), so ~3.1 M at 0.05, with floors and top-ups shrinking. Add *t\_base* 0.05 to the survey before building it.
    - **Leaf table.** At ~2.5–3.3 × 10⁹ leaves the leaf table is ~30–40 GB at ~12 B per row, against tier 1's ~16 GB for 5.3 × 10⁸ posted units. That may push quick past 64 GB. If it does, quick's leaf rows can drop to leaf boundaries and level band (~8 B).
* **Phase 1 begun: a tier builder, checked locally (2026-10-30, the user's request: build the tiers at 1 in 100 on HPC, but do the local development steps first).**
  - **Kernel.** `tree_schemes(..., emit=True)` returns Fitch's changes (hash, node, gain or loss) of the assignment it counts. The Rust test checks, over 3,000 random forests with missing data, that they number exactly Fitch's postings and rebuild every leaf that could hold each k-mer.
  - **Kruskal** now breaks weight ties by node pair rather than input order. `overlap_pairs` emits pairs in thread order, so builds differed slightly from run to run: ±0.01% postings against the survey's HPC run. The build is now deterministic.
  - **`build.py`**, format `tree-v0`:
    - Range postings `(hash, lo, hi, sign)` over leaf positions, sorted by hash: base-band Fitch changes (kind 0), floors on multi-leaf pieces at *τ* (kind 1), and per-leaf top-ups to *n\_1* (kind 2). A leaf holds a k-mer iff the signs covering it sum to 1.
    - Proteins without a leaf get positions after the tree's.
    - Beside the postings: the merge tree as position ranges, `pieces.parquet` and `proteins.parquet`.
    - No single-leaf floors (*n\_single* 0, as every proposed tier).
  - **Test** (`tests/python/test_tree_build.py`, fixture):
    - posting counts per kind equal the survey's for the same configuration;
    - expanding the ranges, every cover sums to 0 or 1;
    - each leaf's base band equals its representative's k-mers;
    - every floor is held across its piece;
    - every protein with a valid k-mer is covered by at least min(*n\_1*, its k-mers).
  - **The three tiers on phase 0's core sample** (1,067,329 proteins; laptop, 4 threads):

    | Tier | Postings | Base / floor / top-up | / tier 1 | Build | Peak RSS | `tree-v0` bytes |
    | --- | --- | --- | --- | --- | --- | --- |
    | quick (0.01, 8, *τ* 0, *n\_1* 1) | 1,069,817 | 624 k / 298 k / 148 k | 1.23× matched (0.62× standard) | 73 s | 4.6 GB | 33 MB |
    | medium (0.02, 16, *τ* 0.3, *n\_1* 2) | 2,093,942 | 1,250 k / 611 k / 233 k | 1.21× | 61 s | 4.9 GB | 52 MB |
    | large (0.05, 16, *τ* 0.3, *n\_1* 4) | 3,775,732 | 3,134 k / 288 k / 353 k | 2.19× standard | 61 s | 4.7 GB | 82 MB |

    - The large tier is now measured: its base band (3.13 M) is as estimated, floors shrink to 288 k, and top-ups are 353 k.
    - At tier 1's all-in compact ~11.6 B per posting, scaled from core: quick ~67 GB, medium ~131 GB, large ~236 GB. These are upper bounds, since sharing grows with scale.
    - `tree-v0` is a development format (18 B per posting plus 9.6 MB of `proteins.parquet` here), not the compact one.
  - *Next:* phase 2 (counting query: lookup, interval sweep, node counts) locally against brute force on these builds. Then a BUILD step for the 1-in-100 subset on HPC (~57 M proteins, ~50× this sample; the exact topology may need the staged build at that size).
* **Phase 2: counting query (2026-10-30, the user's request).**
  - **Kernel.** `_core.contained_sums(lo, hi, w, q_lo, q_hi)` sums, for each query range, the weights of the posting ranges that contain it: one sweep in `lo` order with a Fenwick tree over `hi`, O((P + Q) log P). A Rust test checks it against brute force on random ranges.
  - **Build** now also writes |S\_v|: `pos_size.npy` per leaf position and `node_size.npy` per internal node (Fitch's states), both from that kernel over all postings.
  - **`query.py`:**
    - *Lookup.* Reads (six-frame translated, as today) or proteins are hashed. A hash is kept if it is ≤ *t\_cap* or one of the few keys beyond it, then counted per distinct hash and probed in the sorted key array.
    - *Counting.* Leaves come from one sweep over the hit postings' endpoints, giving piecewise-constant segments, of which the positive ones are expanded. Internal nodes come from `contained_sums`. Each gets distinct hit k-mers, occurrences, |S\_v| and containment.
    - Writes `leaves.parquet`, `nodes.parquet` and `query.json`.
    - *Ceiling:* every internal node is queried (O(N log P)). At full scale, only the touched nodes' ancestors should be.
  - **Test** (fixture, protein sample with duplicates and an unrelated sequence):
    - leaf counts and occurrences equal brute-force holdings joined to the sample;
    - `pos_size` equals each leaf's held k-mers;
    - for every internal node, the k-mers whose containing postings sum to 1 number `node_size`, and their sample hits equal the query's counts.
  - **Sketch beyond *t\_cap*.** Top-ups beyond *t\_cap* are 0.2–1.3% of postings: quick 1,898, medium 8,693, large 50,926, on 1,898 / 6,870 / 27,037 positions. Because of them the largest key is ~2⁶⁴, which would make the query keep every read k-mer. Keeping only hashes ≤ *t\_cap* plus those keys cut sampled hashes 5× (4.4 M to 0.88 M on 125 k read pairs of SRR1438000), with identical results.
  - **Local timings on the core sample's tiers** (4 threads):
    - 125 k read pairs (SRR1438000): 0.4–0.6 s, ~0.9 GB peak RSS. Few hits, as expected at 1 in 5,000: 162 / 332 / 582 hit hashes.
    - 50,000 of the sample's own proteins as the query: counting takes 0.15–0.46 s for 186 k / 363 k / 715 k hit postings. Lookup takes 10–13 s, all of it in the protein path's per-sequence hashing (tests only). Peak RSS 0.44–0.56 GB.
  - **Finding: the leaf minimum covers leaves, not every member.** Of those 50,000 proteins, 1,023 (quick), 296 (medium) and 86 (large) do not hit their own leaf. 99% are partial fragments shorter than their leaf's representative (median 43–66 aa against 134–155). The leaf's base band, floor and top-ups are the representative's k-mers, some outside the fragment.
    - Reads from the fragment's real gene will usually hit the container's k-mers too, so this matters mainly for fragment-only genes.
    - A per-member guarantee would top up each member from its own hashes (posted at its leaf, held by all its members, as floors are). Its cost is not measured. The survey's "every protein covered" holds per leaf.
  - *Next:* phase 3 (reporting, recall first) on these counts, or the HPC BUILD of the 1-in-100 tiers.
* **A leaf per protein, surveyed (2026-10-30, the user's request: every protein, not representatives; survey on HPC first, to size the tiers and choose parameters).**
  - **What was stored.** In phase 0's tree, leaves are distinct sets at the topology rates after the equal-set and containment collapses. Each leaf stores its longest member's k-mers. So contained fragments, and members whose sets differ at a tier's own rates, are not indexed by their own k-mers: phase 2's 1–2% own-leaf misses.
  - **`--leaves proteins`** (survey and build; the builder's default; the user, same day: the topology is built from dense k 8, and proteins share a leaf only when those dense sets match exactly):
    - Each distinct dense 8-mer set (every 8-mer, rate 1) is a leaf, with no containment collapse. Each leaf stores its own k-mers (its longest member's; members have equal 8-mer sets, so nearly always equal sequences).
    - The topology is the storage tree (Kruskal with the virtual empty root) under |A Δ B| of the dense sets, from an exact pair join. `--topology-max-holders` can cap common k-mers, after which pairs are weighed exactly.
    - *First try, superseded:* distinct 11-mer sets up to *t\_cap* as leaves, on the sampled topology. On the quarter sample: 237,875 leaves against 151,474 collapsed, 21–26% more postings. Medium (*t\_base* 0.02, *n\_min* 16, *τ* 0.3, *n\_1* 2): 637,307 against 519,959, i.e. 1.47× tier 1 (432,440) against 1.20×.
  - **Quarter sample, dense 8-mer leaves** (3.3 min, 3.9 GB; the 8-mer pair join took 92 s):
    - 262,391 proteins give 261,954 leaves, so almost every protein has its own. There are 179,062 internal nodes and 82,892 components.
    - Postings are 1.25–1.33× the collapsed tree's at every setting with *n\_1* ≥ 1.
    - Over tier 1 (432,440): quick (0.01/8, *τ* 0, *n\_1* 1) 0.82×; medium (0.02/16, *τ* 0.3, *n\_1* 2) 1.56×; large (0.05/32, *τ* 0.3, *n\_1* 2) 2.96× (tier 1 at 0.05/32 is not built on the quarter).
  - **Set up for HPC** (`hpc/kfp-ablations/tree-index-phase0`, resumed):
    - SURVEY over both leaf modes. Grid: *t\_base* 0.01 / 0.02 / 0.05, *n\_min* 8 / 16 / 32, *τ* 0 / 0.3 / 0.6 / 0.9, *n\_single* 0 / 2 / 4 / 16, *n\_1* 0 / 1 / 2 / 4.
    - TIER1_MATCHED also builds the large tier (`tbase05_all`: 0.05, *n\_min* 32) on the sample, so 0.05/32 rows get `over_matched`.
    - INDEX_SIZES lists every file under the release index, so tier 1's bytes are included.
    - Test profile run locally: 1,152 rows, both modes, matched rows at 0.01/8, 0.02/16 and 0.05/32.
    - *Not included:* tier 1 on core + extra, so `all`'s ratios still compare against core-only tier 1.
* **Survey with a leaf per protein: results and tier parameters (HPC, 2026-10-30; the user's request: choose the tiers from the results).** `tree-index-phase0.tar.gz`, resumed.
  - **Run.** SURVEY took 16 min / 12.5 GB on core and 34 min / 43 GB on all. TIER1_MATCHED built today's large tier (0.05/32) on the sample in 2.5 h: 3,546,868 postings.
  - **INDEX_SIZES gap closed.** The release's tier-1 posting table is the `tier2.*` files (the old name, kept on disk): 90.6 GB of the standard tier's 196 GB. The unit columns are ~105 GB.
  - **Leaves.** Core's 1,067,329 proteins give 1,065,645 leaves (distinct dense 8-mer sets), 736,282 internal nodes and 329,363 components; all's 2,067,329 give 2,062,934. At release scale that is ~5.7 × 10⁹ leaves. Only the 172 proteins with no valid 11-mer are blind (*n\_1* ≥ 1).
  - **Cost of a leaf per protein.** 1.25–1.4× the collapsed tree's postings at every setting. Against tier 1 at matched rates (*τ* 0.3 unless noted):

    | Settings | *n\_1* 1 | *n\_1* 2 | *n\_1* 4 |
    | --- | --- | --- | --- |
    | 0.01 / 8 (*τ* 0) | 1.64× | 2.12× | 3.64× |
    | 0.02 / 16 | 1.47× | 1.58× | 1.99× |
    | 0.05 / 32 | 1.45× | 1.47× | 1.55× |

    - Single-leaf floors stay expensive: at 0.02/16, *n\_single* 2 / 4 / 16 cost +11% / +27% / +79%.
    - *τ* 0.3 costs +1–2% over components; 0.6 costs +6–8%.
  - **Release sizes (estimate).** About 9 B per posting (tier 1's measure on the sample: 15.57 MB for 1,727,926 postings; the release's table is ~9.75 B on core scaling). Per-leaf tables add ~35 GB: |S\_v| and merge level at 2 B each per leaf and per gap of the leaf order, the gap standing for the internal node that splits there. The protein-id map (~30 GB) is read only for reported ranges, and is not counted. Postings are projected from all (low) and core (high):

    | Settings | Release postings | ~Size |
    | --- | --- | --- |
    | 0.01 / 8, *τ* 0, *n\_1* 1 | 6.2–7.7 × 10⁹ | 90–105 GB |
    | 0.01 / 8, *τ* 0, *n\_1* 2 | 7.6–9.9 × 10⁹ | 103–124 GB |
    | 0.02 / 8, *τ* 0.3, *n\_1* 2 | 10.3–11.9 × 10⁹ | 128–142 GB |
    | 0.02 / 16, *τ* 0.3, *n\_1* 2 | 11.9–14.7 × 10⁹ | 142–167 GB |
    | 0.05 / 16, *τ* 0.3, *n\_1* 4 | 23.5–25.7 × 10⁹ | 246–266 GB |
    | 0.05 / 32, *τ* 0.3, *n\_1* 2 | 24.9–28.0 × 10⁹ | 259–287 GB |

    **~64 GB is out of reach with a leaf per protein.** The per-leaf tables alone are ~35 GB, and the cheapest setting with every protein covered needs ~6 × 10⁹ postings.
  - **Tiers chosen** (all with a leaf per dense 8-mer set, the dense 8-mer topology and *n\_single* 0):

    | Tier | *t\_base* | *n\_min* | *τ* | *n\_1* | ~Size | Against today's tier at its rates (format today) |
    | --- | --- | --- | --- | --- | --- | --- |
    | quick | 0.01 | 8 | 0 | 1 | 90–105 GB | 1.64× slim's postings (slim ~91 GB) |
    | medium | 0.02 | 16 | 0.3 | 2 | 140–170 GB | 1.58× the standard tier's (196 GB) |
    | large | 0.05 | 32 | 0.3 | 2 | 260–290 GB | 1.47× the large tier's (~310 GB) |

    - *Quick:* the cheapest rates and floors with every protein covered, for a quick profile. It is ~1.5× the 64 GB target because of the per-protein tables.
    - *Medium:* the standard tier's rates, for a profile followed by the dense tier. *n\_1* 2 lets every protein reach two hits, and *τ* 0.3 gives floors more coherent than tier 1's for +1–2%. If ~128 GB is a hard limit, 0.02/8 (128–142 GB) halves the floors, at the cost in short-protein detection that *n\_min* 16 was chosen to avoid.
    - *Large:* today's large tier's rates (`tbase05_all`, which calibrates `aai`), so the two compare like for like. 0.05/16 with *n\_1* 4 (246–266 GB) fits 256 GB more closely.
    - In today's format, each tree tier is about the size of the tier it replaces or smaller, at per-protein resolution.
* **The three tiers built in parts, set up for HPC (2026-10-30, the user's request).** One process cannot build the 1-in-100 subset (~9.9 × 10⁹ residues, ~53 M proteins). From the quarter sample, a proteins-mode build is 3.4 GB per 262 k proteins (most of it the dense 8-mer join), so ~700 GB at that size. `workflows/tree-index/build.nf` builds it in parts instead:
  - **PARTITION** (`partition.py`). Proteins are linked when they can share a k-11 storage edge at the topology rates, or have the same sequence, so equal dense sets never straddle parts. Components are packed, largest first, into parts of ~2 M proteins.
    - A new kernel, `_core.overlap_components`, unions pairs as `overlap_pairs`' scan finds them (shared scan, batched union-find), so the pairs are never held. On the core sample there would be 72 M of them.
    - Core sample: 311,795 components, the largest 3,454 proteins; 37 s, 2.6 GB (3.8 GB with the pairs kept). That is ~130 GB at 53 M proteins, in a 256 GB job.
  - **BUILD:** `build.py` per part and tier, ~26 GB per 2 M-protein part.
  - **MERGE** (`merge.py`): parts concatenated in name order, positions offset, postings re-sorted by hash. Same format, so `query.py` reads it.
  - **Ceiling:** dense-topology edges between parts are lost, so a k-mer shared across parts is posted in each. Cuts fall between k-11 storage components, which few such edges cross.
  - **Tests:**
    - Rust: `overlap_components` equals union-find over `overlap_pairs` on random sets.
    - `test_partitioned_build_is_sound`: fixture in ≥ 3 parts, built and merged. The merged index passes the single build's brute-force checks: positions one-to-one with dense 8-mer sets, base band, floors, coverage, and leaf and node query counts.
    - `build.nf`'s test profile runs locally: 5 parts × 3 tiers merged.
  - **HPC:** `hpc/kfp-ablations/tree-index-build/run.sh` builds quick, medium and large from `lin100-dense/members.parquet` into `results/index-<tier>/`.
* **Diff schemes on the per-protein design; the spanning tree replaces the merge tree (2026-10-30, the user's request: evaluate the diff schemes while the tiers build; then switch, and estimate the build cost).**
  - **Schemes on the base band** (survey's `schemes`, phase 0's core sample, 1,067,329 proteins), relative to Fitch on Kruskal's merge tree, the same at *t\_base* 0.01, 0.02 and 0.05:

    | Scheme | Collapsed leaves | A leaf per dense 8-mer set |
    | --- | --- | --- |
    | Fitch on the merge tree (phase 0's choice) | 1 | 1 |
    | **Fitch on the spanning tree** | **0.973** | **0.875** |
    | Sequence MST | 1.13–1.15 | 0.94 |
    | Add-only | 1.09–1.10 | 1.26–1.27 |
    | Majority | 1.15 | 1.17 |
    | Flat | 1.63–1.64 | 2.43–2.44 |

    - With a leaf per protein, the sequence MST beats Fitch on the merge tree. Fitch is optimal on its tree, so the merge tree loses something: it joins whole clusters at each Kruskal step and forgets which proteins the edge connected, and chains of fragments and close relatives pay for that.
    - Fitch on the spanning tree itself keeps that adjacency. It can cost no more than the sequence MST (its internal states may take each node's own set), and it is the cheapest scheme everywhere.
  - **The spanning tree as a binary tree** (`survey.spanning_tree`). Each component is rooted at its node with the root edge. A node's subtree is the node as a leaf, joined to its children's subtrees one at a time, lightest edge first; a join's level is that edge's weight.
    - Same leaves, internal-node count and components as the merge tree. Mean node range 33.5 leaves against 37.7 (quarter sample).
    - *Cost:* a single-linkage cluster is still a connected subtree, but no longer always one leaf range, so the neighbour rule's candidates become a few ranges rather than one.
  - **Whole tiers** (quarter sample, survey with both trees):
    - Spanning over merge: quick 0.938, medium 0.947, large 0.917.
    - The base band is −11.6%. At *τ* 0.3 floors rise 0.5–5.7% and top-ups 1–9%, because the spanning tree's nodes are less coherent groupings. At *τ* 0 floors and top-ups are unchanged, since the components are the same.
    - Every surveyed setting is cheaper except 0.01/32 at *τ* 0.3 (+2%), which no tier uses.
  - **Build cost** (quarter sample, medium tier): unchanged. 70.3 s against 70.4 s, peak RSS 4.66 GB against 4.67 GB. The spanning tree takes 0.4 s per 262 k proteins (~3.5 s per 2 M-protein part, a Python loop; a Rust kernel at release scale), and its Fitch pass takes as long as the merge tree's. The index is 3.7% smaller in bytes and has 5.3% fewer postings.
    - *Note:* this build peaked at 4.7 GB, above the 3.4 GB measured earlier for the same configuration. Per 2 M-protein part that is ~36 GB, inside BUILD's 64 GB.
  - **Tier sizes, revised** (release postings × the quarter's ratio; ~9 B per posting + ~35 GB of leaf tables): quick 5.8–7.2 × 10⁹ postings, ~87–100 GB; medium 11.3–13.9 × 10⁹, ~137–160 GB; large 22.8–25.7 × 10⁹, ~240–266 GB.
  - **Switched:** `--topology-tree spanning` is now the default of `survey.py` and `build.py`, so `build.nf` and `hpc/kfp-ablations/tree-index-build` build on it. `--topology-tree merge` keeps phase 0's tree. `survey.json` reports Fitch on both trees per *t\_base*, with timings.
* **Phase 3: reporting, recall first (2026-10-30, the user's request).**
  - **Build** now also writes the spanning tree's edges between positions (`edge_a`, `edge_b`, `edge_w`, in dense 8-mer |A Δ B|) and each position's representative length (`pos_len`). `merge.py` offsets the edges. That is 16 B per position, ~90 GB at ~5.7 × 10⁹ leaves in this development format; a compact one would store each leaf's parent edge only (~8 B).
  - **`query.py`** also writes `hits.parquet`, the sample's k-mers that are keys, for gather.
  - **`report.py`:**
    - *AAI* per hit leaf and internal node: containment^(1/k), with a 95% Clopper-Pearson interval on hit k-mers out of |S\_v|. Naive: low depth lowers containment as divergence does. Phase 4 calibrates it against `presence-calibration`'s truth.
    - *Candidates.* Every leaf with ≥ `--min-hits` hit k-mers (default 1: every hit leaf, as the recall-first pileup takes every hit unit) is a seed. Each adds its single-linkage cluster on the spanning tree at *R\_v* = d(aai\_lo) + d(1 − δ), where d(p) = 2 (len − 7)(1 − p⁸) is the expected dense 8-mer |A Δ B| at identity p. The cluster is walked by Prim from the seed (nearest first in minimax distance), at most `--max-neighbours` (256) leaves beyond it. On the spanning tree a cluster is a connected subtree, not a leaf range (Diff schemes entry above), so this walks edges rather than looking up an ancestor.
    - Candidates are written as position ranges (`candidates.parquet`) and as proteins (`members.parquet`, with a seed flag), the inputs for round 2 and the pileup's `REFERENCE`. *Not wired:* `kfp-pileup reference` takes a unit profile and the unit store; a tree-index store fetched by range is still to build.
    - *Gather* over the seed leaves. The score is remaining hit k-mers × windows / |S\_v| (estimated shared windows, as today's gather divides by *t\_g*). Leaves in the best's component that share a remaining k-mer with it and score within `--tie-z` (2) Poisson standard deviations of it are tied. The group is reported as its smallest common subtree, and its k-mers are removed.
    - *Amended during the run:* ties were first not restricted to the best's component. Groups then spanned components (one covered 260 k of 262 k positions), which is not a subtree.
  - **Test** (`test_tree_build.py`, single and partitioned builds): AAI bounds, and AAI 1 at full containment. Each walk, uncapped, equals the seed's component under edges ≤ radius (scipy's connected components, five radii); capped, it is the uncapped walk's prefix. Gather's groups take disjoint k-mers that together are every hit k-mer a seed holds, and each group's range contains its leaf. Every sampled protein's leaf is a candidate.
  - **Local run** (medium tier on the quarter sample, 262,391 proteins; 4 threads):

    | Query | Seeds | Candidates (+ neighbours) | Capped seeds | Ranges | Gather groups (tied) | Report |
    | --- | --- | --- | --- | --- | --- | --- |
    | 125 k read pairs (SRR1438000) | 384 | 1,017 (+633) | 0 | 143 | 109 (48) | 0.08 s |
    | 50,000 core proteins as the query | 137,437 | 140,670 (+3,233) | 35,766 | 8,289 | 10,131 (5,963) | 37 s |

    - Of the 50,000, 12,162 have a position in the quarter index. All but one hit their own leaf and are candidates. Gather picks the own leaf itself for 38% and reports it inside a group for 97.8%. The other 2.2% are explained away: a group in another component took their k-mers first (k-mers shared across components are posted in each).
    - Tie groups: median 2 leaves, p90 21, p99 158, max 1,642. Common subtrees: median 2 positions, p90 35, p99 317. At the medium tier's ~4 base k-mers per protein, near-identical relatives cannot be told apart, so ties are the expected result, not a defect.
    - The cap binds for a quarter of the protein seeds. A leaf with few k-mers has a wide interval (2 of 2 hit gives aai\_lo ~0.85), so its radius reaches far. That is the rule's intended direction for recall; phase 4 measures what it costs in candidates.
    - *Ceilings:* the walk is a Python loop per seed (12 s for 137 k seeds), and gather recounts all pairs per pick (24 s for 10 k picks). Both need a kernel before HPC-scale samples.
  - *Next:* phase 4 (benchmark against tier 1: δ at 0.05, 0.1 and 0.15, AAI calibration, nearest-match identity scoring), after the HPC build of the 1-in-100 tiers.
* **Candidate recall, measured (2026-10-30, the user's request: recall first, ~100% of proteins at depth > 0.5 into the candidates, since tier 2 is built per sample from them, with few false candidates; explore broadening on the tree, then focus on tier 2).**
  - **Benchmark** (`workflows/tree-index/recall.py`): 2,000 proteins with a position in the quarter sample's index each carry a strain made by `spike.py simulate`, at identity 1 / 0.98 / 0.95 / 0.9 / 0.85 / 0.8 and lognormal depth (μ 0, σ 1.2). The 125 k real read pairs of SRR1438000 are mixed in as background. A strain counts as recalled when its source protein's position is a candidate. *observable* is the number of the source's 11-mers that occur whole in the strain's own reads. Background-only queries give the false-candidate baseline.
  - **Frames.** `stopfree` drops every frame with a stop codon, so reads crossing a gene end are lost, and with them most reads of short genes. At 100% identity and depth ≥ 2, 7 of 7 misses were 28–100 aa proteins whose own k-mers were in the reads. `all` frames raised recall at depth > 0.5 from 65% to 74% (medium) and are used below. *For the main query too:* the recall-first path should not use `stopfree`.
  - **Broadening on the tree is not the lever.** On the medium tier at depth > 0.5: hit leaves 73.6%; gain ranges (a hit k-mer's whole gain range, including leaves that lost it) 74.1%; neighbour rule δ 0.1, cap 256: 77.5%; whole components of the hit leaves, the most any tree walk can add: 77.6%. δ from 0 to 0.3 changed nothing. The cap mattered (16 against 256: +2%).
  - **Why: density.** A medium leaf holds a median 6–8 k-mers. With 0 of them in the strain, recall is 17%; with 1, 77%; with ≥ 3, ≥ 96%. At 80% identity, 70% of strains share none of the source's indexed k-mers.
  - **Window minimizers** (`build.py --window W [--window-min M]`): the base band also takes each leaf's least hash, of any value, in every window of W k-mers. With `--window-min` M, the window shrinks to n\_k // M, so short proteins get about 2M minimizers. Window k-mers are diffed by Fitch with the rest of the base band, so relatives share them. A 50-aa read spans 40 k-mers, so W ≤ 40 puts an indexed k-mer in every whole read of an identical protein. Checked against brute force (`test_window_kmers_match_brute_force`). Minimizers first stopped at *t\_cap*, which left ~11% of w10 windows empty; uncapped: +1.6%. The query already keeps keys beyond *t\_cap*.
  - **Results** (quarter sample, depth > 0.5, 1,538 strains; neighbour rule δ 0.1, cap 256; "short" ≤ 60 aa; bg = leaves hit by the background reads alone):

    | Index | Postings (× medium) | Bytes | bg leaves | Candidates | All | ≥ 95% id | Short, ≥ 95% | > 60 aa, ≥ 95% | observable 8–15 | observable ≥ 16 |
    | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
    | medium | 1 | 20.9 MB | 602 | 65,439 | 0.775 | 0.928 | 0.768 | 0.954 | 0.624 | 0.907 |
    | large | 1.84 | 30.5 MB | 1,126 | 70,432 | 0.842 | 0.958 | 0.786 | 0.987 | 0.734 | 0.969 |
    | medium + w20 | 2.79 | 41.4 MB | 1,361 | 73,783 | 0.873 | 0.968 | 0.804 | 0.996 | 0.815 | 0.995 |
    | medium + w10 | 4.13 | 56.8 MB | 1,705 | 75,297 | 0.901 | 0.977 | 0.839 | 1.0 | 0.925 | 0.998 |
    | medium + w30, min 8 | 3.62 | 51.0 MB | 1,447 | 75,263 | 0.919 | 0.984 | 0.884 | 1.0 | 0.942 | 0.996 |
    | medium + w20, min 8 | 3.99 | 55.2 MB | 1,561 | 75,695 | 0.921 | 0.984 | 0.884 | 1.0 | 0.948 | 0.997 |
    | medium + w10, min 8 | 5.72 | 75.1 MB | 2,336 | 77,160 | 0.932 | 0.984 | 0.884 | 1.0 | 0.988 | 1.0 |

    - Lowering *t\_base* to 0.005 under the windows saved nothing: floors and top-ups grew back (4.16× and 5.9×), with the same recall.
    - `--min-hits` 2 halves to thirds the background leaves (1,561 → 530 on w20 min 8) but loses 5.7% of recall at depth > 0.5, mostly short and divergent strains. Recall first keeps 1.
    - *Candidates* here are 26–29% of the quarter's proteins, for 2,000 strains: ~35 leaves per strain, its near-identical relatives. The background alone seeds 0.6–0.9% (1.4–2.3 k leaves from 125 k pairs).
  - **Reading.**
    - For proteins over 60 aa at ≥ 95% identity, windows reach 100% (w30 min 8 and up). Whenever the reads hold ≥ 16 source 11-mers, w10 min 8 recalls 100% and w20 min 8 99.7%.
    - What is left is not reachable by an exact 11-mer index. Of short strains at depth > 0.5, 20% have no source 11-mer in their reads and 51% fewer than 8. At 80% identity, 18% of all strains have none. Depth > 0.5 on a 30–60 aa gene is often one read pair overlapping it by 15–20 aa.
    - *The user's question (floors for short proteins):* raising them is what `--window-min` does, and it took short proteins with ≥ 8 observable k-mers from 86–90% to 99.1%. Beyond that a higher floor cannot help. A shorter seed for short proteins (k 8–9 or a reduced alphabet) would, at a cost in chance hits; not tried.
    - *Cost:* w20 min 8 is ~4× medium's postings, ~2.6× its bytes here (per-leaf tables are fixed). At release that is roughly 500 GB against medium's 140–170 GB (estimate from the survey's projection). w30 min 8 is ~3.6× for the same recall at ≥ 95%.
* **A k ladder for short proteins (2026-10-30, the user's request: try a shorter seed for short proteins first, as low as k 7; the medium tier should be about half of w30 min 8's size).**
  - **`build.py --ladder k:max_len:w,...`:** leaves shorter than max\_len also hold window minimizers of k-mers of a shorter k, diffed by Fitch with the rest. The query hashes reads at each ladder k and keeps the k-mers that are keys.
    - *Bug found on the way:* k-mer hashes do not depend on k. A 9-mer shared a hash with 11-mers, which put duplicate (leaf, hash) pairs into the Fitch kernel (a panic in `lca`). Ladder hashes are now moved to their own domain per k by a bijection (`query.salted`).
  - **Results** (on medium + w30 min 8, depth > 0.5; neighbour rule δ 0.1, cap 256):

    | Ladder | Postings (× medium) | bg hit hashes | bg leaves | All | ≥ 95% id | ≥ 90% id | Short, ≥ 95% | Short, all | observable 0 |
    | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
    | none | 3.62 | 700 | 1,447 | 0.919 | 0.984 | 0.972 | 0.884 | 0.744 | 0.082 |
    | k 9 < 100 aa (w 4) | 5.45 | 2,089 | 2,613 | 0.949 | 0.991 | 0.987 | 0.938 | 0.843 | 0.342 |
    | k 9 < 100 aa, k 7 < 60 aa (w 2) | 6.63 | 40,915 | 27,472 | 0.973 | 1.0 | 0.995 | 1.0 | 0.973 | 0.795 |
    | k 7 < 100 aa (w 2) | 6.79 | 101,439 | 63,521 | 0.976 | 1.0 | 0.995 | 1.0 | 0.973 | 0.849 |

  - **Reading.**
    - k 7 closes the short-protein gap: 100% at ≥ 95% identity, and 80–85% of strains with no source 11-mer in their reads are found.
    - But the background reads alone then seed 10–24% of the quarter index's leaves. The cause is read-side chance, not the index: a 7-mer occurs by chance about once per 20⁷ ≈ 1.3 × 10⁹ read windows, so in a 10 M-pair sample (~2 × 10⁹ windows) every indexed 7-mer is expected to be hit about once. At release every leaf holding 7-mers would be seeded in every sample.
    - k 9 costs +1.8× postings for +3% recall (short proteins 88% to 94% at ≥ 95% identity), with ~1.8× the background seeds. At release, chance single 9-mer hits are ~2 × 10⁹ × keys / 20⁹, which is large too.
  - **Design proposed (the user's points: the k per leaf from the data at build time, by k-mer count or uniqueness in the database rather than fixed thresholds; promiscuous k-mers removed by an information criterion on the tree rather than a holder cut):**
    - *Evidence of a hit.* A hit on k-mer x is worth b(x) = −log₂ q₀(x) − log₂ g(x) bits. q₀(x) is the chance probability of x in a protein window, from residue composition, so low-complexity k-mers are worth less. g(x) is x's number of independent gains on the tree (Fitch), so a hit's weight is split over the g unrelated clades it seeds. A clade marker (one gain, however many leaves) keeps all its bits. This replaces the 16-unit cut and the plan's fixed gain cut.
    - *Seeds by evidence within one read pair:* a leaf is a seed when its hits from one read pair sum to at least B bits. B is set by an expected number of chance seeds per sample: B = log₂(read windows × leaves / chance seeds). It is the one parameter, in units of false candidates, and it holds for every k: one clean 11-mer hit carries ~46 bits, so 7-mers need two or more hits in the same read to count. Seeds then need per-read hit sets (`FastxHits` gives read ids).
    - *Pruning:* a k-mer is dropped when no read can use it, i.e. when even the most hits of its k that fit in one read alongside it cannot reach B.
    - *The k per leaf:* going down the ladder (11, 9, 8, 7), a leaf takes the next k only while it has fewer window-covering k-mers with positive evidence than a read overlap needs. Evidence comes from the same Fitch pass, so the rule is set by the database, not by protein length.
* **Seeds by evidence, the ladder by information, sized for two budgets (2026-10-30, the user's go-ahead on the design above; budgets: large about w30 min 8's ~425–515 GB, medium about half that).**
  - **Build** writes each key's evidence per hit, b(x) = −log₂ q₀(x) − log₂ g(x) (`keys.npy`, `key_bits.npy`, `key_k.npy`).
    - q₀ is the chance of the k-mer at a protein window, from the database's residue frequencies; g is the key's number of positive postings (Fitch gains, floors, top-ups).
    - `merge.py` adds gains across parts exactly: 2^−b = q₀ g, so merged bits are −log₂ Σ 2^−b.
  - **Ladder by information** (`--ladder k:w,...`, `--min-kmers M`): a leaf takes the next k while it has fewer than M informative k-mers at the k before. Informative means a priori rarer than one chance holder among the N leaves: −log₂ q₀ ≥ log₂ N. Low-complexity leaves therefore descend sooner, at the same length.
  - **Query seeds by evidence** (`query.py`, reads; `seeds.parquet`, `evidence.parquet`):
    - A read pair's hits on a leaf score Σ b − (m − 1) log₂(windows per read pair). Each hit after the first had that many places to fall by chance.
    - Computed with one sweep over (read pair, range endpoints), as the counting sweep. Checked against brute force on the fixture (single and partitioned builds).
    - The bar B is the least score at which the expected chance seeds, Σ over keys with b ≥ B of (read windows of its k) × 2^−b, are at most `--chance-share` (1%) of the seeds scoring ≥ B. The one parameter is a share of false candidates.
    - B is relative, so it should hold as sample and index grow together. *Ceiling:* the chance term counts single hits only.
  - **No pruning by gains was needed.** 7-mers carry 25–35 bits at ~3.5 gains per key: never enough alone (bar ~33–37 bits), but two in one read pair reach ~47 bits. A rule that drops keys unable to contribute removes almost nothing. Size is set by how many minimizers each leaf takes, so that is what the configurations below vary.
  - **Floors are redundant** once every leaf has spread minimizers: without them (`--n-min 0`), recall was unchanged (M2 against M4, L2 against L3) for 0.23 M fewer postings on the quarter sample.
  - **Results** (quarter sample, depth > 0.5, all frames; seeds by evidence at 1% chance share, with the neighbour rule δ 0.1, cap 256; bg = seeds from the background reads alone; release size scaled from medium's projection, 107–132 GB of postings plus ~35 GB of leaf tables, without the ~45 GB of spanning-tree edges the neighbour rule needs):

    | Configuration | Postings (× medium) | ~Release | bg seeds | All | ≥ 95% id | ≥ 90% id | Short, ≥ 95% | > 60 aa, ≥ 95% |
    | --- | --- | --- | --- | --- | --- | --- | --- | --- |
    | medium (hit leaves, `stopfree`) | 1 | 140–170 GB | 602 | 0.775 | 0.928 | 0.888 | 0.768 | 0.954 |
    | **M5:** *t\_base* 0.01, no floors, ~6 spread minimizers (`--window 1000 --window-min 3`), k 7 at w 3 for < 32 informative 11-mers | 1.46 | 190–230 GB | 1,124 | 0.846 | 0.971 | 0.939 | 0.920 | 0.979 |
    | **M4:** *t\_base* 0.02, no floors, ~8 spread minimizers (`--window-min 4`), k 7 at w 2 for < 32 | 2.04 | 255–305 GB | 1,683 | 0.889 | 0.989 | 0.966 | 0.955 | 0.994 |
    | M2: as M4 with floors (*n\_min* 16, *τ* 0.3) | 2.41 | 290–355 GB | 1,829 | 0.895 | 0.989 | 0.969 | 0.955 | 0.994 |
    | E0: w30 min 8, no ladder | 3.62 | 425–515 GB | 1,447 | 0.919 | 0.984 | 0.972 | 0.884 | 1.0 |
    | **L3:** w40 min 8, no floors, k 7 at w 3 for < 48 | 4.09 | 475–575 GB | 2,969 | 0.951 | 1.0 | 0.993 | 1.0 | 1.0 |
    | L2: as L3 with floors | 4.45 | 510–620 GB | 3,081 | 0.952 | 1.0 | 0.993 | 1.0 | 1.0 |
    | E7: w30 min 8, floors, k 7 at w 2 for < 64 | 5.64 | 640–780 GB | 6,027 | 0.962 | 1.0 | 0.994 | 1.0 | 1.0 |

    - Without the evidence bar, the same E7 index's hit leaves are seeded 42.6 k times by the background, against 6.0 k with it, for 0.6% less recall (0.941 against 0.935 before the neighbour rule).
    - The neighbour rule adds 2.5–6% at all-identity recall, mostly divergent strains, for +18–27% candidate members.
  - **Reading.**
    - *Large:* L3 meets the target (100% at ≥ 95% identity, short proteins included; 99.3% at ≥ 90%) at ~475–575 GB, ~10% over w30 min 8's budget.
    - *Medium:* M4 gives 98.9% at ≥ 95% (99.4% over 60 aa, 95.5% for short) at ~255–305 GB. M5 fits ~half the large budget (~190–230 GB) at 97.1%.
    - The remaining misses at ≥ 95% identity are mostly strains whose reads hold few or no source 11-mers.
* **Tiers by depth at L3's recall (2026-10-30, the user's request: L3-level recall with large to depth 0.5, medium to 1–1.5, small to 3–4).**
  - **Second sample** (`sim2`): 3,000 other targets at identity 1 / 0.98 / 0.95 / 0.9, depth lognormal (μ 0.9, σ 0.8), the same background. Pooled with the first sample: 2,991 strains at ≥ 95% identity and 3,998 at ≥ 90% above depth 0.5, and 989 / 1,321 above depth 4.
  - **Minimizer-only builds** (*t\_base* 0.001, no floors, `--window 1000 --window-min M`: about 2M spread minimizers per protein; k 7 at w 4 for < 16 informative 11-mers) are the efficient frontier at ≥ 90% identity.
  - **Recall at ≥ 95% / ≥ 90% identity above each depth** (seeds by evidence at 1% chance share plus the neighbour rule; release size scaled as above, without spanning-tree edges):

    | Index | Postings (× medium) | ~Release | bg seeds | > 0.5 | > 1 | > 1.5 | > 2 | > 3 | > 4 |
    | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
    | S3: *t\_base* 0.01, ~4 minimizers | 0.99 | 140–165 GB | 786 | 0.964 / 0.924 | 0.980 / 0.944 | 0.985 / 0.951 | 0.986 / 0.956 | 0.988 / 0.959 | 0.990 / 0.964 |
    | M5 | 1.46 | 190–230 GB | 1,124 | 0.987 / 0.961 | 0.996 / 0.975 | 0.998 / 0.978 | 0.999 / 0.983 | 0.999 / 0.984 | 1.000 / 0.986 |
    | S4: ~8 minimizers | 1.53 | 200–235 GB | 725 | 0.986 / 0.964 | 0.995 / 0.978 | – | 0.999 / 0.987 | 0.999 / 0.988 | 0.999 / 0.990 |
    | **S5:** ~10 minimizers | 1.94 | 245–290 GB | 887 | 0.991 / 0.975 | 0.998 / 0.987 | 0.999 / 0.990 | 1.000 / 0.994 | 1.000 / 0.996 | 1.000 / 0.998 |
    | M4 | 2.04 | 255–305 GB | 1,683 | 0.994 / 0.978 | 0.999 / 0.988 | 0.999 / 0.990 | 1.000 / 0.992 | 1.000 / 0.994 | 1.000 / 0.995 |
    | **S7:** ~14 minimizers | 2.73 | 330–395 GB | 1,247 | 0.994 / 0.982 | 0.998 / 0.991 | 1.000 / 0.994 | 1.000 / 0.996 | 1.000 / 0.997 | – |
    | **L3** | 4.09 | 475–575 GB | 2,969 | 0.999 / 0.993 | 1.000 / 0.997 | 1.000 / 0.999 | 1.000 / 0.999 | 1.000 / 1.000 | 1.000 / 1.000 |

  - **Reading.**
    - L3's level at depth > 0.5 is 0.999 / 0.993. It is reached by L3 at depth 0.5, by S7 (~14 minimizers) at 1.5, and by S5 (~10) at 2. S4 (~8) reaches 0.999 / 0.988–0.990 at depth 3–4.
    - *Depth stops paying once the reads cover the protein.* What is left at ≥ 90% identity is proteins whose indexed 11-mers all carry a substitution: at identity p a protein with n indexed 11-mers is missed with probability ~(1 − p¹¹)ⁿ, 0.69ⁿ at 90%, at any depth. ~99.3% at 90% therefore needs n ≈ 14, which is S7.
    - So tiers separated by depth compress: below ~1.5× (S3) recall at ≥ 90% stalls near 96% even at depth 4. A small tier at L3's ≥ 90% level is not much cheaper than the medium. At ≥ 95% only (miss ~0.43ⁿ), ~8 minimizers suffice from depth 2.
* **Units above the leaves, the clustered index generalised (2026-10-30, the user's idea: store k-mers only for higher nodes while their members are within ~95% AAI, at higher density; then: pick each node by whether its k-mers can detect it to an AAI, rather than by a radius).**
  - **Radius units with unions** (`build.py --unit-radius r`). Each spanning-tree component is walked from its root, and a leaf joins its parent's unit while its summed identity distance from the unit's first leaf is ≤ r (edge distance from dense 8-mer weights). Every member holds the union of the members' selected k-mers.
    - r 0.025: 156,609 units of 261,954 leaves; 31,912 multi-leaf, holding 137,257 leaves (largest 233). r 0.05: 126,755 units, 167,697 leaves in multi-leaf units.
  - **Detectability units** (`--detect-aai a*`, `--detect-miss ε`). A leaf joins its parent's unit while the unit's members still share N\* = ⌈log ε / log(1 − a\*¹¹)⌉ k-mers (≤ *t\_cap*): N\* is 14 at a\* 0.9 and 8 at 0.95, at ε 0.007. A multi-leaf unit then indexes its 2N\* lowest shared hashes in place of its members' own minimizers; lone leaves keep theirs.
  - **Results** (both samples pooled; S5's settings, ~10 minimizers per leaf, unless noted; ≥ 95% / ≥ 90% identity):

    | Build | Postings (× medium) | ~Release | bg seeds | > 0.5 | > 1 | > 2 | > 4 |
    | --- | --- | --- | --- | --- | --- | --- | --- |
    | S4 (~8 per leaf) | 1.53 | 200–235 GB | 725 | 0.986 / 0.964 | 0.995 / 0.978 | 0.999 / 0.987 | 0.999 / 0.990 |
    | D95 | 1.54 | 200–240 GB | 861 | 0.983 / 0.961 | 0.996 / 0.979 | 1.000 / 0.988 | 1.000 / 0.994 |
    | S5 | 1.94 | 245–290 GB | 887 | 0.991 / 0.975 | 0.998 / 0.987 | 1.000 / 0.994 | 1.000 / 0.998 |
    | U25 (union, r 0.025) | 2.05 | 255–305 GB | 1,088 | 0.991 / 0.975 | 0.998 / 0.987 | 1.000 / 0.994 | 1.000 / 0.998 |
    | U50 (union, r 0.05) | 2.12 | 260–315 GB | 1,277 | 0.991 / 0.975 | 0.998 / 0.987 | 1.000 / 0.994 | 1.000 / 0.998 |
    | D90 | 2.19 | 270–325 GB | 1,050 | 0.991 / 0.976 | 0.998 / 0.987 | 1.000 / 0.994 | 1.000 / 0.998 |
    | D90d (lone leaves ~14) | 2.71 | 325–395 GB | 1,182 | 0.993 / 0.981 | 0.998 / 0.990 | 1.000 / 0.995 | 1.000 / 0.998 |
    | S7 (~14 per leaf) | 2.73 | 325–395 GB | 1,247 | 0.994 / 0.982 | 0.998 / 0.991 | 1.000 / 0.996 | 1.000 / 0.998 |
    | U25d (union, ~20 per leaf) | 4.11 | 475–580 GB | 2,208 | 0.996 / 0.988 | 0.998 / 0.993 | 1.000 / 0.997 | 1.000 / 0.998 |
    | L3 | 4.09 | 475–575 GB | 2,969 | 0.999 / 0.993 | 1.000 / 0.997 | 1.000 / 0.999 | 1.000 / 1.000 |

  - **Reading.**
    - Units neither save postings nor add recall. Unions cost +6–9% for identical recall. Detectability units land on the per-leaf curve (D95 ≈ S4; D90d ≈ S7, D90 between S5 and S7), and their candidate members are about the same as per-leaf builds.
    - *Why:* the tree is already the general clustered index. Fitch posts a k-mer shared by near-identical members once, at their subtree, so per-leaf minimizers of close relatives cost about what a unit's shared set costs. What remains per leaf are k-mers its relatives do not share, and those are exactly the ones a strain diverged from that leaf needs. A union adds relatives' private k-mers, which a strain at 90% to the source does not carry.
    - Recall at ≥ 90% is set by the indexed k-mers each protein has, ~0.69ⁿ misses as above, whether they are its own or its unit's. At equal n the cost is equal, so a unit cannot buy density.
    - *Scale:* sharing grows with the number of relatives, so the per-protein cost of any design should fall at release against this quarter sample. The ratios between designs should hold, since every design shares the same way.
* **Intersections, soft intersections and the tree's own internal states (2026-10-30, the user's follow-ups: intersections instead of unions with denser sampling; then sets weighted toward conserved k-mers; then the internal nodes' diffs themselves).**
  - **Hard intersections** (`--unit-kmers n`): units grow over dense 11-mer sets while their members still share ≥ n. Each unit indexes ~n of the shared k-mers, spread along its first protein (minimizers restricted to shared k-mers), in place of its members' own.
  - **Soft intersections** (`--unit-tau τ`): a unit's set is the k-mers held by ≥ τ of its members; a leaf joins while that set keeps ≥ n and the leaf holds ≥ τ of it.
  - **Cuts on the tree's own states** (`--cut-detect a*`). The cut is the highest node whose Fitch state every leaf below still shares N\* of (as above; the shared count is |state(v)| less the losses on the path, a lower bound), and Fitch's changes below each cut are dropped. Every leaf under a cut then holds the node's state, a parsimony consensus that is a soft intersection taken from the tree with no extra sets. Covers were checked to sum to 0 or 1 (`build.holdings`) on a 20 k subset.
  - **Results** (both samples pooled; *t\_base* 0.001, no floors, k 7 ladder at w 4 for < 16; ≥ 95% / ≥ 90% identity):

    | Build | Postings (× medium) | ~Release | bg seeds | > 0.5 | > 1 | > 2 | > 4 |
    | --- | --- | --- | --- | --- | --- | --- | --- |
    | S4 (~8 per leaf) | 1.53 | 200–235 GB | 725 | 0.986 / 0.964 | 0.995 / 0.978 | 0.999 / 0.987 | 0.999 / 0.990 |
    | I16 (hard, n 16) | 1.60 | 205–245 GB | 897 | 0.983 / 0.962 | 0.993 / 0.977 | 0.998 / 0.987 | 0.999 / 0.992 |
    | S5 (~10 per leaf) | 1.94 | 245–290 GB | 887 | 0.991 / 0.975 | 0.998 / 0.987 | 1.000 / 0.994 | 1.000 / 0.998 |
    | C90 (S5, cuts at 0.9) | 1.94 | 240–290 GB | 887 | identical to S5 | | | |
    | I24 (hard, n 24) | 1.99 | 250–300 GB | 1,082 | 0.989 / 0.971 | 0.997 / 0.985 | 1.000 / 0.993 | 1.000 / 0.997 |
    | T60 (soft τ 0.6, n 24) | 1.99 | 250–300 GB | 1,031 | 0.990 / 0.972 | 0.997 / 0.983 | 0.999 / 0.993 | 1.000 / 0.998 |
    | I32 (hard, n 32) | 2.36 | 290–345 GB | 1,155 | 0.990 / 0.975 | 0.997 / 0.986 | 1.000 / 0.994 | 1.000 / 0.997 |
    | T80 (soft τ 0.8, n 24) | 2.37 | 290–345 GB | 988 | 0.992 / 0.976 | 0.997 / 0.987 | 1.000 / 0.994 | 1.000 / 0.998 |
    | S7 (~14 per leaf) | 2.73 | 325–395 GB | 1,247 | 0.994 / 0.982 | 0.998 / 0.991 | 1.000 / 0.996 | 1.000 / 0.998 |
    | C95d (~20 per leaf, cuts at 0.95) | 2.81 | 335–405 GB | 1,582 | 0.993 / 0.980 | 0.998 / 0.988 | 1.000 / 0.995 | 1.000 / 0.998 |
    | C90d (~20 per leaf, cuts at 0.9) | 3.59 | 420–510 GB | 1,700 | 0.996 / 0.988 | 0.998 / 0.993 | 1.000 / 0.997 | 1.000 / 0.998 |
    | S10 (~20 per leaf, no cuts) | 3.87 | 450–545 GB | 1,722 | 0.996 / 0.988 | 0.998 / 0.993 | 1.000 / 0.997 | 1.000 / 0.998 |
    | L3 | 4.09 | 475–575 GB | 2,969 | 0.999 / 0.993 | 1.000 / 0.997 | 1.000 / 0.999 | 1.000 / 1.000 |

  - **Reading.**
    - Hard and soft intersections sit on or slightly below the per-leaf curve: I24 and T60 match S5 at +3%, and I32 and T80 match it at +22%. Their units are tight groups whose shared k-mers the tree already posts once.
    - *Cuts on the tree's states are the one free saving:* at a\* 0.9 they drop 7% of a dense build's postings (S10 → C90d) with identical recall to three decimals. At S5's density they never bind, since a leaf with ~10 k-mers cannot share 14. Cutting at 0.95 removes 27% but lands back on the curve (C95d ≈ S7).
    - L3 still beats the ~20-per-leaf builds at ≥ 90% identity (0.993 against 0.988) at similar size. Its k 7 ladder reaches more leaves (< 48 informative 11-mers, w 3, against < 16, w 4), so the ladder's reach matters more than 11-mer density for divergent strains.
    - *Next, if pursued:* cuts at 0.9 on L3's settings, for an expected ~7% off the large tier at its recall.
