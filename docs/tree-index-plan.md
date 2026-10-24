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

- **No dense global tier.** Every distinct k-mer needs at least one gain, so the tree's postings are at least the distinct k-mer count, ~5 × 10¹¹ dense 11-mers if about half of MGnify's ~10¹² occurrences are distinct (the main plan's guess). With deletions that is ~5–8 TB, ~2× smaller than a flat dense index but ~50× the standard tier.
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

## Risks

- **Postings above budget.** If MGnify's members are less redundant than assumed (members at 90–94% to their nearest neighbour post more than a fresh set), the tree degenerates into today's per-unit storage plus a node table. Phase 0 decides.
- **Leaf and node tables.** One row per leaf (level, |S\_v|, length, first member) and one per internal node (leaf range, merge level; N − 1 rows). At ~2 × 10⁹ leaves and ~12 B per row each that is ~48 GB, memory-mapped and read only for touched ranges; ~140 GB if partials do not collapse. Postings that carry their leaf range directly (two leaf ids) would remove the internal-node table at ~4 B more per posting. Phase 0 sizes both.
- **Ties at the stored density.** The stored sets define the storage objective exactly (Building the tree), but with ~16 k-mers per set many candidate parents tie, and the tree's topology among relatives is then arbitrary. That does not matter for storage, but it does for reported subtrees, so ties are broken at a denser nested rate.
- **Single-linkage chaining.** A spanning tree chains: a subtree can span distant ends through intermediates. This matters for what a reported subtree means (a function), not for storage. The reporting model has to bound a subtree's diameter, as TreeCluster does, rather than take whole subtrees.
- **Assignment model.** Gather and EM over nested subtrees are new code. Today's component-wise EM does not carry over directly.
- **Build at scale.** The staged k-mer build (Building the tree) is sorts and joins over ~3 × 10¹⁰ digests plus a prefix join on stage 3's representatives. Its unknowns are the size of stage 3's largest signature groups and stage 4's posting lengths after merging, both measured in phase 0. MGnify90 membership is not needed; it is kept only as a check (does the tree recover the clusters?).

## Implementation plan

Laptop work runs one job at a time with capped threads; full-scale runs go to HPC as Nextflow (SLURM + Singularity) for the user to run.

0. **Measure, no index (gate for everything after).** On a sample of MGnify: a few hundred thousand MGnify90 clusters with all their members, plus the representatives that share k-mers with them across clusters (to see cross-cluster sharing). Laptop on a biome subset, HPC for the full-scale sample. At *t\_base* 0.02, *n\_min* 16, k 11, with today's hashing (`crates/core` `protein_kmers`):
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
1. **Build prototype** (Python, Rust kernels as the existing tier, on the laptop scale of phase 0):
   - the staged build: set digests (Δ = 0), pigeonhole signatures (Δ ≤ d), the rarity-ordered prefix join on representatives (Δ > d), Kruskal by integer level with the virtual empty root, and dense tie-breaks; collapses;
   - the merge tree and leaf order; Fitch states per k-mer (stage 6, missing data for leaves below a k-mer's level); range postings with level bands; leaf and node tables; a sequence store reordered by leaf order.
   - Checked against a brute-force reconstruction: every sequence's sampled set from its root path equals its own; and against brute-force Fitch on small trees (equal minimum change counts).
   - *Updates:* a held-out batch (e.g. 10% of phase 0's sample) added through the delta index. Query counts must be identical to a rebuild containing the batch, and identical again after compaction. The batched Kruskal update's total weight must equal a from-scratch tree's, under equal candidate edges. Report the drift (weight) of leaf-only insertion without the batched update.
2. **Query: counting.** Lookup as today (sketch, merge-join), interval sweep per band, node counts and |S\_v|. Checked against brute-force per-sequence containment on fixtures.
3. **Reporting, recall first.** Candidate subtrees above background, members as leaf ranges into the main plan's round-2 dense stage and the pileup's `REFERENCE`. Then greedy gather with LCA ties, and per-node containment AAI.
4. **Benchmark against tier 1**, at matched index size and *t*:
   - `fmh-benchmark` and `presence-calibration` (spiked truth): recall and precision, rolled up to MGnify90 (via membership), Pfam and KOfam;
   - AAI calibration against `presence-calibration`'s truth;
   - *neighbour rule* (Recall-first candidates and neighbours): δ at 0.05, 0.1 and 0.15. Report the share of spiked units with no own hit that reach the candidate set, by depth and identity; the candidate count against "every hit node"; and how often the true divergence exceeds *d̂⁺*. Against the main plan's index-level neighbour graph (recall-first pileup, next step 0), by that step's rule: worth it if it reaches a substantial share of unhit spikes for < ~2× the candidates;
   - cost: index size, build CPU-h, query wall and peak RSS per stage;
   - **gate (proposed):** recall at matched precision ≥ tier 1's at every roll-up level, for ≤ 1.2× query cost and ≤ 1.5× index size. Below that, the tree stays a round-2 candidate generator only.
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
