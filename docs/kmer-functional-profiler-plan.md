# Protein k-mer functional profiler — assessment & implementation plan

Sep 26, 2026 · @Tim

## Verdict

The core idea has been done: fmh-funprofiler already does FracMinHash functional profiling on translated amino-acid k-mers against KEGG orthologs, and kMermaid and UProC are older/parallel k-mer approaches. What is not done is your specific combination: a guaranteed per-gene sampling floor, frame-aware read translation, and a proper abundance model. That is a real but narrow niche, and it targets fmh-funprofiler's documented weakness (poor completeness for low-abundance functions).

Keep, change, drop:

- **Keep: translation to amino acids.** Sound and established. Silent-site tolerance is real; the cost is weaker specificity per residue, so k must rise.
- **Keep, but downgrade: frame detection.** Mostly a speed and noise reduction (about 6 frames to \~1.5), not a sensitivity gain. A stop-codon filter gets most of it for free.
- **Keep, but reformulate: length-scaled density.** Do it as a per-gene (or per-function) FracMinHash threshold with nested thresholds. The post-hoc correction is then exact, not approximate.
- **Drop: start/stop anchoring.** It either breaks sketch consistency (position-based selection that reads cannot reproduce) or adds almost nothing (content-based anchors hit few reads). The per-gene floor already delivers what anchoring was meant to.
- **Drop: stacking minimizers + hash threshold + anchors.** Use one context-free sampler. Three layered samplers make the statistics hard and the reads-vs-reference consistency fragile.

The harder problems are not sparsity. They are sensitivity to divergent homologs (exact amino-acid k-mers fail below roughly 80% identity), k-mers shared across functions (domains, paralogs), and database licensing (KEGG is not freely redistributable). The plan below addresses those first.

## Prior art

The closest tool is fmh-funprofiler; everything in your plan except the density control and frame handling has a precedent.

| Tool | What it does | Relevance to your plan |
| --- | --- | --- |
| [fmh-funprofiler](https://github.com/KoslickiLab/fmh-funprofiler) (Koslicki lab, Bioinformatics 2024) | Sketches each KEGG ortholog group and the translated metagenome with sourmash protein FracMinHash, then uses `sourmash prefetch` overlaps as abundances. Default scaled = 1000. | Direct predecessor. Reported 39–99x faster and up to 40–55x less memory than DIAMOND, with comparable completeness and better purity. Its low completeness on rare KOs comes from downsampling, which is exactly your target. |
| [kMermaid](https://pmc.ncbi.nlm.nih.gov/articles/PMC12507277/) (PLoS Comp Biol 2025) | Assigns reads to taxa-agnostic clusters of homologous proteins using amino-acid k-mer frequency profiles in a nested hash map. | Unsparsified k-mer approach; a speed/accuracy baseline and an example of clustering proteins before indexing to reduce multimapping. |
| [UProC](https://academic.oup.com/bioinformatics/article/31/9/1382/200454) (2015) | "Mosaic matching" of amino-acid words against Pfam or KEGG families, including ORF selection from DNA reads. | Early proof that word-based protein classification beats profile HMMs on 100 bp reads. |
| [Orpheum](https://github.com/czbiohub-sf/orpheum) (CZ Biohub) | Picks the coding frame of each read by six-frame translation and the highest reduced-alphabet k-mer Jaccard against a reference proteome. | Your frame-detection idea, done with database containment rather than codon statistics. Evaluated on bacteria in the [metapangenome paper](https://dib-lab.github.io/2021-paper-metapangenomes/). |
| [Metabuli](https://github.com/steineggerlab/Metabuli) (Nat Methods 2024) and its 2026 update | "Metamers": 8-codon k-mers storing both the amino-acid translation and codon identity; exact amino-acid matching, DNA Hamming distance for specificity. The update adds spaced metamers, reduced alphabets and syncmers. | Taxonomic, not functional, but the best template for your k-mer encoding and sampling. Syncmers halved database size and doubled speed at a small recall cost. |
| [FragGeneScanRs](https://github.com/unipept/FragGeneScanRs) | Rust HMM gene prediction for short, error-prone reads. GPL-3.0. | Option for frame calling, especially with indel errors. The licence matters if you link it. |
| sylph | FracMinHash containment with a statistical model for low coverage (taxonomic). | Template for coverage-aware abundance and containment correction. |
| HUMAnN 3/4, DIAMOND, MMseqs2 | Translated alignment. | Accuracy baselines; HUMAnN is the de facto standard for gene-family profiles. |

Nothing I found anchors k-mer selection to start/stop codons or scales density by gene length. There may be a reason: see the next section.

## Critique of each element

### Sketch consistency is the constraint everything must respect

A sketch works only if a read selects the same k-mers the reference selected, using information the read has. Any rule that depends on position within a gene fails this: a 150 bp read usually does not know where the gene starts. Selection must be a function of the k-mer itself (FracMinHash, syncmers) or of a local window the read also contains (minimizers, with edge losses at read ends).

### Start/stop anchoring

- **Position-based version** ("keep k-mers within x residues of the start") is not reproducible from reads. Broken.
- **Content-based version** ("always keep k-mers containing a start/stop") is consistent, but weak. Internal Met is \~2% of residues, so "contains M" is just a biased subsample. Stops appear only if you append `*` to reference proteins; then only reads crossing the stop produce that anchor, a few percent of reads per gene.
- True starts are ambiguous (GTG/TTG, annotation errors), and termini are usually less conserved than core domains.
- In protein space, the reference is already a set of proteins with defined ends. ORF finding on the reference is unnecessary; it only matters on the read side.

### Translation and frame detection

Translation is the right call. The cost: an exact amino-acid k-mer of length 11 survives at identity *p* with probability *p*^11, so 0.31 at 90%, 0.09 at 80%, 0.02 at 70% identity. Functional orthologs are often below 70%. Reduced alphabets and spaced seeds (Metabuli's route) are how you recover sensitivity.

Frame detection by stop codons is cheap and effective at moderate GC. A random 50-codon frame has no stop with probability (61/64)^50 ≈ 0.09, so a coding read keeps \~1.5 of 6 frames. At 70% GC, stop codons are rarer (\~1.9% per codon), about 0.38 of off-frames survive and you keep \~2.9 frames. Frame filtering is a speed and false-positive optimisation, not a sensitivity gain: off-frame k-mers rarely match a protein database anyway.

### Minimizers vs hash threshold

Random minimizers do nest: with a fixed order, a (w=20)-minimizer is also a (w=10)-minimizer, since it is the minimum of every sub-window containing it. But they depend on neighbouring k-mers, so short reads lose picks at their edges, and the density is not an unbiased fraction. FracMinHash is context-free and gives a clean binomial model. Closed syncmers are the middle ground: context-free with a spacing guarantee (Edgar 2021). Use one sampler, not three.

### Length-scaled density and post-hoc correction

The problem is real. At scaled = 1000, a single 300-aa gene has \~0.29 expected sketched k-mers, so it is detected with probability \~25%. fmh-funprofiler hides this by pooling all genes of a KO, and pays in completeness for rare KOs.

The fix works if thresholds are nested. Give each unit *g* (protein cluster or function) its own threshold *t\_g* = max(*t\_base*, *n\_min* / (*L\_g* − *k* + 1)). The query keeps every k-mer with hash < max(*t\_g*). A hit counts for *g* only if hash < *t\_g*. Because {hash < *t\_g*} ⊆ {hash < *t\_max*}, per-unit containment and coverage are unbiased, with *t\_g* as the known sampling rate. No approximation needed.

The catch: the query side runs at the densest rate any unit needs. So the saving is in index size and lookups, not hashing, which you do for every k-mer regardless (SIMD hashing makes that cheap). Stream the query and count hits; do not store a query sketch.

### What the plan misses

- **Shared k-mers.** Multi-domain proteins and paralogs share k-mers across functions. Summing hits double-counts. Needs EM over read-level equivalence classes (salmon/kallisto style) or a gather-style greedy assignment.
- **Normalisation.** Hits per sketch k-mer ≈ k-mer coverage; converting to copies per genome needs a single-copy marker or average-genome-size normaliser.
- **Database.** KEGG requires a licence for bulk data; eggNOG, Pfam, UniRef and KOfam HMM-derived annotations are open alternatives.

## Revised design

One context-free sampler with nested per-unit thresholds, frame-filtered translated reads, and an EM abundance model.

**Reference (index build)**

1. Input: MGnify Proteins 2026\_07 FASTA plus its cluster membership: 90% clusters (MGnify90) as units, Pfam as the functional label (see Database below). No ORF calling.
2. Optional reduced alphabet (e.g. Murphy-10, Dayhoff) and optional spaced-seed mask; amino-acid k of 8–12, packed into u64.
3. Per 90% cluster *g*: *t\_g* = min(1, max(*t\_base*, *n\_min* / *n\_kmers(g)*)) for non-singleton clusters, *t\_g* = *t\_base* for singletons. Keep k-mers with hash < *t\_g*.
4. Store hash → list of (unit id); store *t\_g* and *n\_kmers(g)* per unit. Sorted array or minimal perfect hash, memory-mapped.

**Query (streaming, no stored sketch)**

1. Read FASTQ; for each read and mate, translate the frames that survive the stop filter (default) or all six (`--frames all`); split at stops; drop segments shorter than k.
2. Hash every k-mer; skip if hash ≥ *t\_max*; look up survivors.
3. Per read (pair), collect the set of units hit, counting a hit for *g* only if hash < *t\_g*. Emit equivalence class counts.

**Abundance**

- Per unit, hits ≈ Poisson(λ\_g × *m\_g*), where *m\_g* is the sketch size and λ\_g the k-mer coverage.
- Two-stage assignment: uniqueness-weighted detection, then EM per connected component of clusters that share kept k-mers (see Assignment below).
- Report: coverage λ\_g, containment (fraction of sketch k-mers seen) with a sylph-style low-coverage correction, a copies-per-genome value using a single-copy marker set, and relative abundance.

**Knobs to ablate:** *t\_base*, *n\_min*, k, alphabet, frame mode, sampler (FracMinHash vs closed syncmers), unit granularity.

A `--sourmash-compat` mode using sourmash's protein encoding and hash lets you check your containment numbers against fmh-funprofiler exactly before diverging.

## Database: MGnify Proteins

Use MGnify's own 90% clusters as units and index the k-mers of all members, not just representatives; skip your own 95% dereplication. MGnify90 clusters are the only grouping level: the floor applies per cluster but only to non-singleton clusters, or the index will not fit in memory.

**What the release provides.** The latest release is 2026\_07: 5.74 billion sequences in 1.66 billion 90% clusters (73% singletons), CC0-licensed, so pre-built indexes can be redistributed. Pfam hits are provided for all proteins (52% have one). There is no 30% membership table (MGnify30 ships only as representative FASTAs of three subsets), so this plan uses MGnify90 clusters only; grouping them into 30% families is future work. Details, file sizes and Parquet layout: [mgnify-2026\_07.md](mgnify-2026_07.md).

**95% vs 90%.**

- No 95% clustering is published. Making one means clustering 5.7 billion sequences yourself.
- The identity threshold only changes k-mer content if you index representatives alone. A member at 90% identity shares only \~0.9^11 ≈ 31% of its 11-mers with the representative, so a rep-only index loses most k-mers of divergent members.
- If member k-mers are indexed and labelled with their cluster, the threshold only sets unit granularity. 90% is already finer than any functional question needs; 95% would add units and shared k-mers with no gain.
- If memory forces a cut, subsample member k-mers (lower *t\_base*), not the identity threshold.

**Index size (rough, to verify on real data).**

| Component | Estimate | Assumption |
| --- | --- | --- |
| Base sketch at *t\_base* = 1/1000 | \~5×10^8 hashes, \~6–10 GB | \~5.7×10^9 proteins × \~200 aa, \~half distinct |
| Floor at 90%-cluster level, *n\_min* = 8 | ≥ 1.3×10^10 hashes, \~160 GB | 1.66×10^9 clusters in 2026\_07 (0.45×10^9 non-singletons) |
| Floor at 90% non-singleton clusters only, *n\_min* = 8 | ≥ 3.6×10^9 hashes, \~45 GB | 0.45×10^9 non-singleton clusters |

So: apply the floor to non-singleton 90% clusters only, and give singletons the base rate. If that is still too large, lower *n\_min* or restrict the floor to clusters with a full-length member.

**Functional labels.** Pfam hits are published for every protein, so each cluster can take the Pfam labels of its members. KO or eggNOG labels would need your own annotation run on cluster representatives (0.45×10^9 non-singletons: expensive).

**Practical options.**

- Development index from a biome subset (e.g. human gut) using MGnify's biome counts per sequence: small enough for a laptop.
- A two-stage screen like sylph's new `.syl2db` format: a very sparse first pass picks candidate clusters, a denser second pass quantifies only those ([sylph releases](https://github.com/bluenote-1577/sylph/releases)).

## Assignment of shared k-mers

Uniqueness-first is a good rule for deciding which units are present, but a poor one for quantifying them; use it as the detection stage and EM for abundances.

| Method | Rule | Weakness here |
| --- | --- | --- |
| sourmash gather | Greedy: take the target with the largest remaining overlap, remove its k-mers, repeat. | Order-dependent; favours big units. |
| sylph profile | Winner-take-all: each shared k-mer goes entirely to the genome with the highest containment ANI, then ANI is recomputed on the reassigned k-mers ([sylph paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12339375/)). | Works at species level (>95% ANI, few close neighbours). Dense 90% protein clusters have many near-equal competitors. |
| Uniqueness-first greedy (your idea) | Rank units by unique or least-shared k-mer support; assign those first. | See below. |
| EM on equivalence classes | Shared hits split in proportion to current abundance estimates; unique hits pin those estimates. | More code; needs candidate pruning at MGnify scale. |

**Why pure uniqueness-first greedy is risky for 90% clusters**

- Unique k-mers are not a random sample of a unit. They sit at the variable positions, which is exactly where a sample strain tends to differ from every reference. Abundance from them is biased low for divergent strains.
- Units in dense neighbourhoods (many near-identical clusters) have few unique k-mers, so estimates are noisy just where competition is highest.
- A sample protein sitting between two references (93% to each) should split roughly evenly; any greedy rule gives it all to one.
- With per-cluster thresholds, raw unique counts must be normalised by each unit's sampling rate *t\_g* before ranking.

**Recommended two stages**

1. **Detection (your idea, formalised).** At build time, store per unit its sketch size *m\_g* and unique sketch k-mers *u\_g* (kept k-mers found in no other cluster). At query time, test unique hits against a Poisson background of false hits, and rank by an IDF-weighted score Σ hits(x) / |units(x)|. Keep units that pass.
2. **Quantification.** Within each connected component, EM over the surviving units: hits(x) \~ Poisson(Σ\_{g ∋ x} λ\_g), zero-inflated for divergent k-mers as in sylph. Components (clusters linked by shared kept k-mers, computed at build time after dropping promiscuous k-mers) are independent, so this runs in parallel.

At the sketch densities planned, most reads carry zero or one hit, so read-level equivalence classes collapse to k-mer-level classes. Run the EM on k-mer hit counts; it is simpler and loses nothing.

Pfam totals do not depend on how hits are split between clusters that carry the same Pfam label, so most errors in the within-component split do not reach the level most users will report.

**Evaluation:** all four rules run on the same hit table, so implement gather, winner-take-all and uniqueness-first as cheap baselines in phase 4 and compare against EM in the phase-7 ablations.

## Homologous groups and k-mer weighting

Treat each 90% group as one unit whose k-mer set is the union of its members, and store how often each k-mer occurs inside the group. That in-group frequency, set against the out-of-group frequency, is both the selection criterion and the likelihood weight.

**How fmh-funprofiler does it.** It builds one FracMinHash sketch per KEGG ortholog from the amino-acid k-mers of all its genes (k = 7, 11, 15; scaled = 1000; abundances tracked). The authors tried per-gene sketches and dropped them because genes are too short to keep enough k-mers. They run `sourmash prefetch`, which reports each KO's overlap independently, and derive abundances from those overlaps ([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11373326/)). So there is no within-group weighting and no reassignment of k-mers shared between KOs: a shared k-mer counts for every KO containing it.

**In-group vs out-of-group frequency.** For k-mer *x* and group *G*:

- *p\_in(x)*: fraction of G's members containing *x*. High means core: a read from an unseen member of G probably carries it. Private k-mers (one member) are weak evidence.
- *p\_out(x)*: number of other groups containing *x*. Low means specific.
- Score: log(*p\_in* / *p\_out*), a naive-Bayes log-likelihood ratio. Precedents: CLARK keeps only target-discriminative k-mers; MetaPhlAn picks markers by in-clade core plus out-of-clade uniqueness.
- Compute *p\_in* over full-length members only, with sequence weights (as in profile HMMs) so over-sampled lineages and partial predictions do not skew it. Quantise it to 2–4 bits per posting.

**Score-based selection is free on the query side.** The query checks every k-mer with hash < *t\_max*. So the index can keep any score-chosen subset of those k-mers, and queries stay consistent. Rules:

1. Among k-mers passing the hash threshold, keep the highest-scoring ones per group until the floor *n\_min* is met.
2. Drop k-mers present in more than N clusters outright (low-complexity runs, common motifs such as Walker A). This also removes the longest posting lists.
3. Record which k-mers were kept, so each group's expected hit count uses the kept set, not a nominal rate.

**Likelihood with *p\_in*.** hits(*x*) \~ Poisson(λ\_G · *p\_in(x)*) turns strain variation into part of the model rather than zero-inflation noise. It also yields a per-group divergence estimate: a containment-AAI analogous to sylph's containment ANI. That is a useful novelty signal per function.

**Tiered lookup ("test the most unique first").**

1. A small tier-1 index of high-score k-mers (high *p\_in*, low *p\_out*), a few per group, stays in RAM. Every query k-mer probes it.
2. Groups with tier-1 support above background become candidates.
3. Only candidates' full k-mer sets (tier 2, memory-mapped) are probed, via a second pass over the reads or buffered query hashes.

This is the same pattern as sylph's two-stage database. At MGnify scale it is probably what makes the index usable on a normal node.

## Index encoding

Store \~12–16 bits per hash, with the high bits implicit in the layout, and spend the real effort on compressing unit IDs, which will dominate the index.

**Fingerprints.** A bare 16-bit key fails: \~5×10^8 hashes over 2^16 values puts thousands of k-mers on every value, and no model recovers the signal. What sets the false-hit rate is the effective hash length:

- Partition by the high bits (e.g. 2^24 buckets addressed directly) or use Elias–Fano over sorted hashes; entries store only the low bits.
- Hashes below the FracMinHash threshold share log₂(1/*t*) leading zero bits for free.
- Stored bits per entry ≈ log₂(1/ε) + 2, for false-hit rate ε per lookup.

| Stored bits per entry | False-hit rate ε per lookup |
| --- | --- |
| 8 | \~1.6×10^-2 |
| 12 | \~10^-3 |
| 16 | \~6×10^-5 |

**Statistics.** A spurious hit lands on a random entry, so group *g* gets an expected Q·ε·*m\_g*/n false hits, where Q is the number of query lookups, *m\_g* the group's kept k-mers and n the index size. Add this as a known Poisson background in detection and EM, and verify it with shuffled or decoy reads.

**The floor that matters more.** Genuine chance matches already occur. With 20 amino acids and k = 11 there are \~2^47.5 possible k-mers; MGnify might hold \~2^39 distinct ones (rough), so \~0.3% of off-target query k-mers hit something by coincidence. A 16-bit fingerprint adds far less. Reduced alphabets shrink the space (Murphy-10 at k = 11 gives \~2^36.5, near saturation), so they need a larger k.

**Unit IDs.** With \~10^9 units an ID costs \~30 bits, and conserved k-mers carry many. Map each k-mer to the ID of its distinct unit set and store each set once, as kallisto, Themisto and Fulgor do; sets repeat heavily. Store component IDs in tier 1 and resolve units only in tier 2.

**Rules.** Compare thresholds on the full 64-bit hash before truncating, so nested per-cluster thresholds stay exact. Keep the fingerprint width a build parameter and measure ε against index size.

## Query cost

Querying at the densest per-unit rate (*t\_max* = *t\_cap* = 0.2) is affordable if lookups are sorted and merged instead of random, and if the query keeps per-(unit, k-mer) counts instead of per-read rows. The density raises the size of the sample sketch, not the cost of hashing, which is paid for every k-mer at any rate.

**Why the query is dense.** A read cannot tell which unit it came from, so it must test every k-mer that any unit could have kept: hash ≤ max *t\_g*. 43% of non-singleton clusters reach *t\_cap* (phase 6, step 3), so *t\_max* = 0.2: 200× fmh-funprofiler's rate (scaled 1000), and 200× what the 73% singletons need. No context-free sampler avoids this; score-based selection keeps the index small but not the query.

**Rough cost per 10 Gbp sample** (6.7×10⁷ reads of 150 bp, ~1.5 stop-free frames, ~60 11-mers per read, full index at the defaults; to be measured):

| Stage | Scales with | Estimate | Now (Python prototype) |
| --- | --- | --- | --- |
| Decompress, translate, hash | reads × frames | ~4×10⁹ hashes, ~1–2 min on 12 cores (kernel ~54 M reads/min); gzip decompression single-threaded unless bgzf | Rust, fine |
| Sampled k-mers | × *t\_max* | ~8×10⁸ at 0.2; distinct perhaps half (error k-mers and low coverage are mostly distinct) | streamed per batch |
| Lookup in tier 2 | sampled k-mers, and index residency | ~31 GB tier 2 + unit thresholds. Resident: 2–4 cache misses per lookup, ~1 min on 12 cores. Not resident (cold mmap on a network file system): a page fault per lookup, hours | numpy `searchsorted` on 8 B per key materialised (~27 GB) |
| Hit rows | hits × set size | hit fraction unknown on MGnify (10–50%?): 10⁸–4×10⁸ rows | one polars row per (unit, hash, read), ~25 B each: 3–10 GB |
| Gather | hit (unit, k-mer) pairs | 10⁶–10⁷ units, ~10⁸ pairs | Python sets and a heap: slow and tens of GB at this size |
| Model fits | detected (unit, k-mer) pairs × iterations | em, zi, zib, zip all fitted, one reported | 4 fits |
| Posterior draws (*D* = 100) | pairs × *D* × sweeps | 5–10× all the rest (fmh benchmark: 23 s → 2m17) | per-read rows needed for the read resampling |
| Dense tier | second read pass + dense table | 15× tables, ~6× query time (phase 4, step 9) | hash-major table, same lookup |

So at full scale the cost is set by three things, in order: per-read hit rows and Python model code (memory), index residency (time), and the posterior (time). Hashing and the sampling rate come after.

**How sylph keeps profiling cheap** (source read at v1.0.0):

- *Sketch the sample first.* Reads become a sample sketch (`.sylsp`): distinct sampled k-mers with counts, at c = 200. Everything downstream works on distinct k-mers, not on read occurrences, and a sketch is reused across databases.
- *Two stages (`.syl2db`, v1.0).* Stage 1 is a pooled minimal perfect hash (boomphf) over all genomes' sparse k-mers (screen c = 3000), with a multi-owner CSR from k-mer to genomes: loaded in full, probed once per distinct sample k-mer, and permissive (85% ANI). Stage 2 stores each genome's dense k-mers as its own Golomb-Rice block at a known offset; only screen survivors are decoded, by positional reads (`pread`) rather than mmap, and decoded blocks are not cached, so RSS stays bounded by the survivors. The screen's k-mer floor is scaled by c / screen c so small genomes are not dropped before the dense stage.
- *Cuckoo filters* are used for read deduplication (a scalable cuckoo filter over (k-mer, read-pair markers), approximate at a false-positive rate of 10⁻⁴), to save memory against an exact hash set. It is an accuracy feature for PCR duplicates, not a lookup structure.

**What transfers, and what does not.**

- *Sketch-first transfers directly*, with one change: the posterior resamples reads, so pass 1 keeps counts per distinct hash only and read ids come back in an optional second pass limited to the detected units' k-mers (as the dense tier's second pass already is). Pass 1 memory is then the sketch, not the hits.
- *Sorted merge instead of random lookup.* Tier 2 is hash-ordered, so a batch of sampled hashes sorted by hash (radix sort) is resolved by one forward walk through the tier-2 arrays: every page read once, in order, whether or not the index is resident. Radix-partitioning the sampled hashes by the build's pack hash ranges (spilled to local disk when a memory budget is set) makes the query an external merge-join: one sequential pass over the index per sample, with memory = one range of the sketch plus one range of tier 2. With the index in the page cache (node-local disk or `/dev/shm`), the same walk is just cache-friendly.
- *Cohorts.* The merge-join takes the k-way merge of many samples' sorted sketches, so one index pass serves a batch of samples. On HPC this is the cheapest mode for large studies; separately, jobs on one node share one page-cached copy of an mmapped index.
- *Two-stage screening transfers only partly.* Sylph's screen works because genomes are long: at c = 3000 a 2 Mbp genome still has ~700 screen k-mers. A floored 90% cluster has ~8 k-mers at *t\_g* up to 0.2, so any sparser screen cannot see it on its own; screening at a coarser unit (connected component, Pfam) would pass a rare short cluster only when relatives are present too. That gives up exactly the low-abundance completeness the floor exists for, so a screen is an opt-in fast mode for ablation, not the default.
- *Sylph's stage-2 layout is the right one for the dense tier.* The dense tier is only read for detected units, so it should be unit-major like `.syl2db`: per-unit blocks of sorted hash fingerprints (Golomb-Rice or Elias-Fano, with *p\_in*), read with `pread` for detected units only, loaded into a small hash map that pass 2's reads probe at the dense rate. Its size then costs disk, not RAM, and its query cost scales with detected units, not the table: the reason the dense tier is off by default (15× tables, 6× time) mostly goes away.
- *Cuckoo or Bloom prefilter.* Tier 2 is already a quotient-style filter (implicit bucket bits plus fingerprints), so a separate cuckoo filter only helps if misses dominate and the index is not resident. The blocked Bloom filter from the partitioned build (`_core.bloom_contains`) can be tried as an in-RAM prefilter (~4 GB at 10 bits per key) in that regime; with sorted merges it is unnecessary. Read deduplication with a cuckoo filter is an accuracy option to ablate, not a cost saving.

**Low memory (a 12 GB laptop).** The partitioned merge-join is the laptop mode. `--memory` sets the partition count, and the cost of a small budget is scratch disk and one sequential index read, not many passes. There is one code path: a node with the index in RAM just runs it with one partition.

- *Shard passes vs disk lookups with caching.* Random lookups against an index on disk lose. After sketch dedup, a 10 Gbp sample leaves ~4×10⁸ distinct lookups, spread uniformly over ~50 GB of index. A page cache of ~8 GB then hits ~15% of them, and an NVMe drive at 10⁵–5×10⁵ random 4 KB reads/s needs 15 min to 1 h+ (macOS and thread count dependent). Caching does little beyond what dedup already does: repeated k-mers of high-coverage genes are merged in the sketch, and what remains has no reuse. Reading the index shard by shard, in order, reads every page once: 50 GB at 2–3 GB/s is under a minute.
- *Spill the sketch, not the reads.* Pass 1 radix-partitions sampled hashes into *R* scratch files by the index's pack hash ranges: 8 B per sampled k-mer, ~6.4 GB for 10 Gbp (16 B, ~13 GB, when read ids are kept for the posterior). Each partition is then sorted, counted and merged with its slice of tier 2, which is mmapped and read sequentially (`MADV_SEQUENTIAL`, pages dropped after), so it sits in the page cache and not in RSS. With *R* = 8, a partition is under 1 GB. Re-reading the FASTQ once per index shard instead needs no scratch space, but costs a decompress-and-hash pass per shard (~2–3 min each on a laptop); keep that only as the fallback when scratch disk is short.
- *Everything after the lookup must also stream.*
  - Hits come out per hash range but the model works per unit, so per-(unit, hash) counts are spilled by unit range.
  - The `max_hash_g` check moves after aggregation and becomes a sorted join with the unit table. That table is mmapped by column and read in unit order, never randomly.
  - Gather, EM and the posterior run one component at a time, reading only that component's units. Peak memory is then set by the largest component, which Q1 has to measure: promiscuous k-mers are cut at 64 clusters, but a giant component is possible.
- *The dense tier needs no second read pass* when *t\_dense* ≤ *t\_max*, because the spilled sketch already holds every dense-rate hash. Detected units' dense blocks are read, re-partitioned by hash range and merged with the same sketch partitions. A second read pass is needed only for *t\_dense* > *t\_max* or for read ids.
- *Rough laptop budget, 10 Gbp:*
  - hashing ~2–3 min (8 cores);
  - spill write and read ~6–13 GB (seconds on NVMe);
  - sorting ~10–20 s;
  - index read ~50 GB, ~0.5–2 min;
  - model: to be measured.

  That is about 1.5–2× the in-RAM time, for ≤ 8 GB RSS. Needed on disk: the index (~50 GB after phase 6's free cuts) plus scratch.
- *Smaller indexes are the other lever.* These trade accuracy for size, unlike the merge, which does not:
  - biome subsets (e.g. human gut), which may fit in RAM whole;
  - a lighter laptop index (*t\_cap* 0.1, *n\_min* 4, singletons out).

**Rules.**

- The shipped query has no per-read state in pass 1; read ids exist only in pass 2, for detected units' k-mers, and only when the posterior or the dense tier needs them.
- One estimator is fitted by default (the one phase 7 picks); the others and the baselines (em, zib, zip, winner-take-all, uniqueness-first, `kmers_out`) are benchmark options.
- Components are independent throughout (gather, EM, posterior), so they run in parallel, and a component with one unit and no shared k-mers needs no Gibbs sweeps: nothing is split, so its interval comes from a 2-D grid over (λ, π) or the read resampling alone.
- Every query records time and peak RSS per stage, so ablations report cost beside accuracy.

## Error vs variation model

Adopt sylph's separation of sequencing error from divergence, with two changes: a mixture-aware likelihood, and a dense second tier so each group has enough k-mers to fit.

**What transfers.** Errors thin every k-mer's count roughly equally and lower the effective coverage λ. Divergence removes a k-mer from all reads and shows up as a structural zero. Estimating λ from the counts of present k-mers gives:

- abundance that stays correct when the sample variant is only \~90% identical to the reference (total hits divided by sketch size understates it);
- coverage-corrected containment, and from it a containment-AAI per group.

In protein space errors are gentler than in DNA: \~25% of random substitutions are synonymous, and a stop-creating error just splits the read. For Illumina this all folds into λ and cancels in relative abundances.

**Model.** For group *G* and kept k-mer *x*:

- *z\_x* \~ Bernoulli(*π\_x*), with *π\_x* informed by *p\_in(x)*;
- count(*x*) | *z\_x* = 1 \~ NegBin(mean λ\_G · *r*, dispersion φ\_G), where *r* is the error-thinning factor; count = 0 if *z\_x* = 0;
- empirical-Bayes priors on λ\_G and on divergence, shared across groups, to stabilise sparse groups.

**Mixtures.** Groups are often covered by several strains at different coverages, so core k-mers get summed coverage and variant k-mers partial coverage. The negative binomial absorbs mild mixing. A dispersion test per group flags groups whose counts don't fit a single Poisson; the flag is reported as "multi-variant". A small latent-variant mixture is a later option if the flag turns out to be common.

**Enough k-mers.** A sketched group has tens of k-mers, too few below \~1x where nearly all counts are 1. Tier 2 is therefore dense (or near-dense) for candidate groups only, giving hundreds of k-mers per fit at bounded cost. Report confidence intervals, not point estimates, at low coverage.

**Absolute coverage (optional).** Estimate *r* from read quality scores, or from near-invariant single-copy marker groups where divergence ≈ 0.

**Out of scope for now.** Long-read indels cause frameshifts, which remove every downstream k-mer in that frame. That is not uniform thinning and needs a frameshift-aware model.

## Development strategy: Python first, Rust at the end

Prototype the algorithm in Python, with stable hot loops in Rust from day one via PyO3, then port the rest to Rust once evaluation has frozen the design. The Python version stays in the repo as the reference oracle for the Rust port.

| Component | Language during development | Why |
| --- | --- | --- |
| FASTQ streaming, translation, stop-filter frames, amino-acid k-mer packing, hashing, threshold filter | Rust from phase 1 (PyO3 module) | Hot, well-specified, unlikely to change. Pure Python would make every experiment on real data too slow. |
| Index lookup (sorted hashes, fingerprints) | Rust once the layout settles (end of phase 3); memory-mapped in phase 6 | Memory and speed matter on MGnify-scale subsets. |
| Index build: clustering tables, p\_in/p\_out scoring, floors, tiers | Python (polars, DuckDB) through phase 5; Rust or a partitioned build in phase 6 | Where most design changes happen, but the Python build cannot hold the full release in memory, so it is the component that blocks full-scale work. |
| Assignment, EM, zero-inflated negative-binomial model, dispersion test | Python (numpy, scipy) | Statistical design still open; easiest to iterate and inspect. |
| Evaluation, plots, ablations | Python + Nextflow | Stays Python permanently. |
| Final CLI and full pipeline | Rust (phase 8) | Release product. |

**Rules for the split**

- The Rust kernels expose batch APIs that return numpy arrays (e.g. read index, frame, hash per kept k-mer), never per-k-mer Python calls.
- Every Rust kernel has a slow pure-Python twin under `reference/`, used only in tests; property tests check they agree.
- Before porting (end of phase 7), write a short algorithm spec: parameters, formulas, file formats. The port implements the spec; the Python outputs on fixtures become golden files.
- A Python component moves to Rust early only if profiling shows it blocks experiments, and its interface has not changed for a phase.

## Implementation plan

Ten phases, each with a go/no-go gate, plus an optional eleventh; phases 1–5 are the Python prototype with Rust kernels, phase 6 makes the method run on all of MGnify Proteins at a reasonable cost, phase 7 runs the ablations against that full-scale method and freezes it, phase 8 ports the rest to Rust, and phase 10 (desirable, not essential) lets users add their own proteins to a released index. Phase 3 is the first point where the prototype must at least match fmh-funprofiler, or the project should stop.

| Phase | Deliverable | Language | Gate |
| --- | --- | --- | --- |
| 0. Skeleton | Mixed repo (uv + maturin + Cargo workspace), CI, pre-commit, stub Python CLI. Fixtures, the parity job, coverage and `bench.yml` move to phase 1, when there are kernels to test. | Both | CI green on Linux and macOS arm64 for both languages. |
| 1. Rust kernels (done) | PyO3 module: FASTQ streaming, codon tables (11, 4), six-frame translation, stop-filter frames, reduced alphabets, amino-acid k-mer packing, hashing, FracMinHash filter; batch numpy outputs. Pure-Python reference twins. | Rust + Python tests | Property tests pass (frame symmetry, threshold nesting, synonymous invariance); Rust matches reference; ≥ 1 M reads/min/thread. |
| 2. Index prototype (built; gate needs a real subset) | Build from a MGnify biome subset (DuckDB/Parquet): member k-mers per 90% group, p\_in/p\_out scores, per-cluster floor (non-singletons), connected components, tier 1/tier 2, fingerprints; stored as Parquet + numpy. | Python | Sizes match the analytical estimates; scoring behaves on hand-checked clusters. |
| 3. Query + naive counts (done; gate passed, see Progress log) | Hit counting and containment; `--sourmash-compat` using the sourmash Python API. Lookup moves to Rust once the layout settles. | Python (+ Rust lookup) | Matches fmh-funprofiler containment in compat mode; ≥ parity in completeness/purity at equal density. **Stop here if not.** |
| 4. Model (in progress: baselines, EM, zero-inflated EM adopted, dense tier, *p\_in*-weighted presence, copies-scaled abundance, posterior intervals, ambiguity groups, presence probability, copies error calibrated on real data; simulation benchmark; see Progress log) | Uniqueness-weighted detection, EM per connected component, zero-inflated negative-binomial model, dispersion flag, genome normalisation, dense tier 2. | Python | Clear completeness gain for low-abundance groups over phase 3 at ≤ 2x index size; calibrated intervals on simulations. |
| 5. Tool benchmarks (in progress: DIAMOND, fmh-funprofiler, kMermaid, HUMAnN 3.9 and 4 built, results pending on HPC; see Progress log) | Benchmarks vs fmh-funprofiler, HUMAnN, kMermaid, DIAMOND on the fmh benchmark, with CPU time and peak memory per tool. Pfam truth and scoring on the same metagenomes (the primary benchmark, see Benchmark labels); KO scoring kept for the direct comparison with fmh-funprofiler. | Python + Nextflow | Every tool scored with the same truth and metrics, on Pfam where it can report Pfam and on KO otherwise; our accuracy and cost relative to each recorded. |
| 6. Full-scale method (in progress: nested builds done and cost model checked; free index cuts, Rust pass-1 kernel and partitioned build (tier 2 packed by hash range, sets deduplicated by set hash) done; whole-release statistics pending on HPC; query cost planned (see Query cost), Q1 ladder under way on HPC; Q3's Rust tier-2 lookup built; unit table memory-mapped by column; see Progress log) | Make the method run on all of MGnify Proteins at a reasonable compute cost, before any ablation. Full-scale cost study: nested all-biome MGnify subsets (1 in 10⁴ to 1 in 10 clusters) and a full-release statistics pass (per-cluster k-mer counts, predicted postings) to fit how storage, build and query cost scale. Query cost: per-stage query time and memory against the nested and full indexes, a Rust sorted-merge lookup, a unit-major dense tier and cohort mode (see Query cost). Scalable build: index build and lookup in Rust (ported ahead of the spec, as the rules above allow for a component that blocks experiments) or a hash-partitioned Nextflow build; memory-mapped lookup. Tune for cost: *t\_base*, *n\_min*, *t\_cap*, dense-tier rate and scope (e.g. non-singletons only), fingerprint width, unit-ID encoding; re-run the Pfam and simulation benchmarks at each candidate to choose defaults on cost vs accuracy (KO alongside, for comparison). Then the full build and a query of the fmh-benchmark metagenomes against it, scored on Pfam truth (the full MGnify index carries Pfam labels only). Expect many iterations and some accuracy given up for cost. | Rust + Python + Nextflow | Full MGnify index builds on one HPC node at a chosen cost/accuracy point (build time, index size, query memory and time recorded), with the accuracy given up versus phases 4–5 recorded. |
| 7. Ablations and freeze | Ablations with the phase-6 method as the baseline, so each measures a change against what will ship: EM vs gather/winner-take-all/uniqueness-first, zero inflation, *p\_in* weighting, dense tier, floors, alphabet and k; divergence ladder. Where an ablation changes the index, it is run on a nested subset whose accuracy phase 6 has tied to the full build. Algorithm spec written; golden outputs recorded. | Python + Nextflow | Defaults chosen on Pfam metrics at full-scale cost; spec reviewed. If the divergence ladder forces a different alphabet or k, phase 6's cost study is repeated for it. |
| 8. Rust port | Query, model and CLI in Rust, implementing the spec (the index build and lookup already ported in phase 6 are brought in line with it). Differential tests against the Python golden outputs (exact for counts, tolerance for EM). | Rust | All golden tests pass; ≥ 10x Python end to end; full-scale results of phase 7 reproduced. |
| 9. Release | Rust binary via cargo-dist, bioconda recipe, Nextflow module for the hybrid profiling pipeline; optional Python wheel of the bindings. | Rust | Tagged release reproduces phase 7 results. |
| 10. Reference extension (desirable, not essential) | `extend`: add a study's own proteins (MAG protein predictions, assembly gene calls with their functional annotations) to a released MGnify index as an overlay, without rebuilding it; query base + overlays together; `compact` folds overlays into a new base. See Extending the reference below. | Rust (+ Nextflow module) | On the divergence ladder, adding the held-out genomes' proteins as an overlay recovers ≥ 90% of the completeness of a full rebuild that includes them, with no loss on units already in the base; extending with one study (~10⁶ proteins) takes minutes on one node, not a rebuild. |

Out of scope initially: long reads (indels break frames; would need FragGeneScan-style frameshift handling), eukaryotic genes, metatranscriptomes.

## Progress log

What each step did, and the choices, results and interpretations behind it, newest phase last. Add an entry with every step; keep superseded choices and say what replaced them.

### Project-wide

* **Decided:** MGnify90 clusters are the only grouping level; no 30% families.
* **Decided:** package name `kmer_functional_profiler` (tool name may still change before release); licence GPL-3.0-or-later, so FragGeneScanRs can be linked.
* **Full-scale runs:** anything over the whole release (subset extraction, index build) ships as a Nextflow pipeline with README and setup scripts for HPC; local work uses samples only.

* **Decided (2026-09-28): full-scale check in phase 6, not phase 4.** No test or benchmark so far uses the full MGnify release: unit tests use fixtures, the simulation 300 synthetic units from 100 MGnify seed proteins, the fmh benchmark KEGG KO indexes (tier 2 72 MB, dense tier 1.09 GB at 1/10), and the only real MGnify build is the pre-fix gut 1-in-10,000 subset (213,645 proteins, about 4×10⁻⁵ of the release). The Python build holds the members table and k-mer tables in memory, so a full-release build needs the Rust port (or a partitioned build) anyway. Phase 6 therefore carries the full-scale cost study and tuning (see the phase table). Expect a lot of tuning there to bring cost down, and some compromises on accuracy; defaults chosen in phase 5 are provisional until phase 6 has priced them, and the phase-4/5 accuracy numbers are the reference that any cost saving is measured against.
* **Decided (2026-09-29): full-scale method before ablations; phases re-ordered.** Supersedes the phase-6 placement above. Ablations measured on the prototype at small scale could pick defaults that the full-scale cost tuning then overturns, so the full-scale method now comes first. Phase 5 keeps only the tool benchmarks; the new phase 6 (full-scale method) takes the cost study, scalable build, tuning and full build from the old phase 6; ablations, divergence ladder, spec and golden outputs move to phase 7 and run against the phase-6 method; the Rust port of query, model and CLI becomes phase 8 and release phase 9. The index build and lookup move to Rust (or a partitioned build) in phase 6, ahead of the spec, because they block full-scale work; phase 8 brings them in line with the spec.
* **Decided (2026-09-29): reference extension is a desirable, non-essential objective (phase 10).** Users should be able to add their own proteins (a study's MAGs and assembly functions) to the MGnify90 index. It does not gate the release, but earlier phases should not rule it out: see the constraints in Extending the reference.
* **Decided (2026-09-29): Pfam is the primary benchmark label; KO is secondary.** The tool labels MGnify90 clusters with Pfam, so defaults and gates from phase 5 on are chosen on Pfam truth. The KO benchmark stays, for direct comparison with fmh-funprofiler (which profiles KOs) and with the phase-3 to phase-5 KO results. The fmh benchmark gains Pfam truth from the same genomes and Pfam scoring for our indexes, DIAMOND and HUMAnN; fmh-funprofiler and kMermaid stay on KO. Details: Benchmark labels. Earlier KO-based decisions (phases 3–4) stand until the Pfam benchmark has rerun them.

### Phase 1

* **Decided (phase 1):** k-mer hash is the splitmix64 finalizer of the bit-packed k-mer (5/4/3 bits per residue for protein/Murphy-10/Dayhoff, so k ≤ 12/16/21), a bijection on u64; threshold rule is keep iff hash ≤ max\_hash, with max\_hash = ⌊*t*·2^64⌋ − 1. The sourmash hash is used only in `--sourmash-compat`.

### Phase 2

* **Decided (phase 2):** the build hashes members three times (distinct k-mers per unit; candidate hashes at *t\_max*; presence of candidates in every unit) so no dense k-mer table is held. *p\_in* is the unweighted fraction of full-length members (all members if a unit has none); sequence weights are deferred until over-sampled lineages are seen to skew it. Score = log2(max(*p\_in*, 0.5/*n*) / *n\_groups*), *n\_groups* counting the unit itself. *p\_in* is quantised to 4 bits in tier 2 (value = unit << 4 | *p\_in* level). Every candidate is kept (no score-based trimming beyond the promiscuity cut, default 64 groups); tier 1 is the top 4 by score per unit, mapped to component ids. Fingerprint layout: leading zero bits implicit, ⌈log2 keys⌉ bucket bits addressing an offsets array, 16 fingerprint bits stored; colliding keys share the union of their unit sets.
* **Phase 2 subsets:** member sequences are spread over the whole 1 TB sequence file, so no real subset fits a laptop download; `workflows/mgnify-subset` extracts one on HPC (default: human gut, 1 in 1000 clusters).
* **Phase 2 gate (gut, 1 in 10,000 clusters: 213,645 proteins, 24,117 clusters):** candidates matched Σ *t\_g*·*n\_kmers* (67,652 vs 67,493); shared k-mers link clusters with a common Pfam (562 of 574 labelled). Four problems, fixed on `index-fixes`:
  1. *Query density.* 90% clusters have few distinct k-mers (median 114), so the floor applied to 8,136 of 8,141 non-singletons and *t\_max* reached 0.57. **Decided:** *t\_g* = max(*t\_base*, min(*t\_cap*, oversample·*n\_min*/*n\_kmers*)) with *t\_cap* = 0.2, so the query samples 20%; clusters with fewer than 40 distinct k-mers get fewer than *n\_min* (430 of 8,141 non-singletons, 4 with none). A 0.05 cap was tried first: 5% sampling but 3,322 clusters below the floor.
  2. *Layout.* 19.4 B per hash, mostly one offsets bucket per key. Now ~4 keys per bucket and the smallest unsigned dtype per array: 7.7 B per hash on the subset (small unit ids; larger ids at full scale need wider values).
  3. *Floor k-mers chosen by hash only.* 30% had *p\_in* = 0 (fragments only). **Decided:** floored clusters sample at 4× (`oversample`) and keep their *n\_min* best-scoring candidates; median *p\_in* rose from 0.13 to 0.87, *p\_in* = 0 fell to 1%.
  4. *Adapter artefacts.* The most shared k-mers were six-frame translations of Illumina TruSeq/Nextera adapters and P5/P7 ends in MGnify proteins (4.1% of members). **Decided:** mask every 6-mer of those translations with X before hashing (shuffled-protein control: 0.08% hit). After masking no k-mer exceeded 64 clusters and the largest component fell from 59 to 6 clusters.

  Rebuilt subset (*t\_cap* = 0.2): 66,337 postings (candidates 210,658 vs 210,940 expected), about 3.7×10⁹ for the full release, so roughly 30–45 GB depending on value widths at full scale.

### Phase 3

* **Decided (phase 3):** the query samples at *t\_max*, looks hashes up in tier 2 and counts a hit for a unit only if hash ≤ *t\_g* of that unit, which also drops most fingerprint false hits. Output per unit: hits, distinct k-mers hit, reads hit, containment (k-mers hit / *m\_g*) and coverage (hits / *m\_g*).
* **Done (phase 3, compat):** `import-sourmash` loads protein signatures (e.g. fmh-funprofiler's KO sketches) as an index; queries on it hash reads with sourmash. On fmh-funprofiler's demo reads (k = 11, scaled = 1000) it reports the same 211 KOs with the same overlaps as `sourmash prefetch`, so the same abundances (fmh-funprofiler's abundance is `f_match_query` = overlap / query hashes, renormalised). Our containment is raw (overlap / *m\_g*); sourmash's `f_query_match` divides it by 1 − (1 − 1/scaled)^(*m\_g*·scaled).
* **Benchmark data (phase 3 gate):** fmh-funprofiler's paper does not publish its simulated metagenomes. It publishes 64-genome-scale inputs to regenerate them with CAMISIM (Zenodo 10055954, CC-BY: KEGG genomes 9.1 GB, `protein_ref_db_giant.faa` 3.2 GB, `present_genes_and_koids.csv` 0.26 GB) and its mean metrics (sourmash k = 11: purity 0.98, completeness 0.61 at 0% error).
* **Benchmark pipeline:** `workflows/fmh-benchmark` builds our KO indexes from the same KEGG proteins as the sketches, simulates metagenomes with InSilicoSeq (64 genomes, ~1 Gbp, lognormal abundances) instead of CAMISIM, derives KO truth by mapping reads back with minimap2 (a KO is present if a read overlaps one of its genes, the paper's rule) and scores purity and completeness for fmh-funprofiler (compat) and our indexes. Runs on HPC; results under Phase 3 below and in `workflows/fmh-benchmark/README.md`.
* **Phase 3 gate: passed (10 InSilicoSeq novaseq metagenomes, 64 genomes, ~1 Gbp, k = 11).** Means over seeds 1..10, sd ≤ 0.015; full table in `workflows/fmh-benchmark/README.md`.

  | Index | Count | min\_hits | Purity | Completeness | Low 25% |
  | --- | --- | --- | --- | --- | --- |
  | fmh\_compat (fmh-funprofiler, scaled 1000) | kmers\_hit | 1 | 0.975 | 0.688 | 0.295 |
  | kfp\_s1000 (ours, same density) | kmers\_hit | 1 | 0.985 | 0.677 | 0.285 |
  | kfp\_s1000\_floor8 | kmers\_hit | 1 | 0.976 | 0.719 | 0.404 |
  | kfp\_s100 (*t\_base* 0.01) | kmers\_hit | 1 | 0.951 | 0.960 | 0.853 |
  | **kfp\_s100** | **kmers\_unique** | **1** | **0.985** | **0.955** | **0.837** |

  *Interpretation.* At equal density we match fmh-funprofiler (parity, as the gate requires); the floor adds 4 points of completeness and 11 on the least-covered quarter. The real gain is density: scaled 100 lifts completeness from 0.69 to 0.96. Its false positives were mostly modular PKS/NRPS KOs whose shared domains carry identical k-mers, so a hit on a shared k-mer is weak evidence.

  **Decided:** add gather-style reassignment (`query.gather`: repeatedly take the unit with the most unassigned hit k-mers, scaled by 1/*t\_g*) as the phase-4 detection baseline; profiles report `kmers_unique` beside `kmers_hit`. Gather removes 70% of false positives (348 → 105 per sample) for 31 true KOs lost, and adds no measurable query time. Recommended setting: `kfp_s100`, `kmers_unique` ≥ 1. `min_hits` > 1 trades too much low-coverage completeness (0.49 at 2 on scaled 1000). The floor does not help at scaled 100 (within 0.005 of *n\_min* 0), so `kfp_s100_floor8` was dropped from the defaults.

  *Caveats.* InSilicoSeq, not CAMISIM, so not directly comparable with the paper's 0.98 / 0.61; KO units from KEGG proteins, not MGnify90 clusters. The phase-2 gate on a real MGnify subset (index size at scale) is still open.

### Phase 4

* **Phase 4, step 1 — EM quantification (`query.em`).** On the units gather keeps, each hit k-mer's count is Poisson with mean Σ λ\_g over the units holding it, and every kept k-mer of a unit (*m\_g*, hit or not) counts in its expectation; the multiplicative EM update gives the MLE. Profiles gain `coverage_em` (0 for units gather drops). Choices: EM on k-mer counts, not read equivalence classes (see Assignment); all components in one sparse product, since they share no k-mers and so update independently (per-component parallelism only matters in the Rust port); stop at max relative change 10⁻⁶ or 1,000 iterations. Tests: a two-unit case against the closed-form MLE, and Σ λ\_g·*m\_g* = observed hits on the fixture reads.

  *Not yet:* zero-inflation and *p\_in* weights, fingerprint false-hit background, uniqueness-first and winner-take-all baselines, dispersion flag, genome normalisation, dense tier 2. *Result pending:* abundance accuracy has not been benchmarked; the fmh benchmark scores detection only (truth has per-KO read and base counts to compare against).

* **Phase 4, step 2 — abundance scoring in the fmh benchmark.** TRUTH adds per-KO `depth` = Σ over the KO's genes of aligned bases / gene length (read depth × copies, the quantity λ should track up to a constant). SCORE pairs each detection count with an abundance estimate: `coverage` (hits / *m\_g*) with `kmers_hit`, `coverage_em` with `kmers_unique`. On the detected KOs it reports `spearman_tp` (rank correlation with depth over true positives, so it measures quantification separately from detection) and `l1` (L1 between relative abundances over true ∪ detected KOs, so misses and false positives count). TRUTH now takes `bench.py` as an input, so `-resume` recomputes the truth tables (one minimap2 pass per sample). Checked on the test profile only; *result pending* on HPC.

  *Expectation to test.* A KO unit is the union of k-mers of every KEGG gene with that KO, and a sample carries only a few of those genes, so λ = hits / *m\_g* is diluted by the share of the unit's k-mers the sample holds, which varies widely between KOs. Plain EM does not correct this; the zero-inflated model (λ from present k-mers, as in sylph) should. If `spearman_tp` is low for both estimates, that is the case for doing zero-inflation next.

* **Phase 4, step 3 — winner-take-all and uniqueness-first baselines (`query.assign_best`).** With a fixed score per unit, both rules reduce to one pass: each hit k-mer goes to the holding unit with the highest score (ties to the lowest id). Gather differs only in that its scores fall as k-mers are taken.
  - *Winner-take-all* (sylph): score = containment, k-mers hit / *m\_g*. One pass; sylph's recomputation of containment on the reassigned k-mers is not repeated.
  - *Uniqueness-first*: score = Σ over the unit's hit k-mers of 1 / (units holding the k-mer), an IDF weight, divided by *t\_g* as this plan requires. The plan's Poisson test of unique hits against a false-hit background waits for that background model.

  Profiles gain `kmers_wta`, `coverage_wta`, `kmers_ufirst`, `coverage_ufirst` (coverage = hits on the assigned k-mers / *m\_g*). The benchmark scores each count at every `min_hits` for detection and pairs it with its coverage for abundance, so all four rules (gather + EM, winner-take-all, uniqueness-first, none) run on the same hit table as planned. Tests: hand-worked cases where the two rules disagree (a fully contained small unit wins under winner-take-all, loses under uniqueness-first unless it is sampled more sparsely), and on fixture reads every hit k-mer and hit is assigned exactly once. *Result pending* on HPC.

* **Phase 4, step 4 — simulated abundance benchmark (`workflows/sim-benchmark/sim.py`).** The fmh benchmark needs HPC, so a local simulation with exact truth guides model work between HPC runs (about 3 s per seed). 100 families from the MGnify sample, each with 3 paralogous units at 85–95% identity to the seed (shared k-mers), 4 members per unit at 97%; half the units present, each as one strain at 100/95/90/85% identity to its centroid with lognormal depth; random-codon CDSs, 150 bp single-end reads, 0.2% substitutions, 20% random decoy reads. Index configs: `dense` (*t\_base* 1), `s10` (0.1), `floor` (MGnify defaults). Metrics: unit purity and completeness; Spearman over true positives; L1 over units and over families (the Pfam-level view); sd of log(estimate / depth); `bias_<identity>` = median estimate / depth at that strain identity over the median at 100%.

  ```bash
  uv run python workflows/sim-benchmark/sim.py --seeds 5
  ```

* **Phase 4, step 5 — zero-inflated EM (`em(zero_inflated=True)`, `coverage_zi`, `present_zi`).** Each unit gets a present fraction π\_g beside λ\_g: only π\_g·*m\_g* of its kept k-mers occur in the sample, the rest are structural zeros. Hits are shared in proportion to λ·π; π\_g = (unit's share of hit k-mers) / (*m\_g*·(1 − e^−λ)), capped at 1. With no shared k-mers the fixed point is sylph's zero-truncated Poisson MLE; with no excess zeros it reduces to plain EM (both tested).

  Simulation, 5 seeds (min\_hits 1):

  | Config | Rule | Purity | Spearman | L1 | L1 family | log-ratio sd | bias 95% | bias 90% | bias 85% |
  | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
  | dense | none (all hits) | 0.58 | 0.75 | 0.55 | 0.41 | 0.74 | 0.58 | 0.34 | 0.22 |
  | dense | gather + EM | 0.91 | 0.79 | 0.50 | 0.41 | 0.78 | 0.56 | 0.30 | 0.16 |
  | dense | winner-take-all | 0.90 | 0.76 | 0.52 | 0.41 | 0.84 | 0.55 | 0.28 | 0.14 |
  | dense | uniqueness-first | 0.90 | 0.76 | 0.52 | 0.41 | 0.84 | 0.55 | 0.28 | 0.14 |
  | **dense** | **gather + ZI EM** | 0.91 | **0.95** | **0.23** | **0.18** | **0.52** | **0.94** | **0.94** | **0.91** |
  | s10 | gather + EM | 0.97 | 0.75 | 0.54 | 0.45 | 0.83 | 0.56 | 0.29 | 0.16 |
  | s10 | gather + ZI EM | 0.97 | 0.92 | 0.21 | 0.17 | 0.65 | 0.95 | 0.98 | 0.89 |
  | floor | gather + EM | 0.96 | 0.71 | 0.60 | 0.49 | 0.77 | 0.54 | 0.32 | 0.25 |
  | floor | gather + ZI EM | 0.96 | 0.84 | 0.39 | 0.32 | 0.72 | 0.94 | 0.92 | 0.84 |

  *Interpretation.* Divergence, not shared k-mers, dominates abundance error. Without zero-inflation, estimates fall roughly as identity^k (0.95^11 = 0.57, 0.90^11 = 0.31, 0.85^11 = 0.17), because a divergent strain carries only that share of its unit's k-mers; members' own variation dilutes further (median π ≈ 0.25 even at 100% identity). Zero-inflation removes the bias at every config. The three simple rules and plain EM differ little from each other (EM best by 0.03 Spearman); all gain purity over counting every hit (0.58 → 0.91 dense), and gather's detection is what does that. Residuals by depth (dense): median log-ratio within ±0.07 above depth 1 (−0.14 at ≤ 1), with sd 1.15 at depth ≤ 1, 0.45 at 1–3, 0.14 at 3–10, 0.09 above 10. With the floor index units have 2–4 hit k-mers and sd stays 0.25–1.0 at every depth. So the remaining error is thin data: low coverage, few k-mers.
  **Decided:** `coverage_zi` is the recommended abundance estimate; gather stays the detection rule. `coverage_em`, `coverage_wta` and `coverage_ufirst` remain as baselines for the fmh benchmark.

* **Phase 4, step 6 — empirical-Bayes prior on the present fraction (`fit_present_prior`, `coverage_zib`, experimental).** A Beta prior on π, matched by moments to units with λ ≥ 3 (where π is measured), with π updated by EM over unhit k-mers' presence. Result (5 seeds): dense log-ratio sd 0.52 → 0.42 and L1 0.23 → 0.21 at a small bias cost (85%: 0.91 → 0.88); s10 sd 0.65 → 0.46 but bias at 85% 0.89 → 0.76; floor worse everywhere (Spearman 0.84 → 0.74, bias at 85% 0.84 → 0.39), because with ~8 kept k-mers the prior outweighs the data and pulls divergent strains' π up, so λ down. **Decided:** not the default. It helps only where units have many k-mers, which is what the dense tier 2 will give, so revisit it there; `coverage_zib` stays in profiles so the HPC fmh run measures it on real data.

  *Next, by this evidence:* (1) the fmh benchmark on HPC, to confirm ZI on real genomes (KO units are far larger unions than MGnify90 clusters, so π is smaller); (2) dense tier 2 for candidate units, since k-mer count per unit limits precision more than the model does; (3) *p\_in*-weighted presence (core k-mers more likely present); (4) confidence intervals, which low-coverage units need.

* **Phase 4, step 7 — dense tier (`IndexParams.t_dense`, `--t-dense`).** The plan's two-stage lookup: the sparse tier 2 detects, a denser table quantifies.
  - *Build:* a fourth hashing pass keeps each unit's k-mers with hash ≤ max(*t\_dense*, *t\_g*) (`max_hash_dense`), minus promiscuous ones (same `max_groups` cut, counted over all units), with quantised *p\_in*, in a `dense` table of the tier-2 layout; `m_dense` per unit. The per-unit max makes each unit's dense set contain its tier-2 set. A uniform *t\_dense* was tried first: at 1 in 50 it gave floored units (tier 2 up to 0.2) fewer k-mers than tier 2, and ZI Spearman fell below the floor alone (0.83 vs 0.84), with completeness 0.78.
  - *Query:* after gather on tier-2 hits, the reads are streamed again at the dense rate and only the detected units' hits are kept; `coverage_em`, `_zi`, `_zib` are then fitted on those over `m_dense`, and `kmers_dense` reports the dense k-mers hit. Detection and the other baselines are unchanged. A second pass rather than buffering hashes, since the kernel runs at ~54 M reads/min and buffered hashes at dense rates would not fit in memory.
  - *EM fix found on the way:* with more shared k-mers visible, EM explains some gather-detected units away entirely; their weight decays towards 0, and zero-inflated EM divided 0 by 0 (NaN). Divisions are now guarded and units left with < 10⁻³ expected hits (`EXPLAINED_AWAY`) get coverage 0. Up to 1.5% of true units are explained away this way (at 1 in 50), so they drop out of true positives (floor 0.879 → 0.867–0.879 completeness).
  - *Test:* sparse tier 2 plus a fully dense tier gives exactly the EM and ZI estimates of an index that is dense throughout.

  Simulation, 5 seeds, zero-inflated EM (`gather_zi`); tables = size of the lookup tables relative to `floor`:

  | Config | Tables | Spearman | L1 | L1 family | log-ratio sd | bias 85% |
  | --- | --- | --- | --- | --- | --- | --- |
  | floor | 1× | 0.835 | 0.391 | 0.317 | 0.719 | 0.84 |
  | floor + dense 1/50 | 3.5× | 0.891 | 0.325 | 0.256 | 0.715 | 0.90 |
  | floor + dense 1/20 | 3.9× | 0.897 | 0.313 | 0.241 | 0.698 | 0.89 |
  | floor + dense 1/5 | 11× | 0.922 | 0.251 | 0.198 | 0.519 | 0.92 |
  | floor + dense 1/1 | 55× | 0.940 | 0.231 | 0.185 | 0.443 | 0.93 |
  | dense index (reference) | 54× | 0.948 | 0.228 | 0.178 | 0.517 | 0.91 |

  *Interpretation.* The dense tier recovers most of the gap between the floor index and a fully dense one: at 1 in 5, Spearman 0.84 → 0.92 and L1 0.39 → 0.25, with detection unchanged. With the dense tier, the empirical-Bayes prior no longer hurts (`gather_zib` at 1 in 5: Spearman 0.925, sd 0.42, but bias at 85% 0.81), as step 6 predicted; it still trades divergence bias for precision. The size ratios here overstate the cost at MGnify scale in one way and understate it in another: the simulated units are all non-singleton and floored (8 tier-2 k-mers each), while MGnify's 73% singletons sit at *t\_base* = 0.001 in tier 2 and would be raised to *t\_dense* in the dense table. The dense table is only probed for detected units and can stay on disk (memory-mapped), so its size bounds disk, not RAM, but the phase-4 gate (≤ 2× index size) is not met at any useful rate here.
  **Decided:** keep the dense tier as an option (default off) and measure it on real data before choosing a rate: `kfp_s100_d10` (tier 2 at 0.01, dense at 0.1) added to the fmh benchmark. Open: restrict the dense tier to non-singleton clusters, or to a rate that keeps it ≤ 1× tier 2 at MGnify scale; decide once the MGnify subset build reports `dense_bytes`.

* **Fix (fmh benchmark):** `index_*/meta.json` (index sizes, incl. `dense_bytes`) was never published: Nextflow folds a file inside a declared directory output into that directory, which did not match the publish pattern. INDEX and IMPORT\_SKETCHES now copy `meta.json` out and publish it to `<outdir>/index_<name>/meta.json`.

* **Phase 4, step 8 — *p\_in*-weighted presence (`em_pin`, `coverage_zip`, `present_zip`, experimental).** Zero-inflated EM with a per-k-mer presence probability from *p\_in(x)* (the plan's *z\_x* ~ Bernoulli(*π\_x*), *π\_x* informed by *p\_in*). The index stores per unit a histogram of kept k-mers over the 16 quantised *p\_in* levels (`pin_hist`, `pin_hist_dense`), so unhit k-mers are summed by level; the query now carries each hit's `pin_q` from tier 2 / the dense table. *p\_in* is the level midpoint, so no k-mer is certainly absent. Hit k-mers count as present by their share, as in `em`. Tests: all k-mers at one level reproduces zero-inflated EM; a k-mer shared between a unit where it is core and one where it is private goes mostly to the first.
  - *First link, logistic:* presence = σ(α\_g + logit *p\_in*), α\_g fitted by a Newton step per iteration. Worse than plain ZI everywhere (5 seeds): dense Spearman 0.948 → 0.936, bias at 85% 0.91 → 0.82, sd 0.52 → 0.55. Divergence scales presence multiplicatively (each k-mer survives with probability ≈ identity^k whatever its *p\_in*), and a logit shift cannot do that: bringing core k-mers down drives low-*p\_in* ones to 0.
  - *Adopted link, multiplicative:* presence = *s\_g*·*p\_in(x)* (*s\_g* = 1: the strain carries k-mers like a random member; < 1: divergent), capped so presence ≤ 1. Closed-form M-step: *s\_g* = expected present k-mers / Σ *p\_in*.

  Simulation, 5 seeds, `gather_zi` → `gather_zip`:

  | Config | Spearman | L1 | log-ratio sd | bias 85% |
  | --- | --- | --- | --- | --- |
  | dense | 0.948 → 0.947 | 0.228 → 0.229 | 0.517 → 0.411 | 0.91 → 0.88 |
  | s10 | 0.919 → 0.918 | 0.213 → 0.209 | 0.653 → 0.481 | 0.89 → 0.84 |
  | floor | 0.835 → 0.829 | 0.391 → 0.394 | 0.719 → 0.672 | 0.84 → 0.82 |
  | floor + dense 1/50 | 0.891 → 0.892 | 0.325 → 0.315 | 0.715 → 0.496 | 0.90 → 0.86 |
  | floor + dense 1/5 | 0.922 → 0.923 | 0.251 → 0.244 | 0.519 → 0.405 | 0.92 → 0.87 |
  | floor + dense 1/1 | 0.940 → 0.940 | 0.231 → 0.227 | 0.443 → 0.348 | 0.93 → 0.89 |

  *Interpretation.* *p\_in* tells the model which zeros are expected (private k-mers) and which are informative (core k-mers missed), so per-unit estimates scatter 20–30% less; ranks and L1 barely move, because they are dominated by the large depth range. The residual bias (≈ 4% at 85% identity) fits the simulation's cluster model, where presence given *p\_in* is not proportional to *p\_in*: members are star-like mutants of a centroid, so a k-mer at *p\_in* = 1/4 is almost never in the centroid, while the model expects it in a quarter of strains. Against the empirical-Bayes prior (step 6), it gives a similar precision gain at a quarter of the bias cost and does not break the floor index. **Decided:** not the default yet: `coverage_zip` joins the fmh benchmark, where KO units are unions over many species and *p\_in* is low for most k-mers, a very different regime. Adopt it if it holds there.

* **Fix (fmh benchmark):** INDEX and IMPORT\_SKETCHES now take the package sources as an input, like PROFILE, so `-resume` rebuilds indexes when the build code changes; cached indexes without `pin_hist` made every PROFILE fail. Indexes from runs before this change are rebuilt once.

* **Phase 4, step 9 — fmh benchmark with abundance (HPC, PR #13 code, 10 InSilicoSeq metagenomes).** Detection reproduced the phase-3 numbers exactly. Abundance against truth depth (Σ over the KO's genes of aligned bases / gene length), true positives, `kmers_unique` ≥ 1:

  | Index | Abundance | Spearman | L1 |
  | --- | --- | --- | --- |
  | fmh\_compat | coverage (hits / *m\_g*) | 0.08 | 1.41 |
  | fmh\_compat | ZI EM | 0.40 | 0.97 |
  | kfp\_s100 | coverage | 0.32 | 1.35 |
  | kfp\_s100 | EM / winner-take-all / uniqueness-first | 0.35 / 0.34 / 0.37 | 1.25 / 1.32 / 1.23 |
  | kfp\_s100 | ZI EM | 0.58 | 0.94 |
  | kfp\_s100 | ZI EM + prior | 0.53 | 0.98 |
  | kfp\_s100\_d10 | ZI EM (dense tier 1/10) | 0.63 | 0.93 |

  As in the simulation, zero-inflation is the largest gain and the three simple rules and plain EM are indistinguishable; the prior does not help. But the level is far below the simulation's 0.95. *Why:* a KO unit is the union of its genes from many species, and a sample holds several of them at different depths; truth sums depth over those copies, while ZI coverage is depth per present k-mer, about one copy's. Scaling by copies present (offline, from the profiles, with copies ≈ `present_zi` × `n_members`, valid when members share few k-mers as KO members do) gives Spearman **0.935** (kfp\_s100) and **0.964** (d10), log-ratio sd 0.63 and 0.45. Hits / *t\_g* / (k-mers per member), with no model at all, gives 0.92: most of the gap was the estimand, not the model.
  **Decided:** profiles report `copies_zi` = present k-mers / `pin_sum` and `abundance_zi` = `coverage_zi` × `copies_zi`. `pin_sum` (`pin_sum_dense`) is Σ *p\_in* over a unit's kept k-mers, computed at build from unquantised *p\_in*: how many kept k-mers an average member holds, so present / `pin_sum` counts member-equivalents present. The 4-bit *p\_in* levels cannot give this for KOs (most KO k-mers have *p\_in* < 1/15). The two estimates answer different questions: `coverage_zi` is depth per copy and stays the per-cluster estimate for MGnify90 units (one strain each; the simulation shows `abundance_zi` there brings back the divergence bias, 0.16 at 85%, because a divergent strain looks like fewer copies); `abundance_zi` is total depth over copies, which is what the KO-unit benchmark measures. For the real tool, function-level abundance is Σ over clusters of `coverage_zi`, per the plan (EM at cluster level, then aggregate). The fmh benchmark now scores `abundance_zi` to confirm the exact form (needs a rerun; these indexes predate `pin_sum`).

  *Dense tier cost on real data* (kfp\_s100\_d10 vs kfp\_s100): tables 1.09 GB vs 72 MB (15×; 7.9 B vs 5.2 B per posting), index build 36.5 GB vs 7.8 GB peak memory, query 2m16–2m46 vs 23 s and 5 GB vs 0.58 GB (the prototype materialises 8 B per key for lookup). For +0.05 Spearman (ZI) or +0.03 (copies-scaled), and a 30% lower log-ratio sd. **Decided:** the dense tier stays off by default; revisit after the Rust lookup (memory-mapped, no materialised keys) and with a rate or singleton restriction chosen from MGnify-subset `dense_bytes`.

  *Also:* `em_pin` read quantised *p\_in* levels as (*l* + 0.5)/16, but the build quantises as round(15·*p\_in*); now *l*/15 with level 0 at 1/30. Simulation after the fix (`gather_zip` vs `gather_zi`): dense sd 0.52 → 0.42 as before, floor Spearman 0.835 → 0.806 (was 0.829). The `gut-lin10000` archive from the same HPC batch is the pre-fix phase-2 build (*t\_max* 0.57, 19.4 B per hash), already recorded under Phase 2.

* **Phase 4, step 10 — confidence intervals (`bootstrap_zi`, `query --bootstrap B`).** 95% percentile intervals for `coverage_zi` and `abundance_zi` (`*_lo`, `*_hi`) from a Poisson bootstrap over reads: each replicate weights every read pair by a Poisson(1) draw, rebuilds the k-mer counts and refits zero-inflated EM (sharing between units included). Reads, not k-mers, are resampled because one 150 bp read hits ~40 consecutive amino-acid k-mers of a unit, so k-mer counts are correlated and a Poisson/Fisher interval would be too narrow. To do this the query now keeps hits per (unit, k-mer, read) for detected units, in both passes. Off by default (B = 0); tests: intervals bracket the estimate, leave every other column unchanged, and are reproducible (fixed seed).
  - *Simulation calibration* (5 seeds, B = 100, `coverage_zi` against true k-mer coverage = depth × a per-sample scale, the median estimate / depth over strains at 100% identity and depth > 5):

  | Config | Coverage (95% nominal) | Coverage, depth ≤ 2 | Median width, log(hi / lo) |
  | --- | --- | --- | --- |
  | dense | 0.964 | 0.962 | 1.17 (×3.2) |
  | s10 | 0.973 | 0.979 | 1.42 (×4.1) |
  | floor | 0.968 | 0.983 | 2.18 (×8.8) |
  | floor + dense 1/50 | 0.958 | 0.966 | 1.63 (×5.1) |
  | floor + dense 1/5 | 0.950 | 0.925 | 1.15 (×3.2) |

  *Interpretation.* Calibrated to slightly conservative everywhere, including at low depth where the point estimate is least reliable; widths track the information per unit (the floor index's ~8 k-mers per unit give ×9 intervals). This meets the phase-4 gate's "calibrated intervals on simulations" for `coverage_zi`. Caveat: the scale is estimated from the same sample; the ±1–2% conservatism is within what that allows.
  - *fmh benchmark:* PROFILE runs with `--bootstrap` (`params.bootstrap`, default 100) and SCORE reports `ci_cover` and `ci_width` for estimates with intervals, against truth depth on the estimate's scale (median estimate / depth over true positives). Result pending on HPC; the dense-tier profile will be the slowest (100 EM refits over its larger hit table).

* **Phase 4, step 11 — fmh benchmark with copies-scaled abundance (HPC, PR #14 code).** Detection and every earlier estimate reproduced the previous run exactly. `abundance_zi` against truth depth (`kmers_unique` ≥ 1, true positives, 10 seeds):

  | Index | `coverage_zi` Spearman / L1 | `abundance_zi` Spearman (sd) / L1 |
  | --- | --- | --- |
  | kfp\_s1000 | 0.373 / 0.97 | 0.807 (0.011) / 0.57 |
  | kfp\_s1000\_floor8 | 0.369 / 1.02 | 0.827 (0.014) / 0.56 |
  | kfp\_s100 | 0.576 / 0.94 | **0.957 (0.003) / 0.22** |
  | kfp\_s100\_d10 | 0.627 / 0.93 | **0.988 (0.002) / 0.14** |

  *Interpretation.* The exact form (present k-mers / Σ *p\_in*) beats the offline approximation (0.935), and on KO units abundance is now as good as detection: kfp\_s100 detects 95.5% of KOs at 98.5% purity and ranks their abundance at 0.96. The dense tier's gain is larger here than for per-copy coverage (L1 0.22 → 0.14, 36% lower), but its cost is unchanged (15× tables, ~6× query time); still off by default, and a candidate for the phase-6 cost study. At scaled 1000 the floor helps abundance a little (0.81 → 0.83).
  - `fmh_compat`: `abundance_zi` equals plain EM (0.09), as expected: imported sketches carry no members, so *p\_in* = 1 and `pin_sum` = *m\_g*, and copies reduce to the present fraction. fmh-funprofiler's own sketches cannot give copy-aware abundance.
  - `coverage_zip` equals `coverage_zi` on our KO indexes (0.576 vs 0.576): nearly all KO k-mers have *p\_in* < 1/15 and fall in quantised level 0, so presence cannot vary between them. On `fmh_compat` it equalled plain EM, which exposed a bug: with level 15 read as exactly 1, a k-mer present with probability 1 stays so, because EM then reads every unhit copy as present. **Fixed:** level probabilities are clipped to [1/30, 1 − 1/30]; a test covers all-core units. The simulation is unchanged by the fix (its `coverage_zip` bias is structural, not this bug). `coverage_zip` stays experimental: it cannot help KO-like units at 4-bit *p\_in*, and its intended case, MGnify90 clusters with varied *p\_in*, has no truth-bearing benchmark until phase 6.
  **Decided:** recommended setting for KO-like units is `kfp_s100`, detection `kmers_unique` ≥ 1, abundance `abundance_zi`; for MGnify90 clusters `coverage_zi`, summed over clusters for function-level abundance. Intervals (PR #15) have not been run on HPC yet.

* **Phase 4, step 12 — bootstrap intervals on real data (HPC, PR #15 code, B = 100).** Detection and all point estimates reproduced the previous run; the `em_pin` fix shows as expected (`fmh_compat` `coverage_zip` 0.09 → 0.395, in line with `coverage_zi`). Interval coverage of truth depth on the estimate's scale (nominal 95%, true positives):

  | Index | `abundance_zi` coverage | median width, log(hi / lo) | `coverage_zi` coverage |
  | --- | --- | --- | --- |
  | kfp\_s1000 | 0.556 | 1.28 | 0.529 |
  | kfp\_s100 | 0.607 | 0.62 | 0.333 |
  | kfp\_s100\_d10 | 0.717 | 0.38 | 0.163 |

  (`coverage_zi` estimates depth per copy, not the KO total, so its low coverage is expected.) *Why `abundance_zi` intervals are too narrow:* the estimate is unbiased at every depth (median log residual within ±0.03 per depth quartile on kfp\_s100), but the bootstrap half-width shrinks with depth (0.88 → 0.12 from the lowest to the highest quartile) while the residual sd does not (0.69 → 0.34). Coverage therefore falls with depth, 0.80 → 0.43 on kfp\_s100 and 0.88 → 0.49 with the dense tier. The bootstrap measures read-sampling noise only; what it misses is a per-unit model error. Adding a constant σ on the log scale in quadrature restores calibration offline: σ = 0.2 gives 0.95 (d10) and 0.88 (s100), σ = 0.3 gives 0.97 and 0.93 (s1000 needs more: 0.82 at 0.3). The simulation had no such error (single-strain units, so copies are exact), which is why it was calibrated there. The likeliest source is the copies estimate: it divides present k-mers by an *average* member's kept k-mers, but the KO genes in a sample are specific members whose lengths vary between species.
  *Cost:* the bootstrap multiplies PROFILE time by 4–6 on kfp\_s100 (23 s → 1m31–2m09) and 8–11 with the dense tier (2m16–2m46 → 18–30 min; peak memory 6 GB).
  **Decided:** intervals are not calibrated on real data yet; `--bootstrap` stays off by default in the tool, and results must not be read as 95% intervals. Next: add a per-unit model-error term, from the spread of member lengths (kept k-mers per member, measurable at build) divided by the square root of the copies present, plus a small calibrated floor, and test it on the next HPC run; reduce the benchmark's default replicates (or bootstrap only non-dense indexes) to cut cost.

* **Phase 4, step 13 — posterior intervals and ambiguity groups (`posterior_zi`, `query --draws D`; replaces the read bootstrap).** *Intent (clarified):* intervals should carry the uncertainty of k-mer matches, i.e. of splitting hits between similar units, including off-target matches from sequencing errors, which mostly land on close relatives. An LCA over a similarity tree was considered and not needed: the likelihood already shows which units the data cannot separate. Aggregating to functions is out of scope for the tool (downstream analysis); only the benchmarks do it.
  - *Method:* each of D draws reweights reads (Poisson(1) per read pair), refits zero-inflated EM (read-sampling uncertainty), then runs 10 Gibbs sweeps on those counts and keeps the last state (split uncertainty). A sweep splits each hit k-mer's count among its holders multinomially ∝ coverage × present fraction; updates coverage from the counts on the k-mers each unit was given, a zero-truncated Poisson likelihood free of the present fraction (Metropolis on log coverage, 5 steps, Gamma(1, 0.01) prior); and draws the present fraction exactly on a 512-point grid from hit k-mers ~ Binomial(*m*, π(1 − e^−λ)). A unit given no hits in a draw is absent in it. *Ambiguity groups:* holders of a shared hit k-mer whose coverage draws correlate below −0.5 are linked; connected units form a group with `ambiguity_group`, `group_size` and interval on the group total (`group_coverage_zi_*`, `group_abundance_zi_*`).
  - *Tried first, and why they failed:* (1) Gibbs with data augmentation of presence and conjugate Gamma/Beta draws, reads reweighted every sweep: intervals too narrow (dense coverage 0.85), because coverage and present fraction trade off along a ridge the chain crawls along. (2) The same with the zero-truncated coverage update: still too narrow (0.81), because a single chain only partly follows each sweep's reweighted counts and so averages over them. Nesting (refit per draw, short chain from it) fixed both.
  - *Simulation with near-identical twins* (5 seeds, D = 100, `--twins 0.3`: 30% of families gain a unit that is a 99% copy of another, each present independently):

  | Config | Coverage | depth ≤ 2 | twins | width | TP grouped | group cover | twin pairs grouped |
  | --- | --- | --- | --- | --- | --- | --- | --- |
  | dense | 0.933 | 0.964 | 0.682 | 1.21 | 0.015 | 0.53 | 0.10 |
  | s10 | 0.959 | 0.982 | 0.682 | 1.71 | 0.014 | 0.88 | 0.22 |
  | floor | 0.969 | 0.966 | 0.883 | 2.66 | 0.000 | – | 0.00 |
  | floor + dense 1/5 | 0.946 | 0.954 | 0.863 | 1.35 | 0.037 | 0.80 | 0.21 |

  *Interpretation.* Intervals are calibrated overall and at low depth. Most twins are resolved by gather (the unsupported twin gets no gather k-mers and is not reported), which is correct. The hard case is both twins present at different depths: the split is biased, the true twin's interval covers 68–88%, and only 10–22% of such pairs form a group, because their draws are not anti-correlated enough. The likely cause is the split step's approximation (weights ∝ coverage × present fraction, ignoring the holders' joint presence on each shared k-mer).
  - *Removed:* `bootstrap_zi` and `--bootstrap` (`--draws` replaces them; the fmh benchmark's `params.bootstrap` is now `params.draws`). *fmh benchmark:* SCORE adds `fp_grouped` (share of false-positive KOs grouped with a true KO) and `group_cover` for `abundance_zi`; result pending on HPC. The copies model error found in step 12 is separate and not addressed here.
  - *Next:* model the holders' joint presence exactly in the split (most shared k-mers have two holders, so enumeration is cheap), and check whether twin coverage and grouping improve; tune the anti-correlation threshold on the twin simulation.

* **Fix (fmh benchmark):** SUMMARY failed when a metric was empty in one sample's score file (e.g. `group_cover` when no ambiguity group formed): Polars reads an all-empty column as String and cannot stack it with Float64. Metrics are now cast to Float64 on read; a regression test stacks a file with an empty metric.

* **Phase 4, step 14 — posterior intervals and ambiguity groups on real data (HPC, PR #17 + #18 code, D = 100).** Point estimates unchanged. `abundance_zi` interval coverage of truth depth (nominal 95%; the step-12 bootstrap in brackets): kfp\_s1000 0.69 (0.56), kfp\_s100 0.72 (0.61), kfp\_s100\_d10 0.75 (0.72), fmh\_compat 0.43 (0.36); widths grew by 10–30%. Better, still not calibrated: the per-KO model error of the copies estimate (step 12) is untouched by sampling. PROFILE time: kfp\_s100 2m17–2m39 (23 s without draws), dense tier 20–25 min, kfp\_s1000 39–44 s.
  - *Ambiguity groups did not form:* none on kfp\_s100 or kfp\_s100\_d10 in any of 10 samples; groups in 2–4 samples on the scaled-1000 and compat indexes; no false-positive KO was ever grouped with a true one (`fp_grouped` empty). *Why:* the false positives left after gather are small units with a median of 1 gather-assigned k-mer (true KOs: 15), borrowing hits from a well-supported true KO. The true KO's coverage barely moves when the small unit takes or loses those hits, so their draws are not anti-correlated: correlation is the wrong signal for lopsided pairs.
  - *Absent in some draws* (lower bound 0) as a false-positive flag: catches 52% of false positives on kfp\_s100 (21% with the dense tier) but also 8% of true KOs, so only 9% of flagged KOs are false. Weak.
  - *Next (to decide):* group by shared evidence instead of correlation: link a unit to another when most of its allocated hits (across draws) lie on k-mers it shares with that unit, which catches lopsided pairs; and report each unit's share of evidence that is its own. Separately, the copies model error still keeps `abundance_zi` intervals from calibrating.

* **Phase 4, step 15 — ambiguity groups by shared evidence (replaces anti-correlation).** Over all posterior draws, each unit's allocated hits are accumulated per k-mer; a unit is linked to another holder when at least half (`shared_evidence` = 0.5) of its allocated hits lie on k-mers the two share, and linked units form a group. A unit never given hits (explained away in every draw) links to every co-holder of its hit k-mers. Profiles gain `own_evidence`: a unit's share of allocated hits on k-mers no other detected unit holds. Test: a small unit with one stray hit of its own and three k-mers shared with a well-supported unit (the lopsided pair correlation missed) is grouped with it, and so is one with no own k-mers; an independent unit is not.
  - *Simulation* (5 seeds, D = 100, `--twins 0.3`), anti-correlation (step 13) → shared evidence: twin pairs grouped dense 0.10 → **1.0**, s10 0.22 → **1.0**, floor + dense 1/5 0.21 → **1.0**, floor 0.00 → 0.17 (≈8 k-mers per unit leave little shared evidence to count); group intervals holding the true total 0.53 → 0.83 (dense), 0.88 → 0.78 (s10), 0.80 → 0.90 (floor + dense 1/5), – → 0.85 (floor); true positives in a group 3–23%. Per-unit interval coverage unchanged (0.93–0.97 overall; twins both present still 0.70–0.88).
  - *Interpretation:* grouping now finds the pairs the data cannot separate, and group totals are closer to calibrated than members, as intended; group intervals still undercover a little (0.78–0.90), for the same split-approximation reason as twins' member intervals. The 0.5 threshold is a default, not tuned; the fmh run (`fp_grouped`) will show how many real false positives it catches. Result pending on HPC.

* **Phase 4, step 16 — shared-evidence groups on real data (HPC, PR #19 code, D = 100).** Point estimates and intervals unchanged from step 14. Per index, over 10 samples (true positives TP, false positives FP, `kmers_unique` ≥ 1):

  | Index | FPs grouped | TPs grouped | groups / sample | groups of true KOs only | FP median `own_evidence` | group cover |
  | --- | --- | --- | --- | --- | --- | --- |
  | kfp\_s1000 | 1.0% | 2.0% | 27 | 99% | 1.0 | 0.58 |
  | kfp\_s100 | 5.6% | 2.3% | 53 | 91% | 1.0 | 0.57 |
  | kfp\_s100\_d10 | 5.0% | 1.3% | 34 | 87% | 0.98 | 0.57 |
  | fmh\_compat | 1.9% | 2.1% | 38 | 98% | 1.0 | 0.36 |

  Of the FPs that are grouped, 84–100% share the group with a true KO (the benchmark's `fp_grouped` reported this conditional share by mistake; see fix below). Groups have a median of 2 KOs (up to 60 on kfp\_s100).
  - *Interpretation.* The FPs left after gather are mostly not sequence-similarity confusions: their hits sit on k-mers no other detected KO holds (median `own_evidence` 1.0). Gather already removed the similarity-driven FPs in phase 3 (348 → 105 per sample), so what remains looks like independent evidence: chance or error matches, or genes the truth rule does not count (collagen VII was the most persistent in step 3). Groups mostly flag pairs of true KOs present together (paralogous KOs), which is the intended output; their totals are not calibrated (0.57 cover), for the same copies-model-error reason as the members' intervals. The per-unit calibration gap on real data is therefore the copies estimate, not ambiguity.
  - *Fix:* `fp_grouped` averaged `is_in` over FPs, but FPs without a group gave null and were skipped, so it reported the share of *grouped* FPs that sit with a true KO (0.84–1.0) instead of the share of all FPs grouped with one (≈ 1–5%). Nulls now count as not grouped; regression test added.
  - *Next options:* model the copies error (spread of member lengths per unit, at build) so KO-level intervals calibrate; examine the ungrouped FPs directly (which KOs, which k-mers, error-rate dependence) to decide whether a background false-hit term is needed; exact joint presence in the posterior split (step 13) only matters for near-identical units, which are rare in this benchmark.

* **Phase 4, step 17 — presence probability, copies error, false-positive diagnostic.** Three changes, from step 16's next options plus a new requirement: uncertainty must also reflect weak evidence, i.e. a unit hit on few k-mers, especially k-mers many units share.
  - *Presence probability (`present_prob`, always computed).* Per unit gather keeps: present vs. detected by background. Evidence is the unit's own k-mers (the tier-2 hit k-mers gather gave it), so hits other units explain do not count. An absent unit is detected by background (off-target homologs, error k-mers) with probability 1 − exp(−μ), μ = `BACKGROUND` × reads × *t\_g*, then has *h* own k-mers ~ Geometric(`CLUMP`). Each k-mer's odds are multiplied by its index *holders*, since conserved motifs shared by many units are what unindexed genes carry too. A present unit's *h* follows a distribution fitted to the sample. With A absent index units, P(present | *h*) = *w\_h* / (*w\_h* + A · P(*h* | absent)), fitted by fixed-point iteration; at the fixed point, *w\_h* ≈ units seen with *h* minus expected background. Posterior draws drop a unit with probability 1 − `present_prob`; group totals keep a dropped member's share. Defaults `BACKGROUND` = 1e-6 per read pair per unit of *t\_g* and `CLUMP` = 0.3 are calibrated on the step-16 fmh run (below).
  - *Tried and dropped (on the step-16 fmh profiles, offline):*
    1. Background fitted per sample (β per posting, presence prior ρ, per-unit hit rate *h*/*m*): **not identifiable**. A free per-unit hit rate explains any few-hit unit as well as background does, so fits went to one extreme or the other. In 6 of 10 s100 samples β ≈ 7e-8 and nothing was flagged; in 4, β ≈ 0.017 and 79–90% of all KOs were flagged. The simulation had not shown this, because its false positives are explained by other detected holders.
    2. A beta-binomial population prior on the hit rate, or a free present-unit *h* distribution with ρ free: ρ → 1 and every unit ≈ 1.
    3. Background ∝ *m\_g* (β per posting): calibratable, but within one-k-mer units it ranks the wrong way (AUC 0.47 at s100, 0.19 at s1000), because false positives are small KOs. ∝ *t\_g* ranks neutrally at fixed *t\_g* and well where *t\_g* varies (floor index: AUC 0.87 within one-k-mer units, where false positives sit on floored KOs).
    4. A Poisson count given background: overconfident at 2 own k-mers (1.00 vs 0.96 true at s100), because background comes in clumps (37% of s100 false positives have ≥ 2 own k-mers). Replaced by the geometric.
    5. (Simulation) a floor of one copy's k-mers at the fitted coverage: flagged 54% of true units, because divergent strains keep 0.85¹¹ = 17% of k-mers at 85% identity.
  - *Calibration (step-16 fmh profiles, offline, 3.3 M read pairs, holders = 1).* P(true | own k-mers = 1, 2, 3), predicted vs. observed: kfp\_s100 0.80/0.85, 0.94/0.96, 0.98/0.98; kfp\_s1000 0.99/0.99, 0.99/1.00, 1.00/1.00; kfp\_s1000\_floor8 0.98/0.95, 0.99/0.99, 0.99/0.99; fmh\_compat 0.98/0.97, 0.99/1.00, 1.00/1.00. Among single-k-mer KOs, no per-unit feature in the profile separates true from false (AUC 0.42–0.53 for hits, reads per k-mer, *m\_g*, hit fraction, coverage, gather rank), so the calibrated probability of the count is what the data support; *t\_g* adds separation only where it varies. Holders cannot be checked on KO indexes (*u\_g* / *m\_g* ≈ 0.99: KO k-mers are almost never shared); it will matter for MGnify90 clusters. False positives are spread over many KOs (694 KOs for 1,048 s100 false positives, 29% of rows from KOs false in ≥ 3 samples), except with the floor, where some are false in all 10 samples (e.g. K03044, an RNA polymerase subunit homologous to bacterial rpoB).
  - *Copies error.* The index stores `len_cv` (`len_cv_dense`): the coefficient of variation of kept k-mers per counting member; its mean is `pin_sum` (tested). Each posterior draw scales copies by exp(N(0, *s*²)), *s*² = log(1 + `len_cv`² / copies) + `copies_error`². *Offline check on the step-16 run:* the excess error in `abundance_zi`, √(residual sd² − posterior sd²), falls with copies (dense tier 0.30 → 0.19, s100 0.54 → 0.33 from the lowest to the fourth copies quintile) but levels off. Fit σ² ≈ a / copies + floor²: a ≈ 0.07 and floor ≈ 0.17 with the dense tier; a ≈ 0.25 and floor ≈ 0.30 at s100. So `len_cv` should cover the low-copies end, but a floor of ≈ 0.2 will be needed. **Decided:** `copies_error` stays 0 until the HPC run measures `len_cv` alone, then the floor is fitted.
  - *Simulation* (5 seeds, D = 100, `--twins 0.3`; step 13/15 values in brackets). With ~10⁴ reads, background is negligible and `present_prob` ≈ 1, so intervals are as before: coverage dense 0.936 (0.933), s10 0.964 (0.959), floor 0.970 (0.969), floor + dense 1/5 0.953 (0.946); twins 0.71, 0.76, 0.88, 0.86; group cover 0.81, 0.84, 0.75, 0.91. Its false positives are paralogs keeping own k-mers after gather; the ambiguity groups handle those, not presence.
  - *False-positive diagnostic (fmh benchmark).* `query --kmers` writes the tier-2 hits per (unit, k-mer) with `hits` and `holders`. A new DETECTED step writes one row per KO gather keeps (`detected.tsv`) with `tp`, counts, `present_prob`, `own_evidence`, group and, over its hit k-mers, median holders, most hits on one k-mer and `in_genome`: the share found in the six-frame translation of the sample genomes. A false positive's k-mers in the genomes are real sequence (another gene or KO); those not in them come from read errors. `--iss_mode perfect` reruns with error-free reads. SCORE adds `prob_tp`, `prob_fp`, `flag_tp`, `flag_fp` for `abundance_zi`. Result pending on HPC (default run, then `--iss_mode perfect`).
  - *Next:* from the HPC run, (1) whether `present_prob` stays calibrated in the tool (own k-mers from gather, reads counted, holders from lookups) and how its interval dropping changes `abundance_zi` coverage; (2) fit `copies_error` with `len_cv` in place; (3) what false-positive k-mers are (in genome or not, error dependence), which decides whether `BACKGROUND` should depend on read errors.

* **Phase 4, step 18 — fmh benchmark with presence, `len_cv` and error-free reads (HPC, PR #22 code, D = 100).** Reads simulated with `--iss_mode perfect`: inferred from fresh SIMULATE tasks and every hit k-mer being in the sample genomes; the run's parameters were not in the results. Detection matches step 16 (read errors) almost exactly: kfp\_s100 completeness 0.956 (0.955), purity 0.985 (0.985), 103 false positives per sample (105).
  - *False positives come from real sequence, not read errors.* Every hit k-mer of every false-positive KO is in the six-frame translation of the sample genomes (`in_genome` = 1.0 on all kfp indexes), and removing read errors removes almost no false positives. They have a median of 1 holder (the k-mers are unique among KOs), so they are sequence in the sample genomes that the truth rule does not count for that KO: another gene, an unannotated gene, or another frame. `BACKGROUND` therefore need not depend on the error rate.
  - *`present_prob` in the tool* reproduces the offline calibration (step 17): P(true | own k-mers = 1, 2, 3), predicted vs. observed, kfp\_s100 0.80/0.85, 0.94/0.96, 0.98/0.98; kfp\_s1000 0.99/0.99, 0.99/1.00, 1.00/1.00; floor8 0.98/0.95, 0.99/0.99, 0.99/0.99; fmh\_compat 0.98/0.97, 0.99/1.00, 1.00/1.00 (corrected: first read from `detected.tsv`, whose `fmh_compat` rows were shifted one column, see step 19). Mean over false positives 0.85 vs. 0.98 over true ones at s100. Within single-k-mer KOs it separates only on the floored index (AUC 0.87; 0.50–0.52 elsewhere), as expected.
  - *`abundance_zi` intervals with `len_cv` (`copies_error` = 0):* coverage of truth depth (step 16 in brackets) kfp\_s100 0.896 (0.72), kfp\_s100\_d10 0.857 (0.75), kfp\_s1000 0.900 (0.69), floor8 0.903 (0.70), fmh\_compat 0.45 (0.43, no members so no `len_cv`); median width log(hi / lo) 1.29 (0.79), 0.58 (0.42), 2.77 (1.64). Group cover 0.80 (0.57), 0.69 (0.57), 0.87 (0.58), 0.82 (0.58). The gain is `len_cv`: coverage among intervals with a positive lower bound is about the same as overall (0.88 vs 0.89 at s100), so presence-dropping only widens already-wide intervals. The lowest two copies quintiles now cover 0.92–0.95; the highest still undercover, with excess error ≈ 0.15 (dense tier) to 0.3–0.45 (s100, s1000), the floor step 17 predicted.
  - **Decided:** `copies_error` = 0.2. Widening the intervals offline by 0.2 in quadrature gives kfp\_s100 0.95, kfp\_s100\_d10 0.97, kfp\_s1000 0.93, floor8 0.92 (0.15: 0.94, 0.96, 0.92, 0.91). Fitted on one intact profile per index (see fix), so the next run must confirm it over 10 samples.
  - *Fix:* PROFILE published `profile.tsv` and `kmers.parquet` under the same name, so 42 of 50 published profiles were k-mer tables. SCORE, DETECTED and the summaries read the task outputs and are unaffected. The tables now go to `profiles/` and `kmers/`.
  - *Next:* rerun the default (read-error) benchmark to confirm `copies_error` over all samples. To find what the false-positive sequence is, locate each false positive's k-mers in the sample genomes against the gene coordinates: inside a gene of another KO, a gene with no KO, or outside genes. That decides whether these are errors of the truth rule or genuine off-target homology.

* **Step 19 — reproducible posterior draws; `detected.tsv` column fix.**
  - `posterior_zi` is seeded, but its draws depended on input row order, which Polars `group_by` does not fix: two identical fmh test runs differed in `own_evidence`. It now sorts its hits by unit, hash and read first; a test checks that shuffled input gives identical output, and two fresh test-profile runs give identical `detected.tsv` and `scores.tsv`.
  - `DETECTED` wrote only the columns a profile has, so `fmh_compat` rows (no `n_members`) had 15 fields under a 16-field header in the stacked `detected.tsv`. Their columns were shifted one place from `n_members` on, or the file failed to parse, depending on which file came first. Every column is now written, null where missing. The step-18 `fmh_compat` presence calibration was recomputed from the shifted rows and corrected; SCORE's `prob_*` metrics read the profiles directly and were right.

* **Phase 4, step 20 — fmh benchmark with `copies_error` = 0.2 and read errors (HPC, PR #24 code, iss novaseq, D = 100; described by the run's `run.json`).** Detection and point estimates as in step 16.
  - *`abundance_zi` intervals are calibrated on real data.* Coverage of truth depth over 10 samples (step 18, `len_cv` only, in brackets): kfp\_s100 **0.950** (sd 0.005; 0.90), kfp\_s100\_d10 0.962 (0.86), kfp\_s1000 0.926 (0.90), floor8 0.926 (0.90), as the offline fit predicted (0.95, 0.97, 0.93, 0.92). By copies quintile, kfp\_s100 covers 0.945–0.956 and kfp\_s100\_d10 0.956–0.971 in every quintile; kfp\_s1000 and floor8 fall to 0.90 in the top two quintiles, so the sparse indexes could take a slightly larger floor. Median width log(hi / lo) 1.50 at s100 (1.29) and 0.95 with the dense tier (0.58). Group cover 0.93, 0.93, 0.88, 0.87 (0.80, 0.69, 0.87, 0.82). fmh\_compat stays at 0.50: its sketches carry no members, so no `len_cv`, and its copies are not meaningful.
  - *Read errors cause few false positives.* 96% of s100 false positives have every hit k-mer in the sample genomes (mean share 0.977); about 4% involve a k-mer from a read error, consistent with 105 false positives per sample here vs. 103 with error-free reads (step 18). 14% of true KOs have some error-derived hit, which does no harm.
  - *`present_prob`* calibration is the same as with error-free reads (kfp\_s100 0.80/0.85, 0.94/0.96, 0.98/0.98 at 1/2/3 own k-mers).
  - The run predates #25 and #26: 8 floor8 DETECTED tasks were again killed at 8 GB before passing on retry, and `detected.tsv` has the column shift (fmh\_compat file first here, so the other indexes' rows carry an extra field). Both are fixed on `main`.
  - **Decided:** `copies_error` = 0.2 stays. With the phase-4 gate's "calibrated intervals on simulations" met and intervals now calibrated on real data at the recommended setting (kfp\_s100), the remaining phase-4 item is where false-positive sequence comes from (step 18, *Next*).

### Phase 5

* **Phase 5, step 1 — other tools in the fmh benchmark (`--tools`).** The fmh benchmark now runs DIAMOND, fmh-funprofiler, kMermaid and HUMAnN on the same 10 metagenomes and scores them with the same truth and metrics as our indexes, so every tool sits in one `summary.tsv`. Each tool's output becomes a profile of `name`, `evidence` (its detection count) and `abundance` (`bench.py tool-profile`); SCORE treats (`evidence`, `abundance`) as one more rule. Choices:
  - *DIAMOND 2.2.8* blastx, best hit per read (`-k 1`), default sensitivity, against the same KEGG proteins our KO indexes are built from: the alignment baseline without a reference handicap. Evidence: read pairs; abundance: Σ aligned / subject length, the quantity truth depth sums.
  - *fmh-funprofiler 1.1.1*, the released `funcprofiler` on its KO sketches (scaled 1000), for its cost; `fmh_compat` already reproduces its calls. `--threshold_bp` = scaled (one shared hash): its default at scaled 1000, and on the test profile (scaled 10) it gives the same 12 KOs as `fmh_compat`, where the default 1000 bp (100 hashes) gave 2.
  - *kMermaid* (commit `edcb4ed`) ships a RefSeq-cluster model with no KO labels, so it is retrained with KOs as clusters through its own `get_kmer_dict_cluster`, as its README allows (`kmermaid_kfp.py`). Clusters are capped at 50 random proteins per KO (`--kmermaid_max_members`), near the shipped model's ~55 per cluster; uncapped KEGG KOs would not fit in memory as its Python dicts. Evidence: read pairs; abundance: reads / mean member length.
  - *HUMAnN 3.9* (MetaPhlAn 4.1.1, `mpa_vJun23_CHOCOPhlAnSGB_202403`) on its own ChocoPhlAn and UniRef90 databases, UniRef90 families regrouped to KOs with its mapping; unstratified RPK. It reports no read counts, so every KO it reports has evidence 1. The sample genomes are in UniRef90, but its KO calls go through UniRef90 → KO rather than KEGG's gene → KO table, so part of any gap will be annotation, not detection.
  - *Cost.* The trace is now raw (ms, bytes) with `%cpu` and `peak_rss`; `bench.py cost` summarises wall hours, CPU hours and peak RSS per step and index or tool. Threads differ (DIAMOND 8, HUMAnN 16, the rest 1), so CPU hours are the comparison. One-off database builds (DIAMOND_DB, KMERMAID_MODEL, HUMANN_DB) are separate rows, like our INDEX.
  - *Containers.* Tool steps run in pinned biocontainers (the kMermaid image is built by `containers/build.sh`, since kMermaid has no package); `-profile docker` locally, `-profile singularity` on HPC. Our own steps still run on the host in `.venv`. HPC kit: `hpc.example.config`, `run_hpc.sh` (Slurm head job), `containers/pull.sh` (pre-pulls images into the Singularity cache).
  - *Checked* on the test profile with Docker (DIAMOND, fmh-funprofiler, kMermaid; HUMAnN only on HPC, its databases are ~45 GB): all steps run and score, and native fmh-funprofiler matches `fmh_compat`. `tests/python/test_fmh_benchmark.py` covers each tool's parsing and `cost`. Also fixed: under `-profile test` the trace was written to `results/`, since `trace.file` resolved before the profile set `outdir`.
  - *Result pending on HPC.* Expected: DIAMOND near-complete detection at a much higher CPU cost (the fmh-funprofiler paper reports 39–99x); kMermaid's k = 5 whole-read scores trade purity for speed; HUMAnN lower completeness on KEGG truth from the UniRef90 → KO mapping. Our query is still the Python prototype, so its cost is an upper bound until the phase-6 port.
  - *Next:* the HPC run; then the divergence ladder and ablations (phase 5's other items).

* **Phase 5, step 2 — HUMAnN 4 (`humann4`); HUMAnN images fixed.** HUMAnN 4 is the current development release, 4.0.0.alpha.2 (GitHub `master`, commit `e07b3a3`, July 2026). It is not on bioconda or PyPI, so `containers/humann4` installs it from GitHub. It runs as distributed: v4 alpha ChocoPhlAn and the EC-filtered UniRef90, the only protein database it distributes, so KOs without an EC number rely on the nucleotide search. MetaPhlAn 4.1.1 uses `mpa_vOct22_CHOCOPhlAnSGB_202403`, the database HUMAnN 4 checks for. KOs come from the v4 `map_ko_uniref90`, scored like `humann` (unstratified adjusted CPM).
  - *Workaround:* HUMAnN 4 alpha.2 checks the MetaPhlAn database by looking for its tag in `metaphlan --version`. MetaPhlAn 4.1.1 prints no tag there, so HUMAnN exits whenever it runs MetaPhlAn itself. HUMANN4 therefore runs MetaPhlAn first (same `-t rel_ab_w_read_stats`) and passes its profile with `--taxonomic-profile`, whose header check accepts it. Its CPU time counts in the same task.
  - *Fix, found while building it:* the HUMAnN 3.9 biocontainer resolves bowtie2 2.2.3, whose `bowtie2-build` rejects `--threads`, so step 1's HUMANN would have failed at its index build on HPC. `containers/humann3` adds bowtie2 2.5.5 to that image, and the HUMAnN 4 image builds on it. Both images ran their bundled demos to gene families, and HUMAnN 4 through `humann_regroup_table`, locally (x86 emulation). A stub run checks the pipeline wiring with both versions. The full databases (~45 GB for 3.9, ~70 GB for 4) are only downloaded on HPC.
  - *Result pending on HPC.* Expected: HUMAnN 4 lower completeness than 3.9 on KOs without an EC number (translated search limited to EC-annotated families). If so, a second run pairing HUMAnN 4 with the 3.9 full UniRef90 (UniRef 2019 on both) would separate the database effect from the algorithm.

### Phase 6

* **Phase 6, step 1 — cost-study pipelines (`workflows/mgnify-subset`, `kmer_functional_profiler.cost`).** Both halves of the cost study run from the subset workflow, so the full-release statistics and the nested builds share one extraction path.
  - *Nested subsets.* `--sample` takes a list (e.g. `100,1000,10000`): MEMBERSHIP and EXTRACT run once at the densest sample, then SUBSET keeps `cluster_rep % N == 0` for each N, so every subset contains the sparser ones, and INDEX builds each under `1inN/index/`. The raw Nextflow trace (as in the fmh benchmark) records each build's time and peak RSS.
  - *Whole-release statistics.* Per-cluster distinct k-mer counts need a cluster's members together, but members are spread over the whole `protein_id` range. EXTRACT therefore tags every member with a bucket (DuckDB `hash(cluster_rep) % --buckets`) and MERGE writes one members table per bucket, so each bucket holds whole clusters: a one-off 1 TB shuffle whose output the full build will reuse. STATS runs pass 1 of the build (`index._n_kmers`, same masking, k and alphabet) per bucket for exact `n_kmers`, plus `n_members`, `n_full_length`, `sum_len`, `max_len`, and samples (hash, cluster) pairs at 1/1000 (`--stats_rate`). COMBINE counts clusters per sampled hash (`groups.tsv`, the promiscuity profile) and predicts sizes (`cost.tsv`).
  - *Cost model (`cost.predict_cost`).* Applies the build's rules to `n_members` and `n_kmers`: *t\_g* from the same expression as the build (`index._t_g`, now shared), candidates *t\_g*·*n\_kmers*·*keep*, where *keep* is the fraction of sampled pairs on k-mers in at most `max_groups` clusters; floored units keep min(candidates, *n\_min*); dense postings max(*t\_dense*, *t\_g*)·*n\_kmers*·*keep*. Bytes use the dtypes `PackedTable.build` picks, with one key per posting and `--sets` value sets per posting (both upper bounds; the build stores each (unit, *p\_in* level) set once, and the nested builds give the real `tier2_sets / postings`). Default grid: *n\_min* 4/8/16 × *t\_cap* 0.05/0.2 × *t\_dense* 0/0.02/0.1, for 1-in-1 to 1-in-10⁴.
  - *Checked:* on random clusters `cluster_stats` gives the build's `n_kmers` exactly, and the predicted postings, dense postings and tier-2 bytes (with the build's set ratio) are within 5% of the built index (`tests/python/test_cost.py`); the workflow's test profile runs nested builds, statistics and COMBINE end to end. Also changed: MEMBERSHIP writes membership sorted by `protein_id` (so each EXTRACT range reads only its row groups instead of the whole table) and checks `--max_protein_id` itself; EXTRACT checks that every member in its range has a sequence, replacing MERGE's global check, which a bucketed MERGE cannot do cheaply.
  - *Expected cost:* pass 1 runs at ~16 M residues/s on a laptop (Rust hashing plus polars distinct), so the 1.07×10¹² residues take ~20 h of one such node, spread over 256 STATS jobs. At full scale the unit id needs 35 bits with the 4 *p\_in* bits, so each set value becomes a u64, and ~4×10⁹ postings put set ids near the u32 limit: the byte estimate will show how much unit-ID encoding (the phase-6 tuning item) has to recover.
  - *Fix after the first HPC run:* MEMBERSHIP on the whole release was killed at 96 GB (third attempt). DuckDB's `memory_limit` was the whole allocation, `--biome root` still joined against all 1.66×10⁹ cluster ids, and the output was sorted globally (5.7×10⁹ rows). Now DuckDB gets 75% of the task's memory, `root` skips the join, and membership is written partitioned by EXTRACT's `protein_id` ranges (`membership/shard=i/`) instead of sorted, so each range still reads only its own files.
  - *Result pending on HPC* (commands in the workflow README). *Next:* check the predictions for the nested subsets against their builds, fit build time and memory against subset size, then choose between the Rust port of the build and a hash-partitioned build.

* **Phase 6, step 2 — nested all-biome builds (`cost-nested`), cost model checked.** The nested run (`--biome root --sample 100,1000,10000`, default build parameters: *t\_base* 0.001, *n\_min* 8, *t\_cap* 0.2, no dense tier, 16-bit fingerprints) completed on HPC. The whole-release statistics run (`cost-full`) is still to come.
  - *Builds* (INDEX rows of the trace; 8 CPUs requested, 3.5–5 used):

    | Subset | Proteins | Residues | Units | Postings | Tier 1 + 2 | Index dir | Build | Peak RSS |
    | --- | --- | --- | --- | --- | --- | --- | --- | --- |
    | 1 in 10⁴ | 0.53 M | 9.8×10⁷ | 0.15 M | 0.34 M | 4.5 MB | 12 MB | 7 s | 4.9 GB |
    | 1 in 10³ | 5.35 M | 1.0×10⁹ | 1.54 M | 3.38 M | 48 MB | 126 MB | 72 s | 6.5 GB |
    | 1 in 10² | 53.2 M | 9.9×10⁹ | 15.4 M | 33.8 M | 514 MB | 1.3 GB | 16.4 min | 39.2 GB |

    Units, postings and table bytes scale linearly (×10.0 per step); ~2.2 postings and ~33 bytes of tier 1 + 2 per unit, tier 2 ~9–10 bytes per posting. Promiscuous candidates dropped (*max\_groups* 64) grow with density: 0, 0.04%, 0.47% of candidates. `tier2_sets / postings` is 0.239, 0.239, 0.243. `postings.parquet` and `units.parquet` make up about 60% of the index directory. The other steps were small: MEMBERSHIP 3.3 min at 3.5 GB, EXTRACT ≤ 73 s per shard, MERGE 6 s at 20 GB.
  - *Build scaling.* Time is ~75 ns per residue at 1 in 10⁴ and 10³ but 99 ns at 1 in 10²; memory is ~4.5 GB of fixed cost plus ~3.5 bytes per residue. Extrapolated to the full release (1.07×10¹² residues, 100× the 1-in-100 subset): ≥ 27 h and ~3.7 TB for one process. The single-process Python build cannot do the full release on one node; the full build has to be partitioned (by k-mer hash, with a global union-find over postings for components) whatever language it is in.
  - *Fix: uint32 overflow in residue sums.* polars' `str.len_bytes()` is UInt32, so `n_residues` in `meta.json` wrapped at 1 in 100 (1.33×10⁹ reported; 9.92×10⁹ true, from 1.33×10⁹ + 2·2³²), and so did the cumulative sum `_batches` splits on: batch ids restarted after each wrap, merging batches into ~2.3× the intended 2×10⁷ residues. Units were never split, so the index is correct, but batch memory was higher than set, which may account for part of the 1-in-100 superlinearity. At a full-release bucket (~4×10⁹ residues) the same would have merged batches too. Both now cast to UInt64 before summing, as does `cluster_stats`' `sum_len` (which cast after).
  - *Cost model against the builds.* STATS and COMBINE run locally on the 1-in-1000 members (47 s, 4.5 GB), predicting 1 in 10³ and 10⁴ with the build's parameters and `--sets 0.24`: postings +0.57% / +0.52%, tier-1 hashes +0.56% / +0.19%, tier-2 bytes +0.63% / +0.52%. Tier-1 bytes were over-predicted (+49% / +99%): the model assumed one tier-1 set per unit, but only units with a posting have one (519 K of 1.54 M), and at 1 in 10⁴ the set ids then fit in u16. `predict_cost` now counts units with a posting as Σ(1 − e^(−*m*)) (Poisson on the expected postings *m*; 525 K predicted vs 519 K built, 52.5 K vs 52.5 K), which puts tier-1 bytes at +0.65% / +0.08%; `tests/python/test_cost.py` now checks tier-1 bytes too. The model is good enough to size the full index from `cost-full` without building it.
  - *Interpretation.* Linear extrapolation gives ~1.5×10⁹ units, ~3.4×10⁹ postings and ~50 GB of tier 1 + 2 for the full release before the u64 value width (35-bit unit ids) and growing promiscuity are counted; `cost-full` will give the real figures. Postings stay under the u32 limit (4.3×10⁹) but not by much, so offsets and set ids need a width check at full scale.
  - *Next:* the `cost-full` run (now with the overflow fix, which matters for its per-bucket batching); then design the hash-partitioned build, and profile where the 1-in-100 build spends its time and memory to decide how much of it moves to Rust.

* **Phase 6, step 3 — feasibility of size and speed reductions (analysis only, no code changed).** Measured on the 1-in-1000 build and its members (local), the cost model (COMBINE grid on the 1-in-1000 statistics, ×1000 for the full release), and the simulation (5 seeds, `gather_zi`). Full-release figures are linear extrapolations until `cost-full` reports.
  - *Where the size is.* At full scale the query-resident index would be ~290 GB, and the tiers are the smaller part. The unit table loads at 154 B per unit (1.54×10⁹ units → ~240 GB): 16 columns, of which `pin_hist` (16 × u32) is 64 B, `name` duplicates `cluster_rep`, and `t_g`/`max_hash_g` are 16 B derivable from `n_kmers`. 66% of units have no posting (90% of singletons), so they can never be hit. Tier 2 is ~31 GB (~34 GB with u64 values), tier 1 ~17 GB. On disk, `postings.parquet` and `units.parquet` add ~75 GB.
  - *Free cuts (no change to query results).* (1) Tier 1 and the build's components are not read by the query, which recomputes components on the hit subgraph (`query.py`); `u_g` is unused too. Dropping them saves ~17 GB and removes the only global graph step (union-find over ~3.4×10⁹ edges) from a partitioned build. Tier 1 comes back with the two-stage lookup, if the Rust query needs it. (2) Keep unit rows only for units with a posting, renumbered (5.3×10⁸ at full scale, 29 bits) and slimmed (`cluster_rep`, `n_kmers`, singleton flag, `m_g`, `pin_sum` f32, `len_cv` f16, `pin_hist` as u8 while `m_g` ≤ 255 (it is ≤ 10 at the defaults); dense-tier columns only with a dense tier): ~36 B per unit, ~19 GB, memory-mapped by column. (3) `postings.parquet` becomes optional (inspection only). Together: ~290 GB → ~53 GB, with no accuracy change.
  - *Precision.* Fingerprints 16 → 8 bits saves ~1 B per posting (11%) but raises ε from ~6×10⁻⁵ to ~1.6×10⁻² per lookup: ~0.3 false hits per 150-bp read at *t\_max* 0.2 (~20 lookups), several times the genuine chance-match background the plan estimates (~0.3% of off-target k-mers, ~0.06 per read). 12 bits saves only if bit-packed (~0.5 B, 5%). Neither is worth it. *p\_in* bits: re-packing the 1-in-1000 postings with 4/3/2/1/0 bits gives 808/783/725/575/522 K sets and 30.9/30.7/30.3/29.1/28.6 MB, so dropping *p\_in* saves ≤ 7% and loses `em_pin`'s 20–30% lower scatter. What *p\_in* bits do change is the value width at full scale: renumbered units (29 bits) plus 4 bits is 33 bits (u64); 3 bits, or *p\_in* in a separate u8 array beside u32 unit values, keeps u32 (~1 B per posting, ~10%).
  - *Density.* *n\_min* sets the size: floored non-singletons hold 96% of postings (412 K × 8), singletons 3.7%. Tier 2 at full scale by *n\_min* 2/4/8/16 (*t\_cap* 0.2): 8.6/16/31/59 GB. Simulation (units only): completeness 0.49/0.73/0.85/0.94, Spearman 0.74/0.77/0.84/0.88. So halving *n\_min* costs 12 points of unit-level completeness; the Pfam-level cost should be smaller (a Pfam present is usually carried by several clusters) and needs the Pfam benchmark. *t\_base* changes only singletons: 0.0001/0.001/0.01 give 30/31/45 GB. The dense tier is the most expensive feature: +120 GB at *t\_dense* 0.02 and +250 GB at 0.1 (*n\_min* 8), against +0.05 Spearman on KO units; it stays off for the first full build. Dropping singletons altogether saves 3.7% of postings and 22% of the (slimmed) unit rows; 90% of singletons already have no posting.
  - *Query rate.* The query looks up every read k-mer with hash ≤ *t\_max* = *t\_cap*, since 43% of non-singletons (fewer than 160 distinct k-mers) reach the cap. *t\_cap* 0.1 halves lookups for 5% fewer postings but leaves 19% of floored units below *n\_min* (3.6% at 0.2; 43% at 0.05); the simulation's units are long enough that 0.1 changes nothing and 0.05 costs 1 point of purity, so short MGnify clusters need the real benchmark. Lookup cost at full scale is not measured yet: the prototype materialises 8 B per key for lookup (~27 GB at full scale), so the memory-mapped Rust lookup comes first anyway.
  - *Build time and memory.* Profile of the 1-in-1000 build (81 s, 4.9 GB locally): pass 1 (exact distinct k-mers per unit) 47 s (58%), almost all in polars `unique` over every k-mer (the Rust hashing is 10.7 s over all passes); `_len_cv` (a fourth hashing pass) 6.3 s; packing and writing 5.7 s. A Rust distinct-count kernel over unit-sorted members would roughly halve build time. Peak memory holds about three copies of the sequences (original, masked, and `partition_by` batches, with `members` kept alive for later joins), which fits 39 GB at 1 in 100 (9.9 GB of sequence); holding one copy is a small Python change.
  - *Partitioned build.* With components gone, everything except *n\_groups* (*p\_out*) is local to a unit and so to a `cost-full` cluster bucket: *n\_kmers*, *t\_g*, candidates, *p\_in*, `len_cv`, the floor's top-*n\_min*. *n\_groups* feeds the promiscuity cut and the score. Counting it from candidates only (a hash-partitioned reduce of each bucket's candidates) matches the exact count for 99.45% of postings at 1 in 1000 and keeps 354 promiscuous candidates above 64 where the exact count drops 4,603. But sharing grows about 4× per decade of density (shared postings 0.25%, 1.1%, 4.0% at 1 in 10⁴, 10³, 10²), so the approximation will be much worse at full scale. The exact alternative broadcasts a Bloom filter of candidate hashes (~10¹⁰ at full scale, ~14 GB at 10 bits each) to every bucket job, which emits (hash, unit) presence for matching k-mers only, then reduces by hash range. `cost-full`'s `groups.tsv` will show how much sharing there is to get right.
  - *Recommended order.* (1) Free cuts: drop tier 1, build components and `u_g`; slim and trim the unit table; make `postings.parquet` optional; hold one copy of the sequences. (2) Rust distinct-count kernel for pass 1. (3) Partitioned build over the `cost-full` buckets, with exact *n\_groups* via the Bloom-filter broadcast. (4) Only then trade accuracy, on the Pfam benchmark: *n\_min* 4 vs 8, *t\_cap* 0.1 vs 0.2, singletons in or out, dense tier rate. Fingerprint and *p\_in* precision stay as they are.

* **Phase 6, step 4 — free cuts (index format 2).** Item (1) of step 3's order, none of which changes what the query computes.
  - *Dropped:* tier 1, the build's connected components and `u_g`, and the `tier1_per_unit` parameter (supersedes phase 2's "tier 1 is the top 4 by score per unit, mapped to component ids"; tier 1 returns with the two-stage lookup if the Rust query needs it). The query already finds components on the hit subgraph. The cost model no longer predicts tier-1 bytes, and sizes set values by units with a posting instead of all clusters.
  - *Unit table:* only units with a posting get a row, renumbered from 0 in `cluster_rep` order (`unit_pfam.parquet` follows). Columns: `unit`, `cluster_rep`, `n_members`, `n_kmers`, `t_g`, `max_hash_g`, `m_g`, `pin_hist` (smallest unsigned type holding the largest `m_g`), `pin_sum` and `len_cv` as f32, plus the dense-tier columns with a dense tier. Dropped `name` (the query derives it from `cluster_rep`; sourmash imports keep theirs), `n_full_length`, `n_counting`, `n_candidates`, `n_promiscuous` (still summed into `promiscuous_dropped`). Stats gain `n_clusters`; `n_units` now counts indexed units. The query reads format-1 indexes too (it ignores `tier1_per_unit`).
  - *Other:* `postings.parquet` is written only with `--postings` (`build_index(..., postings_parquet=True)`); the members table drops its sequences once batched, so the batches hold the only copy.
  - *Result, 1 in 1000 (local):* 525,764 of 1,543,842 units indexed (the cost model predicted 525 K); unit table 238 MB → 34 MB in memory (154 → 64 B per unit); tier 1 (16.8 MB) gone; tier 2 unchanged byte for byte (30.9 MB, same postings and sets); index directory 126 MB → 39 MB; build 81 s / 4.9 GB → 76 s / 4.4 GB peak RSS. Extrapolated ×1000: ~290 GB query-resident → ~31 GB tier 2 + ~34 GB units.
  - *Not done:* deriving `t_g`/`max_hash_g` from `n_kmers` and `n_members` at load (16 B per unit; sourmash imports have no `n_kmers`, so it waits for the Rust loader), f16 `len_cv` (polars has no f16), and memory-mapping the unit table by column (with the Rust lookup). Together they would reach step 3's ~36 B per unit.
  - *Next:* (2) the Rust distinct-count kernel for pass 1.

* **Phase 6, step 5 — Rust distinct-count kernel for pass 1 (`_core.distinct_kmers`).** Item (2) of step 3's order. The kernel takes a batch's sequences and their unit ids (non-decreasing, which `load_members` and `_batches` guarantee, and which it checks), hashes every k-mer of each unit's run of members into one buffer, sorts and dedups it, and returns the count per unit. Runs are split over threads (`std::thread::scope`, slices of about equal residues; no new dependency). It replaces the polars `unique` over every (member, hash) in `index._n_kmers`, which `cost.cluster_stats` (the `cost-full` STATS step) also calls.
  - *Checked:* identical `n_kmers` for all 1,543,842 units at 1 in 1000 (the kernel also returns 0 for the 850 units without a k-mer, which callers filled with 0 before); rebuilt index identical to step 4's (units table, stats and every `.npy`); Rust test on hand-built runs, including unsorted ids and a length mismatch.
  - *Result, 1 in 1000 (local, 12 cores):* pass 1 46.9 s → 1.4 s (~700 M residues/s against ~21 M); whole build 76 s → 30 s, peak RSS unchanged (4.4 GB). STATS on `cost-full` should now be bound by reading its bucket and by the sampled-pairs pass, not by pass 1 (step 1 estimated ~20 node-hours for pass 1 over the release).
  - *Where the time is now (profile, 28 s):* `_kmers` 15.8 s, the passes that keep (unit, member, hash) rows (passes 2 and 3 at *t\_max* 0.2, and `_len_cv`): `hash_proteins` 6.5 s (single-threaded) and polars `unique`/joins the rest; `_len_cv` 5.1 s of it; packing and writing 2.5 s. Peak memory is still set by holding the sequences and the pass-2/3 tables, which the partitioned build (item 3) splits by bucket.
  - *Next:* (3) the partitioned build over the `cost-full` buckets, with exact *n\_groups* via the Bloom-filter broadcast. The HPC environment needs the rebuilt extension before `cost-full` runs.

* **Phase 6, step 6 — partitioned build (`kmer_functional_profiler.partition`, `mgnify-subset --build`).** Item (3) of step 3's order. The build runs over the MERGE buckets (whole clusters each), with units keyed by `cluster_rep` until the last step: CANDIDATES (per bucket: pass 1, unit table, candidate hashes), BLOOM (once: Bloom filter of all candidate hashes, 7 probes by double hashing, `--bloom_bits` per hash; plus 1,025 quantiles of the candidate hashes), PRESENCE (per bucket: every distinct (unit, k-mer) up to the largest candidate hash that passes the filter, with its counting members and whether it is the unit's own candidate, split into `--ranges` hash ranges), GROUPS (per range: exact *n\_groups*, then candidate rows only, split by bucket), POSTINGS (per bucket: score, promiscuity cut, floor, `len_cv`, Pfam) and PACK (once: renumber and write). `build_index` was split into shared helpers (`_prepare`, `_select_postings`, `_unit_pfam`, `finish_index`) that both paths use, so the rules exist once.
  - *Choices.* *n\_groups* is exact, not the candidates-only approximation step 3 rejected: the filter only thins the rows each bucket emits, and GROUPS drops its false hits because they have no candidate row. Hash ranges are cut at candidate-hash quantiles rather than at equal widths: candidates sit below each unit's own *t\_g*, so equal widths put most of them in the lowest range (seen in the test: 1,130 real rows per range against several times that in range 0 before). PRESENCE stops at the largest candidate hash, since no k-mer above it can be a candidate. Stage outputs are new files, never rewritten, because Nextflow stages inputs as links into upstream work directories. No dense tier yet (it stays off for the first full build).
  - *Checked:* on 300 random clusters in 3 buckets at *k* 4 (so k-mers are shared across buckets and the promiscuity cut drops some), with a 2-bit filter so ~40% of presence rows are false hits, the partitioned index equals `build_index`'s file for file (units, postings, unit Pfam, every tier-2 array, stats); the same holds at 1 in 1000 (8 buckets, 8 ranges, 10-bit filter) against step 5's index, and for the workflow's `-profile test,test_build` run (two buckets) against `build_index` on the same members.
  - *Result, 1 in 1000 (local, stages run one after another):* CANDIDATES 11.3 s over 8 buckets (≤ 1.5 s each), BLOOM 0.8 s, PRESENCE 19.1 s (≤ 2.7 s each), GROUPS 0.8 s, POSTINGS 6.4 s, PACK 2.3 s; 42 s in all against 30 s single-process, as every per-bucket stage re-reads, masks and batches its members. 10.8 M candidate keys (counted per bucket), 11.7 M presence rows (1.08 per candidate), 2.6% of them false hits (more in the upper ranges, where fewer k-mers are candidates), ~11.5 B per row on disk, ranges within 10% of each other.
  - *Full-scale expectation (×1000):* Bloom filter ~14 GB, read by every PRESENCE job; ~1.2×10¹⁰ presence rows (~140 GB) over `--ranges` 256, ~47 M rows per GROUPS job; per-bucket steps as STATS plus one pass at *t\_max*. BLOOM's probes are numpy (`bitwise_or.at` over ~7×10¹⁰ positions, likely an hour or more); a Rust kernel if it matters. PACK still holds every posting (~3.4×10⁹) and the string-keyed set dedup of `PackedTable.build` in one process: it is the remaining single-node step, and the next to split (pack tier 2 by hash range, whose keys are already contiguous, with set ids made global afterwards).
  - *Next:* split PACK; then run `--build` on the whole release (it needs `--pfam true`, unlike `cost-full`).

* **Phase 6, step 7 — tier 2 packed by hash range (PACK split into UNITS, PACK\_RANGE, CONCAT).** PACK held every posting and ran `PackedTable.build` over them in one process. Now POSTINGS also computes each unit's index columns (`m_g`, `pin_hist`, `pin_sum`, f32 `len_cv`; `write_index`'s `unit_columns`) and writes its postings sorted by hash. UNITS (once) concatenates the unit tables of units with a posting, numbers them in `cluster_rep` order, casts `pin_hist` to the width the largest `m_g` needs, and fixes the tier-2 layout from the total postings and *t\_max*; pack ranges are the candidate quantiles rounded down to key boundaries, so hashes that share a key (and so pool their value sets) never fall in two ranges. PACK\_RANGE (per range) reads that range's postings from every bucket (row-group pushdown on the hash-sorted files), maps them to unit ids and builds a `PackedPart`: sorted keys, a local set per key, and the range's distinct sets with a 128-bit content hash each. CONCAT (once) deduplicates the sets globally, concatenates fingerprints, set ids and offsets bucket counts, and writes `meta.json` with `write_meta`.
  - *Changed in `PackedTable.build`* (both builds, so they stay identical): value sets are numbered in order of their content hash (two splitmix64 sums over the set's values, 128 bits; a collision would merge two sets, ~10⁻²⁰ at 10⁹ sets) instead of by string rank, which also replaces the string-keyed dedup step 2's ponytail note flagged; and the offsets bucket count comes from the number of (hash, value) rows instead of distinct hashes, because the partitioned layout is fixed before any range has seen its hashes. Postings and distinct hashes differ by 0.2% at 1 in 1000, the same bit width, so the layout did not change: the rebuilt 1-in-1000 tier 2 has the same offsets, fingerprints, bytes (30,920,774) and value set per key as step 5's; only set numbering differs. The cost model already sized buckets by postings. Format stays 2 (the layout fields in `meta.json` mean the same).
  - *Checked:* `PackedTable` built in three ranges with 4-bit fingerprints (many pooled keys, shared sets) and concatenated equals the table built at once, dtypes included; the random-cluster test, the 1-in-1000 run (8 buckets, 8 ranges) and `-profile test,test_build` all equal `build_index` again (`candidates_expected` to float rounding).
  - *Result, 1 in 1000:* UNITS 0.1 s, PACK\_RANGE 1.2 s over 8 ranges, CONCAT 0.7 s (the old PACK took 2.3 s). The 8 ranges hold 2.04 M local sets against 0.81 M distinct ones (2.5×), and 262 K to 750 K postings each: ranges are cut at candidate quantiles, and postings sit lower in hash than candidates.
  - *Full-scale limits:* no step holds every posting now, but CONCAT holds tier 2 (~31 GB) and every range's sets before dedup (more duplicated at 256 ranges than at 8; up to ~2×10⁹ at 16 B of hash each plus contents), and UNITS the ~34 GB unit table. If CONCAT outgrows a 128 GB node, dedup sets in a reduce keyed by set hash first. Pack ranges could be balanced on posting-hash quantiles recorded by POSTINGS.
  - *Next:* run `--build` on the whole release (it needs `--pfam true`, unlike `cost-full`), after `cost-full` has checked the predicted sizes.

* **Phase 6, step 8 — Bloom filter in Rust (`_core.bloom_insert`, `_core.bloom_contains`).** BLOOM's numpy probes (`bitwise_or.at` over 7 positions per key) were the slow part of that step at scale, and PRESENCE's numpy lookups the largest part of its time. Both are now Rust, with a blocked filter: each key's 7 bits lie in one 64-byte block (block from a remix of the hash by multiply-shift, bits from a second remix, 9 bits each), so an insert or lookup touches one cache line instead of up to 7, which matters for a ~14 GB filter that no cache holds. Inserts are threaded without atomics: each thread owns a run of blocks and scans every hash, setting only its own bits. Lookups split the hashes over threads. The filter length is rounded up to whole blocks. PRESENCE memory-maps the filter, so jobs on one node share its pages.
  - *Cost of blocking:* 0.97% false hits at 10 bits per key against 0.82% for the unblocked filter (10.8 M random keys); at 1 in 1000, 3.1% of presence rows are false hits (was 2.6%), all dropped by GROUPS as before.
  - *Checked:* Rust tests (no false negatives, false-hit rate between 0.5% and 2%, threaded insert equal to a sequential one on a block count that does not divide by the threads, partial blocks rejected); the Python filter test; the 1-in-1000 partitioned build and `-profile test,test_build` still equal `build_index`.
  - *Result:* 10.8 M keys insert in 0.03 s against 0.88 s for numpy (29×) and look up in 0.03 s against 0.48 s (16×); 10⁸ keys insert in 0.42 s and look up in 0.25 s (12 cores). At 1 in 1000 BLOOM takes 0.2 s (was 1.0 s; the rest is reading candidates and the quantiles) and PRESENCE 9.9 s (was 19.8 s). At full scale (~10¹⁰ keys) filling the filter should take about a minute, against ~14 min extrapolated for numpy; BLOOM is then bound by reading the buckets' candidate files (~80 GB).

* **Phase 6, step 9 — set dedup reduce before CONCAT (`dedup_sets`, `SetSlice`, `mgnify-subset` DEDUP).** Step 7 left CONCAT holding every pack range's value sets before deduplication (2.5× the distinct sets at 8 ranges, more at 256). Each `PackedPart` keeps its sets in content-hash order, so DEDUP job *s* of `--ranges` takes, from every part, the run of sets whose first hash word lies in [*s*·2⁶⁴/*n*, (*s*+1)·2⁶⁴/*n*), deduplicates them, and writes the distinct sets in hash order with each part's local-to-job id map. Global set ids are then a job's base (the running sum of earlier jobs' counts) plus its job id, so CONCAT only concatenates. `PackedTable.concat` takes the slices; without them it deduplicates in one slice, as `build_index` does, so there is still one code path.
  - *Checked:* `PackedTable` from parts with sets deduplicated in three set-hash ranges equals the table packed at once (dtypes included); the random-cluster test, the 1-in-1000 run (8 buckets, 8 ranges) and `-profile test,test_build` equal `build_index`.
  - *Result, 1 in 1000:* DEDUP 0.7 s over 8 jobs, CONCAT 0.2 s (was 0.6 s); each job holds 100,713–101,380 distinct sets (±0.4%, the content hashes being uniform). CONCAT now holds tier 2 (~31 GB at full scale) and one part's id map at a time; each DEDUP job holds 1/*n* of all local sets.

* **Phase 6, step 10 — query cost plan (analysis only, no code changed).** With the build partitioned, the other half of "reasonable compute" is profiling a metagenome against the full index. Design in Query cost above; sylph v1.0.0's source read for its sample sketch, two-stage `.syl2db` and cuckoo-filter dedup. All figures are estimates until step Q1 measures them; Q1 decides which of Q2–Q7 are built, and in what order.
  - *Where the cost is.* Hashing is fixed per read. At full scale the query is limited by (1) per-read hit rows in polars and the Python gather/model code, (2) random lookups into a ~31 GB tier 2 that may not be resident, (3) posterior draws (5–10× the rest). The dense rate *t\_max* = 0.2 matters through the sketch size (~8×10⁸ sampled k-mers per 10 Gbp), which sorting and spilling bound.
  - *Order.*
    - **Q1. Measure first, on HPC, with a real metagenome.** Development of Q2–Q7 is guided by these numbers, and a step is built only if Q1 shows its stage matters at full scale (see *Rule* below), to avoid optimising what does not cost. Run the current query unchanged, apart from logging, as a grid:
      - *Reads:* one public human gut metagenome of ~40 M read pairs (150 bp, ~12 Gbp; run to be chosen), plus a nested read ladder drawn from it with a fixed seed so each subset contains the smaller ones: 0.4, 1.2, 4, 12 and 40 M pairs.
      - *Indexes:* the nested all-biome builds (1 in 10⁴, 10³, 10²; 1 in 10 if built), default parameters. Lookups do not depend on the index subset and hits grow with it, so the read × index grid gives both scaling exponents and extrapolates to 40 M pairs against the full index.
      - *Settings:* `draws` 0 over the whole grid; `draws` 100 on the 1-in-100 column only. Cells that run out of memory or time are kept as results: the ladder makes sure smaller cells succeed and shows where the prototype breaks.
      - *Recorded per cell:* wall time, CPU time and peak RSS from the Nextflow trace; per-stage time and RSS from timers in `profile` (read/hash, lookup, hit aggregation, gather, each fit, presence, posterior); and the counts the design choices hinge on: sampled k-mers, distinct sampled k-mers (the value of sketch dedup), hits and hit rows (hit fraction), (unit, hash) pairs, detected units, and component sizes (largest component, which bounds the low-memory mode).
      - *Lookup at full size:* a synthetic tier 2 (random hashes, 3.4×10⁹ postings, the real layout), resident vs cold, until the full build exists.
      - *Pipeline:* a query-cost entry in `mgnify-subset` (read subsampling, then query per read subset × index), results in a `query_cost.tsv` next to the build costs; the query prediction in `cost.py` (lookups = reads × k-mers per read × frames × *t\_max*; hits by the fitted exponent) is checked against it, as the build model was.
      - *Rule:* after Q1, a later step goes ahead only if its stage is a large share (tentatively ≥ 20%) of time or peak memory at the extrapolated full-scale cell, or blocks the laptop gate; the ladder is re-run after each step, so priorities are re-read from new data rather than from this plan.
    - **Q2. Free cuts (no change to shipped estimates).** Fit only the shipped estimator by default; em, zib, zip, winner-take-all, uniqueness-first and `kmers_out` behind a benchmark flag. Pass 1 aggregates to per-(unit, hash) counts; per-read rows only when `draws` > 0, and only for detected units (second read pass, in-RAM set of their hashes). Skip Gibbs for one-unit components. Parallelise gather and EM by component.
    - **Q3. Rust lookup and aggregation.** A `_core` query kernel over the memory-mapped tier-2 arrays: sampled hashes in batches, radix-sorted and counted, resolved by one forward walk over `offsets`/`fingerprints`/`set_ids` (no materialised keys), set expansion and the `max_hash_g` check, output per-(unit, hash) counts with `pin_q` and holders as numpy. `madvise` sequential; `--memory` sets the partition count: the sketch is radix-partitioned by the pack hash ranges and spilled to scratch when it exceeds the budget, so a sample needs one pass over the index at any budget (see Low memory). Hit counts are spilled by unit range, the `max_hash_g` check becomes a sorted join with the mmapped unit table, and gather, EM and posterior stream one component at a time. Checked against the Python `unit_hits` path on the fmh fixtures (exact counts).
    - **Q4. Unit-major dense tier.** Rebuild the dense tier as per-unit blocks (sylph's stage 2): sorted hash fingerprints plus *p\_in*, Golomb-Rice or Elias-Fano, read by `pread` for detected units only, probed by pass 2. Then re-run the dense-tier comparison of phase 4, step 9 for cost (disk, RSS, time) and accuracy; the phase-4 gate's "≤ 2× index size" becomes a disk bound.
    - **Q5. Deployment and cohorts.** Index staged once per node (local SSD or `/dev/shm`), several samples per node sharing its pages; a cohort mode that merges many samples' sorted sketches in one index pass; Nextflow module for both. Record per-sample cost at cohort sizes 1, 10, 100.
    - **Q6. Cost ablations on the Pfam benchmark** (with the phase-7 accuracy ablations, same runs; KO alongside), each reporting time, peak RSS and index bytes beside Pfam completeness, purity and abundance error, to draw the cost/accuracy front:

      | Knob | Values | Cost it moves | Expected accuracy risk |
      | --- | --- | --- | --- |
      | *t\_cap* (sets *t\_max*) | 0.2, 0.1, 0.05 | sampled k-mers ∝ *t\_cap*; sketch memory | short clusters below *n\_min* (19% at 0.1, 43% at 0.05) |
      | Frames | stop-free, all six, frame caller (one) | hashing and sampled k-mers ~1.5 : 4 : 1 | high-GC reads, reads crossing stops |
      | Posterior | *D* = 0, 20, 100; Gibbs only on multi-unit components | model time 1–10× | interval calibration (phase 4, step 20) |
      | Estimators | shipped one only vs all | model time ~3–4× | none for the shipped column |
      | Dense tier (unit-major) | off, 0.02, 0.1, full | disk; pass-2 time ∝ detected units | none; accuracy gain measured |
      | Singletons | in at *t\_base*, out | unit table −22%, postings −4% | rare singleton functions |
      | Lookup mode / `--memory` | random resident, sorted resident, partitioned spill at 64 / 12 / 8 GB, re-read reads per shard, cohort | time vs RSS vs scratch disk | none (identical counts) |
      | Two-stage screen | off, component-level sparse screen | lookups for dense stage only | loses rare short clusters: the floor's gain |
      | Read dedup (cuckoo, as sylph) | off, on | sketch memory +, per-read state | abundance on duplicate-heavy libraries |
    - **Q7. Only if Q1–Q6 leave the query too costly:** the two-stage screen as an opt-in fast mode, a Bloom prefilter for non-resident random lookups, 12-bit packed fingerprints.
  - *Gate for this part of phase 6:* a 10 Gbp metagenome profiles against the full index on one HPC node in ≤ 1 h and ≤ 64 GB RSS with the index not preloaded (≤ 16 GB with the partitioned spill), with counts identical across lookup modes; and on a laptop (12 GB RAM, 8 cores, NVMe) at `--memory 8G`, peak RSS ≤ 8 GB and wall time ≤ 2× the HPC node's, with identical results. Targets to confirm after Q1.

* **Phase 6, step 11 — Q1 query-cost workflow and stage timers (`profile(timer=)`, `query --stats`, `mgnify-subset --query`).** Q1's measurement, built; the run is pending on HPC.
  - *Reads:* ERR7746321 (PRJEB49206, Hadza gut, NovaSeq 6000, ~140 bp paired, 29 Gbp, ~106 M pairs by ENA's counts; to be confirmed by the LADDER count). Deeper than the planned ~40 M pairs, so every ladder rung, 40 M included, is a proper random subset, and deeper rungs can be added later without a new sample.
  - *Stage timers.* `Timer` records wall time, CPU time and the process's peak RSS per stage: load, keys (tier 2's materialised keys, built once on first lookup), hash (the Rust read kernel, timed per batch), lookup (`unit_hits`: tier-2 lookup, set expansion, `max_hash_g`), aggregate (per-read rows and per-(unit, hash) counts), components, gather, presence, dense, fit\_em, fit\_zi, fit\_zib, fit\_zip, result, posterior, baselines, total. Peak RSS is the high-water mark at each stage's last exit (`getrusage`), so the stage that raises it shows as a step. Counts: reads, sampled k-mers, distinct sampled k-mers (estimated on 1 in 256 of hash space, so the count costs no memory), hit k-mers (sampled occurrences with ≥ 1 unit), hit rows, per-read rows, (unit, hash) pairs, hit and detected units, and the connected components of hit units linked by shared k-mers (number, largest in units and pairs). The query is otherwise unchanged.
  - *Workflow* (`--query true`): FETCH\_READS (ENA, md5-checked; or `--query_reads` local files), LADDER (nested subsets by rank in one seeded permutation, so each rung holds the smaller ones, exact sizes, mates in step), QUERY per rung × index (`draws` 0 everywhere; `--query_draws` 100 on `--query_draws_on` 1in100 only), QUERY\_COST → `query_cost.tsv` (one row per cell, stage columns `{stage}_wall_s`, `_cpu_s`, `_peak_rss`). QUERY runs at `--query_memory` (128 GB) without retries and failures are ignored, so cells that break the prototype show as FAILED in `trace.tsv` and the rest still run.
  - *Fixed:* `posterior_zi` crashed when no unit was detected (an empty `bincount` bound); it now returns an empty table.
  - *Checked:* the `test,test_query` profile (fixture reads, ladder 10 and 40 pairs, both test indexes, draws 3 on 1in1) and a stub run of the download path; the CLI test checks the stats against the profile (hit rows = total hits, hit units = rows); a test checks that ladder subsets nest and pair.
  - *Not built yet:* the synthetic full-size tier 2 for resident vs cold lookups, and the check of `cost.py`'s query prediction against `query_cost.tsv`; both come after the run.

* **Phase 6, step 12 — first Q1 run: every cell out of memory; stats written per stage.** Against `full-build`, all five rungs (0.4–40 M pairs) were OOM-killed at 128 GB within minutes, the 0.4 M rung included, so the cost does not come from the reads. `--stats` was written only on success, so the failed cells recorded nothing.
  - *Suspect:* the `keys` stage. `stored_keys` materialises tier 2's full keys as uint64 (~27 GB at ~3.4×10⁹ keys) through temporaries of the same size (bucket repeat, shift, fingerprint cast).
  - *Changed:* `Timer` rewrites the stats JSON at every stage start and end, with `running` naming the stage in progress, and logs each stage's start and end with current and peak RSS to stderr (`.command.err`), so a killed cell shows where it died. The ladder gains 10 k and 100 k rungs, where the read-side cost is negligible.
  - *Re-run (10 k to 40 M pairs):* every cell killed in 30–70 s; the 10 k cell's stats show `load` done at 43.5 GB peak RSS (20.6 s wall, the unit table read into memory) and `running: keys`.
  - *Fixed:* `stored_keys` shifts and fills the key array in place (`<<=`, `|=` with a buffered cast) instead of building the bucket, shifted, cast and result arrays side by side. Measured with `tracemalloc` on 2×10⁸ keys at ~4 keys per bucket: peak 4.8 → 2.4 GB, identical keys; ×17 for the full tier 2 (~3.4×10⁹ keys), ~82 → ~41 GB, so `load` + `keys` ≈ 85 GB, under 128 GB.
  - *Re-run:* memory fits (10 k cell: `load` 43.5 GB, `keys` 87.9 GB peak, no later stage higher), but every cell exited 1 after `gather`: `index.units["pin_hist"].to_numpy()` panics with Polars' maximum length, the full index's 16 levels per unit (~8.5×10⁹ values) exceeding 2³². `rt64` stays rejected (see the Pfam join fix), so `em_pin` takes the column as a Series and gathers the hit units' rows only.
  - *Lookup is I/O-bound:* 70.7 s wall against 3.7 s CPU for 239 K sampled k-mers (~0.3 ms per lookup): random page faults into the memory-mapped tier 2 on `/hps`. Q1's resident vs cold question, measured on real data; QUERY now reads tier 2 once into the page cache before the query (`--query_preload`, on by default; its time goes to `.command.err`), so the next run shows the resident case. It also copies the files the query reads (meta, units, tiers) to node-local scratch first (`--query_scratch`, default `$TMPDIR`; in place if the copy fails), so pages evicted by other jobs or by the query's own growth are re-read locally, not from `/hps` (which may be GPFS, whose cache the preload does not fill). `--query_in_memory` (`query --in-memory`, off by default) reads the tiers into RSS instead of memory-mapping them: no page faults, ~31 GB more peak RSS.
  - *First complete cell (10 k pairs, before the scratch copy):* 127 s wall, 98.6 GB peak RSS. `load` 49 s and 43.5 GB, `keys` 42 s and 87.9 GB, `lookup` 28 s wall against 3.4 s CPU; every per-read stage under 0.5 s. `load` and `keys` were ~3× slower than the run before, likely from the ladder's cells reading `/hps` at once. `fit_zip` and `result` raised peak RSS by 10.7 GB at 10 k reads: `em_pin`'s `pin_hist` rows, the unit columns joined into the result and the baselines' `m_g`/`t_g` all went through the full unit table.
  - *Fixed:* a `hit_units` stage filters the unit table to the hit units once (~32 K rows at 10 k pairs); `em_pin` takes those rows (unit, `pin_hist`), and the result and baselines join them. Measured on a 5×10⁷-unit table (50 chunks, 16 × u8 `pin_hist`): +0.78 → +0.02 GB RSS, identical output; ~8 GB at full scale.
  - *100 k pairs (scratch copy, preload):* 68 s wall, 111 GB peak RSS. The copy made `load` 10 s, `keys` 11 s and `lookup` CPU-bound (1.2 s wall = CPU for 2.4 M sampled k-mers); `lookup` raised RSS by 14 GB of mapped tier-2 pages (file-backed, reclaimable). `gather` and `presence` each added ~4.3 GB: both converted the full `t_g` column (~5.3×10⁸ × 8 B) to index it at hit units. The benchmark fits (`fit_em`, `fit_zib`, `fit_zip`, baselines) took 36.7 s of the 68 s, `fit_zip` alone 30.6 s (11× for 10× reads); cutting them waits until the benchmark has settled which estimator ships (Q2).
  - *Fixed:* after `hit_units`, units are renumbered 0.. in the hit units' order, so `t_g`, `m_g`, `pin_hist`, `pin_sum` and `len_cv` (and their dense forms) come from the hit units' rows; the output maps back to the index's ids. `em_pin` takes a plain array again. The only full-column conversions left are `max_hash_g`/`max_hash_dense`, which the lookup indexes per batch (Q3). Profiles on the fixtures (18 of 20 units hit, with and without a dense tier, 5 draws) are byte-identical to before.
  - *Next:* Q3.
  - *Next:* re-run the ladder. Walking the memory-mapped `offsets`/`fingerprints` without materialising keys (Q3) and memory-mapping the unit table remain the full-scale fixes.


* **Phase 6, step 13 — Q3, first part: Rust tier-2 lookup over the memory-mapped arrays (`_core.unit_hits`, `_core.packed_lookup`).** The 10 k and 100 k cells showed materialising tier 2's keys (`keys`: ~44 GB of the 88 GB peak, 11–42 s) as the largest fixed cost of a query, so the lookup moves to Rust first.
  - *Kernel* (`crates/core/src/packed.rs`): per sampled hash, the bucket from its top bits addresses `offsets` directly, the bucket's few fingerprints are scanned, and the set is expanded with the `max_hash_g` check; rows (unit, hash, read, `pin_q`, holders) in the same order as the numpy path. Arrays of any unsigned width are read in place, so a memory-mapped table stays on disk except for the buckets and sets touched. `stored_keys` and the `keys` stage are gone; `PackedTable.lookup` uses the same kernel.
  - *Checked:* equal to the replaced numpy path (materialised keys + `searchsorted`), set ids and hit rows in order, on a random table with 4-bit fingerprints (shared, merged sets), random `max_hash_g`, hashes above `max_hash`, in memory and memory-mapped; fixture profiles (with and without a dense tier, 5 draws) byte-identical.
  - *Result, 5×10⁷ keys (memory-mapped, 407 MB), 2.4 M lookups (a 100 k-pair batch):* 0.62 s and +0.52 GB RSS (the table's touched pages) against 0.85 s and +2.6 GB for the numpy path with its keys. At full scale this removes ~44 GB of anonymous memory and the `keys` stage.
  - *Not built (Q1 rule: no measured cost yet):* batch radix sort and a forward walk (the per-hash lookup is O(1) via `offsets`; sorting only helps a non-resident index, which the scratch copy and preload now avoid); aggregation to per-(unit, hash) counts in Rust (`aggregate` took 0.03 s at 100 k pairs); `madvise`, `--memory` partitioning and spill, component streaming; the `max_hash_g` check as a join with a memory-mapped unit table (the lookup still takes `max_hash_g` as a full ~4 GB array). Threads: the kernel is single-threaded (~0.26 µs per lookup here).
  - *Next:* re-run the ladder; then memory-map the unit table by column (`load`'s 43.5 GB), which also serves `max_hash_g` to the lookup without a copy.
* **Phase 6, step 14 — unit table memory-mapped by column (`UnitTable`, `units.<column>.npy`, `unit-columns`).** `load` read the whole unit table into memory (43.5 GB at full scale, 20–49 s), and the query then used ~32 K of its ~5.3×10⁸ rows plus `max_hash_g` for the lookup.
  - *Format:* beside `units.parquet` (kept, for inspection and string columns), `write_index` and `concat` write each numeric column but `unit` (the row number) as `units.<column>.npy`, fixed-size arrays (`pin_hist`) as 2-D; one column in memory at a time. `Index.units` is a `UnitTable`: `units[column]` memory-maps the column (read into memory with `--in-memory`), `rows(units)` gathers whole rows for the hit units (strings, e.g. sourmash `name`, by a filtered parquet scan), `frame()` reads the whole table for tests. Without the `.npy` files (indexes built before this step, e.g. `full-build`) a column is read from the parquet on first use, so old indexes still query, at one column's memory rather than the table's; `kmer-functional-profiler unit-columns <index>` adds the files in place. QUERY's scratch copy includes them. `max_hash_g`/`max_hash_dense` go to `_core.unit_hits` as the mapped arrays, no copy.
  - *Checked:* rows gathered from the columns equal the parquet's; fixture profiles (sparse + dense tier, 20 draws) identical with the column files and with the parquet fallback; the sourmash import (string `name`) through the scan.
  - *Result, 5×10⁷ units (full-build schema, 16 × u8 `pin_hist`; 1.8 GB parquet, 2.8 GB of `.npy`), 32 K hit units:* the unit columns the query uses plus the hit units' rows in 0.08 s against 1.6–2.1 s for reading the parquet and filtering, with no anonymous copy of the table. Peak RSS still reads ~2.1 GB (vs 3.5–5.0 GB): 32 K random rows touch most pages of each column at this size; those pages are file-backed and reclaimable, and at full scale (10× the units, same hits) a far smaller share. Writing the columns took 3.1 s. At full scale ~32 GB of `.npy` beside the parquet (60 B per unit at the current dtypes); the `load` stage's 43.5 GB should go.
  - *Baseline, 10 k pairs with Q3 and the unit table still read whole (`full-build`, scratch copy, preload):* 15.1 s wall, 53.9 GB peak RSS (was 127 s and 98.6 GB before Q3): the `keys` stage is gone. `load` 10.3 s (80 s CPU) and 43.5 GB, 81% of peak: the unit table this step removes. `lookup` 0.39 s CPU-bound for 239 K sampled k-mers (~1.6 µs each), +9.8 GB of touched tier-2 pages (file-backed). `hit_units` 1.27 s (filtering the full table; ~0.1 s expected from the columns). `fit_zip` 2.1 s is the largest stage after `load` (14% of wall, ~40% once `load` goes): `em_pin` iterates all 17.9 K detected units jointly until the slowest converges, with an `np.add.at` per iteration. Every other stage under 0.3 s. 32 K hit units, 17.9 K detected, largest component 58 units.
  - *Fixed, first run against `full-build` (no `.npy` files yet):* `hit_units` panicked with Polars' maximum length: the parquet fallback read `pin_hist` whole (16 × ~5.3×10⁸ values > 2³²), as `unit-columns` would have. Columns are now read in slices of 2²⁶ rows (`COLUMN_SLICE`), into a preallocated array (fallback) or straight into the `.npy` via `open_memmap` (`unit-columns`, so one slice in memory). Checked: slices of 3 rows on the fixture index give the same arrays both ways; 5×10⁷ units in 12 slices, identical `pin_hist`, 9.4 s against 3.9 s unsliced (~8 slices at full scale).
  - *Next:* `unit-columns` on `full-build`, then re-run the ladder. If `fit_zip` stays a large share: `bincount` over `col * 16 + level` instead of `np.add.at` (same output), or Q2's fit-only-the-shipped-estimator.

* **Phase 6, step 15 — 100 k pairs with memory-mapped unit columns; `em_pin` without `np.add.at`.**
  - *Result, 100 k pairs against `full-build` (scratch copy, preload):* 58.5 s wall, 65.7 GB peak RSS. `load` 0.02 s and 0.19 GB (was 10.3 s and 43.5 GB). `lookup` 1.46 s for 2.4 M sampled k-mers (~0.6 µs each), RSS to 36.6 GB: nearly all of tier 2's mapped pages touched (file-backed). `hit_units` 14.7 s (82.6 s CPU) and +29 GB: the unit columns read whole from the parquet (~60 B × 5.3×10⁸ units), i.e. the queried index had no `units.*.npy`, so the parquet fallback ran; `unit-columns` on the index removes it. `fit_zip` 30.6 s (52% of wall), then `fit_zi` 3.1, `baselines` 2.7, `fit_zib` 2.2, `gather` 1.3, `fit_em` 1.1 s. 292 K hit units, 164 K detected, largest component 71 units, 527 K (unit, hash) pairs.
  - *`fit_zip`:* `em_pin` runs all detected units jointly until the slowest converges, so it hits `max_iter` (1000); `np.add.at` into the (unit, level) table was 37% of it. Replaced by a `bincount` over `unit * 16 + level`: identical output; on a synthetic table at this size (164 K units, 330 K pairs) 18.1 → 12.4 s. Fewer iterations is not a fix: at 50/100/200/1000 iterations, 18.7 K/9.0 K/4.7 K/530 units are > 0.1% off the 5000-iteration coverage.
  - *Next:* `unit-columns` on `full-build`, re-run. Then `em_pin` per component (components are independent; the largest is 71 units): iterate only components not yet converged, which leaves the per-iteration cost with the slow units only. Stopping becomes per component, so estimates change within `tol`; check against the joint fit on the fixtures and the fmh benchmark. `fit_zi`/`em` share the pattern.

## Libraries

Most of the plumbing exists; the amino-acid k-mer hashing and the translation LUT are small enough to write yourself.

| Need | Crate | Note |
| --- | --- | --- |
| FASTA/FASTQ parsing | `paraseq` or `needletail` | paraseq has parallel paired-end record sets; needletail is the mature default. |
| Compression | `niffler` + `flate2` (zlib-rs backend), `zstd` | Transparent gz/zst input. |
| DNA minimizers/hashing | [`simd-minimizers`](https://github.com/rust-seq/simd-minimizers), `packed-seq`, `seq-hash` | Needs AVX2 or NEON; also accepts ASCII text, so could sample amino-acid strings. |
| General hashing | `xxhash-rust` or `wyhash` | For u64-packed amino-acid k-mers. |
| sourmash compatibility | `sourmash` (Rust core) | Protein/Dayhoff/HP encodings and its hash, for the compat mode. |
| Parallelism | `rayon`, `crossbeam-channel` | Reader → worker → reducer pipeline. |
| Index storage | `memmap2`, `rkyv` or `bincode`; `ptr_hash` or sorted arrays | Zero-copy, memory-mapped index. |
| Stats | `statrs` | Poisson/binomial for coverage correction. |
| CLI, errors, logs | `clap`, `anyhow`/`thiserror`, `tracing`, `indicatif` | Standard. |
| Testing | `proptest`, `insta`, `assert_cmd`, `cargo-nextest` | Invariants, golden outputs, CLI tests. |
| Benchmarks | `criterion` or `divan` | Micro-benchmarks of hashing and translation. |
| Frame calling (optional) | [FragGeneScanRs](https://github.com/unipept/FragGeneScanRs) | GPL-3.0: call as a subprocess or keep your crate GPL. |

For binaries using the SIMD crates, pin `target-cpu=x86-64-v3` in release builds; Groot Koerkamp's [post on distributing SIMD binaries](https://curiouscoding.nl/posts/distributing-rust-simd-binaries/) covers bioconda and the `ensure_simd` check.

**Bridge and Python side (phases 0–5)**

| Need | Library | Note |
| --- | --- | --- |
| Rust ↔ Python | `pyo3`, `maturin`, `numpy` (the crate `rust-numpy`) | Zero-copy numpy outputs from batch kernels; release the GIL inside kernels. |
| Environment and builds | `uv` with the maturin build backend | `uv sync` rebuilds the extension; `maturin develop --release` for fast iteration. |
| MGnify subsets | `duckdb` over MGnify's Parquet files | Remote querying of MGnify Proteins via Parquet/DuckDB is documented by MGnify. |
| Tables and storage | `polars`, `pyarrow`, `numpy` | Parquet for index tables during prototyping. |
| Statistics | `scipy` (optimize, stats), `numpy` | EM, zero-inflated negative binomial, dispersion tests. |
| sourmash compat | `sourmash` Python package | Direct comparison with fmh-funprofiler sketches. |
| CLI | `typer` | Prototype CLI mirroring the planned Rust one. |
| Testing | `pytest`, `hypothesis`, `pytest-benchmark` | Property tests mirror the Rust `proptest` ones. |
| Lint, format, types | `ruff`, `mypy` | Run in pre-commit and CI. |
| Plots, notebooks | `matplotlib`, Jupyter or marimo | Evaluation only; kept out of the package. |

## Benchmark labels

Pfam is the primary label for benchmarking; KO is kept as a secondary label.

- **Why Pfam.** The released index labels MGnify90 clusters with MGnify's own Pfam hits (no KO or eggNOG annotation run), so Pfam is what users will report and what defaults should be chosen on. Pfam is also open, so truth sets and indexes can be shared.
- **Why keep KO.** fmh-funprofiler profiles KEGG orthologs and is the direct predecessor, so a KO benchmark is the only like-for-like comparison with it (the phase-3 gate was passed on KO). KO results are reported beside Pfam results, but decisions and gates from phase 5 on use Pfam.
- **Pfam truth.** Annotate the benchmark genomes' proteins with Pfam-A (`hmmsearch --cut_ga`, the Pfam release of MGnify's `mgy_proteins_pfam`). A Pfam is present if a read overlaps one of its domain hits (the domain's coordinates, not the whole gene), and its depth is Σ over those domains of aligned bases / domain length, as for KOs. Multi-domain proteins count towards each of their Pfams.
- **Our profiles.** Two index types: Pfam units built from the benchmark's reference proteins (the analogue of the KO units, isolating the method from the reference), and MGnify90 cluster indexes where each Pfam's abundance is Σ over clusters carrying it of `coverage_zi` (the shipped configuration).
- **Other tools on Pfam.** DIAMOND: a read counts towards the Pfams whose domains its alignment overlaps on the subject. HUMAnN: UniRef90 families regrouped to Pfam (`humann_regroup_table -g uniref90_pfam`). fmh-funprofiler and kMermaid stay on KO only: fmh-funprofiler ships KO sketches, and kMermaid's disjoint clusters do not fit multi-domain proteins.

## Test and benchmark data

Three tiers: tiny hand-built fixtures in the repo, a small simulated community for CI, and large benchmarks run outside CI.

| Tier | Data | Tests | Where |
| --- | --- | --- | --- |
| Unit fixtures (< 1 MB) | \~20 proteins with their CDS; hand-built reads: in-frame, each of 6 frames, reverse strand, spanning a stop, with synonymous and non-synonymous changes, with N, with an indel. | Translation, frame filter, silent-change invariance, threshold nesting. | `tests/data/`, generated by a script in `scripts/` so they are reproducible. |
| CI community (\~5–10 MB) | 5 complete genomes spanning GC: *S. aureus* (low GC), *E. coli* K-12, *B. subtilis* 168, *P. aeruginosa* PAO1 (high GC), one archaeon. Reads from InSilicoSeq with known origins; truth = per-gene coverage from read origins mapped through annotations. | End-to-end golden output (insta snapshot), GC effect on frame filter, abundance error bounds. | Repo via Git LFS, or Zenodo with a checksum-verified download step. |
| fmh-funprofiler datasets | Their simulated metagenomes and KO ground truth, plus Pfam truth from the same genomes. | Head-to-head at equal density; compat-mode parity (KO); primary accuracy (Pfam). | Zenodo record from their repo. |
| Divergence ladder | Hold out the source genomes' proteins from the index, keeping relatives at \~95/90/80/70% amino-acid identity. | Sensitivity decay vs alphabet, k, spaced seeds. The key experiment for any amino-acid k-mer method. | Built by the eval pipeline. |
| CAMI II (marine, strain madness, plant-associated) | Public simulated metagenomes with genome-level truth. | Functional truth derived by annotating the source genomes with Pfam (primary) and KofamScan KOs (for fmh-funprofiler) and projecting read origins. | CAMI data portal. |
| Mock communities | ZymoBIOMICS standards (e.g. the gut standard) with published reference genomes. | Real sequencing error and library bias with known composition. | ENA/SRA runs. |
| Real cohort | A subset of HMP2/IBDMDB with published HUMAnN outputs. | Concordance with HUMAnN; runtime at scale. | IBDMDB portal. |
| Negative controls | Shuffled reads, human reads, intergenic-only simulated reads. | False-positive rate. | Generated. |

**Reference databases:** develop on a human-gut biome subset of MGnify Proteins 2026\_07; release on the full set (CC0, so pre-built indexes can be shared). Ground-truth genomes for the CI community and CAMI must be annotated against the same unit definitions, e.g. by mapping their predicted proteins to MGnify clusters, and with the same Pfam release as MGnify's labels.

## Repo setup

One GitHub repo holding a uv-managed Python package built with maturin and a Cargo workspace; pre-commit runs fast checks for both languages, and CI tests both plus Rust-vs-Python parity.

**Create**

Done in phase 0 (repo `timrozday-mgnify/kmer-functional-profiler`, package `kmer_functional_profiler`, GPL-3.0-or-later). `uv init --build-backend maturin` was not used: it puts a single crate at the root, which conflicts with the workspace layout below. Instead `pyproject.toml` points maturin at `crates/py` (`manifest-path`, `python-source = "python"`, `module-name = "kmer_functional_profiler._core"`). Dependencies are added in the phase that first uses them (polars, pyarrow, duckdb in phase 2; scipy in phase 4; sourmash in phase 3; pytest-benchmark and rust-numpy in phase 1).

**Layout**

```text
pyproject.toml        # maturin backend; module <name>._core
uv.lock
Cargo.toml            # workspace, shared deps, lints, release profile
rust-toolchain.toml   # pinned stable
.cargo/config.toml    # target-cpu=x86-64-v3 for release
crates/core/          # kernels: fastx, translation, k-mers, hashing, samplers, lookup
crates/py/            # PyO3 bindings over core -> <name>._core
crates/cli/           # phase 8: final Rust binary (index, query, model)
python/<name>/        # prototype: index build, query, model, typer CLI
python/<name>/reference/  # slow pure-Python twins of the kernels (test oracles)
tests/python/  tests/data/  tests/golden/   # golden outputs from the prototype
scripts/              # fixture generation, data download
eval/                 # Nextflow benchmarks, notebooks
docs/spec.md          # algorithm spec, written at end of phase 7
```

Keeping `crates/core` free of PyO3 means phases 6 and 8 reuse the kernels unchanged; only `crates/py` knows about Python.

**Pre-commit** (the `pre-commit` framework, or its Rust port `prek`):

- `pre-commit-hooks`: trailing-whitespace, end-of-file-fixer, check-yaml, check-toml, check-added-large-files (limit 1 MB, so data goes to LFS/Zenodo).
- Python: `ruff check --fix` and `ruff format`; `mypy` on `python/` at the pre-push stage.
- Rust: `cargo fmt --all -- --check`; `cargo clippy --workspace --all-targets -- -D warnings` at the pre-push stage.
- `uv lock --check` so the lockfile never drifts from `pyproject.toml`.
- `typos` for spelling, `taplo` for TOML formatting.

**CI (GitHub Actions)**

- `ci.yml` on push and pull request, three jobs:
  - **rust**: fmt, clippy with `-D warnings`, `cargo nextest`, `cargo doc`, MSRV build, `cargo deny`.
  - **python**: `astral-sh/setup-uv`, `uv sync` (builds the extension in release mode), ruff, mypy, pytest with hypothesis, on Python 3.11–3.13.
  - **parity**: Rust kernels vs `reference/` twins on fixtures and the CI community; from phase 8, the Rust CLI vs `tests/golden/`.
  - Matrix: ubuntu-latest and macos-latest (arm64, covers NEON). Coverage from `pytest-cov` and `cargo-llvm-cov` to Codecov.
- Use `dtolnay/rust-toolchain`, `Swatinem/rust-cache`, and `taiki-e/install-action` for nextest, cargo-deny and llvm-cov.
- `bench.yml`: manual or weekly `pytest-benchmark` and criterion runs; not a merge gate.
- `release.yml` (phase 9): `cargo-dist` binaries on tags, optional wheels via `PyO3/maturin-action`; bioconda recipe after the first tagged release.
- Dependabot for Cargo, uv and Actions; branch protection on `main` requiring `ci.yml`.

## Extending the reference (desirable)

Users should be able to add their own proteins to a released index: typically MGnify90 plus a study's MAG protein predictions plus the genes and functional annotations of its assemblies. This is desirable, not essential; it is phase 10, after the release, and nothing earlier depends on it. It matters because a study's own genomes are the closest references its reads will ever have, and much of a new study's protein content is not yet in MGnify.

**Approach: overlays, not rebuilds.** A full build is an HPC job over the whole release, so extension adds a small overlay index next to the base and never rewrites it.

1. *Assign added proteins to units.* Proteins that match an existing MGnify90 cluster at ≥ 90% identity join it; the rest are clustered among themselves at 90% (MMseqs2 linclust, as MGnify does) into new units. Candidate clusters come cheaply from the index itself (dense-tier containment), then an alignment against their representatives confirms them, so no search of all 1.66×10⁹ representatives is needed.
2. *Build the overlay with the same rules.* New units get *t\_g* from their own *n\_kmers*, exactly as in the build; existing units that gain members keep their *t\_g* and add the new members' k-mers under it, so nested thresholds stay consistent and the query needs no change in sampling. The overlay has the base's layout (tier 1, tier 2, dense) and unit ids that continue the base's.
3. *Global quantities.* *n\_groups* (score and promiscuity cut) is the base value plus the overlay's, found by probing the base; components are merged at query time with a union-find over units linked across base and overlay. *p\_in* of existing units is not recomputed (added members are few against the cluster), and this approximation is recorded in the overlay's metadata.
4. *Labels.* Added units carry the study's own annotations (e.g. KO, eggNOG or Pfam from its assembly pipeline) beside the MGnify Pfam labels; existing units that gain members keep their labels and record the added ones.
5. *Query and compaction.* The query probes the base and each overlay and runs detection and EM over the union. `compact` folds overlays into a new base when they grow large or many (a partitioned rebuild of only the affected components).

**Constraints on earlier phases** (keep extension possible; do not build for it):

- Keep per-unit build rules local to the unit (*t\_g*, floor, *p\_in*), with only *n\_groups* and components global, so an overlay can compute its part alone.
- Keep unit ids appendable and the unit table separate from the hash tables; record the base release and parameters in `meta.json` so an overlay can check it matches.
- Let the phase-6 Rust lookup probe more than one table and merge their hits.

**Evaluation.** Reuse the divergence ladder: hold genomes out of the index, add their proteins back as an overlay, and compare against a rebuild that includes them and against the base alone. Also check that adding unrelated proteins changes no base unit's profile.

## Risks and open questions

The biggest risk is that the gain over fmh-funprofiler with a lower scaled value is too small to justify a new tool.

- **MGnify scale.** A floor on all 1.66×10^9 clusters would need \~160 GB of index; on the 0.45×10^9 non-singletons, \~45 GB. Phase 6 (see Progress log): the Python prototype cannot build the full release on one node, so full-scale storage and compute are measured and tuned there, before the phase-7 ablations.
- **Component size.** Promiscuous k-mers can chain clusters into one giant component, which makes the EM serial. Measure component sizes on the development subset in phase 2 and tune the N-clusters cut-off.

* **Marginal novelty.** Setting sourmash to scaled = 100 may recover most of the completeness gap at modest cost. Run that baseline in phase 3 before building phase 4.
* **Query density set by the smallest units.** If *n\_min* forces *t\_max* near 1, the query does almost no sparsification. Measure the distribution of *t\_g* on the real database early.
* **Divergence.** Exact amino-acid k-mers miss distant homologs regardless of sampling. If the ladder shows steep decay below 80% identity, reduced alphabets or spaced seeds become mandatory, not optional.
* **Shared k-mers and hierarchy.** EM at protein-cluster level, then aggregation to function, is likely better than EM directly on functions. Untested.
* **Normalisation.** Which single-copy marker set, and whether to report per-genome copies by default.
* **Frame filter at high GC.** Keeps \~3 frames at 70% GC; acceptable, but check false positives there specifically.
* **Reference extension (phase 10).** Overlay approximations (*p\_in* of existing units not updated, *n\_groups* summed across tables) may bias scores for clusters that gain many members; `compact` bounds the drift. Assigning added proteins to MGnify90 clusters by alignment is the costliest part of `extend`.
* **Future work:** 30% families by mapping MGnify90 representatives onto the 128.7 M MGnify30-C2 representatives, as a coarser level for floors, EM partitions and annotation.
* **Open:** whether KO/eggNOG labels are worth the annotation run for users, or Pfam suffices. Benchmarking uses Pfam (see Benchmark labels), with KO only for the fmh-funprofiler comparison.

## Sources

- [Metagenomic functional profiling: to sketch or not to sketch? (Bioinformatics 2024)](https://academic.oup.com/bioinformatics/article/40/Supplement_2/ii165/7749078)
- [fmh-funprofiler repository](https://github.com/KoslickiLab/fmh-funprofiler)
- [fmh-funprofiler preprint v2 (completeness discussion)](https://www.biorxiv.org/content/10.1101/2023.11.06.565843v2.full.pdf)
- [kMermaid (PLoS Comp Biol)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12507277/)
- [UProC (Bioinformatics 2015)](https://academic.oup.com/bioinformatics/article/31/9/1382/200454)
- [Orpheum repository](https://github.com/czbiohub-sf/orpheum)
- [Protein k-mers enable assembly-free microbial metapangenomics](https://dib-lab.github.io/2021-paper-metapangenomes/)
- [Metabuli repository](https://github.com/steineggerlab/Metabuli)
- [Metabuli update: spaced metamers, reduced alphabets, syncmers (bioRxiv 2026)](https://www.biorxiv.org/content/10.64898/2026.03.13.711249v2)
- [simd-minimizers](https://github.com/rust-seq/simd-minimizers) and [docs](https://docs.rs/simd-minimizers/latest/simd_minimizers/)
- [Distributing Rust SIMD binaries](https://curiouscoding.nl/posts/distributing-rust-simd-binaries/)
- [FragGeneScanRs (BMC Bioinformatics 2022)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9148508/)

* [MGnify Proteins resource (2026\_07)](https://docs.mgnify.org/src/docs/mgnify-proteins.html)
* [MGnify: the microbiome sequence data analysis resource in 2023 (NAR)](https://academic.oup.com/nar/article/51/D1/D753/6880769)
* [sylph (Nature Biotechnology 2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12339375/)
* [sylph releases (two-stage .syl2db)](https://github.com/bluenote-1577/sylph/releases)
