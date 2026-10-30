//! Tree-index kernels (docs/tree-index-plan.md, phase 0): candidate pairs by shared stored
//! k-mers, Kruskal spanning forests and their merge trees, and the postings each diff scheme
//! writes on a merge forest.
//!
//! Sets are in CSR form: set `i` is `hashes[offsets[i]..offsets[i + 1]]`, sorted and distinct.

use crate::{Error, threads};

/// No node: the parent of a top, or an unreachable node's range.
pub const NONE: u32 = u32::MAX;

fn check_sets(offsets: &[u64], hashes: &[u64]) -> Result<usize, Error> {
    let ok = !offsets.is_empty()
        && offsets[0] == 0
        && offsets.windows(2).all(|w| w[0] <= w[1])
        && offsets[offsets.len() - 1] == hashes.len() as u64;
    if !ok {
        return Err(Error::BadOffsets);
    }
    Ok(offsets.len() - 1)
}

fn set<'a>(offsets: &[u64], hashes: &'a [u64], i: usize) -> &'a [u64] {
    &hashes[offsets[i] as usize..offsets[i + 1] as usize]
}

/// Candidate pairs `a < b` and their edge weights.
#[derive(Debug, Default, PartialEq, Eq)]
pub struct Pairs {
    pub a: Vec<u32>,
    pub b: Vec<u32>,
    /// Stored k-mers both hold.
    pub shared: Vec<u32>,
    /// |S\_a Δ S\_b| over hashes both could hold (≤ the smaller of their thresholds).
    pub weight: Vec<u32>,
    /// Counter increments of the join: its work, Σ over k-mers of holders².
    pub probes: u64,
    /// K-mers skipped for more than `max_holders` holders.
    pub skipped: u64,
}

/// Every pair of sets sharing at least `min_shared` stored k-mers, from an inverted index.
///
/// With `storage`, only pairs whose edge can be in the storage MST are kept: on the cycle
/// a–b–(empty root) an edge heavier than both root edges (|S\_a|, |S\_b|) never is. K-mers held
/// by more than `max_holders` sets are skipped (exact when `usize::MAX`).
pub fn overlap_pairs(
    offsets: &[u64],
    hashes: &[u64],
    max_hash: &[u64],
    min_shared: u32,
    storage: bool,
    max_holders: usize,
) -> Result<Pairs, Error> {
    let n = check_sets(offsets, hashes)?;
    if max_hash.len() != n {
        return Err(Error::LengthMismatch);
    }
    let size = u32::try_from(n).map_err(|_| Error::NodeOutOfRange(n as u64))?;
    let mut inverted: Vec<(u64, u32)> = Vec::with_capacity(hashes.len());
    for i in 0..size {
        inverted.extend(set(offsets, hashes, i as usize).iter().map(|&h| (h, i)));
    }
    inverted.sort_unstable();
    let posting = |h: u64| {
        let lo = inverted.partition_point(|&(x, _)| x < h);
        let hi = inverted.partition_point(|&(x, _)| x <= h);
        &inverted[lo..hi]
    };
    let below = |i: usize, m: u64| set(offsets, hashes, i).partition_point(|&h| h <= m) as u32;
    let n_threads = threads().min(n.max(1));
    let parts: Vec<Pairs> = std::thread::scope(|scope| {
        let workers: Vec<_> = (0..n_threads)
            .map(|t| {
                scope.spawn(move || {
                    let mut out = Pairs::default();
                    let mut count = vec![0u32; n];
                    let mut touched: Vec<u32> = Vec::new();
                    for a in (t..n).step_by(n_threads) {
                        for &h in set(offsets, hashes, a) {
                            let list = posting(h);
                            if list.len() > max_holders {
                                out.skipped += 1;
                                continue;
                            }
                            let after = list.partition_point(|&(_, b)| b as usize <= a);
                            for &(_, b) in &list[after..] {
                                if count[b as usize] == 0 {
                                    touched.push(b);
                                }
                                count[b as usize] += 1;
                            }
                            out.probes += (list.len() - after) as u64;
                        }
                        let size_a = (offsets[a + 1] - offsets[a]) as u32;
                        for &b in &touched {
                            let (bi, shared) = (b as usize, count[b as usize]);
                            count[bi] = 0;
                            if shared < min_shared {
                                continue;
                            }
                            let m = max_hash[a].min(max_hash[bi]);
                            let weight = below(a, m) + below(bi, m) - 2 * shared;
                            let size_b = (offsets[bi + 1] - offsets[bi]) as u32;
                            if storage && weight > size_a.max(size_b) {
                                continue;
                            }
                            out.a.push(a as u32);
                            out.b.push(b);
                            out.shared.push(shared);
                            out.weight.push(weight);
                        }
                        touched.clear();
                    }
                    out
                })
            })
            .collect();
        workers.into_iter().map(|w| w.join().unwrap()).collect()
    });
    let mut out = Pairs::default();
    for p in parts {
        out.a.extend(p.a);
        out.b.extend(p.b);
        out.shared.extend(p.shared);
        out.weight.extend(p.weight);
        out.probes += p.probes;
        out.skipped += p.skipped;
    }
    Ok(out)
}

/// Shared stored k-mers and weight (as in [`overlap_pairs`]) of given pairs (`a[i]`, `b[i]`).
pub fn pair_weights(
    offsets: &[u64],
    hashes: &[u64],
    max_hash: &[u64],
    a: &[u32],
    b: &[u32],
) -> Result<(Vec<u32>, Vec<u32>), Error> {
    let n = check_sets(offsets, hashes)?;
    if max_hash.len() != n || a.len() != b.len() {
        return Err(Error::LengthMismatch);
    }
    if let Some(&x) = a.iter().chain(b).find(|&&x| x as usize >= n) {
        return Err(Error::NodeOutOfRange(u64::from(x)));
    }
    let mut shared = vec![0u32; a.len()];
    let mut weight = vec![0u32; a.len()];
    let per_thread = a.len().div_ceil(threads()).max(1);
    std::thread::scope(|scope| {
        for (t, (sh, we)) in shared
            .chunks_mut(per_thread)
            .zip(weight.chunks_mut(per_thread))
            .enumerate()
        {
            scope.spawn(move || {
                for k in 0..sh.len() {
                    let (i, j) = (
                        a[t * per_thread + k] as usize,
                        b[t * per_thread + k] as usize,
                    );
                    let m = max_hash[i].min(max_hash[j]);
                    let (x, y) = (set(offsets, hashes, i), set(offsets, hashes, j));
                    let (x, y) = (
                        &x[..x.partition_point(|&h| h <= m)],
                        &y[..y.partition_point(|&h| h <= m)],
                    );
                    let (mut p, mut q, mut c) = (0, 0, 0u32);
                    while p < x.len() && q < y.len() {
                        match x[p].cmp(&y[q]) {
                            std::cmp::Ordering::Less => p += 1,
                            std::cmp::Ordering::Greater => q += 1,
                            std::cmp::Ordering::Equal => {
                                c += 1;
                                p += 1;
                                q += 1;
                            }
                        }
                    }
                    sh[k] = c;
                    we[k] = (x.len() + y.len()) as u32 - 2 * c;
                }
            });
        }
    });
    Ok((shared, weight))
}

fn find(parent: &mut [u32], mut u: u32) -> u32 {
    while parent[u as usize] != u {
        let up = parent[parent[u as usize] as usize];
        parent[u as usize] = up; // path halving
        u = up;
    }
    u
}

/// A Kruskal spanning forest and its single-linkage merge forest.
///
/// Leaves are nodes `0..n`; internal node `n + i` merges `left[i]` and `right[i]` at
/// `level[i]`, so a parent's id is above its children's.
#[derive(Debug, Default, PartialEq, Eq)]
pub struct Kruskal {
    /// Per input edge: in the spanning forest.
    pub mst: Vec<bool>,
    pub left: Vec<u32>,
    pub right: Vec<u32>,
    pub level: Vec<u32>,
    /// Roots of the merge forest, in order of their smallest leaf.
    pub tops: Vec<u32>,
}

/// Minimum spanning forest of nodes `0..n` under edges (`u`, `v`, `w`), lightest first and
/// ties by node pair (so input order does not matter), with its merge forest.
pub fn kruskal(n: usize, u: &[u32], v: &[u32], w: &[u32]) -> Result<Kruskal, Error> {
    if u.len() != v.len() || u.len() != w.len() {
        return Err(Error::LengthMismatch);
    }
    if let Some(&x) = u.iter().chain(v).find(|&&x| x as usize >= n) {
        return Err(Error::NodeOutOfRange(u64::from(x)));
    }
    let size = u32::try_from(n).map_err(|_| Error::NodeOutOfRange(n as u64))?;
    let mut by_weight: Vec<u32> = (0..u.len() as u32).collect();
    by_weight.sort_by_key(|&e| {
        let e = e as usize;
        (w[e], u[e].min(v[e]), u[e].max(v[e]))
    });
    let mut parent: Vec<u32> = (0..size).collect();
    let mut node: Vec<u32> = (0..size).collect(); // merge-tree node of each union-find root
    let mut out = Kruskal {
        mst: vec![false; u.len()],
        ..Kruskal::default()
    };
    for e in by_weight {
        let e = e as usize;
        let (a, b) = (find(&mut parent, u[e]), find(&mut parent, v[e]));
        if a == b {
            continue;
        }
        out.mst[e] = true;
        let id = size + out.left.len() as u32;
        out.left.push(node[a as usize]);
        out.right.push(node[b as usize]);
        out.level.push(w[e]);
        let root = a.min(b);
        parent[a.max(b) as usize] = root;
        node[root as usize] = id;
    }
    for x in 0..size {
        if find(&mut parent, x) == x {
            out.tops.push(node[x as usize]);
        }
    }
    Ok(out)
}

/// Leaf order of a merge forest: depth first from each top, left child first.
#[derive(Debug, Default, PartialEq, Eq)]
pub struct Layout {
    /// The leaf at each position.
    pub order: Vec<u32>,
    /// Per node, its leaf range `lo..hi` in positions (`NONE` if no top reaches it).
    pub lo: Vec<u32>,
    pub hi: Vec<u32>,
    pub parent: Vec<u32>,
}

fn check_forest(n: usize, left: &[u32], right: &[u32], tops: &[u32]) -> Result<(), Error> {
    if left.len() != right.len() {
        return Err(Error::LengthMismatch);
    }
    let total = n + left.len();
    for (i, (&l, &r)) in left.iter().zip(right).enumerate() {
        // children are created before their parent
        if l as usize >= n + i || r as usize >= n + i {
            return Err(Error::NodeOutOfRange(u64::from(l.max(r))));
        }
    }
    match tops.iter().find(|&&t| t as usize >= total) {
        Some(&t) => Err(Error::NodeOutOfRange(u64::from(t))),
        None => Ok(()),
    }
}

pub fn layout(n: usize, left: &[u32], right: &[u32], tops: &[u32]) -> Result<Layout, Error> {
    check_forest(n, left, right, tops)?;
    let total = n + left.len();
    let mut out = Layout {
        order: Vec::new(),
        lo: vec![NONE; total],
        hi: vec![NONE; total],
        parent: vec![NONE; total],
    };
    let mut seen = vec![false; total];
    let mut stack: Vec<u32> = Vec::new();
    for &top in tops {
        stack.push(top);
        while let Some(x) = stack.pop() {
            if std::mem::replace(&mut seen[x as usize], true) {
                return Err(Error::NodeOutOfRange(u64::from(x))); // a node reached twice
            }
            if (x as usize) < n {
                let p = out.order.len() as u32;
                out.order.push(x);
                (out.lo[x as usize], out.hi[x as usize]) = (p, p + 1);
            } else {
                let i = x as usize - n;
                out.parent[left[i] as usize] = x;
                out.parent[right[i] as usize] = x;
                stack.push(right[i]);
                stack.push(left[i]);
            }
        }
    }
    for i in 0..left.len() {
        let x = n + i;
        if seen[x] {
            out.lo[x] = out.lo[left[i] as usize];
            out.hi[x] = out.hi[right[i] as usize];
        }
    }
    Ok(out)
}

/// Postings each diff scheme writes for each k-mer on a merge forest.
///
/// Each top hangs below a fixed empty root, so a k-mer present at a top costs a gain there.
#[derive(Debug, Default, PartialEq, Eq)]
pub struct Schemes {
    pub hash: Vec<u64>,
    pub holders: Vec<u32>,
    /// Fitch parsimony: the minimum number of state changes (the design's postings).
    pub fitch: Vec<u32>,
    /// Of those, gains (absent to present, root to leaves) in one most-parsimonious
    /// assignment: Fitch's top-down pass, keeping the parent's state where it can.
    pub gains: Vec<u32>,
    /// Majority: a node is present when most of its leaves that could hold the k-mer do.
    pub majority: Vec<u32>,
    /// Add-only (ProPhyle): one posting per largest subtree whose leaves all hold it.
    pub add_only: Vec<u32>,
    /// With `emit`: Fitch's changes in that assignment, one row each (the k-mer's hash, the
    /// node, gain or loss). A leaf holds a k-mer iff the nearest change on its root path is a
    /// gain, so as leaf ranges, gains minus losses covering a leaf is 1 or 0.
    pub change_hash: Vec<u64>,
    pub change_node: Vec<u32>,
    pub change_gain: Vec<bool>,
}

// Fitch state sets: bit 0 "absent", bit 1 "present"; 0 is missing data ("?").
const ABSENT: u8 = 1;
const PRESENT: u8 = 2;

fn combine(x: u8, y: u8) -> (u8, u32) {
    if x == 0 {
        return (y, 0);
    }
    if y == 0 {
        return (x, 0);
    }
    match x & y {
        0 => (x | y, 1),
        both => (both, 0),
    }
}

/// The state Fitch's top-down pass gives a node with set `s` under a parent in `state`.
fn pick(state: u8, s: u8) -> u8 {
    if s == 0 || s & state != 0 { state } else { s }
}

/// The forest's lookups, shared by every k-mer.
struct Tree<'a> {
    n: usize,
    left: &'a [u32],
    right: &'a [u32],
    layout: Layout,
    /// Per node, the top above it.
    top: Vec<u32>,
    /// Range maximum over `between` (the lowest common ancestor of positions p and p + 1,
    /// `NONE` across tops): level j holds the maximum of 2^j entries.
    sparse: Vec<Vec<u32>>,
    /// Per level l, leaves at positions `..p` with level ≥ l: `nm[l * (n_pos + 1) + p]`.
    nm: Vec<u32>,
    n_pos: usize,
}

impl Tree<'_> {
    fn children(&self, x: u32) -> (u32, u32) {
        let i = x as usize - self.n;
        (self.left[i], self.right[i])
    }

    fn sibling(&self, x: u32) -> u32 {
        let (l, r) = self.children(self.layout.parent[x as usize]);
        if l == x { r } else { l }
    }

    fn range(&self, x: u32) -> (u32, u32) {
        (self.layout.lo[x as usize], self.layout.hi[x as usize])
    }

    /// Leaves below `x` that could hold a k-mer of level `l` (not missing data).
    fn nonmissing(&self, x: u32, l: usize) -> u32 {
        let (lo, hi) = self.range(x);
        let row = &self.nm[l * (self.n_pos + 1)..];
        row[hi as usize] - row[lo as usize]
    }

    /// Lowest common ancestor of the leaves at positions `i < j`, or `NONE` across tops.
    fn lca(&self, i: u32, j: u32) -> u32 {
        let (lo, hi) = (i as usize, j as usize); // between[lo..hi]
        let k = (hi - lo).ilog2() as usize;
        self.sparse[k][lo].max(self.sparse[k][hi - (1 << k)])
    }

    fn is_ancestor(&self, a: u32, b: u32) -> bool {
        let (a_lo, a_hi) = self.range(a);
        let (b_lo, b_hi) = self.range(b);
        a_lo <= b_lo && b_hi <= a_hi
    }
}

/// One virtual-tree node of a k-mer: a holder leaf or a lowest common ancestor of holders.
#[derive(Clone, Copy)]
struct Virtual {
    node: u32,
    vparent: u32, // index into the k-mer's virtual nodes, NONE for a root
    holders: u32,
    nm: u32,
    fitch: u8,
    state: u8,
    /// This node's chain: `chain[start..end]` holds (node, set) for the nodes above it up to
    /// the virtual parent's child on its side (or its top), lowest first, stopping where the
    /// set becomes "absent" (the rest of the chain is absent too).
    start: u32,
    end: u32,
    /// The set at the top of the chain (the virtual parent's child, or the top).
    head: u8,
    /// That node: the chain's top.
    stop: u32,
}

/// Postings of one k-mer held at the sorted leaf positions `holders`, at level `l`; Fitch's
/// changes (node, gain) are appended to `emitted`.
fn kmer_postings(
    tree: &Tree<'_>,
    holders: &[u32],
    l: usize,
    emitted: &mut Vec<(u32, bool)>,
) -> [u32; 4] {
    let leaf = |p: u32| tree.layout.order[p as usize];
    let mut nodes: Vec<u32> = holders.iter().map(|&p| leaf(p)).collect();
    nodes.extend(
        holders
            .windows(2)
            .map(|w| tree.lca(w[0], w[1]))
            .filter(|&x| x != NONE),
    );
    nodes.sort_unstable();
    nodes.dedup();
    // virtual parents: preorder (lo ascending, ancestors first) with a stack of ancestors
    let mut pre: Vec<u32> = (0..nodes.len() as u32).collect();
    pre.sort_by_key(|&i| {
        let (lo, hi) = tree.range(nodes[i as usize]);
        (lo, std::cmp::Reverse(hi))
    });
    let mut vparent = vec![NONE; nodes.len()];
    let mut stack: Vec<u32> = Vec::new();
    for &i in &pre {
        while let Some(&s) = stack.last() {
            if tree.is_ancestor(nodes[s as usize], nodes[i as usize]) {
                break;
            }
            stack.pop();
        }
        vparent[i as usize] = stack.last().copied().unwrap_or(NONE);
        stack.push(i);
    }
    let count = |x: u32| {
        let (lo, hi) = tree.range(x);
        (holders.partition_point(|&p| p < hi) - holders.partition_point(|&p| p < lo)) as u32
    };
    let mut v: Vec<Virtual> = nodes
        .iter()
        .zip(&vparent)
        .map(|(&node, &vparent)| Virtual {
            node,
            vparent,
            holders: count(node),
            nm: tree.nonmissing(node, l),
            fitch: 0,
            state: 0,
            start: 0,
            end: 0,
            head: 0,
            stop: 0,
        })
        .collect();
    let mut kids: Vec<Vec<u32>> = vec![Vec::new(); v.len()];
    for (i, x) in v.iter().enumerate() {
        if x.vparent != NONE {
            kids[x.vparent as usize].push(i as u32);
        }
    }
    let off = |x: u32| if tree.nonmissing(x, l) > 0 { ABSENT } else { 0 };
    // chain from `v[i]` up to `stop`: its records and the set at `stop`; returns the cost
    let walk = |i: usize, stop: u32, v: &mut [Virtual], chain: &mut Vec<(u32, u8)>| {
        let (mut x, mut s, mut cost) = (v[i].node, v[i].fitch, 0);
        v[i].start = chain.len() as u32;
        while x != stop && s != ABSENT {
            let (next, c) = combine(s, off(tree.sibling(x)));
            x = tree.layout.parent[x as usize];
            s = next;
            cost += c;
            chain.push((x, s));
        }
        v[i].end = chain.len() as u32;
        v[i].head = s;
        v[i].stop = stop;
        cost
    };
    let mut chain: Vec<(u32, u8)> = Vec::new();
    let mut fitch = 0u32;
    // bottom-up (ids ascending: children first)
    for i in 0..v.len() {
        let x = v[i].node;
        if (x as usize) < tree.n {
            v[i].fitch = PRESENT;
        } else {
            let (cl, cr) = tree.children(x);
            let mut side = [off(cl), off(cr)];
            for &k in &kids[i] {
                let k = k as usize;
                let c = if tree.is_ancestor(cl, v[k].node) {
                    0
                } else {
                    1
                };
                fitch += walk(k, [cl, cr][c], &mut v, &mut chain);
                side[c] = v[k].head;
            }
            let (s, cost) = combine(side[0], side[1]);
            fitch += cost;
            v[i].fitch = s;
        }
        if v[i].vparent == NONE {
            fitch += walk(i, tree.top[x as usize], &mut v, &mut chain);
            fitch += combine(ABSENT, v[i].head).1;
        }
    }
    // top-down (ids descending: parents first), from an absent root above each top
    let (mut changes, mut gains) = (0u32, 0u32);
    let mut change = |from: u8, to: u8, node: u32| {
        if from != to {
            changes += 1;
            gains += u32::from(to == PRESENT);
            emitted.push((node, to == PRESENT));
        }
    };
    for i in (0..v.len()).rev() {
        let above = match v[i].vparent {
            NONE => ABSENT,
            p => v[p as usize].state,
        };
        let records = &chain[v[i].start as usize..v[i].end as usize];
        let mut state = above;
        if v[i].head == ABSENT {
            // the chain's top (or the part above its last record) is absent
            change(state, ABSENT, v[i].stop);
            state = ABSENT;
        }
        for (r, &(_, s)) in records.iter().enumerate().rev() {
            if r + 1 == records.len() && s == ABSENT {
                continue; // the head, counted above
            }
            let next = pick(state, s);
            change(state, next, records[r].0);
            state = next;
            let below = if r == 0 { v[i].node } else { records[r - 1].0 };
            if state == PRESENT && off(tree.sibling(below)) == ABSENT {
                change(PRESENT, ABSENT, tree.sibling(below));
            }
        }
        let next = pick(state, v[i].fitch);
        change(state, next, v[i].node);
        v[i].state = next;
        if next == PRESENT && (v[i].node as usize) >= tree.n {
            let (cl, cr) = tree.children(v[i].node);
            for c in [cl, cr] {
                let has_kid = kids[i]
                    .iter()
                    .any(|&k| tree.is_ancestor(c, v[k as usize].node));
                if !has_kid && off(c) == ABSENT {
                    change(PRESENT, ABSENT, c);
                }
            }
        }
    }
    debug_assert_eq!(changes, fitch);
    // majority: present where 2 × holders > leaves that could hold it
    let present = |holders: u32, nm: u32| 2 * holders > nm;
    let mut majority = 0u32;
    for i in 0..v.len() {
        let m = present(v[i].holders, v[i].nm);
        let mut x = v[i].node;
        let mut prev = m;
        let stop = match v[i].vparent {
            NONE => tree.top[x as usize],
            p => {
                let (cl, cr) = tree.children(v[p as usize].node);
                if tree.is_ancestor(cl, x) { cl } else { cr }
            }
        };
        while x != stop && prev {
            let s = tree.sibling(x);
            x = tree.layout.parent[x as usize];
            let now = present(v[i].holders, tree.nonmissing(x, l));
            majority += u32::from(now != prev);
            if now && tree.nonmissing(s, l) > 0 {
                majority += 1; // an absent sibling under a present node
            }
            prev = now;
        }
        let above = match v[i].vparent {
            NONE => false,
            p => present(v[p as usize].holders, v[p as usize].nm),
        };
        majority += u32::from(prev != above);
        if m && (v[i].node as usize) >= tree.n {
            let (cl, cr) = tree.children(v[i].node);
            for c in [cl, cr] {
                let has_kid = kids[i]
                    .iter()
                    .any(|&k| tree.is_ancestor(c, v[k as usize].node));
                if !has_kid && tree.nonmissing(c, l) > 0 {
                    majority += 1;
                }
            }
        }
    }
    // add-only: one posting per virtual node that all its could-hold leaves hold, unless its
    // virtual parent does too (the largest such subtree lies on its chain)
    let all = |x: &Virtual| x.holders == x.nm;
    let add_only = v
        .iter()
        .filter(|x| all(x) && (x.vparent == NONE || !all(&v[x.vparent as usize])))
        .count() as u32;
    [fitch, gains, majority, add_only]
}

/// Postings of every k-mer of the leaves' sets (`offsets`, `hashes`) on the merge forest
/// (`left`, `right`, `tops`) under each diff scheme.
///
/// Leaf `i` has level `leaf_level[i]`: it could hold a hash `h` iff
/// `h <= level_hash[leaf_level[i]]`; otherwise it is missing data for that k-mer.
/// `level_hash` is ascending. Leaves no top reaches are ignored. With `emit`, Fitch's changes
/// are returned too.
#[allow(clippy::too_many_arguments)]
pub fn schemes(
    n: usize,
    left: &[u32],
    right: &[u32],
    tops: &[u32],
    offsets: &[u64],
    hashes: &[u64],
    leaf_level: &[u8],
    level_hash: &[u64],
    emit: bool,
) -> Result<Schemes, Error> {
    if check_sets(offsets, hashes)? != n || leaf_level.len() != n || level_hash.is_empty() {
        return Err(Error::LengthMismatch);
    }
    if let Some(&l) = leaf_level.iter().find(|&&l| l as usize >= level_hash.len()) {
        return Err(Error::NodeOutOfRange(u64::from(l)));
    }
    let layout = layout(n, left, right, tops)?;
    let n_pos = layout.order.len();
    let total = n + left.len();
    let mut top = vec![NONE; total];
    for &t in tops {
        let (lo, hi) = (layout.lo[t as usize], layout.hi[t as usize]);
        for &leaf in &layout.order[lo as usize..hi as usize] {
            top[leaf as usize] = t;
        }
    }
    for i in 0..left.len() {
        top[n + i] = top[left[i] as usize];
    }
    let mut between = vec![NONE; n_pos.saturating_sub(1)];
    for i in 0..left.len() {
        let x = n + i;
        if layout.lo[x] != NONE {
            between[layout.hi[left[i] as usize] as usize - 1] = x as u32;
        }
    }
    let mut sparse = vec![between];
    while 1 << sparse.len() <= sparse[0].len() {
        let (prev, half) = (&sparse[sparse.len() - 1], 1 << (sparse.len() - 1));
        let next = (0..prev.len() - half)
            .map(|i| prev[i].max(prev[i + half]))
            .collect();
        sparse.push(next);
    }
    let levels = level_hash.len();
    let mut nm = vec![0u32; levels * (n_pos + 1)];
    for l in 0..levels {
        let row = &mut nm[l * (n_pos + 1)..(l + 1) * (n_pos + 1)];
        for (p, &leaf) in layout.order.iter().enumerate() {
            row[p + 1] = row[p] + u32::from(leaf_level[leaf as usize] as usize >= l);
        }
    }
    let mut inverted: Vec<(u64, u32)> = Vec::with_capacity(hashes.len());
    for (p, &leaf) in layout.order.iter().enumerate() {
        inverted.extend(
            set(offsets, hashes, leaf as usize)
                .iter()
                .map(|&h| (h, p as u32)),
        );
    }
    inverted.sort_unstable();
    let tree = Tree {
        n,
        left,
        right,
        layout,
        top,
        sparse,
        nm,
        n_pos,
    };
    // runs of one hash, split into about equal parts per thread
    let mut starts: Vec<usize> = (0..inverted.len())
        .filter(|&i| i == 0 || inverted[i].0 != inverted[i - 1].0)
        .collect();
    let n_kmers = starts.len();
    starts.push(inverted.len());
    let per_thread = n_kmers.div_ceil(threads()).max(1);
    let parts: Vec<Schemes> = std::thread::scope(|scope| {
        let workers: Vec<_> = (0..n_kmers)
            .step_by(per_thread)
            .map(|first| {
                let (tree, inverted, starts) = (&tree, &inverted, &starts);
                scope.spawn(move || {
                    let mut out = Schemes::default();
                    let mut holders: Vec<u32> = Vec::new();
                    let mut emitted: Vec<(u32, bool)> = Vec::new();
                    for k in first..(first + per_thread).min(n_kmers) {
                        let run = &inverted[starts[k]..starts[k + 1]];
                        let h = run[0].0;
                        holders.clear();
                        holders.extend(run.iter().map(|&(_, p)| p));
                        let l = level_hash.partition_point(|&x| x < h).min(levels - 1);
                        emitted.clear();
                        let [fitch, gains, majority, add_only] =
                            kmer_postings(tree, &holders, l, &mut emitted);
                        if emit {
                            for &(node, gain) in &emitted {
                                out.change_hash.push(h);
                                out.change_node.push(node);
                                out.change_gain.push(gain);
                            }
                        }
                        out.hash.push(h);
                        out.holders.push(holders.len() as u32);
                        out.fitch.push(fitch);
                        out.gains.push(gains);
                        out.majority.push(majority);
                        out.add_only.push(add_only);
                    }
                    out
                })
            })
            .collect();
        workers.into_iter().map(|w| w.join().unwrap()).collect()
    });
    let mut out = Schemes::default();
    for p in parts {
        out.hash.extend(p.hash);
        out.holders.extend(p.holders);
        out.fitch.extend(p.fitch);
        out.gains.extend(p.gains);
        out.majority.extend(p.majority);
        out.add_only.extend(p.add_only);
        out.change_hash.extend(p.change_hash);
        out.change_node.extend(p.change_node);
        out.change_gain.extend(p.change_gain);
    }
    Ok(out)
}

/// For each query range `q_lo[j]..q_hi[j]`, the summed weights of the ranges `lo[i]..hi[i]`
/// that contain it (`lo <= q_lo` and `q_hi <= hi`): a range posting's hits at a leaf (a
/// one-position range) or at a node (its leaf range). One sweep in `lo` order with a Fenwick
/// tree over the ranges' `hi`: O((P + Q) log P).
pub fn contained_sums(
    lo: &[u32],
    hi: &[u32],
    w: &[i64],
    q_lo: &[u32],
    q_hi: &[u32],
) -> Result<Vec<i64>, Error> {
    if lo.len() != hi.len() || lo.len() != w.len() || q_lo.len() != q_hi.len() {
        return Err(Error::LengthMismatch);
    }
    let mut his: Vec<u32> = hi.to_vec();
    his.sort_unstable();
    his.dedup();
    // Fenwick over his, indexed in descending hi order so a prefix is "hi >= x"
    let m = his.len();
    let mut fen = vec![0i64; m + 1];
    let slot = |h: u32| m - his.partition_point(|&x| x < h); // 1-based, descending
    let mut by_lo: Vec<u32> = (0..lo.len() as u32).collect();
    by_lo.sort_unstable_by_key(|&i| lo[i as usize]);
    let mut by_q: Vec<u32> = (0..q_lo.len() as u32).collect();
    by_q.sort_unstable_by_key(|&j| q_lo[j as usize]);
    let mut out = vec![0i64; q_lo.len()];
    let mut next = 0;
    for j in by_q {
        let j = j as usize;
        while next < by_lo.len() && lo[by_lo[next] as usize] <= q_lo[j] {
            let i = by_lo[next] as usize;
            let mut x = slot(hi[i]);
            while x <= m {
                fen[x] += w[i];
                x += x & x.wrapping_neg();
            }
            next += 1;
        }
        // ranges with hi >= q_hi: the first slots, down to the smallest such hi
        let mut x = m - his.partition_point(|&h| h < q_hi[j]);
        let mut sum = 0;
        while x > 0 {
            sum += fen[x];
            x &= x - 1;
        }
        out[j] = sum;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A small deterministic generator (tests only).
    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self, m: u64) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            (self.0 >> 33) % m
        }
    }

    /// Full-tree Fitch, majority and add-only on node `x` for a k-mer: (set, cost) and the
    /// subtree's (holders, could-hold leaves).
    fn brute(
        n: usize,
        left: &[u32],
        right: &[u32],
        held: &[bool],
        could: &[bool],
        x: u32,
    ) -> (u8, u32, u32, u32) {
        if (x as usize) < n {
            let s = if held[x as usize] {
                PRESENT
            } else if could[x as usize] {
                ABSENT
            } else {
                0
            };
            return (
                s,
                0,
                u32::from(held[x as usize]),
                u32::from(could[x as usize]),
            );
        }
        let i = x as usize - n;
        let a = brute(n, left, right, held, could, left[i]);
        let b = brute(n, left, right, held, could, right[i]);
        let (s, c) = combine(a.0, b.0);
        (s, a.1 + b.1 + c, a.2 + b.2, a.3 + b.3)
    }

    /// Majority and add-only postings below `x` given the parent's state.
    fn brute_other(
        n: usize,
        left: &[u32],
        right: &[u32],
        held: &[bool],
        could: &[bool],
        x: u32,
        above: (bool, bool),
    ) -> (u32, u32) {
        let (_, _, h, nm) = brute(n, left, right, held, could, x);
        if nm == 0 {
            return (0, 0); // all missing data: takes its parent's state
        }
        let m = 2 * h > nm;
        let all = h > 0 && h == nm;
        let mut out = (u32::from(m != above.0), u32::from(all && !above.1));
        if (x as usize) >= n {
            let i = x as usize - n;
            for c in [left[i], right[i]] {
                let r = brute_other(n, left, right, held, could, c, (m, all));
                out = (out.0 + r.0, out.1 + r.1);
            }
        }
        out
    }

    #[test]
    fn schemes_match_brute_force() {
        let mut rng = Lcg(7);
        for _ in 0..3000 {
            let n = 2 + rng.next(40) as usize;
            // random forest: merge random clusters, sometimes leave several tops
            let mut roots: Vec<u32> = (0..n as u32).collect();
            let (mut left, mut right) = (Vec::new(), Vec::new());
            let merges = n - 1 - rng.next(3).min(n as u64 - 1) as usize;
            for _ in 0..merges {
                let a = roots.swap_remove(rng.next(roots.len() as u64) as usize);
                let b = roots.swap_remove(rng.next(roots.len() as u64) as usize);
                left.push(a);
                right.push(b);
                roots.push((n + left.len() - 1) as u32);
            }
            let level_hash = [10u64, 20, 30];
            let leaf_level: Vec<u8> = (0..n).map(|_| rng.next(3) as u8).collect();
            let mut sets: Vec<Vec<u64>> = vec![Vec::new(); n];
            for h in [5u64, 15, 25, 7] {
                for (i, s) in sets.iter_mut().enumerate() {
                    if h <= level_hash[leaf_level[i] as usize] && rng.next(3) > 0 {
                        s.push(h);
                    }
                }
            }
            for s in &mut sets {
                s.sort_unstable();
            }
            let mut offsets = vec![0u64];
            let mut hashes = Vec::new();
            for s in &sets {
                hashes.extend(s);
                offsets.push(hashes.len() as u64);
            }
            let got = schemes(
                n,
                &left,
                &right,
                &roots,
                &offsets,
                &hashes,
                &leaf_level,
                &level_hash,
                true,
            )
            .unwrap();
            let lay = layout(n, &left, &right, &roots).unwrap();
            for (k, &h) in got.hash.iter().enumerate() {
                let held: Vec<bool> = sets.iter().map(|s| s.contains(&h)).collect();
                let could: Vec<bool> = leaf_level
                    .iter()
                    .map(|&l| h <= level_hash[l as usize])
                    .collect();
                let (mut fitch, mut majority, mut add_only) = (0, 0, 0);
                for &t in &roots {
                    let (s, c, _, _) = brute(n, &left, &right, &held, &could, t);
                    fitch += c + combine(ABSENT, s).1;
                    let (m, a) = brute_other(n, &left, &right, &held, &could, t, (false, false));
                    majority += m;
                    add_only += a;
                }
                assert_eq!(got.fitch[k], fitch, "fitch, hash {h}");
                assert_eq!(got.majority[k], majority, "majority, hash {h}");
                assert_eq!(got.add_only[k], add_only, "add-only, hash {h}");
                assert!(got.gains[k] >= 1 && got.gains[k] <= got.fitch[k]);
                // the emitted changes: as many as Fitch's, and they rebuild every leaf that
                // could hold the k-mer (nearest change above it, or absent)
                let mine: Vec<(u32, bool)> = (0..got.change_hash.len())
                    .filter(|&c| got.change_hash[c] == h)
                    .map(|c| (got.change_node[c], got.change_gain[c]))
                    .collect();
                assert_eq!(mine.len() as u32, got.fitch[k], "changes, hash {h}");
                for leaf in (0..n).filter(|&i| could[i]) {
                    let mut x = leaf as u32;
                    let mut state = false;
                    loop {
                        if let Some(&(_, g)) = mine.iter().find(|&&(node, _)| node == x) {
                            state = g;
                            break;
                        }
                        match lay.parent[x as usize] {
                            NONE => break,
                            p => x = p,
                        }
                    }
                    assert_eq!(state, held[leaf], "leaf {leaf}, hash {h}");
                }
            }
        }
    }

    #[test]
    fn contained_sums_match_brute_force() {
        let mut rng = Lcg(11);
        for _ in 0..500 {
            let (p, q) = (rng.next(30) as usize, rng.next(30) as usize);
            let span = |rng: &mut Lcg| {
                let a = rng.next(20) as u32;
                (a, a + 1 + rng.next(8) as u32)
            };
            let (mut lo, mut hi, mut w) = (Vec::new(), Vec::new(), Vec::new());
            for _ in 0..p {
                let (a, b) = span(&mut rng);
                lo.push(a);
                hi.push(b);
                w.push(rng.next(7) as i64 - 3);
            }
            let (mut q_lo, mut q_hi) = (Vec::new(), Vec::new());
            for _ in 0..q {
                let (a, b) = span(&mut rng);
                q_lo.push(a);
                q_hi.push(b);
            }
            let got = contained_sums(&lo, &hi, &w, &q_lo, &q_hi).unwrap();
            for j in 0..q {
                let want: i64 = (0..p)
                    .filter(|&i| lo[i] <= q_lo[j] && q_hi[j] <= hi[i])
                    .map(|i| w[i])
                    .sum();
                assert_eq!(got[j], want);
            }
        }
    }

    #[test]
    fn kruskal_builds_the_merge_tree() {
        // 0-1 (1), 1-2 (3), 0-2 (2): MST {0-1, 0-2}; 3 alone
        let got = kruskal(4, &[0, 1, 0], &[1, 2, 2], &[1, 3, 2]).unwrap();
        assert_eq!(got.mst, [true, false, true]);
        assert_eq!(
            (got.left, got.right, got.level),
            (vec![0, 4], vec![1, 2], vec![1, 2])
        );
        assert_eq!(got.tops, [5, 3]);
        let lay = layout(4, &[0, 4], &[1, 2], &[5, 3]).unwrap();
        assert_eq!(lay.order, [0, 1, 2, 3]);
        assert_eq!((lay.lo[5], lay.hi[5], lay.lo[4], lay.hi[4]), (0, 3, 0, 2));
    }

    #[test]
    fn overlap_pairs_weights() {
        // sets: {1,2,3,4}, {1,2,3,5}, {9}, {1,2,3,4,6,7,8,10} ; thresholds all high
        let offsets = [0u64, 4, 8, 9, 17];
        let hashes = [1, 2, 3, 4, 1, 2, 3, 5, 9, 1, 2, 3, 4, 6, 7, 8, 10];
        let mh = [u64::MAX; 4];
        let got = overlap_pairs(&offsets, &hashes, &mh, 1, true, usize::MAX).unwrap();
        let mut rows: Vec<_> = (0..got.a.len())
            .map(|i| (got.a[i], got.b[i], got.shared[i], got.weight[i]))
            .collect();
        rows.sort_unstable();
        // 0-1: Δ 2; 0-3: Δ 4 ≤ 8; 1-3: shared 3, Δ 1 + 5 = 6 ≤ 8
        assert_eq!(rows, [(0, 1, 3, 2), (0, 3, 4, 4), (1, 3, 3, 6)]);
        let (sh, w) = pair_weights(&offsets, &hashes, &mh, &[0, 1, 2], &[1, 3, 3]).unwrap();
        assert_eq!((sh, w), (vec![3, 3, 0], vec![2, 6, 9]));
        // a lower threshold on set 3 restricts its side of the difference
        let mh = [u64::MAX, u64::MAX, u64::MAX, 4];
        let got = overlap_pairs(&offsets, &hashes, &mh, 1, true, usize::MAX).unwrap();
        let w: Vec<_> = (0..got.a.len())
            .map(|i| (got.a[i], got.b[i], got.weight[i]))
            .collect();
        assert!(w.contains(&(0, 3, 0)) && w.contains(&(1, 3, 1)));
    }
}
