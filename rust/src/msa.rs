//! SeqAn 2.4's `globalMsaAlignment(align, SimpleScore(5, -3, -1, -3))` exactly as
//! Porechop_ABI's `msa_consensus` uses it, followed by its profile consensus.
//!
//! The pipeline (graph_align_tcoffee_msa.h):
//!   1. a T-Coffee library of segment matches from every pair: a global Gotoh alignment
//!      (which also gives the pair distance) and local alignments (four Smith-Waterman
//!      "clumped" alignments, or one DP local alignment for 50+ sequences),
//!   2. every match rescored as sum(10 + score) and refined into an alignment graph,
//!   3. a neighbour-joining guide tree, triplet library extension (restricted to subtree
//!      groups for 10+ sequences),
//!   4. progressive alignment with heaviest common subsequences along the guide tree,
//!   5. the alignment graph converted to a gapped matrix (components in topological order).
//!
//! For 50+ sequences SeqAn takes the banded code path with an uninitialised
//! `MsaOptions::bandWidth`. The Bioconda (Linux, gcc) build behaves exactly like the unbanded
//! path there, which is what is implemented (the macOS clang build differs).

use crate::align::{traceback_with, Op, DIAG, HOR, HOR_OPEN, MAX_H, MAX_V, VER, VER_OPEN};
use std::collections::{BTreeSet, HashMap};

const MATCH: i32 = 5;
const MISMATCH: i32 = -3;
const GAP_OPEN: i32 = -3;
const GAP_EXT: i32 = -1;
const NEG: i32 = i32::MIN / 4;
const GAP: u8 = 4;

/// A gapless segment match between two sequences (SeqAn `Fragment<>`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Frag {
    pub s1: usize,
    pub b1: usize,
    pub s2: usize,
    pub b2: usize,
    pub len: usize,
}

#[inline(always)]
fn sub(a: u8, b: u8) -> i32 {
    if a == b { MATCH } else { MISMATCH }
}

/// SeqAn's `_maxScore`: with `CompleteTrace` ties keep both traces, with `SingleTrace` the left one.
#[inline(always)]
fn max_trace(left: i32, right: i32, tl: u8, tr: u8, single: bool) -> (i32, u8) {
    if left == right {
        (left, if single { tl } else { tl | tr })
    } else if left < right {
        (right, tr)
    } else {
        (left, tl)
    }
}

/// SeqAn's DP framework with affine gaps, complete trace and gaps placed left: global without
/// free end gaps (`AlignConfig<>`), or local. Returns the diagonal segments of the traceback.
fn dp_align(h: &[u8], v: &[u8], local: bool, s1: usize, s2: usize, out: &mut Vec<Frag>) {
    let n = h.len();
    let m = v.len();
    let rows = m + 1;
    let mut trace = vec![0u8; (n + 1) * rows];
    let mut s = vec![0i32; rows];
    let mut hh = vec![NEG; rows];
    if !local {
        let mut vert = NEG;
        for i in 1..=m {
            let (vv, tv) = max_trace(vert + GAP_EXT, s[i - 1] + GAP_OPEN, VER, VER_OPEN, false);
            vert = vv;
            s[i] = vv;
            trace[i] = tv | MAX_V;
        }
    }
    let (mut best, mut best_col, mut best_row) = (0i32, 0usize, 0usize);
    let mut best_vh = (NEG, NEG);
    for j in 1..=n {
        let col = &mut trace[j * rows..(j + 1) * rows];
        let mut diag = s[0];
        if local {
            s[0] = 0;
            hh[0] = NEG;
            col[0] = 0;
        } else {
            let (hv, th) = max_trace(hh[0] + GAP_EXT, s[0] + GAP_OPEN, HOR, HOR_OPEN, false);
            hh[0] = hv;
            s[0] = hv;
            col[0] = th | MAX_H;
        }
        let mut vert = NEG;
        let hc = h[j - 1];
        for i in 1..=m {
            let inter = diag + sub(hc, v[i - 1]);
            diag = s[i];
            let (hv, th) = max_trace(hh[i] + GAP_EXT, s[i] + GAP_OPEN, HOR, HOR_OPEN, local);
            hh[i] = hv;
            let (vv, tv) = max_trace(vert + GAP_EXT, s[i - 1] + GAP_OPEN, VER, VER_OPEN, local);
            vert = vv;
            let tmp = th | tv;
            let (g, t2) = max_trace(vv, hv, MAX_V, MAX_H, local);
            let (mut sv, mut t) = max_trace(inter, g, DIAG | tmp, t2 | tmp, local);
            if local && sv <= 0 {
                // _maxScore(score, NONE, score, NONE, trace): a tie keeps NONE
                sv = 0;
                t = 0;
            }
            s[i] = sv;
            col[i] = t;
            if local && sv > best {
                best = sv;
                best_col = j;
                best_row = i;
                best_vh = (vv, hv);
            }
        }
    }
    if !local {
        best_col = n;
        best_row = m;
    } else if best > 0 {
        // _correctTraceValue (single trace): start in a gap matrix when it holds the maximum
        let t = &mut trace[best_col * rows + best_row];
        if best_vh.0 == best {
            *t = (*t & !DIAG) | MAX_V;
        } else if best_vh.1 == best {
            *t = (*t & !DIAG) | MAX_H;
        }
    }
    let (ops, _, _) = traceback_with(|c, r| trace[c * rows + r], 0, n, m, best_col, best_row, true);
    let (mut i, mut j) = (0usize, 0usize);
    let mut run: Option<(usize, usize, usize)> = None;
    for op in ops {
        match op {
            Op::Diag => {
                run = Some(match run {
                    Some((a, b, l)) => (a, b, l + 1),
                    None => (i, j, 1),
                });
                i += 1;
                j += 1;
            }
            Op::Hor | Op::Vert => {
                if let Some((a, b, l)) = run.take() {
                    out.push(Frag { s1, b1: a, s2, b2: b, len: l });
                }
                if op == Op::Hor {
                    i += 1;
                } else {
                    j += 1;
                }
            }
        }
    }
    if let Some((a, b, l)) = run {
        out.push(Frag { s1, b1: a, s2, b2: b, len: l });
    }
}

// Smith-Waterman trace directions (graph_align_tcoffee_library.h).
const SW_DIAG: u8 = 0;
const SW_HOR: u8 = 1;
const SW_VER: u8 = 2;
const SW_STOP: u8 = 3;

/// One `_localAlignment(..., SmithWatermanClump)`: best local alignment avoiding forbidden
/// cells, whose cells then become forbidden. Returns the score.
fn sw_once(a: &[u8], b: &[u8], forbidden: &mut [bool], s1: usize, s2: usize, out: &mut Vec<Frag>) -> i32 {
    let len1 = a.len();
    let len2 = b.len();
    // Per cell: bit 3 vertical-open, bit 2 horizontal-open, bits 0-1 direction.
    let mut trace = vec![0u8; len1 * len2];
    let mut mat = vec![0i32; len2 + 1];
    let mut horizontal = vec![GAP_OPEN - GAP_EXT; len2 + 1];
    let mut score_max = 0i32;
    let (mut best_row, mut best_col) = (0usize, 0usize);
    let mut initial_dir = SW_STOP;
    let mut tv_mat = 0u8;
    for col in 1..=len1 {
        let mut diag_val = mat[0];
        mat[0] = 0;
        let mut vert = GAP_OPEN - GAP_EXT;
        for row in 1..=len2 {
            let cell = (col - 1) * len2 + (row - 1);
            let max_val;
            if forbidden[cell] {
                trace[cell] = SW_STOP;
                max_val = 0;
                vert = 0;
                horizontal[row] = 0;
            } else {
                let mut t = 0u8;
                let x = mat[row - 1] + GAP_OPEN;
                let y = vert + GAP_EXT;
                if x > y {
                    vert = x;
                    t |= 8;
                } else {
                    vert = y;
                }
                let x = mat[row] + GAP_OPEN;
                let y = horizontal[row] + GAP_EXT;
                if x > y {
                    horizontal[row] = x;
                    t |= 4;
                } else {
                    horizontal[row] = y;
                }
                let mut mv = diag_val + sub(a[col - 1], b[row - 1]);
                tv_mat = SW_DIAG;
                if vert > mv {
                    mv = vert;
                    tv_mat = SW_VER;
                }
                if horizontal[row] > mv {
                    mv = horizontal[row];
                    tv_mat = SW_HOR;
                }
                if 0 >= mv {
                    mv = 0;
                    tv_mat = SW_STOP;
                }
                trace[cell] = t | tv_mat;
                max_val = mv;
            }
            diag_val = mat[row];
            mat[row] = max_val;
            if max_val > score_max {
                best_row = row;
                best_col = col;
                score_max = max_val;
                initial_dir = tv_mat;
            }
        }
    }
    // _alignSmithWatermanTrace
    let (mut l1, mut l2) = (best_col, best_row);
    if l1 == 0 || l2 == 0 {
        return score_max;
    }
    let at = |c: usize, r: usize| trace[(c - 1) * len2 + (r - 1)];
    let mut print = |p1: usize, p2: usize, seg: usize, tv: u8| {
        if seg != 0 && tv == SW_DIAG {
            out.push(Frag { s1, b1: p1, s2, b2: p2, len: seg });
        }
    };
    let mut next = at(l1, l2);
    let mut tv = SW_DIAG;
    if initial_dir == SW_DIAG {
        tv = next & 3;
    } else if initial_dir == SW_HOR {
        if (next >> 2) & 1 != 0 {
            l1 -= 1;
        } else {
            tv = SW_HOR;
        }
    } else if initial_dir == SW_VER {
        if (next >> 3) & 1 != 0 {
            l2 -= 1;
        } else {
            tv = SW_VER;
        }
    }
    let mut seg = 0usize;
    let mut tv_old = tv;
    loop {
        next = at(l1, l2);
        if next & 3 == SW_STOP {
            break;
        }
        forbidden[(l1 - 1) * len2 + (l2 - 1)] = true;
        if tv == SW_DIAG {
            tv = next & 3;
        } else if tv == SW_HOR {
            tv = if (next >> 2) & 1 != 0 { SW_DIAG } else { SW_HOR };
        } else if tv == SW_VER {
            tv = if (next >> 3) & 1 != 0 { SW_DIAG } else { SW_VER };
        }
        if tv == SW_DIAG {
            if tv != tv_old {
                if tv_old == SW_VER {
                    l2 -= 1;
                } else {
                    l1 -= 1;
                }
                seg += 1;
                print(l1, l2, seg, tv_old);
                tv_old = tv;
                seg = 0;
            } else {
                seg += 1;
                l1 -= 1;
                l2 -= 1;
            }
        } else if tv == SW_HOR {
            if tv != tv_old {
                print(l1, l2, seg, tv_old);
                if (next >> 2) & 1 != 0 {
                    l1 -= 1;
                    tv = SW_DIAG;
                    seg = 0;
                } else {
                    tv_old = tv;
                    seg = 1;
                    l1 -= 1;
                }
            } else {
                seg += 1;
                l1 -= 1;
            }
        } else if tv == SW_VER {
            if tv != tv_old {
                print(l1, l2, seg, tv_old);
                if (next >> 3) & 1 != 0 {
                    l2 -= 1;
                    tv = SW_DIAG;
                    seg = 0;
                } else {
                    tv_old = tv;
                    seg = 1;
                    l2 -= 1;
                }
            } else {
                seg += 1;
                l2 -= 1;
            }
        }
        if l1 == 0 || l2 == 0 {
            break;
        }
    }
    if seg != 0 {
        print(l1, l2, seg, tv_old);
    }
    score_max
}

/// `_localAlignment(str, matches, scores, sc, 4, SmithWatermanClump())`.
fn sw_clump(a: &[u8], b: &[u8], s1: usize, s2: usize, out: &mut Vec<Frag>) {
    let mut forbidden = vec![false; a.len() * b.len()];
    let mut last = 0i32;
    for _ in 0..4 {
        let from = out.len();
        let score = sw_once(a, b, &mut forbidden, s1, s2, out);
        if 2 * score < last {
            out.truncate(from);
            break;
        }
        last = score;
    }
}

/// The library fragments of one pair, as `appendSegmentMatches` produces them.
pub(crate) fn pair_fragments(a: &[u8], b: &[u8], s1: usize, s2: usize, deep: bool) -> (Vec<Frag>, Vec<Frag>) {
    let mut global = Vec::new();
    dp_align(a, b, false, s1, s2, &mut global);
    let mut local = Vec::new();
    if deep {
        dp_align(a, b, true, s1, s2, &mut local);
    } else {
        sw_clump(a, b, s1, s2, &mut local);
    }
    (global, local)
}

/// Pair distance from the global alignment fragments (`_setDistanceValue`).
fn pair_distance(a: &[u8], b: &[u8], frags: &[Frag]) -> f64 {
    let mut matches = 0.0f64;
    let mut aligned = 0usize;
    for f in frags {
        aligned += f.len;
        for k in 0..f.len {
            if a[f.b1 + k] == b[f.b2 + k] {
                matches += 1.0;
            }
        }
    }
    let align_len = (aligned + (a.len() - aligned) + (b.len() - aligned)) as f64;
    1.0 - matches / align_len
}

/// The refined alignment graph: per sequence a partition into segments, and weighted edges.
struct AlignGraph {
    /// vertex -> (sequence, begin, length)
    vertices: Vec<(usize, usize, usize)>,
    /// per sequence: segment begins (ascending) and their vertices
    seg_begin: Vec<Vec<usize>>,
    seg_vertex: Vec<Vec<usize>>,
    /// adjacency: (neighbour, weight)
    adj: Vec<Vec<(usize, u64)>>,
}

impl AlignGraph {
    fn find_vertex(&self, seq: usize, pos: usize) -> usize {
        let b = &self.seg_begin[seq];
        let k = b.partition_point(|&x| x <= pos) - 1;
        self.seg_vertex[seq][k]
    }
}

/// `buildAlignmentGraph(matches, scores, g, score, ReScore())`: exact match refinement, then
/// every edge weighs 1 plus, for each match through it, (segment length * match score) / match
/// length, the match score being sum(10 + score) over its aligned pairs.
fn build_graph(seqs: &[Vec<u8>], frags: &[Frag]) -> (AlignGraph, HashMap<(usize, usize), u64>) {
    let nseq = seqs.len();
    // fragment sides per sequence: (fragment, side)
    let mut sides: Vec<Vec<(usize, usize)>> = vec![Vec::new(); nseq];
    for (fi, f) in frags.iter().enumerate() {
        sides[f.s1].push((fi, 0));
        sides[f.s2].push((fi, 1));
    }
    let side = |f: &Frag, sd: usize| if sd == 0 { (f.s1, f.b1, f.s2, f.b2) } else { (f.s2, f.b2, f.s1, f.b1) };
    let mut cuts: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); nseq];
    let mut stack: Vec<(usize, usize)> = Vec::new();
    for f in frags {
        for sd in 0..2 {
            let (s, b, _, _) = side(f, sd);
            for p in [b, b + f.len] {
                if cuts[s].insert(p) {
                    stack.push((s, p));
                }
            }
        }
    }
    while let Some((s, p)) = stack.pop() {
        for &(fi, sd) in &sides[s] {
            let f = &frags[fi];
            let (_, b, o, ob) = side(f, sd);
            if b < p && p < b + f.len {
                let q = ob + (p - b);
                if cuts[o].insert(q) {
                    stack.push((o, q));
                }
            }
        }
    }
    let mut g = AlignGraph { vertices: Vec::new(), seg_begin: vec![Vec::new(); nseq], seg_vertex: vec![Vec::new(); nseq], adj: Vec::new() };
    for s in 0..nseq {
        let len = seqs[s].len();
        let mut bounds: Vec<usize> = cuts[s].iter().copied().filter(|&p| p < len).collect();
        if bounds.first() != Some(&0) {
            bounds.insert(0, 0);
        }
        for (k, &b) in bounds.iter().enumerate() {
            let e = if k + 1 < bounds.len() { bounds[k + 1] } else { len };
            g.seg_begin[s].push(b);
            g.seg_vertex[s].push(g.vertices.len());
            g.vertices.push((s, b, e - b));
        }
    }
    let mut edges: HashMap<(usize, usize), u64> = HashMap::new();
    for f in frags {
        let (a, b) = (&seqs[f.s1], &seqs[f.s2]);
        let score: u64 = (0..f.len).map(|k| (10 + sub(a[f.b1 + k], b[f.b2 + k])) as u64).sum();
        let (mut p1, mut p2) = (f.b1, f.b2);
        let end = f.b1 + f.len;
        while p1 < end {
            let v1 = g.find_vertex(f.s1, p1);
            let v2 = g.find_vertex(f.s2, p2);
            let vlen = g.vertices[v1].2 as u64;
            let key = (v1.min(v2), v1.max(v2));
            *edges.entry(key).or_insert(1) += vlen * score / f.len as u64;
            p1 += vlen as usize;
            p2 += vlen as usize;
        }
    }
    (g, edges)
}

/// Guide tree: leaves 0..n are the sequences; `children` in SeqAn adjacency order.
struct Tree {
    children: Vec<Vec<usize>>,
    parent: Vec<usize>,
    root: usize,
}

impl Tree {
    fn add_vertex(&mut self) -> usize {
        self.children.push(Vec::new());
        self.parent.push(usize::MAX);
        self.children.len() - 1
    }
    /// SeqAn trees prepend new children.
    fn add_edge(&mut self, parent: usize, child: usize) {
        self.children[parent].insert(0, child);
        self.parent[child] = parent;
    }
    fn bfs(&self) -> Vec<usize> {
        let mut order = vec![self.root];
        let mut k = 0;
        while k < order.len() {
            let v = order[k];
            order.extend_from_slice(&self.children[v]);
            k += 1;
        }
        order
    }
    fn collect_leaves(&self, v: usize, out: &mut Vec<usize>) {
        if self.children[v].is_empty() {
            out.push(v);
        } else {
            for &c in &self.children[v] {
                self.collect_leaves(c, out);
            }
        }
    }
}

/// `njTree` (topology only: branch lengths are not used downstream).
fn nj_tree(mat_in: &[f64], nseq: usize) -> Tree {
    let mut t = Tree { children: Vec::new(), parent: Vec::new(), root: 0 };
    if nseq == 1 {
        t.root = t.add_vertex();
        return t;
    }
    if nseq == 2 {
        let v1 = t.add_vertex();
        let v2 = t.add_vertex();
        let iv = t.add_vertex();
        t.add_edge(iv, v1);
        t.add_edge(iv, v2);
        t.root = iv;
        return t;
    }
    let mut norm = 0.0f64;
    for &x in mat_in {
        norm = norm.max(x);
    }
    if norm == 0.0 {
        norm = 1.0;
    }
    let mut mat: Vec<i64> = mat_in.iter().map(|&x| (10000000.0 * (x / norm)) as i64).collect();
    const NIL: usize = usize::MAX;
    let mut connector: Vec<usize> = (0..nseq).collect();
    for i in 0..nseq {
        t.add_vertex();
        mat[i * nseq + i] = 0;
    }
    let mut fnseqs = nseq as i64;
    for _ in 0..nseq - 3 {
        let mut sum_of_branches = 0i64;
        for col in 1..nseq {
            for row in 0..col {
                mat[col * nseq + row] = mat[row * nseq + col];
                sum_of_branches += mat[row * nseq + col];
            }
        }
        let mut not_found = true;
        let mut tmin = 0i64;
        let (mut mini, mut minj) = (0usize, 0usize);
        for col in 1..nseq {
            if connector[col] == NIL {
                continue;
            }
            for row in 0..col {
                if connector[row] == NIL {
                    continue;
                }
                let mut di = 0i64;
                let mut dj = 0i64;
                for i in 0..nseq {
                    di += mat[i * nseq + row];
                    dj += mat[i * nseq + col];
                }
                let mut total = di + dj + (fnseqs - 2) * mat[row * nseq + col] + 2 * (sum_of_branches - di - dj);
                total /= 2 * (fnseqs - 2);
                if not_found || total < tmin {
                    not_found = false;
                    tmin = total;
                    mini = row;
                    minj = col;
                }
            }
        }
        let iv = t.add_vertex();
        t.add_edge(iv, connector[mini]);
        t.add_edge(iv, connector[minj]);
        fnseqs -= 1;
        connector[minj] = NIL;
        connector[mini] = iv;
        for j in 0..nseq {
            if connector[j] != NIL {
                if mini < j {
                    mat[mini * nseq + j] = (mat[mini * nseq + j] + mat[minj * nseq + j]) / 2;
                }
                if mini > j {
                    mat[j * nseq + mini] = (mat[mini * nseq + j] + mat[minj * nseq + j]) / 2;
                }
            }
        }
        for j in 0..nseq {
            mat[j * nseq + minj] = 0;
            mat[minj * nseq + j] = 0;
        }
    }
    let l: Vec<usize> = (0..nseq).filter(|&i| connector[i] != NIL).collect();
    let iv = t.add_vertex();
    t.add_edge(iv, connector[l[0]]);
    t.add_edge(iv, connector[l[1]]);
    let root = t.add_vertex();
    t.add_edge(root, connector[l[2]]);
    t.add_edge(root, iv);
    t.root = root;
    t
}

/// `_subTreeSearch`: groups of at least `min_members` leaves, plus the leftovers attached to the
/// closest group. Returns the group label of every sequence.
fn subtree_groups(t: &Tree, nseq: usize, min_members: usize) -> Vec<usize> {
    let root = t.root;
    let mut vertices = t.bfs();
    vertices.reverse();
    let mut num_leaves = vec![0usize; t.children.len()];
    for &v in &vertices {
        num_leaves[v] = if t.children[v].is_empty() { 1 } else { t.children[v].iter().map(|&c| num_leaves[c]).sum() };
    }
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut group_root: Vec<usize> = Vec::new();
    for &v in &vertices {
        if num_leaves[v] >= min_members {
            let mut g = Vec::new();
            t.collect_leaves(v, &mut g);
            groups.push(g);
            group_root.push(v);
            num_leaves[v] = 0;
            if v != root {
                let mut p = t.parent[v];
                while p != root {
                    num_leaves[p] = 0;
                    p = t.parent[p];
                }
                num_leaves[p] = 0;
            }
        }
    }
    if groups.is_empty() {
        let mut g = Vec::new();
        t.collect_leaves(root, &mut g);
        groups.push(g);
        group_root.push(root);
    }
    let mut closest = vec![(0usize, 0usize); t.children.len()];
    for (i, &gr) in group_root.iter().enumerate() {
        let mut v = gr;
        let mut dist = 0;
        while v != root {
            dist += 1;
            v = t.parent[v];
            if closest[v].0 == 0 || closest[v].0 > dist {
                closest[v] = (dist, i);
            }
        }
    }
    let mut grouped = vec![false; nseq];
    for g in &groups {
        for &s in g {
            grouped[s] = true;
        }
    }
    for s in (0..nseq).filter(|&s| !grouped[s]) {
        let mut v = s;
        while v != root {
            v = t.parent[v];
            if closest[v].0 != 0 {
                groups[closest[v].1].push(s);
                break;
            }
        }
    }
    let mut labels = vec![0usize; nseq];
    for (i, g) in groups.iter().enumerate() {
        for &s in g {
            labels[s] = i;
        }
    }
    labels
}

/// `tripletLibraryExtension`: for every vertex and every pair of its neighbours in different
/// sequences, the pair's edge gains the smaller of the two weights (computed from the original
/// weights). With `labels`, only edges within one group take part. Returns the adjacency lists.
fn triplet_extension(g: &AlignGraph, edges: &HashMap<(usize, usize), u64>, labels: Option<&[usize]>) -> Vec<Vec<(usize, u64)>> {
    let nv = g.vertices.len();
    let seq_of = |v: usize| g.vertices[v].0;
    let mut all: Vec<Vec<(usize, u64)>> = vec![Vec::new(); nv];
    // Neighbours taking part in the extension: every group only extends within itself, and a
    // vertex's restricted neighbours all share its group, so one pass covers all groups.
    let mut nb: Vec<Vec<(usize, u64)>> = vec![Vec::new(); nv];
    for (&(a, b), &c) in edges {
        all[a].push((b, c));
        if labels.map(|l| l[seq_of(a)] == l[seq_of(b)]).unwrap_or(true) {
            nb[a].push((b, c));
            nb[b].push((a, c));
        }
    }
    // Each new weight of (u, w), u < w, accumulated in a dense row while u is fixed.
    let mut row = vec![0u64; nv];
    let mut touched: Vec<usize> = Vec::new();
    let mut adj: Vec<Vec<(usize, u64)>> = vec![Vec::new(); nv];
    for u in 0..nv {
        for &(w, c) in &all[u] {
            if row[w] == 0 {
                touched.push(w);
            }
            row[w] += c;
        }
        let su = seq_of(u);
        for &(v, cuv) in &nb[u] {
            for &(w, cvw) in &nb[v] {
                if w > u && seq_of(w) != su {
                    if row[w] == 0 {
                        touched.push(w);
                    }
                    row[w] += cuv.min(cvw);
                }
            }
        }
        for &w in &touched {
            adj[u].push((w, row[w]));
            adj[w].push((u, row[w]));
            row[w] = 0;
        }
        touched.clear();
    }
    adj
}

/// `heaviestIncreasingSubsequence`: returns the chosen positions, last first.
fn heaviest_increasing_subsequence(seq: &[usize], weights: &[u64]) -> Vec<usize> {
    use std::ops::Bound::{Excluded, Unbounded};
    let mut list: BTreeSet<(usize, u64, usize)> = BTreeSet::new();
    let mut pred: Vec<Option<usize>> = vec![None; seq.len()];
    for (pos, &val) in seq.iter().enumerate() {
        let mut w = weights[pos];
        if w == 0 {
            continue;
        }
        let a_k = list.range(..(val, 0, 0)).next_back().copied();
        let mut b_l = match a_k {
            None => list.iter().next().copied(),
            Some(k) => list.range((Excluded(k), Unbounded)).next().copied(),
        };
        if let Some(k) = a_k {
            w += k.1;
        }
        while let Some(b) = b_l {
            if w < b.1 {
                break;
            }
            b_l = list.range((Excluded(b), Unbounded)).next().copied();
            list.remove(&b);
        }
        if b_l.map(|b| val < b.0).unwrap_or(true) {
            list.insert((val, w, pos));
        }
        if let Some(k) = a_k {
            pred[pos] = Some(k.2);
        }
    }
    let mut positions = Vec::new();
    if let Some(&(_, _, mut v)) = list.iter().next_back() {
        loop {
            positions.push(v);
            match pred[v] {
                Some(p) => v = p,
                None => break,
            }
        }
    }
    positions
}

type Profile = Vec<Vec<usize>>;

/// `heaviestCommonSubsequence(g, str1, str2, align)`.
fn heaviest_common_subsequence(g: &AlignGraph, str1: &Profile, str2: &Profile) -> Profile {
    let m = str1.len();
    let n = str2.len();
    let mut map = vec![usize::MAX; g.vertices.len()];
    for (pos, col) in str1.iter().enumerate() {
        for &v in col {
            map[v] = pos;
        }
    }
    let mut occupied: Vec<usize> = Vec::new();
    for (p2, col) in str2.iter().enumerate() {
        for &v in col {
            for &(u, _) in &g.adj[v] {
                if map[u] != usize::MAX {
                    occupied.push(map[u] * n + (n - p2 - 1));
                }
            }
        }
    }
    occupied.sort_unstable();
    occupied.dedup();
    let slot_to_pos = occupied;
    let mut weights = vec![0u64; slot_to_pos.len()];
    for (p2, col) in str2.iter().enumerate() {
        for &v in col {
            for &(u, c) in &g.adj[v] {
                if map[u] != usize::MAX {
                    let key = map[u] * n + (n - p2 - 1);
                    weights[slot_to_pos.partition_point(|&x| x < key)] += c;
                }
            }
        }
    }
    let seq: Vec<usize> = slot_to_pos.iter().map(|&s| n - 1 - (s % n)).collect();
    let positions = heaviest_increasing_subsequence(&seq, &weights);

    // _heaviestCommonSubsequence
    let num_matches = positions.len();
    let align_len = num_matches + (n - num_matches) + (m - num_matches);
    let mut align: Profile = Vec::with_capacity(align_len);
    let frag_len = |col: &Vec<usize>| g.vertices[col[0]].2;
    let (mut pos1, mut pos2) = (0usize, 0usize);
    let mut p = positions.len() as isize - 1;
    while align.len() < align_len {
        let (mut i, mut j) = (m, n);
        if p >= 0 {
            let s = slot_to_pos[positions[p as usize]];
            i = s / n;
            j = n - 1 - (s % n);
        }
        let mut first_i = true;
        if i != pos1 && j != pos2 {
            if pos1 == 0 && pos2 == 0 {
                let len1: usize = (pos1..i).map(|k| frag_len(&str1[k])).sum();
                let len2: usize = (pos2..j).map(|k| frag_len(&str2[k])).sum();
                if len1 > len2 {
                    first_i = false;
                }
            } else if i == m && i == n {
                // SeqAn compares i with n here (not j)
                let len1: usize = (pos1..i).map(|k| frag_len(&str1[k])).sum();
                let len2: usize = (pos2..j).map(|k| frag_len(&str2[k])).sum();
                if len1 < len2 {
                    first_i = false;
                }
            }
        }
        if first_i {
            while i != pos1 {
                align.push(str1[pos1].clone());
                pos1 += 1;
            }
            while j != pos2 {
                align.push(str2[pos2].clone());
                pos2 += 1;
            }
        } else {
            while j != pos2 {
                align.push(str2[pos2].clone());
                pos2 += 1;
            }
            while i != pos1 {
                align.push(str1[pos1].clone());
                pos1 += 1;
            }
        }
        if p >= 0 {
            let mut col = str1[pos1].clone();
            col.extend_from_slice(&str2[pos2]);
            align.push(col);
            pos1 += 1;
            pos2 += 1;
            p -= 1;
        }
    }
    align
}

/// `progressiveAlignment` along the guide tree; returns the alignment columns (vertex lists).
fn progressive(g: &AlignGraph, t: &Tree) -> Profile {
    let mut order = t.bfs();
    order.reverse();
    let mut seg: Vec<Profile> = vec![Vec::new(); t.children.len()];
    for &v in &order {
        if t.children[v].is_empty() {
            seg[v] = g.seg_vertex[v].iter().map(|&x| vec![x]).collect();
        } else {
            let (c1, c2) = (t.children[v][0], t.children[v][1]);
            let a = std::mem::take(&mut seg[c1]);
            let b = std::mem::take(&mut seg[c2]);
            seg[v] = heaviest_common_subsequence(g, &a, &b);
        }
    }
    std::mem::take(&mut seg[t.root])
}

/// `_createAlignmentGraph` + `convertAlignment(graph, mat)`: rows of base codes or GAP.
fn alignment_matrix(g: &AlignGraph, seqs: &[Vec<u8>], cols: &Profile) -> Vec<Vec<u8>> {
    let nseq = seqs.len();
    let ncomp = cols.len();
    // every column is one connected component, numbered in column order
    let mut comp_of = vec![usize::MAX; g.vertices.len()];
    let mut comp_len = vec![0usize; ncomp];
    for (c, col) in cols.iter().enumerate() {
        for &v in col {
            comp_of[v] = c;
            comp_len[c] = comp_len[c].max(g.vertices[v].2);
        }
    }
    // component graph: consecutive components along each sequence; SeqAn prepends adjacency
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); ncomp];
    for s in 0..nseq {
        let cs: Vec<usize> = g.seg_vertex[s].iter().map(|&v| comp_of[v]).collect();
        for w in cs.windows(2) {
            adj[w[0]].insert(0, w[1]);
        }
    }
    // topologicalSort: iterative DFS, then decreasing finishing time
    let mut visited = vec![false; ncomp];
    let mut finish = vec![0usize; ncomp];
    let mut time = 0usize;
    for r in 0..ncomp {
        if visited[r] {
            continue;
        }
        let mut stack: Vec<(usize, bool)> = vec![(r, false)];
        while let Some((v, fin)) = stack.pop() {
            if fin {
                time += 1;
                finish[v] = time;
            } else if !visited[v] {
                visited[v] = true;
                time += 1;
                stack.push((v, true));
                for &u in adj[v].iter().rev() {
                    if !visited[u] {
                        stack.push((u, false));
                    }
                }
            }
        }
    }
    let mut order: Vec<usize> = (0..ncomp).collect();
    order.sort_by(|&a, &b| finish[b].cmp(&finish[a]));
    let total: usize = comp_len.iter().sum();
    let mut rows = vec![vec![GAP; total]; nseq];
    for s in 0..nseq {
        let mut col = 0usize;
        let mut ci = 0usize;
        for &v in &g.seg_vertex[s] {
            let c = comp_of[v];
            while ci < order.len() && order[ci] != c {
                col += comp_len[order[ci]];
                ci += 1;
            }
            let (_, b, l) = g.vertices[v];
            rows[s][col..col + l].copy_from_slice(&seqs[s][b..b + l]);
            col += comp_len[c];
            ci += 1;
        }
    }
    rows
}

/// Global fragments, local fragments and distance of one pair of sequences.
type PairLibrary = (Vec<Frag>, Vec<Frag>, f64);

/// SeqAn's multiple alignment of Dna sequences (codes 0..=3); rows use 4 for gaps.
pub(crate) fn global_msa(seqs: &[Vec<u8>]) -> Vec<Vec<u8>> {
    let nseq = seqs.len();
    if nseq == 1 {
        return vec![seqs[0].clone()];
    }
    let deep = nseq >= 50;
    let threshold = if deep { 30 } else { 10 };
    // Groups are mostly copies of a few sequences: align each pair of distinct sequences once.
    let mut uniq: Vec<&Vec<u8>> = Vec::new();
    let uid: Vec<usize> = seqs
        .iter()
        .map(|s| match uniq.iter().position(|u| *u == s) {
            Some(k) => k,
            None => {
                uniq.push(s);
                uniq.len() - 1
            }
        })
        .collect();
    let mut cache: HashMap<(usize, usize), PairLibrary> = HashMap::new();
    let mut frags: Vec<Frag> = Vec::new();
    let mut local: Vec<Frag> = Vec::new();
    let mut dist = vec![0.0f64; nseq * nseq];
    for i in 0..nseq {
        for j in i + 1..nseq {
            let (gf, lf, d) = cache.entry((uid[i], uid[j])).or_insert_with(|| {
                let (gf, lf) = pair_fragments(&seqs[i], &seqs[j], 0, 1, deep);
                let d = pair_distance(&seqs[i], &seqs[j], &gf);
                (gf, lf, d)
            });
            dist[i * nseq + j] = *d;
            frags.extend(gf.iter().map(|f| Frag { s1: i, s2: j, ..*f }));
            local.extend(lf.iter().map(|f| Frag { s1: i, s2: j, ..*f }));
        }
    }
    frags.extend(local);
    let (mut g, edges) = build_graph(seqs, &frags);
    for d in dist.iter_mut() {
        *d = ((*d * 1e10) as i64) as f64 / 1e10;
    }
    let tree = nj_tree(&dist, nseq);
    g.adj = if nseq < threshold {
        triplet_extension(&g, &edges, None)
    } else {
        let labels = subtree_groups(&tree, nseq, threshold / 2);
        triplet_extension(&g, &edges, Some(&labels))
    };
    let cols = progressive(&g, &tree);
    alignment_matrix(&g, seqs, &cols)
}

/// `msa_consensus`: each sequence repeated by its count, aligned, then per column the most
/// frequent of A, C, G, T, gap (first one on ties); gap columns are dropped.
pub fn msa_consensus(seqs: &[(Vec<u8>, u64)]) -> String {
    let expanded: Vec<Vec<u8>> = seqs.iter().flat_map(|(s, c)| std::iter::repeat_n(s.clone(), *c as usize)).collect();
    if expanded.is_empty() {
        return String::new();
    }
    let rows = global_msa(&expanded);
    let width = rows[0].len();
    let mut out = String::new();
    for c in 0..width {
        let mut count = [0u64; 5];
        for r in &rows {
            count[r[c] as usize] += 1;
        }
        let mut bi = 0;
        for b in 1..5 {
            if count[b] > count[bi] {
                bi = b;
            }
        }
        if bi < 4 {
            out.push(b"ACGT"[bi] as char);
        }
    }
    out
}
