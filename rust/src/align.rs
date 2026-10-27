//! Pairwise alignment that reproduces SeqAn 2's `globalAlignment` with free end gaps
//! (`AlignConfig<true, true, true, true>`) bit-for-bit, including its DP evaluation order,
//! complete-trace tie bits, max-cell scouting order and GapsLeft traceback preferences.
//!
//! Two flavours are needed:
//!   * affine gaps (Gotoh) with Dna5 sequences: Porechop's adapter alignment
//!   * linear gaps with Dna sequences: the ab-initio "compatibility" check
//!
//! Sequence H (horizontal, columns) is the first SeqAn row; sequence V (vertical, rows) the second.

use std::cell::RefCell;

pub const GAP: u8 = 255;

// Trace bits.
pub(crate) const DIAG: u8 = 1;
pub(crate) const HOR: u8 = 2;
pub(crate) const VER: u8 = 4;
pub(crate) const HOR_OPEN: u8 = 8;
pub(crate) const VER_OPEN: u8 = 16;
pub(crate) const MAX_H: u8 = 32;
pub(crate) const MAX_V: u8 = 64;

const NEG: i32 = i32::MIN / 2;

/// SeqAn's Char -> Dna5 table: ACGT (any case) -> 0..3, U -> T, everything else -> N (4).
pub static DNA5: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t[b'U' as usize] = 3;
    t[b'u' as usize] = 3;
    t
};

/// SeqAn's Char -> Dna table: like Dna5 but unknown characters become A.
pub static DNA4: [u8; 256] = {
    let mut t = [0u8; 256];
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t[b'U' as usize] = 3;
    t[b'u' as usize] = 3;
    t
};

pub fn encode_dna5(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&c| DNA5[c as usize]).collect()
}

#[derive(Clone, Copy, Debug)]
pub struct Scoring {
    pub matching: i32,
    pub mismatch: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}

struct Workspace {
    trace: Vec<u8>,
    s: Vec<i32>,
    h: Vec<i32>,
}

thread_local! {
    static WS: RefCell<Workspace> = const { RefCell::new(Workspace { trace: Vec::new(), s: Vec::new(), h: Vec::new() }) };
}

/// The gapped alignment: per column, the code of each row or GAP.
pub struct Alignment {
    pub row_h: Vec<u8>,
    pub row_v: Vec<u8>,
    pub score: i32,
}

#[inline(always)]
fn max_trace(left: i32, right: i32, tl: u8, tr: u8) -> (i32, u8) {
    if left == right {
        (left, tl | tr)
    } else if left < right {
        (right, tr)
    } else {
        (left, tl)
    }
}

#[inline(always)]
fn max_trace_bl(left: i32, right: i32, tl: u8, tr: u8) -> (i32, u8) {
    let v = if left >= right { left } else { right };
    let t = (if left >= right { tl } else { 0 }) | (if right >= left { tr } else { 0 });
    (v, t)
}

/// Substitution scores of every adapter position against each Dna5 code.
fn profile(v: &[u8], sc: &Scoring) -> [Vec<i32>; 5] {
    let mut p: [Vec<i32>; 5] = Default::default();
    for (c, row) in p.iter_mut().enumerate() {
        *row = v.iter().map(|&x| if x == c as u8 { sc.matching } else { sc.mismatch }).collect();
    }
    p
}

/// One DP column (rows 1..=m). `s`/`hh` hold the previous column on entry and this column on exit.
#[inline(always)]
fn affine_column(s: &mut [i32], hh: &mut [i32], col: &mut [u8], prof: &[i32], sc: &Scoring) {
    let mut diag_s = 0i32; // cacheDiag score (row 0 is a free end gap: 0)
    let mut vert_s = 0i32; // cacheVert score
    let mut vert_v = NEG; // cacheVert vertical score (fresh per column)
    let (go, ge) = (sc.gap_open, sc.gap_extend);
    let it = s[1..].iter_mut().zip(hh[1..].iter_mut()).zip(col[1..].iter_mut()).zip(prof.iter());
    for (((si, hi), ci), &sub) in it {
        let inter = diag_s + sub;
        diag_s = *si;
        let (hv, th) = max_trace_bl(*hi + ge, *si + go, HOR, HOR_OPEN);
        *hi = hv;
        let (vv, tv) = max_trace_bl(vert_v + ge, vert_s + go, VER, VER_OPEN);
        vert_v = vv;
        let tmp = th | tv;
        let (gap_best, t2) = max_trace_bl(vv, hv, MAX_V, MAX_H);
        let (sv, t) = max_trace_bl(inter, gap_best, DIAG | tmp, t2 | tmp);
        *si = sv;
        vert_s = sv;
        *ci = t;
    }
}

/// Same recurrence without trace output (used for the score-only pass).
#[inline(always)]
fn affine_column_score(s: &mut [i32], hh: &mut [i32], prof: &[i32], sc: &Scoring) {
    let mut diag_s = 0i32;
    let mut vert_s = 0i32;
    let mut vert_v = NEG;
    let (go, ge) = (sc.gap_open, sc.gap_extend);
    for ((si, hi), &sub) in s[1..].iter_mut().zip(hh[1..].iter_mut()).zip(prof.iter()) {
        let inter = diag_s + sub;
        diag_s = *si;
        let hv = (*hi + ge).max(*si + go);
        *hi = hv;
        let vv = (vert_v + ge).max(vert_s + go);
        vert_v = vv;
        let sv = inter.max(vv.max(hv));
        *si = sv;
        vert_s = sv;
    }
}

/// Fills the trace matrix (column-major, (m+1) rows) and returns (score, best col, best row).
/// Mirrors SeqAn's evaluation order: column by column, the initial column fully zero (free
/// leading gaps), row 0 of every column zero, max tracked (strictly greater) on the last row of
/// every column and on every cell of the final column.
fn fill_affine(h: &[u8], v: &[u8], sc: &Scoring, ws: &mut Workspace) -> (i32, usize, usize) {
    fill_affine_range(h, v, sc, ws, 0, h.len(), None)
}

/// Fills columns (first, last] starting from the given column state (or the initial column).
fn fill_affine_range(
    h: &[u8],
    v: &[u8],
    sc: &Scoring,
    ws: &mut Workspace,
    first: usize,
    last: usize,
    state: Option<(&[i32], &[i32])>,
) -> (i32, usize, usize) {
    let n = h.len();
    let m = v.len();
    let rows = m + 1;
    let prof = profile(v, sc);
    ws.trace.clear();
    ws.trace.resize(rows * (last - first + 1), 0);
    ws.s.clear();
    ws.h.clear();
    match state {
        Some((s0, h0)) => {
            ws.s.extend_from_slice(s0);
            ws.h.extend_from_slice(h0);
        }
        None => {
            ws.s.resize(rows, 0);
            ws.h.resize(rows, NEG);
        }
    }
    let mut best = NEG;
    let (mut best_col, mut best_row) = (0usize, 0usize);
    if first == 0 {
        best = 0;
        best_col = 0;
        best_row = m;
    }
    for j in first + 1..=last {
        let col = &mut ws.trace[(j - first) * rows..(j - first + 1) * rows];
        ws.s[0] = 0;
        col[0] = 0;
        affine_column(&mut ws.s, &mut ws.h, col, &prof[h[j - 1] as usize], sc);
        if j == n {
            for i in 0..=m {
                if ws.s[i] > best {
                    best = ws.s[i];
                    best_col = j;
                    best_row = i;
                }
            }
        } else if ws.s[m] > best {
            best = ws.s[m];
            best_col = j;
            best_row = m;
        }
    }
    (best, best_col, best_row)
}

fn fill_linear(h: &[u8], v: &[u8], sc: &Scoring, ws: &mut Workspace) -> (i32, usize, usize) {
    let n = h.len();
    let m = v.len();
    let rows = m + 1;
    let gap = sc.gap_extend;
    ws.trace.clear();
    ws.trace.resize(rows * (n + 1), 0);
    ws.s.clear();
    ws.s.resize(rows, 0);
    let s = &mut ws.s;
    let trace = &mut ws.trace;
    let mut best = 0;
    let (mut best_col, mut best_row) = (0usize, m);
    for j in 1..=n {
        let hc = h[j - 1];
        let is_final = j == n;
        let col = &mut trace[j * rows..(j + 1) * rows];
        s[0] = 0;
        let mut diag_s = 0i32;
        let mut vert_s = 0i32;
        col[0] = 0;
        for i in 1..=m {
            let vc = v[i - 1];
            let sub = if hc == vc { sc.matching } else { sc.mismatch };
            let inter = diag_s + sub;
            diag_s = s[i];
            let (g, tv) = max_trace(vert_s + gap, s[i] + gap, VER | MAX_V, HOR | MAX_H);
            let (sv, t) = max_trace(inter, g, DIAG, tv);
            s[i] = sv;
            vert_s = sv;
            col[i] = t;
            if (is_final || i == m) && sv > best {
                best = sv;
                best_col = j;
                best_row = i;
            }
        }
    }
    (best, best_col, best_row)
}

#[derive(Clone, Copy, PartialEq)]
pub(crate) enum Op {
    Diag,
    Vert,
    Hor,
}

/// Returns the ops and the (column, row) where the traceback stopped. `first` is the column that
/// trace column 0 corresponds to.
fn traceback(trace: &[u8], rows: usize, first: usize, n: usize, m: usize, col: usize, row: usize, affine: bool) -> (Vec<Op>, usize, usize) {
    traceback_with(|c, r| trace[(c - first) * rows + r], first, n, m, col, row, affine)
}

pub(crate) fn traceback_with<F: Fn(usize, usize) -> u8>(at: F, first: usize, n: usize, m: usize, mut col: usize, mut row: usize, affine: bool) -> (Vec<Op>, usize, usize) {
    // Ops are collected from the end of the alignment backwards.
    let mut ops: Vec<Op> = Vec::with_capacity(n + m);
    // Trailing gaps.
    for _ in row..m {
        ops.push(Op::Vert);
    }
    for _ in col..n {
        ops.push(Op::Hor);
    }
    let mut tv = at(col, row);
    if affine {
        // PreferGapsAtEnd_ for affine gaps.
        if tv & MAX_V != 0 {
            tv &= VER | VER_OPEN | MAX_V;
        } else if tv & MAX_H != 0 {
            tv &= HOR | HOR_OPEN | MAX_H;
        }
    }
    'main: while col > first && row > 0 && tv != 0 {
        if tv & DIAG != 0 {
            col -= 1;
            row -= 1;
            tv = at(col, row);
            ops.push(Op::Diag);
        } else if tv & MAX_V != 0 && tv & VER != 0 {
            if affine {
                while (tv & VER_OPEN == 0 || tv & VER != 0) && row != 1 {
                    row -= 1;
                    tv = at(col, row);
                    ops.push(Op::Vert);
                }
            }
            row -= 1;
            tv = at(col, row);
            ops.push(Op::Vert);
        } else if tv & MAX_V != 0 && tv & VER_OPEN != 0 {
            row -= 1;
            tv = at(col, row);
            ops.push(Op::Vert);
        } else if tv & MAX_H != 0 && tv & HOR != 0 {
            if affine {
                while (tv & HOR_OPEN == 0 || tv & HOR != 0) && col != 1 && col > first {
                    col -= 1;
                    tv = at(col, row);
                    ops.push(Op::Hor);
                }
                if col == first && first > 0 {
                    break 'main; // reached the window boundary: caller widens the window
                }
            }
            col -= 1;
            tv = at(col, row);
            ops.push(Op::Hor);
        } else if tv & MAX_H != 0 && tv & HOR_OPEN != 0 {
            col -= 1;
            tv = at(col, row);
            ops.push(Op::Hor);
        } else {
            break;
        }
    }
    let (stop_col, stop_row) = (col, row);
    // Leading gaps.
    for _ in 0..row {
        ops.push(Op::Vert);
    }
    for _ in 0..col {
        ops.push(Op::Hor);
    }
    ops.reverse();
    (ops, stop_col, stop_row)
}

pub(crate) fn build_rows(h: &[u8], v: &[u8], ops: &[Op], score: i32) -> Alignment {
    let mut row_h = Vec::with_capacity(ops.len());
    let mut row_v = Vec::with_capacity(ops.len());
    let (mut i, mut j) = (0usize, 0usize);
    for op in ops {
        match op {
            Op::Diag => {
                row_h.push(h[j]);
                row_v.push(v[i]);
                i += 1;
                j += 1;
            }
            Op::Vert => {
                row_h.push(GAP);
                row_v.push(v[i]);
                i += 1;
            }
            Op::Hor => {
                row_h.push(h[j]);
                row_v.push(GAP);
                j += 1;
            }
        }
    }
    Alignment { row_h, row_v, score }
}

/// Affine-gap free-end-gap global alignment (SeqAn Gotoh). Sequences are pre-encoded.
/// SeqAn dispatches to the linear algorithm when gap open == gap extend.
pub fn align_free_ends(h: &[u8], v: &[u8], sc: &Scoring) -> Alignment {
    let affine = sc.gap_open != sc.gap_extend;
    if affine && h.len() > WINDOW_MIN_LEN {
        return align_affine_windowed(h, v, sc);
    }
    WS.with(|ws| {
        let mut ws = ws.borrow_mut();
        let (score, col, row) = if affine { fill_affine(h, v, sc, &mut ws) } else { fill_linear(h, v, sc, &mut ws) };
        let (ops, _, _) = traceback(&ws.trace, v.len() + 1, 0, h.len(), v.len(), col, row, affine);
        build_rows(h, v, &ops, score)
    })
}

const WINDOW_MIN_LEN: usize = 600;
const CHECKPOINT: usize = 64;

/// Exact equivalent of the full affine alignment for long sequences: a score-only pass over the
/// whole matrix (with periodic column checkpoints) finds the best cell, then the trace matrix is
/// only computed for a window of columns ending there, widened until the traceback fits.
fn align_affine_windowed(h: &[u8], v: &[u8], sc: &Scoring) -> Alignment {
    let n = h.len();
    let m = v.len();
    let rows = m + 1;
    let prof = profile(v, sc);
    // pass 1
    let mut s = vec![0i32; rows];
    let mut hh = vec![NEG; rows];
    let mut checkpoints: Vec<(Vec<i32>, Vec<i32>)> = Vec::with_capacity(n / CHECKPOINT + 1);
    checkpoints.push((s.clone(), hh.clone()));
    let mut best = 0;
    let (mut best_col, mut best_row) = (0usize, m);
    for j in 1..=n {
        s[0] = 0;
        affine_column_score(&mut s, &mut hh, &prof[h[j - 1] as usize], sc);
        if j == n {
            for i in 0..=m {
                if s[i] > best {
                    best = s[i];
                    best_col = j;
                    best_row = i;
                }
            }
        } else if s[m] > best {
            best = s[m];
            best_col = j;
            best_row = m;
        }
        if j % CHECKPOINT == 0 {
            checkpoints.push((s.clone(), hh.clone()));
        }
    }
    // pass 2
    let mut span = 2 * m + 32;
    WS.with(|ws| {
        let mut ws = ws.borrow_mut();
        loop {
            let want = best_col.saturating_sub(span);
            let cp = want / CHECKPOINT;
            let first = cp * CHECKPOINT;
            let state = if first == 0 { None } else { Some((checkpoints[cp].0.as_slice(), checkpoints[cp].1.as_slice())) };
            fill_affine_range(h, v, sc, &mut ws, first, best_col, state);
            let (ops, stop_col, stop_row) = traceback(&ws.trace, rows, first, n, m, best_col, best_row, true);
            if first == 0 || !(stop_col == first && stop_row > 0) {
                return build_rows(h, v, &ops, best);
            }
            span *= 2;
        }
    })
}

/// Result of Porechop's ScoredAlignment, as seen by the Python side after the
/// "%f" string round trip.
#[derive(Clone, Copy, Debug)]
pub struct AdapterHit {
    pub read_start: i64, // -1 when the alignment failed
    pub read_end: i64,   // exclusive (Python adds 1); 0 when failed
    pub aligned_region_identity: f64,
    pub full_adapter_identity: f64,
}

/// Emulates C++ std::to_string(double) ("%f") followed by Python float().
#[inline]
pub fn six_decimals(x: f64) -> f64 {
    let s = format!("{:.6}", x);
    s.parse().unwrap()
}

/// Porechop's adapterAlignment + ScoredAlignment + Python's align_adapter().
/// `read` and `adapter` must be Dna5-encoded (codes 0..=4).
pub fn align_adapter(read: &[u8], adapter: &[u8], sc: &Scoring) -> AdapterHit {
    let failed = AdapterHit { read_start: -1, read_end: 0, aligned_region_identity: 0.0, full_adapter_identity: 0.0 };
    if read.is_empty() || adapter.is_empty() {
        return failed;
    }
    let aln = align_free_ends(read, adapter, sc);
    scored_alignment(&aln).unwrap_or(failed)
}

pub(crate) fn scored_alignment(aln: &Alignment) -> Option<AdapterHit> {
    let r = &aln.row_h;
    let a = &aln.row_v;
    let len = r.len();
    if len == 0 {
        return None;
    }
    let mut start = None;
    let (mut rs, mut as_) = (false, false);
    for i in 0..len {
        rs |= r[i] != GAP;
        as_ |= a[i] != GAP;
        if rs && as_ {
            start = Some(i);
            break;
        }
    }
    let mut end = None;
    let (mut re, mut ae) = (false, false);
    for i in (0..len).rev() {
        re |= r[i] != GAP;
        ae |= a[i] != GAP;
        if re && ae {
            end = Some(i);
            break;
        }
    }
    let (start, end) = (start?, end?);
    let ad_start = a.iter().position(|&c| c != GAP)?;
    let ad_end = a.iter().rposition(|&c| c != GAP)?;
    let matches = |lo: usize, hi: usize| (lo..=hi).filter(|&i| a[i] == r[i]).count();
    let region = 100.0 * matches(start, end) as f64 / (end - start + 1) as f64;
    let full = 100.0 * matches(ad_start, ad_end) as f64 / (ad_end - ad_start + 1) as f64;
    let read_start = r[..start].iter().filter(|&&c| c != GAP).count() as i64;
    let read_end_pos = r[..end].iter().filter(|&&c| c != GAP).count() as i64;
    Some(AdapterHit {
        read_start,
        read_end: read_end_pos + 1,
        aligned_region_identity: six_decimals(region),
        full_adapter_identity: six_decimals(full),
    })
}

/// Minimum edit distance between the whole `pattern` and any substring of `text`
/// (Myers' bit-vector algorithm, Hyyrö formulation). Requires pattern.len() <= 64.
pub fn myers_semiglobal(pattern: &[u8], text: &[u8]) -> usize {
    let m = pattern.len();
    debug_assert!(m > 0 && m <= 64);
    let mut peq = [0u64; 5];
    for (i, &c) in pattern.iter().enumerate() {
        peq[c as usize] |= 1u64 << i;
    }
    let high = 1u64 << (m - 1);
    let mut pv: u64 = if m == 64 { !0 } else { (1u64 << m) - 1 };
    let mut mv: u64 = 0;
    let mut score = m;
    let mut best = m;
    for &c in text {
        let eq = peq[c as usize];
        let xv = eq | mv;
        let xh = (((eq & pv).wrapping_add(pv)) ^ pv) | eq;
        let ph = mv | !(xh | pv);
        let mh = pv & xh;
        score += ((ph & high) != 0) as usize;
        score -= ((mh & high) != 0) as usize;
        let ph = ph << 1;
        let mh = mh << 1;
        pv = mh | !(xv | ph);
        mv = ph & xv;
        best = best.min(score);
    }
    best
}

/// Upper bound on the full-adapter identity Porechop can report for this pair.
pub fn identity_upper_bound(adapter: &[u8], text: &[u8]) -> f64 {
    let m = adapter.len();
    if m == 0 || m > 64 {
        return 100.0;
    }
    let ed = myers_semiglobal(adapter, text);
    100.0 * m as f64 / (m + ed) as f64
}

/// The ab-initio compatibility flag (compatibility.cpp).
/// 0: not compatible, 1: compatible, 2: compatible and the shorter one is included in the longer one.
pub fn check_compatibility(raw1: &[u8], raw2: &[u8]) -> i32 {
    let e1: Vec<u8> = raw1.iter().map(|&c| DNA4[c as usize]).collect();
    let e2: Vec<u8> = raw2.iter().map(|&c| DNA4[c as usize]).collect();
    let (h, v) = if e1.len() < e2.len() { (e2, e1) } else { (e1, e2) };
    if h.is_empty() || v.is_empty() {
        return 0;
    }
    let sc = Scoring { matching: 2, mismatch: -1, gap_open: -1, gap_extend: -1 };
    let aln = align_free_ends(&h, &v, &sc);
    let val = |c: u8| c;
    let r1 = &aln.row_h;
    let r2 = &aln.row_v;
    let len = r1.len() as i64;
    let mut start: i64 = -1;
    let (mut s1, mut s2) = (false, false);
    for i in 0..len {
        s1 |= val(r1[i as usize]) != GAP;
        s2 |= val(r2[i as usize]) != GAP;
        if s1 && s2 {
            start = i;
            break;
        }
    }
    let mut end: i64 = -1;
    let (mut e1, mut e2) = (false, false);
    for i in (0..len).rev() {
        e1 |= val(r1[i as usize]) != GAP;
        e2 |= val(r2[i as usize]) != GAP;
        if e1 && e2 {
            end = i;
            break;
        }
    }
    let mapped = end - start + 1;
    let mut lev = 0i64;
    let mut i = start;
    while i < end {
        if val(r1[i as usize]) != val(r2[i as usize]) {
            lev += 1;
        }
        i += 1;
    }
    // toSourcePosition(row1, pos): number of source characters before view position pos.
    let to_source = |pos: i64| -> i64 {
        if pos < 0 {
            return 0;
        }
        r1[..(pos as usize).min(r1.len())].iter().filter(|&&c| c != GAP).count() as i64
    };
    let s1_start = to_source(start);
    let s1_end = to_source(end);
    let s1_len = h.len() as i64;
    let included = s1_start > 3 || (s1_len - s1_end) > 3;
    if mapped <= 0 {
        return 0;
    }
    let identity = ((mapped - lev) * 100 / mapped) as f32;
    let mut flag = 0;
    if identity >= 87.5 {
        flag = 1;
    }
    if included && flag == 1 {
        flag = 2;
    }
    flag
}
