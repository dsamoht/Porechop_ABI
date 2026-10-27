//! AVX2 inter-sequence batch version of the affine free-end-gap alignment: one adapter against
//! up to 16 equally long read regions at once (one 16-bit lane per region). It computes exactly
//! the same scores and trace bits as `align.rs` (same comparisons, same tie handling), then does
//! the traceback per lane with the shared scalar code.

use crate::align::*;

/// Aligns `adapter` against every region; returns one hit per region.
pub fn align_adapter_many(regions: &[&[u8]], adapter: &[u8], sc: &Scoring) -> Vec<AdapterHit> {
    let mut out = vec![None; regions.len()];
    #[cfg(target_arch = "x86_64")]
    {
        if simd_enabled() && sc.gap_open != sc.gap_extend && !adapter.is_empty() && adapter.len() <= 512 {
            // group equal-length regions
            let mut by_len: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
            for (i, r) in regions.iter().enumerate() {
                if !r.is_empty() && r.len() + adapter.len() <= 1500 && scores_fit(sc) {
                    by_len.entry(r.len()).or_default().push(i);
                }
            }
            for (_, idx) in by_len {
                if idx.len() < 4 {
                    continue; // not worth a vector pass
                }
                for group in idx.chunks(16) {
                    let regs: Vec<&[u8]> = group.iter().map(|&i| regions[i]).collect();
                    // SAFETY: AVX2 availability checked by simd_enabled().
                    let hits = unsafe { avx2::batch(&regs, adapter, sc) };
                    for (k, &i) in group.iter().enumerate() {
                        out[i] = Some(hits[k]);
                    }
                }
            }
        }
    }
    out.into_iter()
        .enumerate()
        .map(|(i, h)| h.unwrap_or_else(|| align_adapter(regions[i], adapter, sc)))
        .collect()
}

fn scores_fit(sc: &Scoring) -> bool {
    [sc.matching, sc.mismatch, sc.gap_open, sc.gap_extend].iter().all(|v| v.abs() <= 8)
}

pub fn simd_enabled() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        use std::sync::OnceLock;
        static ON: OnceLock<bool> = OnceLock::new();
        *ON.get_or_init(|| std::env::var_os("PORECHOP_RS_NO_SIMD").is_none() && is_x86_feature_detected!("avx2"))
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

#[cfg(target_arch = "x86_64")]
mod avx2 {
    use super::*;
    use std::arch::x86_64::*;

    const NEG: i16 = -30000;

    #[target_feature(enable = "avx2")]
    #[allow(unsafe_op_in_unsafe_fn)]
    pub unsafe fn batch(regs: &[&[u8]], adapter: &[u8], sc: &Scoring) -> Vec<AdapterHit> {
        let lanes = regs.len();
        let n = regs[0].len();
        let m = adapter.len();
        let rows = m + 1;
        // column-major trace, 16 lanes of u16 per cell
        let mut trace: Vec<__m256i> = vec![_mm256_setzero_si256(); (n + 1) * rows];
        let mut s: Vec<__m256i> = vec![_mm256_setzero_si256(); rows];
        let mut hh: Vec<__m256i> = vec![_mm256_set1_epi16(NEG); rows];
        let ge = _mm256_set1_epi16(sc.gap_extend as i16);
        let go = _mm256_set1_epi16(sc.gap_open as i16);
        let mat = _mm256_set1_epi16(sc.matching as i16);
        let mis = _mm256_set1_epi16(sc.mismatch as i16);
        let b = |x: u8| _mm256_set1_epi16(x as i16);
        let (diag_b, hor_b, ver_b, hopen_b, vopen_b, maxh_b, maxv_b) =
            (b(DIAG), b(HOR), b(VER), b(HOR_OPEN), b(VER_OPEN), b(MAX_H), b(MAX_V));
        let adapter_v: Vec<__m256i> = adapter.iter().map(|&c| _mm256_set1_epi16(c as i16)).collect();
        // ge(a, b): lanes where a >= b, as all-ones masks
        #[inline(always)]
        unsafe fn ge_mask(a: __m256i, b: __m256i) -> __m256i {
            _mm256_or_si256(_mm256_cmpgt_epi16(a, b), _mm256_cmpeq_epi16(a, b))
        }
        #[inline(always)]
        unsafe fn sel(mask: __m256i, bits: __m256i) -> __m256i {
            _mm256_and_si256(mask, bits)
        }
        // scout state
        let mut best = _mm256_setzero_si256(); // initial column, last row: score 0
        let mut best_col = _mm256_setzero_si256();
        let mut best_row = _mm256_set1_epi16(m as i16);
        let mut codes = [0i16; 16];
        for j in 1..=n {
            for l in 0..16 {
                codes[l] = regs[l.min(lanes - 1)][j - 1] as i16;
            }
            let cv = _mm256_loadu_si256(codes.as_ptr() as *const __m256i);
            let base = j * rows;
            let mut diag_s = _mm256_setzero_si256();
            let mut vert_s = _mm256_setzero_si256();
            let mut vert_v = _mm256_set1_epi16(NEG);
            *s.get_unchecked_mut(0) = _mm256_setzero_si256();
            for i in 1..=m {
                let eq = _mm256_cmpeq_epi16(cv, *adapter_v.get_unchecked(i - 1));
                let sub = _mm256_blendv_epi8(mis, mat, eq);
                let inter = _mm256_adds_epi16(diag_s, sub);
                let si = *s.get_unchecked(i);
                diag_s = si;
                let hl = _mm256_adds_epi16(*hh.get_unchecked(i), ge);
                let hr = _mm256_adds_epi16(si, go);
                let hv = _mm256_max_epi16(hl, hr);
                let th = _mm256_or_si256(sel(ge_mask(hl, hr), hor_b), sel(ge_mask(hr, hl), hopen_b));
                *hh.get_unchecked_mut(i) = hv;
                let vl = _mm256_adds_epi16(vert_v, ge);
                let vr = _mm256_adds_epi16(vert_s, go);
                let vv = _mm256_max_epi16(vl, vr);
                let tv = _mm256_or_si256(sel(ge_mask(vl, vr), ver_b), sel(ge_mask(vr, vl), vopen_b));
                vert_v = vv;
                let tmp = _mm256_or_si256(th, tv);
                let gap = _mm256_max_epi16(vv, hv);
                let t2 = _mm256_or_si256(sel(ge_mask(vv, hv), maxv_b), sel(ge_mask(hv, vv), maxh_b));
                let sv = _mm256_max_epi16(inter, gap);
                let t = _mm256_or_si256(
                    sel(ge_mask(inter, gap), _mm256_or_si256(diag_b, tmp)),
                    sel(ge_mask(gap, inter), _mm256_or_si256(t2, tmp)),
                );
                *s.get_unchecked_mut(i) = sv;
                vert_s = sv;
                *trace.get_unchecked_mut(base + i) = t;
            }
            let jv = _mm256_set1_epi16(j as i16);
            if j == n {
                for i in 0..=m {
                    let v = *s.get_unchecked(i);
                    let gt = _mm256_cmpgt_epi16(v, best);
                    best = _mm256_blendv_epi8(best, v, gt);
                    best_col = _mm256_blendv_epi8(best_col, jv, gt);
                    best_row = _mm256_blendv_epi8(best_row, _mm256_set1_epi16(i as i16), gt);
                }
            } else {
                let v = *s.get_unchecked(m);
                let gt = _mm256_cmpgt_epi16(v, best);
                best = _mm256_blendv_epi8(best, v, gt);
                best_col = _mm256_blendv_epi8(best_col, jv, gt);
            }
        }
        let mut bv = [0i16; 16];
        let mut bc = [0i16; 16];
        let mut br = [0i16; 16];
        _mm256_storeu_si256(bv.as_mut_ptr() as *mut __m256i, best);
        _mm256_storeu_si256(bc.as_mut_ptr() as *mut __m256i, best_col);
        _mm256_storeu_si256(br.as_mut_ptr() as *mut __m256i, best_row);
        let tr: &[u16] = std::slice::from_raw_parts(trace.as_ptr() as *const u16, trace.len() * 16);
        let failed = AdapterHit { read_start: -1, read_end: 0, aligned_region_identity: 0.0, full_adapter_identity: 0.0 };
        (0..lanes)
            .map(|l| {
                let at = |c: usize, r: usize| tr[(c * rows + r) * 16 + l] as u8;
                let (ops, _, _) = traceback_with(at, 0, n, m, bc[l] as usize, br[l] as usize, true);
                let aln = build_rows(regs[l], adapter, &ops, bv[l] as i32);
                scored_alignment(&aln).unwrap_or(failed)
            })
            .collect()
    }
}

/// Minimum semi-global edit distance of each pattern (<= 64 bp) against `text`
/// (same as `myers_semiglobal`), four patterns per AVX2 pass.
pub fn myers_many(patterns: &[&[u8]], text: &[u8]) -> Vec<usize> {
    #[cfg(target_arch = "x86_64")]
    {
        if simd_enabled() {
            let mut out = Vec::with_capacity(patterns.len());
            for g in patterns.chunks(4) {
                if g.iter().all(|p| !p.is_empty() && p.len() <= 64) {
                    // SAFETY: AVX2 checked above.
                    out.extend_from_slice(&unsafe { avx2_myers::four(g, text) }[..g.len()]);
                } else {
                    out.extend(g.iter().map(|p| if p.is_empty() || p.len() > 64 { 0 } else { myers_semiglobal(p, text) }));
                }
            }
            return out;
        }
    }
    patterns.iter().map(|p| if p.is_empty() || p.len() > 64 { 0 } else { myers_semiglobal(p, text) }).collect()
}

#[cfg(target_arch = "x86_64")]
mod avx2_myers {
    use std::arch::x86_64::*;

    #[target_feature(enable = "avx2")]
    #[allow(unsafe_op_in_unsafe_fn)]
    pub unsafe fn four(pats: &[&[u8]], text: &[u8]) -> [usize; 4] {
        let mut peq = [[0u64; 4]; 5];
        let mut pv0 = [0u64; 4];
        let mut high = [0u64; 4];
        let mut score0 = [0i64; 4];
        for l in 0..4 {
            let p = pats[l.min(pats.len() - 1)];
            let m = p.len();
            for (i, &c) in p.iter().enumerate() {
                peq[c as usize][l] |= 1u64 << i;
            }
            pv0[l] = if m == 64 { !0 } else { (1u64 << m) - 1 };
            high[l] = 1u64 << (m - 1);
            score0[l] = m as i64;
        }
        let peqv: [__m256i; 5] = std::array::from_fn(|c| _mm256_loadu_si256(peq[c].as_ptr() as *const __m256i));
        let highv = _mm256_loadu_si256(high.as_ptr() as *const __m256i);
        let mut pv = _mm256_loadu_si256(pv0.as_ptr() as *const __m256i);
        let mut mv = _mm256_setzero_si256();
        let mut score = _mm256_loadu_si256(score0.as_ptr() as *const __m256i);
        let mut best = score;
        let ones = _mm256_set1_epi64x(-1);
        let zero = _mm256_setzero_si256();
        for &c in text {
            let eq = *peqv.get_unchecked(c as usize);
            let xv = _mm256_or_si256(eq, mv);
            let sum = _mm256_add_epi64(_mm256_and_si256(eq, pv), pv);
            let xh = _mm256_or_si256(_mm256_xor_si256(sum, pv), eq);
            let ph = _mm256_or_si256(mv, _mm256_xor_si256(_mm256_or_si256(xh, pv), ones));
            let mh = _mm256_and_si256(pv, xh);
            // +1 where ph has the high bit, -1 where mh has it (cmpeq gives -1 when bit set)
            let inc = _mm256_xor_si256(_mm256_cmpeq_epi64(_mm256_and_si256(ph, highv), zero), ones);
            let dec = _mm256_xor_si256(_mm256_cmpeq_epi64(_mm256_and_si256(mh, highv), zero), ones);
            score = _mm256_add_epi64(_mm256_sub_epi64(score, inc), dec);
            let ph2 = _mm256_slli_epi64(ph, 1);
            let mh2 = _mm256_slli_epi64(mh, 1);
            pv = _mm256_or_si256(mh2, _mm256_xor_si256(_mm256_or_si256(xv, ph2), ones));
            mv = _mm256_and_si256(ph2, xv);
            let lt = _mm256_cmpgt_epi64(best, score);
            best = _mm256_blendv_epi8(best, score, lt);
        }
        let mut b = [0i64; 4];
        _mm256_storeu_si256(b.as_mut_ptr() as *mut __m256i, best);
        [b[0] as usize, b[1] as usize, b[2] as usize, b[3] as usize]
    }
}
