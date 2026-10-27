//! Port of approx_counter.cpp: exact k-mer count on sampled read ends, then approximate
//! (<= 2 edits) recount of the most frequent k-mers.
//!
//! The original uses a bidirectional FM index and SeqAn's optimum search schemes. For each
//! k-mer it sums, over e in {0, 1, 2}, the number of reads in which the search reported an
//! occurrence with exactly e errors. The set of strings reported with e errors is exactly the
//! set of strings at edit distance e (e = 0, 1) or <= 2 (e = 2) from the k-mer (verified by
//! exhaustive enumeration in `scheme.rs`), so a read contributes 3 - d when its best
//! semi-global edit distance d to the k-mer is <= 2, and nothing otherwise. That is computed
//! here with Myers' bit-parallel algorithm, several k-mers per machine word.

use rayon::prelude::*;
use std::collections::HashMap;

pub type Kmer = u64; // 2-bit packed, first base in the most significant position

pub fn kmer_to_string(km: Kmer, k: usize) -> String {
    let mut s = vec![0u8; k];
    let mut v = km;
    for i in (0..k).rev() {
        s[i] = b"ACGT"[(v & 3) as usize];
        v >>= 2;
    }
    String::from_utf8(s).unwrap()
}

/// DUST-like dimer score (approx_counter's getComplexity), in f32 like the original.
pub fn complexity(kmer: Kmer, k: usize) -> f32 {
    let mut counts = [0u64; 16];
    let mut km = kmer;
    for _ in 0..k.saturating_sub(1) {
        counts[(km & 15) as usize] += 1;
        km >>= 2;
    }
    let sum: u64 = counts.iter().map(|&v| v * v.saturating_sub(1)).sum();
    sum as f32 / (2 * (k as i64 - 2)) as f32
}

pub fn adjust_threshold(c_old: f32, k_old: usize, k_new: usize) -> f32 {
    let num = ((k_new as f64 - 2.0 + 1.0).powi(2)) as f32;
    let den = ((k_old as f64 - 2.0 + 1.0).powi(2)) as f32;
    c_old * (num / den)
}

/// CompareCount: count descending, then complexity ascending, then k-mer value descending.
/// This is a total order, so partial selection followed by sorting gives the same prefix.
pub fn sort_counts(v: &mut Vec<(Kmer, u64)>, k: usize, limit: Option<usize>) {
    let mut keyed: Vec<(u64, f32, Kmer)> = v.iter().map(|&(km, c)| (c, complexity(km, k), km)).collect();
    let cmp = |a: &(u64, f32, Kmer), b: &(u64, f32, Kmer)| {
        b.0.cmp(&a.0).then_with(|| if a.1 == b.1 { b.2.cmp(&a.2) } else { a.1.partial_cmp(&b.1).unwrap() })
    };
    if let Some(l) = limit {
        if keyed.len() > l && l > 0 {
            keyed.select_nth_unstable_by(l - 1, cmp);
            keyed.truncate(l);
        } else if l == 0 {
            keyed.clear();
        }
    }
    keyed.sort_unstable_by(cmp);
    *v = keyed.into_iter().map(|(c, _, km)| (km, c)).collect();
}

/// Exact count of all k-mers without N, skipping low-complexity and forbidden ones.
/// Returns (counts, number of k-mers skipped because of N).
pub fn count_kmers(samples: &[&[u8]], k: usize, lc: f32, forbidden: &std::collections::HashSet<Kmer>) -> (HashMap<Kmer, u64>, u64) {
    let mask: u64 = if k == 32 { !0 } else { (1u64 << (2 * k)) - 1 };
    let parts: Vec<(Vec<Kmer>, u64)> = samples
        .par_chunks(1024)
        .map(|chunk| {
            let mut all: Vec<Kmer> = Vec::with_capacity(chunk.len() * 100);
            let mut had_n = 0u64;
            for s in chunk {
                if s.len() < k {
                    continue;
                }
                let mut val: u64 = 0;
                let mut since_n = 0usize; // number of valid bases since the last N
                for (i, &c) in s.iter().enumerate() {
                    if c >= 4 {
                        since_n = 0;
                        val = 0;
                    } else {
                        val = ((val << 2) | c as u64) & mask;
                        since_n += 1;
                    }
                    if i + 1 >= k {
                        if since_n >= k {
                            all.push(val);
                        } else {
                            had_n += 1;
                        }
                    }
                }
            }
            (all, had_n)
        })
        .collect();
    let had_n = parts.iter().map(|p| p.1).sum();
    let mut all: Vec<Kmer> = Vec::with_capacity(parts.iter().map(|p| p.0.len()).sum());
    for (p, _) in parts {
        all.extend(p);
    }
    all.par_sort_unstable();
    let mut runs: Vec<(Kmer, u64)> = Vec::new();
    let mut i = 0;
    while i < all.len() {
        let mut j = i + 1;
        while j < all.len() && all[j] == all[i] {
            j += 1;
        }
        runs.push((all[i], (j - i) as u64));
        i = j;
    }
    let total: HashMap<Kmer, u64> = runs
        .into_par_iter()
        .filter(|(km, _)| complexity(*km, k) < lc && !forbidden.contains(km))
        .collect();
    (total, had_n)
}

/// For each k-mer, sum over reads of max(0, 3 - min edit distance) (see module docs).
pub fn approximate_counts(samples: &[&[u8]], kmers: &[Kmer], k: usize) -> Vec<u64> {
    #[cfg(target_arch = "x86_64")]
    {
        if k <= 16 && std::env::var_os("PORECHOP_RS_NO_SIMD").is_none() && is_x86_feature_detected!("avx2") {
            return approximate_counts_avx2(samples, kmers, k);
        }
    }
    if k <= 16 {
        approximate_counts_swar(samples, kmers, k)
    } else {
        approximate_counts_scalar(samples, kmers, k)
    }
}

fn kmer_codes(km: Kmer, k: usize) -> Vec<u8> {
    (0..k).map(|i| ((km >> (2 * (k - 1 - i))) & 3) as u8).collect()
}

fn approximate_counts_scalar(samples: &[&[u8]], kmers: &[Kmer], k: usize) -> Vec<u64> {
    kmers
        .par_iter()
        .map(|&km| {
            let pat = kmer_codes(km, k);
            let mut total = 0u64;
            for s in samples {
                let d = crate::align::myers_semiglobal(&pat, s);
                if d <= 2 {
                    total += 3 - d as u64;
                }
            }
            total
        })
        .collect()
}

/// Four k-mers (k <= 16) per u64, one 16-bit lane each, Myers' algorithm in SWAR form.
fn approximate_counts_swar(samples: &[&[u8]], kmers: &[Kmer], k: usize) -> Vec<u64> {
    const LANES: usize = 4;
    const LO: u64 = 0x0001_0001_0001_0001;
    const HI15: u64 = 0x8000_8000_8000_8000;
    let lane_mask: u64 = if k == 16 { 0xFFFF } else { (1u64 << k) - 1 };
    let full_mask = lane_mask * LO;
    let high_bit = 1u64 << (k - 1);
    let high = high_bit * LO;
    let not_bit0 = !LO;
    // lane-wise add without carries between lanes (k <= 16 bits per lane)
    #[inline(always)]
    fn add16(x: u64, y: u64) -> u64 {
        let low = (x & !HI15).wrapping_add(y & !HI15);
        low ^ ((x ^ y) & HI15)
    }
    let groups: Vec<&[Kmer]> = kmers.chunks(LANES).collect();
    let per_group: Vec<[u64; LANES]> = groups
        .par_iter()
        .map(|g| {
            let mut peq = [0u64; 5];
            for (lane, &km) in g.iter().enumerate() {
                for (i, c) in kmer_codes(km, k).into_iter().enumerate() {
                    peq[c as usize] |= 1u64 << (i + 16 * lane);
                }
            }
            let mut totals = [0u64; LANES];
            for s in samples {
                let mut pv = full_mask;
                let mut mv = 0u64;
                // scores per lane, start at k
                let mut score = (k as u64) * LO;
                let mut le0 = 0u64;
                let mut le1 = 0u64;
                let mut le2 = 0u64;
                // offsets so that bit 15 of each lane is set when score > t (scores < 0x8000)
                let gt0 = (0x7FFF - 0) * LO;
                let gt1 = (0x7FFF - 1) * LO;
                let gt2 = (0x7FFF - 2) * LO;
                for &c in s.iter() {
                    let eq = peq[c as usize];
                    let xv = eq | mv;
                    let xh = (add16(eq & pv, pv) ^ pv) | eq;
                    let ph = (mv | !(xh | pv)) & full_mask;
                    let mh = pv & xh;
                    let inc = (ph & high) >> (k - 1);
                    let dec = (mh & high) >> (k - 1);
                    score = score + inc - dec;
                    let ph2 = (ph << 1) & not_bit0 & full_mask;
                    let mh2 = (mh << 1) & not_bit0 & full_mask;
                    pv = (mh2 | !(xv | ph2)) & full_mask;
                    mv = ph2 & xv;
                    le0 |= !(score + gt0) & HI15;
                    le1 |= !(score + gt1) & HI15;
                    le2 |= !(score + gt2) & HI15;
                }
                for (lane, t) in totals.iter_mut().enumerate().take(g.len()) {
                    let b = 15 + 16 * lane;
                    *t += ((le0 >> b) & 1) + ((le1 >> b) & 1) + ((le2 >> b) & 1);
                }
            }
            totals
        })
        .collect();
    let mut out = Vec::with_capacity(kmers.len());
    for (g, t) in groups.iter().zip(per_group) {
        out.extend_from_slice(&t[..g.len()]);
    }
    out
}

pub fn encode_kmer(s: &[u8]) -> Option<Kmer> {
    let mut v = 0u64;
    for &c in s {
        let code = crate::align::DNA5[c as usize];
        if code >= 4 {
            return None;
        }
        v = (v << 2) | code as u64;
    }
    Some(v)
}

/// Sixteen k-mers (k <= 16) per AVX2 register, one 16-bit lane each.
#[cfg(target_arch = "x86_64")]
fn approximate_counts_avx2(samples: &[&[u8]], kmers: &[Kmer], k: usize) -> Vec<u64> {
    let groups: Vec<&[Kmer]> = kmers.chunks(16).collect();
    let per_group: Vec<[u64; 16]> = groups
        .par_iter()
        .map(|g| {
            let mut peq = [[0u16; 16]; 5];
            for (lane, &km) in g.iter().enumerate() {
                for (i, c) in kmer_codes(km, k).into_iter().enumerate() {
                    peq[c as usize][lane] |= 1u16 << i;
                }
            }
            // SAFETY: only called after AVX2 detection.
            unsafe { avx2_group(&peq, samples, k) }
        })
        .collect();
    let mut out = Vec::with_capacity(kmers.len());
    for (g, t) in groups.iter().zip(per_group) {
        out.extend_from_slice(&t[..g.len()]);
    }
    out
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn avx2_group(peq_lanes: &[[u16; 16]; 5], samples: &[&[u8]], k: usize) -> [u64; 16] {
    unsafe { avx2_group_inner(peq_lanes, samples, k) }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(unsafe_op_in_unsafe_fn)]
unsafe fn avx2_group_inner(peq_lanes: &[[u16; 16]; 5], samples: &[&[u8]], k: usize) -> [u64; 16] {
    use std::arch::x86_64::*;
    let mut peq = [_mm256_setzero_si256(); 5];
    for c in 0..5 {
        peq[c] = _mm256_loadu_si256(peq_lanes[c].as_ptr() as *const __m256i);
    }
    let lane_mask: i16 = if k == 16 { -1 } else { ((1u32 << k) - 1) as i16 };
    let full = _mm256_set1_epi16(lane_mask);
    let high = _mm256_set1_epi16((1u32 << (k - 1)) as u16 as i16);
    let shift = _mm_cvtsi32_si128((k - 1) as i32);
    let ones = _mm256_set1_epi16(-1);
    let t0 = _mm256_set1_epi16(1); // score < 1  <=> score <= 0
    let t1 = _mm256_set1_epi16(2);
    let t2 = _mm256_set1_epi16(3);
    let mut totals = [0u64; 16];
    let mut acc = [_mm256_setzero_si256(); 3];
    let mut pending = 0u32;
    for s in samples {
        let mut pv = full;
        let mut mv = _mm256_setzero_si256();
        let mut score = _mm256_set1_epi16(k as i16);
        let mut le0 = _mm256_setzero_si256();
        let mut le1 = _mm256_setzero_si256();
        let mut le2 = _mm256_setzero_si256();
        for &c in s.iter() {
            let eq = *peq.get_unchecked(c as usize);
            let xv = _mm256_or_si256(eq, mv);
            let sum = _mm256_add_epi16(_mm256_and_si256(eq, pv), pv);
            let xh = _mm256_or_si256(_mm256_xor_si256(sum, pv), eq);
            let ph = _mm256_and_si256(_mm256_or_si256(mv, _mm256_xor_si256(_mm256_or_si256(xh, pv), ones)), full);
            let mh = _mm256_and_si256(pv, xh);
            let inc = _mm256_srl_epi16(_mm256_and_si256(ph, high), shift);
            let dec = _mm256_srl_epi16(_mm256_and_si256(mh, high), shift);
            score = _mm256_sub_epi16(_mm256_add_epi16(score, inc), dec);
            let ph2 = _mm256_and_si256(_mm256_slli_epi16(ph, 1), full);
            let mh2 = _mm256_and_si256(_mm256_slli_epi16(mh, 1), full);
            pv = _mm256_and_si256(_mm256_or_si256(mh2, _mm256_xor_si256(_mm256_or_si256(xv, ph2), ones)), full);
            mv = _mm256_and_si256(ph2, xv);
            le0 = _mm256_or_si256(le0, _mm256_cmpgt_epi16(t0, score));
            le1 = _mm256_or_si256(le1, _mm256_cmpgt_epi16(t1, score));
            le2 = _mm256_or_si256(le2, _mm256_cmpgt_epi16(t2, score));
        }
        // each le lane is 0 or -1: subtracting accumulates +1 per read
        acc[0] = _mm256_sub_epi16(acc[0], le0);
        acc[1] = _mm256_sub_epi16(acc[1], le1);
        acc[2] = _mm256_sub_epi16(acc[2], le2);
        pending += 1;
        if pending == 30000 {
            flush(&mut acc, &mut totals);
            pending = 0;
        }
    }
    flush(&mut acc, &mut totals);
    #[target_feature(enable = "avx2")]
    #[allow(unsafe_op_in_unsafe_fn)]
    unsafe fn flush(acc: &mut [__m256i; 3], totals: &mut [u64; 16]) {
        for a in acc.iter_mut() {
            let mut lanes = [0u16; 16];
            _mm256_storeu_si256(lanes.as_mut_ptr() as *mut __m256i, *a);
            for i in 0..16 {
                totals[i] += lanes[i] as u64;
            }
            *a = _mm256_setzero_si256();
        }
    }
    totals
}
