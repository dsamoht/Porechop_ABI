//! `porechop_abi._native`: replacements for Porechop_ABI's former C++ parts
//! (`cpp_functions.so`, `compatibility.so`, `approx_counter`, `msa_consensus`), taken from
//! porechop_rs. Every function returns exactly what the Python code got from the original
//! component; batch functions release the GIL and run on `threads` threads.

mod align;
mod align_simd;
mod counter;
mod kmer;
mod msa;
mod rng;

use align::{align_adapter as align_one, check_compatibility, encode_dna5, identity_upper_bound, AdapterHit, Scoring, DNA4};
use align_simd::{align_adapter_many, myers_many};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rayon::prelude::*;
use std::io::Write;

/// (full adapter identity, aligned region identity, read start, read end), as returned by
/// nanopore_read.align_adapter().
type Hit = (f64, f64, i64, i64);

/// (adapter index, full adapter identity, read start, read end) of one middle adapter hit.
type MiddleHit = (usize, f64, i64, i64);

fn hit(h: AdapterHit) -> Hit {
    (h.full_adapter_identity, h.aligned_region_identity, h.read_start, h.read_end)
}

fn scoring(v: &[i32]) -> PyResult<Scoring> {
    match v {
        [a, b, c, d] => Ok(Scoring { matching: *a, mismatch: *b, gap_open: *c, gap_extend: *d }),
        _ => Err(PyValueError::new_err("scoring scheme must have 4 values")),
    }
}

fn pool(threads: usize) -> PyResult<rayon::ThreadPool> {
    rayon::ThreadPoolBuilder::new().num_threads(threads.max(1)).build().map_err(|e| PyValueError::new_err(e.to_string()))
}

/// nanopore_read.align_adapter(read_seq, adapter_seq, scoring_scheme_vals).
#[pyfunction]
fn align_adapter(read_seq: &str, adapter_seq: &str, scoring_scheme_vals: Vec<i32>) -> PyResult<Hit> {
    let sc = scoring(&scoring_scheme_vals)?;
    Ok(hit(align_one(&encode_dna5(read_seq.as_bytes()), &encode_dna5(adapter_seq.as_bytes()), &sc)))
}

/// align_adapter for every (region, adapter) pair: result[r][a], None where adapters[a] is None.
#[pyfunction]
fn align_regions(py: Python<'_>, regions: Vec<String>, adapters: Vec<Option<String>>, scoring_scheme_vals: Vec<i32>, threads: usize) -> PyResult<Vec<Vec<Option<Hit>>>> {
    let sc = scoring(&scoring_scheme_vals)?;
    let pool = pool(threads)?;
    Ok(py.detach(|| {
        let regs: Vec<Vec<u8>> = regions.iter().map(|r| encode_dna5(r.as_bytes())).collect();
        let ads: Vec<Option<Vec<u8>>> = adapters.iter().map(|a| a.as_ref().map(|s| encode_dna5(s.as_bytes()))).collect();
        pool.install(|| {
            regs.par_chunks(256)
                .flat_map_iter(|chunk| {
                    let refs: Vec<&[u8]> = chunk.iter().map(|r| r.as_slice()).collect();
                    let mut out: Vec<Vec<Option<Hit>>> = vec![vec![None; ads.len()]; chunk.len()];
                    for (ai, a) in ads.iter().enumerate() {
                        if let Some(a) = a {
                            for (ri, h) in align_adapter_many(&refs, a, &sc).into_iter().enumerate() {
                                out[ri][ai] = Some(hit(h));
                            }
                        }
                    }
                    out
                })
                .collect()
        })
    }))
}

/// The alignment loop of NanoporeRead.find_middle_adapters for each sequence: the successive
/// hits (adapter index, full identity, read start, read end) with full identity >= threshold,
/// masking each hit before the next alignment. Alignments that provably cannot reach the
/// threshold are skipped.
#[pyfunction]
fn find_middle_hits(py: Python<'_>, seqs: Vec<String>, adapters: Vec<String>, scoring_scheme_vals: Vec<i32>, threshold: f64, threads: usize) -> PyResult<Vec<Vec<MiddleHit>>> {
    let sc = scoring(&scoring_scheme_vals)?;
    let pool = pool(threads)?;
    Ok(py.detach(|| {
        let ads: Vec<Vec<u8>> = adapters.iter().map(|a| encode_dna5(a.as_bytes())).collect();
        let pats: Vec<&[u8]> = ads.iter().map(|a| a.as_slice()).collect();
        pool.install(|| {
            seqs.par_iter()
                .map(|s| {
                    let mut masked = encode_dna5(s.as_bytes());
                    // Masking only turns bases into N, which cannot lower the edit distance of N-free adapters.
                    let eds = myers_many(&pats, &masked);
                    let mut hits = Vec::new();
                    for (ai, a) in ads.iter().enumerate() {
                        let m = a.len();
                        if m > 0 && m <= 64 && !a.contains(&4) && (100.0 * m as f64 / (m + eds[ai]) as f64) < threshold - 1e-6 {
                            continue;
                        }
                        let mut first = true;
                        loop {
                            if !first && identity_upper_bound(a, &masked) < threshold - 1e-6 {
                                break;
                            }
                            first = false;
                            let h = align_one(&masked, a, &sc);
                            if h.full_adapter_identity < threshold || h.read_end <= h.read_start {
                                break;
                            }
                            masked[h.read_start as usize..h.read_end as usize].fill(4); // '-' is N in Dna5
                            hits.push((ai, h.full_adapter_identity, h.read_start, h.read_end));
                        }
                    }
                    hits
                })
                .collect()
        })
    }))
}

/// consensus.all_vs_all_matrix: compatibility flags, -1 on the diagonal.
#[pyfunction]
fn compatibility_matrix(py: Python<'_>, sequences: Vec<String>) -> Vec<Vec<i32>> {
    py.detach(|| {
        let n = sequences.len();
        let mut mat = vec![vec![-1; n]; n];
        for i in 0..n {
            for j in i + 1..n {
                let f = check_compatibility(sequences[i].as_bytes(), sequences[j].as_bytes());
                mat[i][j] = f;
                mat[j][i] = f;
            }
        }
        mat
    })
}

/// stdout of `msa_consensus` on a FASTA file holding each sequence `count` times, in order.
#[pyfunction]
fn msa_consensus(py: Python<'_>, sequences: Vec<(String, u64)>) -> String {
    py.detach(|| {
        // msa_consensus reads a DnaString: anything not ACGT becomes A.
        let seqs: Vec<(Vec<u8>, u64)> = sequences.iter().map(|(s, c)| (s.bytes().map(|b| DNA4[b as usize]).collect(), *c)).collect();
        format!("{}\n", msa::msa_consensus(&seqs))
    })
}

/// One `approx_counter` call (see execApproxCounter): writes `<out_file>_<i>.start/.end`.
#[pyfunction]
#[pyo3(signature = (input, config, out_file, runs, verbosity, threads, seed=None))]
#[allow(clippy::too_many_arguments)]
fn approx_counter(py: Python<'_>, input: String, config: String, out_file: String, runs: usize, verbosity: i32, threads: usize, seed: Option<u64>) -> PyResult<()> {
    let pool = pool(threads)?;
    py.detach(|| {
        pool.install(|| {
            let mut out = std::io::stdout().lock();
            counter::approx_counter(&input, &config, &out_file, runs, verbosity, threads, seed, &mut out);
            let _ = out.flush();
        })
    });
    Ok(())
}

#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(align_adapter, m)?)?;
    m.add_function(wrap_pyfunction!(align_regions, m)?)?;
    m.add_function(wrap_pyfunction!(find_middle_hits, m)?)?;
    m.add_function(wrap_pyfunction!(compatibility_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(msa_consensus, m)?)?;
    m.add_function(wrap_pyfunction!(approx_counter, m)?)?;
    Ok(())
}
