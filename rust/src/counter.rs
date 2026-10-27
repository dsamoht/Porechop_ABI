//! In-process replacement for the `approx_counter` program: same parameters, messages and
//! output files.

use crate::align::DNA5;
use crate::kmer::*;
use crate::rng::Rng;
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};

fn die(msg: &str) -> ! {
    eprintln!("{msg}");
    std::process::exit(1)
}

struct Params {
    k: usize,
    sl: usize,
    sn: usize,
    limit: usize,
    lc: f32,
    solid: u64,
    forbidden: Option<String>,
    exact_out: Option<String>,
}

fn parse_config(path: &str) -> Params {
    let text = std::fs::read_to_string(path).unwrap_or_else(|_| die(&format!("/!\\\t Invalide configuration file.\n/!\\\t {path}")));
    let mut map: HashMap<String, String> = HashMap::new();
    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }
        let line: String = line.chars().filter(|&c| c != ' ').collect();
        if let Some((a, b)) = line.split_once('=') {
            map.insert(a.to_string(), b.to_string());
        }
    }
    let gi = |k: &str, d: usize| map.get(k).and_then(|v| v.trim().parse::<f64>().ok()).map(|v| v as usize).unwrap_or(d);
    let p = Params {
        k: gi("k", 16),
        sl: gi("sl", 100),
        sn: gi("sn", 40000),
        limit: gi("lim", 500),
        lc: map.get("lc").and_then(|v| v.trim().parse::<f32>().ok()).unwrap_or(1.0),
        solid: gi("sk", 0) as u64,
        forbidden: map.get("fk").filter(|s| !s.is_empty()).cloned(),
        exact_out: map.get("e").filter(|s| !s.is_empty()).cloned(),
    };
    if p.k < 2 || p.k > 32 {
        die("/!\\ ERROR: kmer size must be between 2 and 32 (included)");
    }
    if p.k > p.sl {
        die("/!\\ ERROR: kmer size must be smaller than the sampling length (k <= sl)");
    }
    p
}

fn strip(b: &[u8]) -> &[u8] {
    let is_ws = |c: &u8| matches!(c, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c);
    let start = b.iter().position(|c| !is_ws(c)).unwrap_or(b.len());
    let end = b.iter().rposition(|c| !is_ws(c)).map(|p| p + 1).unwrap_or(start);
    &b[start..end.max(start)]
}

/// Calls `f` on each sequence of a FASTA or FASTQ file, plain or gzipped (Porechop's
/// parsing rules: the format is given by the first character of the file).
fn for_each_sequence(path: &str, mut f: impl FnMut(&[u8])) {
    let fail = |msg: String| -> ! { die(&msg) };
    let bad = || -> ! { die(&format!("\nError: {path} could not be parsed - is it formatted correctly?")) };
    let open = || File::open(path).unwrap_or_else(|e| fail(format!("Error: {e}")));
    let mut magic = [0u8; 3];
    let gz = open().read(&mut magic).unwrap_or(0) == 3 && magic == [0x1f, 0x8b, 0x08];
    let mut r: Box<dyn BufRead> = if gz { Box::new(BufReader::new(MultiGzDecoder::new(BufReader::new(open())))) } else { Box::new(BufReader::new(open())) };
    let mut line = Vec::new();
    let mut next = |line: &mut Vec<u8>| {
        line.clear();
        r.read_until(b'\n', line).unwrap_or_else(|e| fail(format!("Error reading {path}: {e}"))) > 0
    };
    if !next(&mut line) {
        bad();
    }
    match line[0] {
        b'@' => loop {
            if !strip(&line).is_empty() {
                // header, sequence, '+', qualities
                if !next(&mut line) {
                    bad();
                }
                f(strip(&line));
                if !next(&mut line) || !next(&mut line) {
                    bad();
                }
            }
            if !next(&mut line) {
                break;
            }
        },
        b'>' => {
            let mut seq: Option<Vec<u8>> = None;
            loop {
                let l = strip(&line);
                if l.first() == Some(&b'>') {
                    if let Some(s) = seq.take() {
                        f(&s);
                    }
                    seq = Some(Vec::new());
                } else if let Some(s) = seq.as_mut() {
                    s.extend_from_slice(l);
                }
                if !next(&mut line) {
                    break;
                }
            }
            if let Some(s) = seq {
                f(&s);
            }
        }
        _ => bad(),
    }
}

/// Dna5-encoded read starts and ends eligible for sampling (reads at least 2*sl long),
/// and the total number of reads.
fn load_read_ends(input: &str, sl: usize) -> (Vec<Vec<u8>>, Vec<Vec<u8>>, usize) {
    let enc = |s: &[u8]| s.iter().map(|&c| DNA5[c as usize]).collect::<Vec<u8>>();
    let (mut starts, mut ends, mut total) = (Vec::new(), Vec::new(), 0);
    for_each_sequence(input, |s| {
        total += 1;
        if s.len() >= 2 * sl {
            starts.push(enc(&s[..sl]));
            // approx_counter takes suffix(seq, len - 1 - cut): the last cut + 1 bases.
            ends.push(enc(&s[s.len() - 1 - sl..]));
        }
    });
    (starts, ends, total)
}

/// One approximate count of one read end: (k-mer, count) sorted like approx_counter's output.
fn count_run(p: &Params, sn: usize, pool: &[Vec<u8>], rng: &mut Rng, forbidden: &HashSet<Kmer>, label: &str, suffix: &str) -> Vec<(String, u64)> {
    let mut idx: Vec<usize> = (0..pool.len()).collect();
    rng.shuffle(&mut idx);
    idx.truncate(sn);
    let samples: Vec<&[u8]> = idx.iter().map(|&i| pool[i].as_slice()).collect();
    let (counts, had_n) = count_kmers(&samples, p.k, adjust_threshold(p.lc, 16, p.k), forbidden);
    if had_n > 0 {
        eprintln!(
            "/!\\ WARNING: This dataset contained sequences with 'N' symbols. /!\\ WARNING: Current implementation ignores k-mers containing 'N'./!\\ WARNING: A total of {had_n} k-mers were ignored."
        );
    }
    let mut first: Vec<(Kmer, u64)> = counts.into_iter().collect();
    if p.solid != 0 {
        first.retain(|x| x.1 >= p.solid);
    } else {
        sort_counts(&mut first, p.k, Some(p.limit));
    }
    if let Some(e) = &p.exact_out {
        let path = format!("{e}{suffix}.{label}");
        let text: String = first.iter().map(|(km, c)| format!("{}\t{}\n", kmer_to_string(*km, p.k), c)).collect();
        if std::fs::write(&path, text).is_err() {
            die(&format!("/!\\ ERROR: Failed to export exact k-mer count\nPath: {path}"));
        }
    }
    let kmers: Vec<Kmer> = first.iter().map(|x| x.0).collect();
    let approx = approximate_counts(&samples, &kmers, p.k);
    let mut res: Vec<(Kmer, u64)> = kmers.into_iter().zip(approx).collect();
    sort_counts(&mut res, p.k, Some(p.limit));
    res.into_iter().map(|(km, c)| (kmer_to_string(km, p.k), c)).collect()
}

/// `approx_counter <input> --config <config> -o <out_prefix> -mr <runs> -nt <threads>`:
/// `runs` independent samplings, written to `<out_prefix>_<i>.start` / `.end`. `threads` is
/// only reported; the caller sets up the rayon pool.
#[allow(clippy::too_many_arguments)]
pub fn approx_counter(input: &str, config: &str, out_prefix: &str, runs: usize, v: i32, threads: usize, seed: Option<u64>, out: &mut dyn Write) {
    let p = parse_config(config);
    let forbidden: HashSet<Kmer> = match &p.forbidden {
        Some(f) => {
            let text = std::fs::read_to_string(f).unwrap_or_else(|_| die("/!\\ ERROR: COULD NOT OPEN EXCLUDED KMER FILE, must quit"));
            text.lines().filter_map(|l| encode_kmer(l.as_bytes())).collect()
        }
        None => HashSet::new(),
    };
    let (starts, ends, total_reads) = load_read_ends(input, p.sl);
    if v > 0 {
        let _ = writeln!(out, "Kmer size:             {}", p.k);
        let _ = writeln!(out, "Sampled sequences:     {}", p.sn);
        let _ = writeln!(out, "Sampling length        {}", p.sl);
        let _ = writeln!(out, "LC filter threshold:   {}", p.lc);
        let _ = writeln!(out, "Adjusted LC threshold: {}", adjust_threshold(p.lc, 16, p.k));
        let _ = writeln!(out, "Nb thread:             {threads}");
        if p.solid != 0 {
            let _ = writeln!(out, "Solid kmers:           {}", p.solid);
        } else {
            let _ = writeln!(out, "Number of kept kmer:   {}", p.limit);
        }
        let _ = writeln!(out, "Number of runs:        {runs}");
        let _ = writeln!(out, "Verbosity level:       {v}");
        if runs > 1 {
            let _ = writeln!(out, "\nA total of {runs} runs will be performed.");
        }
    }
    for r in 0..runs {
        if v > 0 && runs > 1 {
            let _ = writeln!(out, "Starting run number {}", r + 1);
        }
        if r == 0 && p.sn > total_reads {
            let _ = out.flush();
            eprintln!("/!\\ WARNING: Sequence set too small for the requested sample size");
            eprintln!("/!\\ WARNING: The whole set will be used.");
        }
    }
    let _ = out.flush();
    let sn = p.sn.min(total_reads);
    let mut master = match seed {
        Some(s) => Rng::new(s),
        None => Rng::from_entropy(),
    };
    let seeds: Vec<u64> = (0..runs).map(|_| master.next_u64()).collect();
    seeds.par_iter().enumerate().for_each(|(i, &s)| {
        let mut rng = Rng::new(s);
        let suffix = format!("_{i}");
        let st = count_run(&p, sn, &starts, &mut rng, &forbidden, "start", &suffix);
        let en = count_run(&p, sn, &ends, &mut rng, &forbidden, "end", &suffix);
        for (counts, label) in [(st, "start"), (en, "end")] {
            let path = format!("{out_prefix}{suffix}.{label}");
            let text: String = counts.iter().map(|(k, c)| format!("{k}\t{c}\n")).collect();
            if std::fs::write(&path, text).is_err() {
                die(&format!("/!\\ ERROR: Failed to export approximate k-mer count\nPath: {path}"));
            }
        }
    });
}
