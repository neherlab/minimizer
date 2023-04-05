// #![warn(dead_code, unused_imports, unused_variables, deprecated_in_future)]

use bio::io::fasta;
use bio::io::fasta::Records;
use std::cmp::min;
use std::fs::File;
use std::io::BufReader;

const REFERENCE_PATH: &str = "data/sc2_reference.fasta";
const QUERY_PATH: &str = "data/sc2_query_long.fasta";

fn read_fasta(filename: &str) -> Records<BufReader<File>> {
    let path = std::path::Path::new(filename);
    let reader = fasta::Reader::new(File::open(path).unwrap());

    reader.records()
}

fn base_to_int(base: u8) -> Option<u8> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn hash_kmer(kmer: &[u8]) -> u64 {
    let mut hash: u64 = 0;
    for &base in kmer {
        if let Some(normalized_base) = base_to_int(base) {
            hash = hash.wrapping_mul(4).wrapping_add(normalized_base as u64);
        } else {
            return u64::MAX;
        }
    }
    hash
}

struct MinimizerQueue {
    queue: Vec<(u64, usize)>,
    w: usize,
    start: usize, // Start of queu
    end: usize,   // End of queue
}

impl MinimizerQueue {
    fn new(w: usize) -> Self {
        Self {
            queue: Vec::with_capacity(w),
            w,
            start: 0,
            end: 0,
        }
    }
}

fn get_minimizers(seq: &[u8], k: usize, w: usize) -> Vec<(u64, usize)> {
    // Set of minimizers
    let mut minimizers = Vec::new();

    let mut queue = MinimizerQueue::new(w);

    // Iterate over windows
    for window_start in 0..seq.len() - k - w {
        let mut minimum = (u64::MAX, 0);
        // Use dynamic programming
        // Step 1: Can cache the hash of each kmer start
        // Step 2: Cache minimum of each window
        // [0, 1, 1, 1, 3, 3, 3, 3, 5]
        // Only the preceding w kmer hashes are relevant
        // A hash minimum only extends up to w down
        // Keep all the minimums in a queue
        // Keep hash and position
        // Pop if position is outside of window
        // New minimum only added when smaller than miminum in queue
        // Queue is sorted by hash
        // But for every hash value, find position in queue and kick all larger hashes
        // 1. Calculate new kmer hash
        // 2. Place at appropriate position in queue (according to first item in tuple)
        // 3. Wipe everything after it (using second item in tuple)
        // 4. Check if beginning needs to be wiped
        // 5. Check if new hash is minimum
        // 6. If it is -> add to minimizers
        for kmer_start in window_start..window_start + w {
            let kmer = &seq[kmer_start..kmer_start + k];
            let hash = hash_kmer(kmer);
            if hash < minimum.0 {
                minimum = (hash, kmer_start);
            }
        }
        minimizers.push(minimum);
    }
    // Deduplicate minimizers
    minimizers.sort();
    minimizers.dedup();
    // minimizers.reverse();
    minimizers
}

fn main() {
    // Read reference and query sequences
    // Create minimizer index from reference
    // Create minimizers from query
    // Search minimizers in index
    // Report results
    let reference = read_fasta(REFERENCE_PATH).next().unwrap().unwrap();

    // Convert reference to bytes
    let reference = reference.seq();

    let minimizers = get_minimizers(reference, 10, 300);
    print!("{:?}", minimizers);

    let query = read_fasta(QUERY_PATH);

    for record in query {
        let query_minimizers = get_minimizers(record.unwrap().seq(), 10, 300);

        for (hash, pos) in query_minimizers {
            if let Ok(index) = minimizers.binary_search(&(hash, pos)) {
                println!("Found match at {}", index);
            }
        }
    }
}
