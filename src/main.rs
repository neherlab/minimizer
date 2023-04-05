// #![warn(dead_code, unused_imports, unused_variables, deprecated_in_future)]

use bio::io::fasta;
use bio::io::fasta::Records;
use std::fs::File;
use std::io::BufReader;

const REFERENCE_PATH: &str = "data/sc2_reference.fasta";
// const QUERY_PATH: &str = "data/sc2_query_long.fasta";
const QUERY_PATH: &str = "/Users/corneliusromer/Downloads/gisaid_hcov-19_2023_04_05_18.fasta";

fn read_fasta(filename: &str) -> Records<BufReader<File>> {
    let path = std::path::Path::new(filename);
    let reader = fasta::Reader::new(File::open(path).unwrap());

    reader.records()
}

// Don't inline
#[inline(never)]
fn invertible_hash(x: u32, p: u32) -> u32 {
    let m: u32 = 2u32.pow(p) - 1;
    let mut x: u32 = (!x + (x << 21)) & m;
    x = x ^ (x >> 24);
    x = (x + (x << 3) + (x << 8)) & m;
    x = x ^ (x >> 14);
    x = (x + (x << 2) + (x << 4)) & m;
    x = x ^ (x >> 28);
    x = (x + (x << 31)) & m;
    x
}

// Don't inline
#[inline(never)]
fn hash_kmer(kmer: &[u8]) -> u32 {
    let mut hash: u32 = 0;
    for &base in kmer {
        hash <<= 2;
        hash += match base {
            b'G' => 0,
            b'C' => 1,
            b'T' => 2,
            b'A' => 3,
            _ => return u32::MAX,
        };
    }
    invertible_hash(hash, 32)
}

struct MinimizerQueue {
    queue: Vec<(u32, usize)>,
    w: usize,
    start: usize, // Start of queu
    end: usize,   // End of queue
    wrap: bool,
}

impl MinimizerQueue {
    fn new(w: usize) -> Self {
        Self {
            queue: vec![(u32::MAX, 0); w],
            w,
            start: 0,
            end: 1,
            wrap: false,
        }
    }

    // Don't inline
    #[inline(never)]
    fn insert(&mut self, hash: u32, pos: usize) -> Option<(u32, usize)> {
        // Insert hash at appropriate position
        // Search linearly for first hash that is larger
        // Delete everything after it
        // Insert it there
        // Need to wrap around
        let mut new_min = false;
        // Insert new hash and wipe everything after it
        let mut early_exit = false;
        if !self.wrap {
            for i in self.start..self.end {
                if hash <= self.queue[i].0 {
                    self.queue[i] = (hash, pos);
                    self.end = i + 1;
                    if self.end == self.w {
                        self.wrap = true;
                        self.end = 0;
                    }
                    if i == self.start {
                        new_min = true;
                    }
                    early_exit = true;
                    break;
                }
            }
        } else {
            for i in self.start..self.end + self.w {
                let i = i % self.w;
                if hash <= self.queue[i].0 {
                    self.queue[i] = (hash, pos);
                    self.end = i + 1;
                    if self.end == self.w {
                        self.wrap = true;
                        self.end = 0;
                    }
                    if i == self.start {
                        new_min = true;
                    }
                    early_exit = true;
                    break;
                }
            }
        }
        if !early_exit {
            self.queue[self.end] = (hash, pos);
            self.end += 1;
            if self.end == self.w {
                self.wrap = true;
                self.end = 0;
            }
        }
        // Check if beginning needs to be wiped
        if self.queue[self.start].1 < pos.saturating_sub(self.w) {
            self.start += 1;
            if self.start == self.w {
                self.start = 0;
                self.wrap = false;
            }
            new_min = true;
        }
        if new_min {
            Some(self.queue[self.start])
        } else {
            None
        }
    }
}

fn get_minimizers(seq: &[u8], k: usize, w: usize) -> Vec<(u32, usize)> {
    // Set of minimizers
    let mut minimizers = Vec::new();

    let mut queue = MinimizerQueue::new(w);

    // Iterate over k-mers
    for window_start in 0..seq.len() - k {
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
        let new_hash = hash_kmer(&seq[window_start..window_start + k]);
        if let Some(new_minimizer) = queue.insert(new_hash, window_start) {
            minimizers.push(new_minimizer);
        }
        //
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

    let minimizers = get_minimizers(reference, 13, 1000);
    // println!("{:?}", minimizers);
    // println!("{} minimizers", minimizers.len());

    let query = read_fasta(QUERY_PATH);

    for record in query {
        let query_minimizers = get_minimizers(record.unwrap().seq(), 13, 1000);
        // println!("{:?}", query_minimizers);

        let mut matches = Vec::new();

        for (hash, pos) in query_minimizers {
            if let Ok(index) = minimizers.binary_search_by_key(&hash, |(hash, _)| *hash) {
                matches.push((pos, minimizers[index].1));
            }
        }
        matches.sort();
        // println!("{:?}", matches);
        // println!("{} matches", matches.len());
    }
}
