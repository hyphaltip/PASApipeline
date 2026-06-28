/// Tier 2 Optimization: Lobject with bitset-based containment tracking.
///
/// The original C++ implementation uses `vector<bool>` which is a packed bit
/// representation but with poor cache locality and no SIMD exploitation. This
/// Rust implementation uses a custom bitset backed by `Vec<u64>` for:
///   - O(1) set/clear per bit
///   - O(n/64) intersection counting (64 bits per word)
///   - Better cache locality (contiguous u64 words)
///
/// The key algorithmic improvement is in `num_unique_contained`: instead of
/// iterating bit-by-bit O(n), we XOR the two bitsets and popcount the result,
/// giving O(n/64) with hardware popcount instructions.
pub struct Lobject {
    pub index: usize,
    pub num_alignments: usize,

    /// Bitset of which alignments are contained within this Lobject.
    /// Stored as a Vec<u64> for efficient bitwise operations.
    contained_cdna_indices: Vec<u64>,

    pub lscore_f: i32,
    pub lscore_r: i32,
    pub combined_score: i32,

    /// Index into the Lobjects vector (replaces raw pointers in C++).
    pub to_lptr: Option<usize>,
    pub from_lptr: Option<usize>,

    pub trace_indices: Vec<usize>,
    pub num_contained_indices: i32,
}

impl Lobject {
    pub fn new(index: usize, num_alignments: usize) -> Self {
        let num_words = (num_alignments + 63) / 64;
        Self {
            index,
            num_alignments,
            contained_cdna_indices: vec![0u64; num_words],
            lscore_f: 0,
            lscore_r: 0,
            combined_score: 0,
            to_lptr: None,
            from_lptr: None,
            trace_indices: Vec::new(),
            num_contained_indices: 0,
        }
    }

    pub fn set_contained_indices(&mut self, indices: &[usize]) {
        self.contained_cdna_indices.fill(0);
        self.num_contained_indices = 0;
        self.lscore_f = 0;
        self.lscore_r = 0;

        for &i in indices {
            let word = i / 64;
            let bit = 1u64 << (i % 64);
            self.contained_cdna_indices[word] |= bit;
            self.lscore_f += 1;
            self.lscore_r += 1;
            self.num_contained_indices += 1;
        }
    }

    /// Returns the number of alignments contained in `self` but NOT in `other`.
    ///
    /// Tier 2 optimization: Uses bitwise operations to count in O(n/64) instead
    /// of O(n) per-element iteration. Each 64-bit word is processed with:
    ///   `self_word & !other_word` → bits in self but not in other
    ///   `.count_ones()` → hardware popcount
    pub fn num_unique_contained(&self, other: &Lobject) -> i32 {
        let mut num = 0i32;
        let min_words = self.contained_cdna_indices.len().min(other.contained_cdna_indices.len());

        for i in 0..min_words {
            let diff = self.contained_cdna_indices[i] & !other.contained_cdna_indices[i];
            num += diff.count_ones() as i32;
        }

        // Handle case where self has more words than other
        for i in min_words..self.contained_cdna_indices.len() {
            num += self.contained_cdna_indices[i].count_ones() as i32;
        }

        num
    }

    pub fn get_contained_bit(&self, i: usize) -> bool {
        let word = i / 64;
        let bit = 1u64 << (i % 64);
        if word < self.contained_cdna_indices.len() {
            self.contained_cdna_indices[word] & bit != 0
        } else {
            false
        }
    }

    pub fn set_trace_indices(&mut self, traces: Vec<usize>) {
        self.trace_indices = traces;
    }

    pub fn get_trace_indices(&self) -> &[usize] {
        &self.trace_indices
    }

    pub fn to_string(&self) -> String {
        format!(
            "Lobject index: [{}] has LscoreF: {}, LscoreR: {}, combinedScore: {}\n",
            self.index, self.lscore_f, self.lscore_r, self.combined_score
        )
    }
}
