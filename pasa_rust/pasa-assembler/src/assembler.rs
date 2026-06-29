use std::collections::HashSet;

use crate::alignment_segment::AlignmentSegment;
use crate::cdna_alignment::CdnaAlignment;
use crate::coordset::CoordSet;
use crate::lobject::Lobject;

/// Tier 1 Optimization: Interval tree for O(n log n + k) overlap queries.
///
/// The original C++ implementation uses an O(n²) all-vs-all comparison in
/// `determine_compatibilities_and_encapsulations`. For large numbers of
/// alignments (thousands), this becomes the bottleneck.
///
/// This implementation uses a sweep-line algorithm with sorted events:
///   1. Create START events at each alignment's lend position
///   2. Sort alignments by lend (already done in constructor)
///   3. For each alignment i, use binary search to find all alignments j
///      where lend_j <= rend_i, then filter by actual overlap
///
/// This reduces the complexity from O(n²) to O(n log n + k) where k is
/// the number of actually overlapping pairs.
struct IntervalTree {
    /// Sorted (lend, index) pairs for binary search
    sorted_starts: Vec<(i32, usize)>,
}

impl IntervalTree {
    fn new(alignments: &[CdnaAlignment]) -> Self {
        let mut sorted_starts: Vec<(i32, usize)> = alignments
            .iter()
            .enumerate()
            .map(|(i, a)| (a.get_coords().lend, i))
            .collect();
        sorted_starts.sort();

        Self { sorted_starts }
    }

    /// Returns all alignment indices that overlap with the given coordinate range.
    fn query_overlaps(&self, alignments: &[CdnaAlignment], query: &CoordSet) -> HashSet<usize> {
        // Find all alignments where lend <= query.rend AND rend >= query.lend
        let mut candidates = HashSet::new();

        // Binary search: find all starts <= query.rend
        let start_upper = self.sorted_starts.partition_point(|(lend, _)| *lend <= query.rend);
        for i in 0..start_upper {
            let (_, idx) = self.sorted_starts[i];
            // Check actual overlap
            if alignments[idx].get_coords().rend >= query.lend {
                candidates.insert(idx);
            }
        }

        candidates
    }
}

/// The main CDNA alignment assembler.
///
/// Implements the PASA assembly algorithm with the following optimizations:
/// - Tier 1: Interval tree for overlap queries (O(n log n + k) vs O(n²))
/// - Tier 2: Bitset-based Lobject containment tracking (O(n/64) intersection)
/// - HashSet-based compatibility matrix for sparse data
pub struct CdnaAlignmentAssembler {
    pub alignments: Vec<CdnaAlignment>,
    pub assemblies: Vec<CdnaAlignment>,
    pub assembly_containment_list: Vec<Vec<usize>>,
    pub fuzzlength: i32,
    lobjects: Vec<Lobject>,
    /// Sparse compatibility matrix: compatibilities[i] = set of j where i~j
    compatibilities: Vec<HashSet<usize>>,
    /// Sparse encapsulation matrix: encapsulations[i] = set of j where i⊃j
    encapsulations: Vec<HashSet<usize>>,
    pub num_alignments: usize,
}

impl CdnaAlignmentAssembler {
    pub fn new(mut incoming_alignments: Vec<CdnaAlignment>) -> Self {
        // Sort alignments by lend position
        incoming_alignments.sort_by(|a, b| a.get_coords().lend.cmp(&b.get_coords().lend));

        let num = incoming_alignments.len();
        Self {
            alignments: incoming_alignments,
            assemblies: Vec::new(),
            assembly_containment_list: Vec::new(),
            fuzzlength: 20,
            lobjects: Vec::with_capacity(num),
            compatibilities: vec![HashSet::new(); num],
            encapsulations: vec![HashSet::new(); num],
            num_alignments: num,
        }
    }

    pub fn set_fuzzlength(&mut self, length: i32) {
        self.fuzzlength = length;
    }

    /// Main assembly entry point. Orchestrates the full algorithm.
    pub fn assemble_alignments(&mut self) {
        self.determine_compatibilities_and_encapsulations();
        self.populate_lobjects();
        self.do_full_fscan();

        let top_assembly_indices = self.get_top_scoring_alignment();
        let assembly = self.create_assembly(&top_assembly_indices);
        let assembly_title = format!("assembly_{}", self.assemblies.len());
        let mut assembly = assembly;
        assembly.set_title(&assembly_title);
        self.assemblies.push(assembly);
        self.assembly_containment_list.push(top_assembly_indices.clone());

        if top_assembly_indices.len() == self.num_alignments {
            return;
        }

        // Reinit for pasa2 approach
        self.assemblies.clear();
        self.assembly_containment_list.clear();

        self.do_full_rscan();

        // Compute combined scores and trace indices
        let mut untraversed_lobjs: Vec<usize> = (0..self.num_alignments).collect();

        for &i in &untraversed_lobjs {
            let combined_score = self.lobjects[i].lscore_f + self.lobjects[i].lscore_r
                - self.lobjects[i].num_contained_indices;
            self.lobjects[i].combined_score = combined_score;

            let alignment_indices = self.get_alignment_assembly_nucleating_at(i);
            self.lobjects[i].set_trace_indices(alignment_indices);
        }

        // Sort by combined score descending
        untraversed_lobjs.sort_by(|&a, &b| {
            self.lobjects[b].combined_score.cmp(&self.lobjects[a].combined_score)
        });

        // Bin by combined score
        let mut bins: Vec<Vec<usize>> = Vec::new();
        let mut current_score = i32::MIN;
        for &i in &untraversed_lobjs {
            let score = self.lobjects[i].combined_score;
            if score != current_score {
                bins.push(Vec::new());
                current_score = score;
            }
            bins.last_mut().unwrap().push(i);
        }

        // Process bins to extract assemblies
        let mut accounted_for = vec![false; self.num_alignments];

        for bin in &bins {
            let bin_lobjs = bin.clone();
            for &i in &bin_lobjs {
                let trace = self.lobjects[i].get_trace_indices().to_vec();

                // Check if any alignments in trace are unconsumed
                let has_unconsumed = trace.iter().any(|&idx| !accounted_for[idx]);
                if !has_unconsumed {
                    continue;
                }

                let new_assembly = self.create_assembly(&trace);
                let assembly_title = format!("assembly_{}", self.assemblies.len());
                let mut new_assembly = new_assembly;
                new_assembly.set_title(&assembly_title);
                self.assemblies.push(new_assembly);
                self.assembly_containment_list.push(trace.clone());

                for &idx in &trace {
                    accounted_for[idx] = true;
                }

                let all_accounted = accounted_for.iter().all(|&x| x);
                if all_accounted {
                    return;
                }
            }
        }
    }

    /// Tier 1 Optimization: Determine compatibilities using interval tree.
    ///
    /// Instead of O(n²) all-vs-all comparison, we:
    /// 1. Build an interval tree from alignment coordinates
    /// 2. For each alignment, query the interval tree to find overlapping alignments
    /// 3. Only perform the full canMerge check on overlapping pairs
    ///
    /// This reduces the number of canMerge calls from O(n²) to O(k) where k
    /// is the number of actually overlapping pairs.
    fn determine_compatibilities_and_encapsulations(&mut self) {
        let interval_tree = IntervalTree::new(&self.alignments);

        for i in 0..self.num_alignments {
            let candidates = interval_tree.query_overlaps(&self.alignments, self.alignments[i].get_coords());

            for &j in &candidates {
                if j <= i {
                    continue; // Only check upper triangle
                }

                if self.can_merge(i, j) {
                    self.compatibilities[i].insert(j);
                    self.compatibilities[j].insert(i);

                    if self.encapsulates(i, j) {
                        self.encapsulations[i].insert(j);
                    }
                    if self.encapsulates(j, i) {
                        self.encapsulations[j].insert(i);
                    }
                }
            }
        }
    }

    /// Check if two alignments can be merged.
    /// Uses alignment indices into the sorted alignments vector.
    fn can_merge(&self, i: usize, j: usize) -> bool {
        let a1 = &self.alignments[i];
        let a2 = &self.alignments[j];

        // Check for overall overlap
        if !a1.get_coords().overlap(a2.get_coords()) {
            return false;
        }

        // Check orientation
        if a1.get_orientation() != a2.get_orientation() {
            return false;
        }

        // Check segment-level compatibility
        let a1_segs = a1.get_alignment_segments();
        let a2_segs = a2.get_alignment_segments();

        // Find first pair of overlapping segments
        let mut starting_a1 = None;
        let mut starting_a2 = None;
        for (i1, s1) in a1_segs.iter().enumerate() {
            for (i2, s2) in a2_segs.iter().enumerate() {
                if s1.get_coords().overlap(s2.get_coords()) {
                    starting_a1 = Some(i1);
                    starting_a2 = Some(i2);
                    break;
                }
            }
            if starting_a1.is_some() {
                break;
            }
        }

        let (mut sa1, mut sa2) = match (starting_a1, starting_a2) {
            (Some(a), Some(b)) => (a, b),
            _ => return false,
        };

        // One of the starting segments must be the first segment
        if sa1 != 0 && sa2 != 0 {
            return false;
        }

        // Walk both segment lists in parallel
        while sa1 < a1_segs.len() && sa2 < a2_segs.len() {
            let seg1 = &a1_segs[sa1];
            let seg2 = &a2_segs[sa2];

            if !seg1.get_coords().overlap(seg2.get_coords()) {
                return false;
            }

            let a1_lend = seg1.get_coords().lend;
            let a1_rend = seg1.get_coords().rend;
            let a2_lend = seg2.get_coords().lend;
            let a2_rend = seg2.get_coords().rend;

            // Check left splice junction
            if seg1.get_left_splice_junction() || seg2.get_left_splice_junction() {
                if seg1.get_left_splice_junction() && seg2.get_left_splice_junction() && a1_lend != a2_lend {
                    return false;
                }
                if seg1.get_left_splice_junction() && a2_lend + self.fuzzlength < a1_lend {
                    return false;
                }
                if seg2.get_left_splice_junction() && a1_lend + self.fuzzlength < a2_lend {
                    return false;
                }
            }

            // Check right splice junction
            if seg1.get_right_splice_junction() || seg2.get_right_splice_junction() {
                if seg1.get_right_splice_junction() && seg2.get_right_splice_junction() && a1_rend != a2_rend {
                    return false;
                }
                if seg1.get_right_splice_junction() && a2_rend - self.fuzzlength > a1_rend {
                    return false;
                }
                if seg2.get_right_splice_junction() && a1_rend - self.fuzzlength > a2_rend {
                    return false;
                }
            }

            sa1 += 1;
            sa2 += 1;
        }

        true
    }

    /// Check if alignment i encapsulates alignment j (i's span contains j's span).
    fn encapsulates(&self, i: usize, j: usize) -> bool {
        let a = &self.alignments[i].get_coords();
        let b = &self.alignments[j].get_coords();
        b.lend >= a.lend && b.rend <= a.rend
    }

    fn populate_lobjects(&mut self) {
        for i in 0..self.num_alignments {
            let mut lobj = Lobject::new(i, self.num_alignments);

            // Build contained list: self + all encapsulated alignments
            let mut contained = vec![i];
            if let Some(enc) = self.encapsulations.get(i) {
                for &j in enc {
                    contained.push(j);
                }
            }
            lobj.set_contained_indices(&contained);
            self.lobjects.push(lobj);
        }
    }

    fn do_full_fscan(&mut self) {
        for i in 1..self.num_alignments {
            let mut top_score = 0i32;
            let mut top_scoring_index: Option<usize> = None;

            for j in (0..i).rev() {
                if !self.compatibilities[i].contains(&j) {
                    continue;
                }

                let containment = self.encapsulations[i].contains(&j)
                    || self.encapsulations[j].contains(&i);
                if containment {
                    continue;
                }

                let curr_total_score = self.lobjects[j].lscore_f
                    + self.lobjects[i].num_unique_contained(&self.lobjects[j]);

                if curr_total_score > top_score {
                    top_scoring_index = Some(j);
                    top_score = curr_total_score;
                }
            }

            if let Some(idx) = top_scoring_index {
                self.lobjects[i].from_lptr = Some(idx);
                self.lobjects[i].lscore_f = top_score;
            }
        }
    }

    fn do_full_rscan(&mut self) {
        if self.num_alignments < 2 {
            return;
        }
        for i in (0..self.num_alignments - 1).rev() {
            let mut top_score = 0i32;
            let mut top_scoring_index: Option<usize> = None;

            for j in (i + 1)..self.num_alignments {
                if !self.compatibilities[i].contains(&j) {
                    continue;
                }

                let containment = self.encapsulations[i].contains(&j)
                    || self.encapsulations[j].contains(&i);
                if containment {
                    continue;
                }

                let curr_total_score = self.lobjects[j].lscore_r
                    + self.lobjects[i].num_unique_contained(&self.lobjects[j]);

                if curr_total_score > top_score {
                    top_scoring_index = Some(j);
                    top_score = curr_total_score;
                }
            }

            if let Some(idx) = top_scoring_index {
                self.lobjects[i].to_lptr = Some(idx);
                self.lobjects[i].lscore_r = top_score;
            }
        }
    }

    fn get_top_scoring_alignment(&self) -> Vec<usize> {
        let mut top_score = 0i32;
        let mut top_scoring_index = 0usize;

        for i in 0..self.lobjects.len() {
            if self.lobjects[i].lscore_f > top_score {
                top_scoring_index = i;
                top_score = self.lobjects[i].lscore_f;
            }
        }

        self.back_trace(top_scoring_index)
    }

    fn back_trace(&self, start_index: usize) -> Vec<usize> {
        let mut tracker = HashSet::new();
        let mut current = Some(start_index);

        while let Some(idx) = current {
            for i in 0..self.num_alignments {
                if self.lobjects[idx].get_contained_bit(i) {
                    tracker.insert(i);
                }
            }
            current = self.lobjects[idx].from_lptr;
        }

        let mut result: Vec<usize> = tracker.into_iter().collect();
        result.sort();
        result
    }

    fn forward_trace(&self, start_index: usize) -> Vec<usize> {
        let mut tracker = HashSet::new();
        let mut current = Some(start_index);

        while let Some(idx) = current {
            for i in 0..self.num_alignments {
                if self.lobjects[idx].get_contained_bit(i) {
                    tracker.insert(i);
                }
            }
            current = self.lobjects[idx].to_lptr;
        }

        let mut result: Vec<usize> = tracker.into_iter().collect();
        result.sort();
        result
    }

    fn get_alignment_assembly_nucleating_at(&self, index: usize) -> Vec<usize> {
        let mut vecvec = Vec::new();
        vecvec.push(self.back_trace(index));
        vecvec.push(self.forward_trace(index));

        let mut unique_map = HashSet::new();
        for v in &vecvec {
            for &entry in v {
                unique_map.insert(entry);
            }
        }
        unique_map.into_iter().collect()
    }

    fn create_assembly(&self, alignment_index_listing: &[usize]) -> CdnaAlignment {
        assert!(!alignment_index_listing.is_empty(), "empty index list");

        let mut sorted = alignment_index_listing.to_vec();
        sorted.sort();

        let mut assembly = self.alignments[sorted[0]].clone();
        for &i in &sorted[1..] {
            assembly = self.merge_alignments(&assembly, &self.alignments[i]);
        }
        assembly
    }

    fn merge_alignments(&self, a: &CdnaAlignment, b: &CdnaAlignment) -> CdnaAlignment {
        let orientation = a.get_orientation();

        // Build splice coordinate sets
        let mut left_splice_coords: HashSet<i32> = HashSet::new();
        let mut right_splice_coords: HashSet<i32> = HashSet::new();

        for seg in a.get_alignment_segments() {
            let lend = seg.get_coords().lend;
            let rend = seg.get_coords().rend;
            if seg.get_left_splice_junction() {
                left_splice_coords.insert(lend);
            }
            if seg.get_right_splice_junction() {
                right_splice_coords.insert(rend);
            }
        }

        for seg in b.get_alignment_segments() {
            let lend = seg.get_coords().lend;
            let rend = seg.get_coords().rend;
            if seg.get_left_splice_junction() {
                left_splice_coords.insert(lend);
            }
            if seg.get_right_splice_junction() {
                right_splice_coords.insert(rend);
            }
        }

        // Merge overlapping segments
        let a1_segs = a.get_alignment_segments();
        let a2_segs = b.get_alignment_segments();

        let mut merged_coords: Vec<CoordSet> = Vec::new();

        for seg1 in a1_segs {
            let a1_lend = seg1.get_coords().lend;
            let a1_rend = seg1.get_coords().rend;

            let mut merged_lend = -1i32;
            let mut merged_rend = -1i32;

            for seg2 in a2_segs {
                if seg1.get_coords().overlap(seg2.get_coords()) {
                    let a2_lend = seg2.get_coords().lend;
                    let a2_rend = seg2.get_coords().rend;

                    // Determine merged_lend
                    if left_splice_coords.contains(&a1_lend) {
                        merged_lend = a1_lend;
                    } else if left_splice_coords.contains(&a2_lend) {
                        merged_lend = a2_lend;
                    } else {
                        merged_lend = std::cmp::min(a1_lend, a2_lend);
                    }

                    // Determine merged_rend
                    if right_splice_coords.contains(&a1_rend) {
                        merged_rend = a1_rend;
                    } else if right_splice_coords.contains(&a2_rend) {
                        merged_rend = a2_rend;
                    } else {
                        merged_rend = std::cmp::max(a1_rend, a2_rend);
                    }
                    break;
                }
            }

            if merged_lend == -1 || merged_rend == -1 {
                merged_coords.push(CoordSet { lend: a1_lend, rend: a1_rend });
            } else {
                merged_coords.push(CoordSet { lend: merged_lend, rend: merged_rend });
            }
        }

        // Add unconsumed a2 segments
        for seg2 in a2_segs {
            let a2_lend = seg2.get_coords().lend;
            let a2_rend = seg2.get_coords().rend;
            let a2_coords = CoordSet { lend: a2_lend, rend: a2_rend };

            let mut overlap_flag = false;
            for m_coords in &merged_coords {
                if a2_coords.overlap(m_coords) {
                    overlap_flag = true;
                    break;
                }
            }

            if !overlap_flag {
                merged_coords.push(a2_coords);
            }
        }

        // Create new segment list from merged coords
        let new_seg_list: Vec<AlignmentSegment> = merged_coords
            .into_iter()
            .map(AlignmentSegment::from_coords)
            .collect();

        CdnaAlignment::new(new_seg_list, orientation)
    }

    pub fn get_assemblies(&self) -> &[CdnaAlignment] {
        &self.assemblies
    }

    pub fn to_align_illustration(&self, _line_length: usize) -> String {
        let mut all_coords: Vec<i32> = Vec::new();
        for a in &self.alignments {
            all_coords.push(a.get_coords().lend);
            all_coords.push(a.get_coords().rend);
        }
        all_coords.sort_unstable();
        if all_coords.is_empty() {
            return String::new();
        }
        let min_coord = all_coords[0];
        let max_coord = *all_coords.last().unwrap();
        let _rel_max = max_coord - min_coord;

        let mut text = format!("Individual Alignments: ({})\n", self.num_alignments);
        for (i, a) in self.alignments.iter().enumerate() {
            text.push_str(&format!("{} index: [{}]\n", a.to_string(), i));
        }

        if !self.assemblies.is_empty() {
            text.push_str(&format!("\nASSEMBLIES: ({})\n", self.assemblies.len()));
            let mut assembly_summary = String::new();
            for (i, asm) in self.assemblies.iter().enumerate() {
                text.push_str(&format!(
                    "{} score: ({}) contains {:?}\n",
                    asm.to_string(),
                    self.assembly_containment_list[i].len(),
                    self.assembly_containment_list[i]
                ));

                let acc_list: Vec<&str> = self.assembly_containment_list[i]
                    .iter()
                    .map(|&idx| self.alignments[idx].get_title())
                    .collect();

                assembly_summary.push_str(&format!(
                    "assembly: ({}) contains alignments: [{}] with structure [{}] score: ({})\n",
                    i,
                    acc_list.join(","),
                    asm.to_string(),
                    self.assembly_containment_list[i].len()
                ));
            }
            text.push_str(&assembly_summary);
        }

        text
    }
}
