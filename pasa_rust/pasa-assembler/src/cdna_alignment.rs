use crate::alignment_segment::{AlignmentSegment, SegmentType};
use crate::coordset::CoordSet;

#[derive(Clone, Debug)]
pub struct CdnaAlignment {
    pub title: String,
    pub orient: char,
    pub alignment_segs: Vec<AlignmentSegment>,
    pub coords: CoordSet,
    pub num_segments: usize,
}

impl CdnaAlignment {
    pub fn new(mut seglist: Vec<AlignmentSegment>, orient: char) -> Self {
        assert!(!seglist.is_empty(), "empty segment list");
        Self {
            title: String::new(),
            orient,
            alignment_segs: Vec::new(),
            coords: CoordSet::default(),
            num_segments: 0,
        }
        .init(&mut seglist)
    }

    fn init(mut self, seglist: &mut Vec<AlignmentSegment>) -> Self {
        self.alignment_segs = std::mem::take(seglist);
        self.refine_alignment();
        self
    }

    fn refine_alignment(&mut self) {
        self.alignment_segs
            .sort_by_key(|s| s.coords.lend);

        if self.alignment_segs.is_empty() {
            return;
        }

        let min_lend = self.alignment_segs[0].coords.lend;
        let max_rend = self.alignment_segs.last().unwrap().coords.rend;
        self.coords = CoordSet { lend: min_lend, rend: max_rend };
        self.num_segments = self.alignment_segs.len();

        let n = self.alignment_segs.len();
        for i in 0..n {
            let seg_type = if n == 1 {
                SegmentType::Single
            } else if i == 0 {
                SegmentType::First
            } else if i == n - 1 {
                SegmentType::Last
            } else {
                SegmentType::Internal
            };
            self.alignment_segs[i].set_type(seg_type);
            match seg_type {
                SegmentType::First => self.alignment_segs[i].set_right_splice_junction(true),
                SegmentType::Last => self.alignment_segs[i].set_left_splice_junction(true),
                SegmentType::Internal => {
                    self.alignment_segs[i].set_left_splice_junction(true);
                    self.alignment_segs[i].set_right_splice_junction(true);
                }
                SegmentType::Single => {}
                SegmentType::Unknown => {}
            }
        }
    }

    pub fn set_title(&mut self, t: &str) {
        self.title = t.to_string();
    }

    pub fn get_title(&self) -> &str {
        &self.title
    }

    pub fn get_orientation(&self) -> char {
        self.orient
    }

    pub fn get_coords(&self) -> &CoordSet {
        &self.coords
    }

    pub fn get_alignment_segments(&self) -> &[AlignmentSegment] {
        &self.alignment_segs
    }

    pub fn get_alignment_segments_mut(&mut self) -> &mut Vec<AlignmentSegment> {
        &mut self.alignment_segs
    }

    pub fn to_string(&self) -> String {
        let segs: Vec<String> = self.alignment_segs.iter().map(|s| s.to_string()).collect();
        format!("{},{},{}", self.title, self.orient, segs.join(","))
    }
}
