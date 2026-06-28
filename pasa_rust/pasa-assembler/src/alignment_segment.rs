use crate::coordset::CoordSet;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SegmentType {
    First,
    Last,
    Internal,
    Single,
    Unknown,
}

#[derive(Clone, Debug)]
pub struct AlignmentSegment {
    pub coords: CoordSet,
    pub has_left_splice_junction: bool,
    pub has_right_splice_junction: bool,
    pub seg_type: SegmentType,
}

impl AlignmentSegment {
    pub fn new(lend: i32, rend: i32) -> Self {
        Self {
            coords: CoordSet::new(lend, rend),
            has_left_splice_junction: false,
            has_right_splice_junction: false,
            seg_type: SegmentType::Unknown,
        }
    }

    pub fn from_coords(coords: CoordSet) -> Self {
        Self {
            coords,
            has_left_splice_junction: false,
            has_right_splice_junction: false,
            seg_type: SegmentType::Unknown,
        }
    }

    pub fn get_coords(&self) -> &CoordSet {
        &self.coords
    }

    pub fn set_coords(&mut self, lend: i32, rend: i32) {
        self.coords = CoordSet::new(lend, rend);
    }

    pub fn set_type(&mut self, t: SegmentType) {
        self.seg_type = t;
    }

    pub fn get_type(&self) -> &SegmentType {
        &self.seg_type
    }

    pub fn set_left_splice_junction(&mut self, val: bool) {
        self.has_left_splice_junction = val;
    }

    pub fn get_left_splice_junction(&self) -> bool {
        self.has_left_splice_junction
    }

    pub fn set_right_splice_junction(&mut self, val: bool) {
        self.has_right_splice_junction = val;
    }

    pub fn get_right_splice_junction(&self) -> bool {
        self.has_right_splice_junction
    }

    pub fn to_string(&self) -> String {
        format!("{}-{}", self.coords.lend, self.coords.rend)
    }
}
