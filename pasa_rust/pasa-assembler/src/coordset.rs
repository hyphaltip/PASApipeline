use std::cmp::{max, min};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct CoordSet {
    pub lend: i32,
    pub rend: i32,
}

impl CoordSet {
    pub fn new(lend: i32, rend: i32) -> Self {
        if lend > rend {
            Self { lend: rend, rend: lend }
        } else {
            Self { lend, rend }
        }
    }

    pub fn overlap(&self, other: &CoordSet) -> bool {
        self.lend <= other.rend && self.rend >= other.lend
    }

    pub fn contains(&self, other: &CoordSet) -> bool {
        self.lend <= other.lend && self.rend >= other.rend
    }

    pub fn span(&self) -> i32 {
        self.rend - self.lend
    }

    pub fn merged(&self, other: &CoordSet) -> CoordSet {
        CoordSet {
            lend: min(self.lend, other.lend),
            rend: max(self.rend, other.rend),
        }
    }
}
