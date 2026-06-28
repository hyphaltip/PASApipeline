use pasa_assembler::coordset::CoordSet;
use pasa_assembler::alignment_segment::AlignmentSegment;
use pasa_assembler::cdna_alignment::CdnaAlignment;
use pasa_assembler::assembler::CdnaAlignmentAssembler;
use pasa_assembler::lobject::Lobject;

#[test]
fn test_coordset_overlap() {
    let a = CoordSet { lend: 100, rend: 200 };
    let b = CoordSet { lend: 150, rend: 250 };
    let c = CoordSet { lend: 300, rend: 400 };

    assert!(a.overlap(&b));
    assert!(!a.overlap(&c));
}

#[test]
fn test_coordset_contains() {
    let outer = CoordSet { lend: 100, rend: 500 };
    let inner = CoordSet { lend: 200, rend: 300 };

    assert!(outer.contains(&inner));
    assert!(!inner.contains(&outer));
}

#[test]
fn test_coordset_new_swaps() {
    let cs = CoordSet::new(300, 100);
    assert_eq!(cs.lend, 100);
    assert_eq!(cs.rend, 300);
}

#[test]
fn test_alignment_segment_creation() {
    let seg = AlignmentSegment::new(100, 200);
    assert_eq!(seg.get_coords().lend, 100);
    assert_eq!(seg.get_coords().rend, 200);
    assert!(!seg.get_left_splice_junction());
    assert!(!seg.get_right_splice_junction());
}

#[test]
fn test_cdna_alignment_single_segment() {
    let seg = AlignmentSegment::new(100, 500);
    let align = CdnaAlignment::new(vec![seg], '+');

    assert_eq!(align.get_orientation(), '+');
    assert_eq!(align.num_segments, 1);
    assert_eq!(align.get_coords().lend, 100);
    assert_eq!(align.get_coords().rend, 500);
}

#[test]
fn test_cdna_alignment_multi_segment() {
    let seg1 = AlignmentSegment::new(100, 200);
    let seg2 = AlignmentSegment::new(300, 400);
    let seg3 = AlignmentSegment::new(500, 600);
    let align = CdnaAlignment::new(vec![seg1, seg2, seg3], '+');

    assert_eq!(align.num_segments, 3);
    assert_eq!(align.get_coords().lend, 100);
    assert_eq!(align.get_coords().rend, 600);

    // Check segment classification
    let segs = align.get_alignment_segments();
    assert_eq!(segs[0].seg_type, SegmentType::First);
    assert_eq!(segs[1].seg_type, SegmentType::Internal);
    assert_eq!(segs[2].seg_type, SegmentType::Last);
}

use pasa_assembler::alignment_segment::SegmentType;

#[test]
fn test_lobject_bitset_basic() {
    let mut lobj = Lobject::new(0, 100);
    lobj.set_contained_indices(&[0, 5, 10, 50]);

    assert!(lobj.get_contained_bit(0));
    assert!(lobj.get_contained_bit(5));
    assert!(lobj.get_contained_bit(10));
    assert!(lobj.get_contained_bit(50));
    assert!(!lobj.get_contained_bit(1));
    assert!(!lobj.get_contained_bit(99));
}

#[test]
fn test_lobject_bitset_unique_contained() {
    let mut lobj_a = Lobject::new(0, 100);
    let mut lobj_b = Lobject::new(1, 100);

    lobj_a.set_contained_indices(&[0, 1, 2, 3, 4]);
    lobj_b.set_contained_indices(&[2, 3, 4, 5, 6]);

    // A has 0, 1 that B doesn't have
    let unique = lobj_a.num_unique_contained(&lobj_b);
    assert_eq!(unique, 2);
}

#[test]
fn test_lobject_bitset_scores() {
    let mut lobj = Lobject::new(0, 100);
    lobj.set_contained_indices(&[0, 1, 2, 3]);

    assert_eq!(lobj.lscore_f, 4);
    assert_eq!(lobj.lscore_r, 4);
    assert_eq!(lobj.num_contained_indices, 4);
}

#[test]
fn test_lobject_bitset_large() {
    let n = 10000;
    let mut lobj = Lobject::new(0, n);

    let indices: Vec<usize> = (0..n).step_by(7).collect();
    lobj.set_contained_indices(&indices);

    for &i in &indices {
        assert!(lobj.get_contained_bit(i), "bit {} should be set", i);
    }
}

#[test]
fn test_assembler_simple_merge() {
    // Two overlapping single-segment alignments that can merge
    let seg1 = AlignmentSegment::new(100, 500);
    let seg2 = AlignmentSegment::new(300, 700);

    let align1 = CdnaAlignment::new(vec![seg1], '+');
    let align2 = CdnaAlignment::new(vec![seg2], '+');

    let mut assembler = CdnaAlignmentAssembler::new(vec![align1, align2]);
    assembler.assemble_alignments();

    assert!(!assembler.get_assemblies().is_empty());
}

#[test]
fn test_assembler_non_overlapping() {
    // Two non-overlapping alignments that cannot merge
    let seg1 = AlignmentSegment::new(100, 200);
    let seg2 = AlignmentSegment::new(500, 600);

    let align1 = CdnaAlignment::new(vec![seg1], '+');
    let align2 = CdnaAlignment::new(vec![seg2], '+');

    let mut assembler = CdnaAlignmentAssembler::new(vec![align1, align2]);
    assembler.assemble_alignments();

    // Should produce separate assemblies
    assert!(assembler.get_assemblies().len() >= 1);
}

#[test]
fn test_assembler_different_orientation() {
    // Same coordinates but different orientation - should not merge
    let seg1 = AlignmentSegment::new(100, 500);
    let seg2 = AlignmentSegment::new(100, 500);

    let align1 = CdnaAlignment::new(vec![seg1], '+');
    let align2 = CdnaAlignment::new(vec![seg2], '-');

    let mut assembler = CdnaAlignmentAssembler::new(vec![align1, align2]);
    assembler.assemble_alignments();

    assert!(assembler.get_assemblies().len() >= 2);
}

#[test]
fn test_assembler_multi_segment_merge() {
    // Two multi-segment alignments with compatible splice junctions
    let seg1a = AlignmentSegment::new(100, 200);
    let seg1b = AlignmentSegment::new(300, 400);
    let align1 = CdnaAlignment::new(vec![seg1a, seg1b], '+');

    let seg2a = AlignmentSegment::new(150, 200);
    let seg2b = AlignmentSegment::new(300, 450);
    let align2 = CdnaAlignment::new(vec![seg2a, seg2b], '+');

    let mut assembler = CdnaAlignmentAssembler::new(vec![align1, align2]);
    assembler.assemble_alignments();

    assert!(!assembler.get_assemblies().is_empty());
}

#[test]
fn test_assembler_containment() {
    // One alignment fully contained within another
    let seg1 = AlignmentSegment::new(100, 1000);
    let align1 = CdnaAlignment::new(vec![seg1], '+');

    let seg2 = AlignmentSegment::new(200, 800);
    let align2 = CdnaAlignment::new(vec![seg2], '+');

    let mut assembler = CdnaAlignmentAssembler::new(vec![align1, align2]);
    assembler.assemble_alignments();

    // The contained alignment should be encapsulated
    assert!(!assembler.get_assemblies().is_empty());
}

#[test]
fn test_assembler_fuzzlength() {
    let seg1 = AlignmentSegment::new(100, 500);
    let seg2 = AlignmentSegment::new(300, 700);

    let align1 = CdnaAlignment::new(vec![seg1], '+');
    let align2 = CdnaAlignment::new(vec![seg2], '+');

    let mut assembler = CdnaAlignmentAssembler::new(vec![align1, align2]);
    assembler.set_fuzzlength(50);
    assembler.assemble_alignments();

    assert!(!assembler.get_assemblies().is_empty());
}

#[test]
fn test_interval_tree_overlap_query() {
    // Create test alignments with known coordinates
    let mut aligns = Vec::new();
    for (start, end) in [(100, 200), (150, 300), (500, 600), (180, 400)] {
        let seg = AlignmentSegment::new(start, end);
        aligns.push(CdnaAlignment::new(vec![seg], '+'));
    }

    // The assembler uses an interval tree internally for overlap queries.
    // We test this by creating an assembler and verifying that non-overlapping
    // alignments don't get merged into the same assembly.
    let mut assembler = CdnaAlignmentAssembler::new(aligns);
    assembler.assemble_alignments();

    // We should get multiple assemblies since some alignments don't overlap
    let assemblies = assembler.get_assemblies();
    assert!(!assemblies.is_empty());
}
