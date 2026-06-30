use cdbtools::faidx_reader::{FaidxEntry, FaidxReader, FastaDatabase};
use cdbtools::{cdb_hash, read_u32_le, IdxData};
use std::io::Write;

#[test]
fn test_cdb_hash_known_values() {
    // Test the djb2 hash function
    // h = 5381, then for each byte: h = h + (h << 5); h ^= c
    assert_eq!(cdb_hash(b""), 5381);

    // For 'a': h = 5381 + (5381 << 5) = 5381 + 172192 = 177573; then h ^= 97 = 177573 ^ 97
    let expected: u32 = 5381u32.wrapping_add(5381 << 5) ^ (b'a' as u32);
    assert_eq!(cdb_hash(b"a"), expected);
}

#[test]
fn test_cdb_hash_consistency() {
    // Same input should produce same output
    assert_eq!(cdb_hash(b"hello"), cdb_hash(b"hello"));
    assert_ne!(cdb_hash(b"hello"), cdb_hash(b"world"));
}

#[test]
fn test_read_u32_le() {
    assert_eq!(read_u32_le(&[0, 0, 0, 0], 0), 0);
    assert_eq!(read_u32_le(&[1, 0, 0, 0], 0), 1);
    assert_eq!(read_u32_le(&[0xFF, 0xFF, 0xFF, 0xFF], 0), 0xFFFFFFFF);
    assert_eq!(read_u32_le(&[0x78, 0x56, 0x34, 0x12], 0), 0x12345678);
}

#[test]
fn test_idx_data_32bit() {
    let data = [100, 0, 0, 0, 200, 0, 0, 0]; // fpos=100, reclen=200
    let idx = IdxData::from_32bit(&data);
    assert_eq!(idx.fpos, 100);
    assert_eq!(idx.reclen, 200);
}

#[test]
fn test_idx_data_64bit() {
    let mut data = [0u8; 12];
    data[0..8].copy_from_slice(&500u64.to_le_bytes());
    data[8..12].copy_from_slice(&300u32.to_le_bytes());
    let idx = IdxData::from_64bit(&data);
    assert_eq!(idx.fpos, 500);
    assert_eq!(idx.reclen, 300);
}

#[test]
fn test_faidx_entry_parsing() {
    let entry = FaidxEntry {
        name: "chr1".to_string(),
        length: 1000,
        offset: 50,
        line_bases: 60,
        line_width: 61,
    };
    assert_eq!(entry.name, "chr1");
    assert_eq!(entry.length, 1000);
}

#[test]
fn test_faidx_reader_open() {
    // Create a temporary .fai file
    let fai_content = "chr1\t1000\t6\t60\t61\nchr2\t2000\t1010\t60\t61\n";
    let temp_path = "/tmp/test_faidx.fai";
    let mut file = std::fs::File::create(temp_path).unwrap();
    file.write_all(fai_content.as_bytes()).unwrap();

    let reader = FaidxReader::open(temp_path).unwrap();

    assert_eq!(reader.len(), 2);
    assert!(reader.get("chr1").is_some());
    assert!(reader.get("chr2").is_some());
    assert!(reader.get("chr3").is_none());

    let chr1 = reader.get("chr1").unwrap();
    assert_eq!(chr1.length, 1000);
    assert_eq!(chr1.offset, 6);
    assert_eq!(chr1.line_bases, 60);
    assert_eq!(chr1.line_width, 61);

    std::fs::remove_file(temp_path).ok();
}

#[test]
fn test_fasta_database_fetch() {
    // Create a temporary FASTA file and .fai index
    let fasta_content = ">chr1\nACGTACGTACGTACGT\n>chr2\nGGGGCCCCAAAATTTT\n";
    let fai_content = "chr1\t16\t6\t16\t17\nchr2\t16\t29\t16\t17\n";

    let fasta_path = "/tmp/test_fetch.fa";
    let fai_path = "/tmp/test_fetch.fa.fai";

    let mut file = std::fs::File::create(fasta_path).unwrap();
    file.write_all(fasta_content.as_bytes()).unwrap();

    let mut fai_file = std::fs::File::create(fai_path).unwrap();
    fai_file.write_all(fai_content.as_bytes()).unwrap();

    let mut db = FastaDatabase::open(fasta_path).unwrap();

    // Fetch chr1
    let result = db.fetch("chr1").unwrap().unwrap();
    assert!(result.contains("ACGTACGTACGTACGT"));
    assert!(result.contains(">chr1"));

    // Fetch chr2
    let result = db.fetch("chr2").unwrap().unwrap();
    assert!(result.contains("GGGGCCCCAAAATTTT"));
    assert!(result.contains(">chr2"));

    // Fetch non-existent
    let result = db.fetch("chr3").unwrap();
    assert!(result.is_none());

    std::fs::remove_file(fasta_path).ok();
    std::fs::remove_file(fai_path).ok();
}

#[test]
fn test_fasta_database_fetch_range() {
    let fasta_content = ">chr1\nACGTACGTACGTACGT\n";
    let fai_content = "chr1\t16\t6\t16\t17\n";

    let fasta_path = "/tmp/test_range.fa";
    let fai_path = "/tmp/test_range.fa.fai";

    let mut file = std::fs::File::create(fasta_path).unwrap();
    file.write_all(fasta_content.as_bytes()).unwrap();

    let mut fai_file = std::fs::File::create(fai_path).unwrap();
    fai_file.write_all(fai_content.as_bytes()).unwrap();

    let mut db = FastaDatabase::open(fasta_path).unwrap();

    // Fetch range chr1:5-10 (should be ACGTAC)
    let result = db.fetch_range("chr1", 5, 10).unwrap().unwrap();
    assert!(result.contains("ACGTAC"));
    assert!(result.contains(">chr1:5-10"));

    std::fs::remove_file(fasta_path).ok();
    std::fs::remove_file(fai_path).ok();
}
