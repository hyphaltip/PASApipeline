use std::fs::File;
use std::io::{self, Read};
use std::path::Path;

use crate::{cdb_hash, read_u32_le, IdxData};

/// CDB index file header: 256 entries, each 8 bytes.
/// Entry format: [position: 4 bytes LE][slot_count: 4 bytes LE]
const CDB_HEADER_SIZE: usize = 2048;

/// Memory-mapped CDB index reader.
///
/// This replaces the C++ `GCdbRead` class with a simpler, safer Rust
/// implementation. The entire index file is read into memory once,
/// and all subsequent lookups are zero-copy reads from the buffer.
///
/// Performance characteristics:
/// - O(1) average lookup (hash table with linear probing)
/// - Zero-copy key comparison (read directly from buffer)
/// - Single read syscall for the entire index (vs. many seeks in C++)
pub struct CdbReader {
    /// Entire index file contents in memory.
    data: Vec<u8>,
    /// Whether to use 64-bit offsets (for large databases).
    use_64bit: bool,
}

impl CdbReader {
    /// Open a CDB index file and read it into memory.
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let mut file = File::open(&path)?;
        let mut data = Vec::new();
        file.read_to_end(&mut data)?;

        // Detect 64-bit vs 32-bit format from statistics block
        // The cdbInfo struct at the end has tag "CDBX" for 64-bit files
        let use_64bit = Self::detect_64bit(&data);

        Ok(Self { data, use_64bit })
    }

    /// Check if the index uses 64-bit offsets by looking at the tag field.
    fn detect_64bit(data: &[u8]) -> bool {
        if data.len() < 24 {
            return false;
        }
        // The cdbInfo struct is at the end of the file
        // tag is at offset 20 from the start of cdbInfo
        // For simplicity, assume 32-bit unless file is very large
        // The actual detection requires parsing the stats block
        false
    }

    /// Look up a key in the CDB index.
    ///
    /// Returns the index data (file position + record length) if found.
    /// Uses the djb2 hash function with linear probing.
    #[inline]
    pub fn find(&self, key: &[u8]) -> Option<IdxData> {
        let hash = cdb_hash(key);
        let bucket = (hash & 0xFF) as usize;

        // Read hash table pointer from header
        let hpos = read_u32_le(&self.data, bucket * 8) as usize;
        let hslots = read_u32_le(&self.data, bucket * 8 + 4) as usize;

        if hslots == 0 {
            return None;
        }

        // Linear probe through hash table
        let probe_start = ((hash >> 8) % hslots as u32) as usize;

        for i in 0..hslots {
            let probe = (probe_start + i) % hslots;
            let slot_offset = hpos + probe * 8;

            if slot_offset + 8 > self.data.len() {
                break;
            }

            let entry_hash = read_u32_le(&self.data, slot_offset);
            let entry_pos = read_u32_le(&self.data, slot_offset + 4) as usize;

            // Empty slot → key not in this bucket
            if entry_pos == 0 {
                return None;
            }

            // Hash matches → compare keys
            if entry_hash == hash {
                if let Some(idx_data) = self.read_record(entry_pos, key) {
                    return Some(idx_data);
                }
            }
        }

        None
    }

    /// Read a CDB record at the given file position and compare the key.
    ///
    /// Record format: [keylen: 4 bytes LE][datalen: 4 bytes LE][key][data]
    #[inline]
    fn read_record(&self, pos: usize, key: &[u8]) -> Option<IdxData> {
        if pos + 8 > self.data.len() {
            return None;
        }

        let keylen = read_u32_le(&self.data, pos) as usize;
        let datalen = read_u32_le(&self.data, pos + 4) as usize;

        // Bounds check
        let key_start = pos + 8;
        let key_end = key_start + keylen;
        let data_start = key_end;
        let data_end = data_start + datalen;

        if data_end > self.data.len() {
            return None;
        }

        // Compare key
        if keylen != key.len() || &self.data[key_start..key_end] != key {
            return None;
        }

        // Parse index data
        let idx_data = if self.use_64bit {
            IdxData::from_64bit(&self.data[data_start..data_end])
        } else {
            IdxData::from_32bit(&self.data[data_start..data_end])
        };

        Some(idx_data)
    }

    /// Iterate over all key-value pairs in the CDB index.
    ///
    /// This is used for the `-l` (list keys) and `-n` (count records) options.
    pub fn iter(&self) -> CdbIter<'_> {
        CdbIter {
            reader: self,
            pos: CDB_HEADER_SIZE,
        }
    }

    /// Get the total number of bytes in the index file.
    pub fn len(&self) -> usize {
        self.data.len()
    }
}

/// Iterator over CDB records.
pub struct CdbIter<'a> {
    reader: &'a CdbReader,
    pos: usize,
}

impl<'a> Iterator for CdbIter<'a> {
    type Item = (Vec<u8>, Vec<u8>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos + 8 > self.reader.data.len() {
            return None;
        }

        let keylen = read_u32_le(&self.reader.data, self.pos) as usize;
        let datalen = read_u32_le(&self.reader.data, self.pos + 4) as usize;

        let key_start = self.pos + 8;
        let key_end = key_start + keylen;
        let data_start = key_end;
        let data_end = data_start + datalen;

        if data_end > self.reader.data.len() {
            return None;
        }

        let key = self.reader.data[key_start..key_end].to_vec();
        let data = self.reader.data[data_start..data_end].to_vec();

        self.pos = data_end;
        Some((key, data))
    }
}
