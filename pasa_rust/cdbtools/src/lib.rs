pub mod cdb_reader;
pub mod faidx_reader;

/// FASTA record retrieved from the database.
pub struct FastaRecord {
    pub header: String,
    pub sequence: String,
}

/// CDB hash function (djb2 variant).
///
/// `h = 5381; for each byte c: h = h + (h << 5); h ^= c`
///
/// This must exactly match the C++ implementation in `gcdb.cpp:682-695`
/// for backward compatibility with existing `.cidx` index files.
pub fn cdb_hash(data: &[u8]) -> u32 {
    let mut h: u32 = 5381;
    for &c in data {
        h = h.wrapping_add(h << 5);
        h ^= c as u32;
    }
    h
}

/// Read a little-endian u32 from a byte slice at the given offset.
#[inline]
pub fn read_u32_le(data: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes([
        data[offset],
        data[offset + 1],
        data[offset + 2],
        data[offset + 3],
    ])
}

/// Index data stored in CDB records.
/// Maps accession names to (file_position, record_length) pairs.
#[derive(Clone, Copy, Debug)]
pub struct IdxData {
    pub fpos: u64,
    pub reclen: u32,
}

impl IdxData {
    /// Parse 8-byte 32-bit offset format.
    pub fn from_32bit(data: &[u8]) -> Self {
        Self {
            fpos: read_u32_le(data, 0) as u64,
            reclen: read_u32_le(data, 4),
        }
    }

    /// Parse 12-byte 64-bit offset format.
    pub fn from_64bit(data: &[u8]) -> Self {
        Self {
            fpos: u64::from_le_bytes([
                data[0], data[1], data[2], data[3],
                data[4], data[5], data[6], data[7],
            ]),
            reclen: read_u32_le(data, 8),
        }
    }
}

// read_u32_le is already pub, no re-export needed
