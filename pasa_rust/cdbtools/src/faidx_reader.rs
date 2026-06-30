use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

/// A single entry in the samtools faidx `.fai` index.
///
/// The `.fai` format is tab-separated with 5 columns:
/// ```text
/// NAME   LENGTH   OFFSET   LINEBASES   LINEWIDTH
/// chr1    248956422   56      60          61
/// ```
///
/// - **NAME**: Sequence identifier (first whitespace-delimited token of the `>` defline)
/// - **LENGTH**: Total sequence length in bases
/// - **OFFSET**: Byte offset of the first base of this sequence within the FASTA file
/// - **LINEBASES**: Number of bases per line (excluding newline)
/// - **LINEWIDTH**: Number of bytes per line (including newline character)
#[derive(Clone, Debug)]
pub struct FaidxEntry {
    pub name: String,
    pub length: u64,
    pub offset: u64,
    pub line_bases: u64,
    pub line_width: u64,
}

/// Samtools faidx index reader.
///
/// This reads the `.fai` index file produced by `samtools faidx` and
/// provides O(1) lookup of sequence entries by name. Combined with
/// the FASTA database file, it enables efficient random-access retrieval
/// of sequences and sub-ranges.
///
/// This is the preferred alternative to `cdbyank` because:
/// 1. `samtools faidx` is already a PASA dependency
/// 2. The `.fai` format is plain text (human-readable, easy to debug)
/// 3. The index is created by a well-maintained external tool
/// 4. The lookup is a simple HashMap get (O(1) average)
pub struct FaidxReader {
    entries: HashMap<String, FaidxEntry>,
}

impl FaidxReader {
    /// Open and parse a `.fai` index file.
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(&path)?;
        let reader = BufReader::new(file);
        let mut entries = HashMap::new();

        for line in reader.lines() {
            let line = line?;
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 5 {
                continue;
            }

            let name = fields[0].to_string();
            let length: u64 = fields[1].parse().unwrap_or(0);
            let offset: u64 = fields[2].parse().unwrap_or(0);
            let line_bases: u64 = fields[3].parse().unwrap_or(0);
            let line_width: u64 = fields[4].parse().unwrap_or(0);

            entries.insert(name.clone(), FaidxEntry {
                name,
                length,
                offset,
                line_bases,
                line_width,
            });
        }

        Ok(Self { entries })
    }

    /// Look up a sequence entry by name.
    pub fn get(&self, name: &str) -> Option<&FaidxEntry> {
        self.entries.get(name)
    }

    /// Get all sequence names in the index.
    pub fn names(&self) -> impl Iterator<Item = &String> {
        self.entries.keys()
    }

    /// Get the total number of sequences in the index.
    pub fn len(&self) -> usize {
        self.entries.len()
    }
}

/// A FASTA database file with random-access sequence retrieval.
///
/// Uses the `.fai` index for byte-offset lookups, then reads the
/// sequence directly from the FASTA file using `seek()`.
pub struct FastaDatabase {
    file: File,
    index: FaidxReader,
}

impl FastaDatabase {
    /// Open a FASTA database file and its `.fai` index.
    ///
    /// The `.fai` file must exist at `<fasta_path>.fai`. If it doesn't,
    /// the user should run `samtools faidx <fasta_path>` first.
    pub fn open<P: AsRef<Path>>(fasta_path: P) -> io::Result<Self> {
        let path = fasta_path.as_ref();
        let fai_path = path.with_extension("fai");

        // If .fai doesn't exist, try <path>.fai (common convention)
        let fai_path = if fai_path.exists() {
            fai_path
        } else {
            let mut p = path.as_os_str().to_owned();
            p.push(".fai");
            p.into()
        };

        let file = File::open(path)?;
        let index = FaidxReader::open(fai_path)?;

        Ok(Self { file, index })
    }

    /// Retrieve a full sequence by name.
    ///
    /// Returns the FASTA record (header + sequence) as a string.
    pub fn fetch(&mut self, name: &str) -> io::Result<Option<String>> {
        let entry = match self.index.get(name) {
            Some(e) => e,
            None => return Ok(None),
        };

        // Seek to the start of this sequence in the FASTA file
        self.file.seek(SeekFrom::Start(entry.offset))?;

        // Read the entire sequence (length bytes, handling newlines)
        let mut seq = String::with_capacity(entry.length as usize + 1);
        let mut buf = vec![0u8; 65536]; // 64KB buffer
        let mut bytes_read = 0u64;

        while bytes_read < entry.length {
            let remaining = entry.length - bytes_read;
            let to_read = (remaining as usize).min(buf.len());
            let n = self.file.read(&mut buf[..to_read])?;
            if n == 0 { break; }

            // Convert to string, stripping newlines
            let chunk = std::str::from_utf8(&buf[..n]).unwrap_or("");
            for c in chunk.chars() {
                if c != '\n' && c != '\r' {
                    seq.push(c);
                }
            }
            bytes_read += n as u64;
        }

        Ok(Some(format!(">{}\n{}", name, seq)))
    }

    /// Retrieve a sub-range of a sequence by name.
    ///
    /// `start` and `end` are 1-based, inclusive coordinates.
    /// Returns the FASTA record (header + sub-sequence) as a string.
    pub fn fetch_range(&mut self, name: &str, start: u64, end: u64) -> io::Result<Option<String>> {
        let entry = match self.index.get(name) {
            Some(e) => e,
            None => return Ok(None),
        };

        if start < 1 || end < start || end > entry.length {
            return Ok(None);
        }

        // Compute byte offset for the start position
        let seq_offset = start - 1; // 0-based offset within sequence
        let byte_offset = entry.offset
            + (seq_offset / entry.line_bases) * entry.line_width
            + (seq_offset % entry.line_bases);

        self.file.seek(SeekFrom::Start(byte_offset))?;

        // Read (end - start + 1) bases, handling newlines
        let num_bases = end - start + 1;
        let mut seq = String::with_capacity(num_bases as usize);
        let mut buf = vec![0u8; 65536]; // 64KB buffer
        let mut bases_read = 0u64;

        while bases_read < num_bases {
            let remaining = num_bases - bases_read;
            let to_read = (remaining as usize).min(buf.len());
            let n = self.file.read(&mut buf[..to_read])?;
            if n == 0 { break; }

            let chunk = std::str::from_utf8(&buf[..n]).unwrap_or("");
            for c in chunk.chars() {
                if bases_read >= num_bases { break; }
                if c != '\n' && c != '\r' {
                    seq.push(c);
                    bases_read += 1;
                }
            }
        }

        Ok(Some(format!(">{}:{}-{}\n{}", name, start, end, seq)))
    }
}
