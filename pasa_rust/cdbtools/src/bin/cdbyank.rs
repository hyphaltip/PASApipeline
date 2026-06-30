use std::env;
use std::fs::File;
use std::io::{self, BufRead, Read, Seek, SeekFrom, Write};
use std::process;

use cdbtools::cdb_reader::CdbReader;
use cdbtools::IdxData;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        print_usage();
        process::exit(1);
    }

    let mut index_file: Option<String> = None;
    let mut db_file: Option<String> = None;
    let mut output_file: Option<String> = None;
    let mut accession: Option<String> = None;
    let mut list_keys = false;
    let mut num_records = false;
    let mut defline_only = false;
    let mut position_only = false;
    let mut range_mode = false;
    let mut summary = false;
    let mut case_insensitive = false;
    let mut _many = false;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-a" => { accession = args.get(i + 1).cloned(); i += 2; }
            "-d" => { db_file = args.get(i + 1).cloned(); i += 2; }
            "-o" => { output_file = args.get(i + 1).cloned(); i += 2; }
            "-R" => { range_mode = true; i += 1; }
            "-l" => { list_keys = true; i += 1; }
            "-n" => { num_records = true; i += 1; }
            "-F" => { defline_only = true; i += 1; }
            "-P" => { position_only = true; i += 1; }
            "-s" => { summary = true; i += 1; }
            "-i" => { case_insensitive = true; i += 1; }
            "-x" => { _many = true; i += 1; }
            "-h" | "--h" => { print_usage(); process::exit(0); }
            _ => {
                if index_file.is_none() {
                    index_file = Some(args[i].clone());
                }
                i += 1;
            }
        }
    }

    let index_file = match index_file {
        Some(f) => f,
        None => { eprintln!("Error: no index file specified"); process::exit(1); }
    };

    let reader = CdbReader::open(&index_file).unwrap_or_else(|e| {
        eprintln!("Error opening index file {}: {}", index_file, e);
        process::exit(1);
    });

    // Handle informational queries first
    if summary {
        println!("Index file: {}", index_file);
        println!("Index size: {} bytes", reader.len());
        return;
    }

    if list_keys {
        for (key, _data) in reader.iter() {
            println!("{}", String::from_utf8_lossy(&key));
        }
        return;
    }

    if num_records {
        let mut count = 0;
        for _ in reader.iter() {
            count += 1;
        }
        println!("{}", count);
        return;
    }

    // Determine database file
    let db_file = db_file.unwrap_or_else(|| {
        // Default: strip .cidx from index filename
        index_file.trim_end_matches(".cidx").to_string()
    });

    // Open database file
    let mut db = File::open(&db_file).unwrap_or_else(|e| {
        eprintln!("Error opening database file {}: {}", db_file, e);
        process::exit(1);
    });

    // Determine output destination
    let mut out: Box<dyn Write> = match &output_file {
        Some(f) => Box::new(File::create(f).unwrap_or_else(|e| {
            eprintln!("Error creating output file {}: {}", f, e);
            process::exit(1);
        })),
        None => Box::new(io::stdout()),
    };

    // Process queries
    if let Some(acc) = &accession {
        // Single accession query
        let search_key = if case_insensitive {
            acc.to_lowercase()
        } else {
            acc.clone()
        };

        if let Some(idx_data) = reader.find(search_key.as_bytes()) {
            fetch_record(&mut db, &idx_data, &mut out, defline_only, position_only);
        }
    } else if range_mode {
        // Range extraction mode: read "key start end" from stdin
        let stdin = io::stdin();
        for line in stdin.lock().lines() {
            let line = line.unwrap_or_default();
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 3 { continue; }

            let key = parts[0];
            let start: u64 = parts[1].parse().unwrap_or(0);
            let end: u64 = parts[2].parse().unwrap_or(0);

            if let Some(idx_data) = reader.find(key.as_bytes()) {
                fetch_range(&mut db, &idx_data, start, end, &mut out);
            }
        }
    } else {
        // Default: read keys from stdin, one per line
        let stdin = io::stdin();
        for line in stdin.lock().lines() {
            let key = line.unwrap_or_default();
            if key.is_empty() { continue; }

            let search_key = if case_insensitive {
                key.to_lowercase()
            } else {
                key.clone()
            };

            if let Some(idx_data) = reader.find(search_key.as_bytes()) {
                fetch_record(&mut db, &idx_data, &mut out, defline_only, position_only);
            }
        }
    }
}

fn fetch_record(
    db: &mut File,
    idx_data: &IdxData,
    out: &mut impl Write,
    defline_only: bool,
    position_only: bool,
) {
    if position_only {
        writeln!(out, "{}", idx_data.fpos).ok();
        return;
    }

    db.seek(SeekFrom::Start(idx_data.fpos)).ok();
    let mut buf = vec![0u8; idx_data.reclen as usize];
    db.read_exact(&mut buf).ok();

    if defline_only {
        // Output only the first line (the defline)
        let text = String::from_utf8_lossy(&buf);
        if let Some(line) = text.lines().next() {
            writeln!(out, "{}", line).ok();
        }
    } else {
        out.write_all(&buf).ok();
        // Ensure trailing newline to match C++ cdbyank behavior
        if !buf.ends_with(b"\n") {
            writeln!(out).ok();
        }
    }
}

fn fetch_range(
    db: &mut File,
    idx_data: &IdxData,
    start: u64,
    end: u64,
    out: &mut impl Write,
) {
    db.seek(SeekFrom::Start(idx_data.fpos)).ok();
    let mut buf = vec![0u8; idx_data.reclen as usize];
    db.read_exact(&mut buf).ok();

    // Parse the FASTA record to extract the sub-range
    let text = String::from_utf8_lossy(&buf);
    let mut lines = text.lines();
    let header = lines.next().unwrap_or("");

    // Collect sequence, stripping newlines
    let seq: String = lines.collect::<Vec<&str>>().join("");

    // Extract sub-range (1-based, inclusive)
    let start_idx = ((start - 1) as usize).min(seq.len());
    let end_idx = (end as usize).min(seq.len());
    let sub_seq = &seq[start_idx..end_idx];

    // Output with 60 chars per line
    writeln!(out, "{}:{}-{}", header.trim_start_matches('>'), start, end).ok();
    for chunk in sub_seq.as_bytes().chunks(60) {
        out.write_all(chunk).ok();
        writeln!(out).ok();
    }
}

fn print_usage() {
    eprintln!("Usage: cdbyank_rust <index_file> [options]");
    eprintln!();
    eprintln!("Options:");
    eprintln!("  -a <key>       Accession name to retrieve");
    eprintln!("  -d <file>      Database file (default: index without .cidx)");
    eprintln!("  -o <file>      Output file (default: stdout)");
    eprintln!("  -R             Range extraction mode (stdin: key start end)");
    eprintln!("  -l             List all keys in index");
    eprintln!("  -n             Display number of records");
    eprintln!("  -F             Pull only defline");
    eprintln!("  -P             Display file offset positions only");
    eprintln!("  -s             Display indexing summary");
    eprintln!("  -i             Case-insensitive query");
    eprintln!("  -x             Allow multiple records per key");
    eprintln!("  -h             Show this help message");
}
