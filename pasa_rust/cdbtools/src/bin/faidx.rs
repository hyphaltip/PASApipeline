use std::env;
use std::fs::File;
use std::io::{self, BufRead, Write};
use std::process;

use cdbtools::faidx_reader::FastaDatabase;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        print_usage();
        process::exit(1);
    }

    let mut fasta_file: Option<String> = None;
    let mut output_file: Option<String> = None;
    let mut regions: Vec<String> = Vec::new();
    let mut list_keys = false;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-o" => { output_file = args.get(i + 1).cloned(); i += 2; }
            "-l" => { list_keys = true; i += 1; }
            "-h" | "--help" => { print_usage(); process::exit(0); }
            _ => {
                if fasta_file.is_none() {
                    fasta_file = Some(args[i].clone());
                } else {
                    regions.push(args[i].clone());
                }
                i += 1;
            }
        }
    }

    let fasta_file = match fasta_file {
        Some(f) => f,
        None => { eprintln!("Error: no FASTA file specified"); process::exit(1); }
    };

    let mut db = FastaDatabase::open(&fasta_file).unwrap_or_else(|e| {
        eprintln!("Error opening FASTA database {}: {}", fasta_file, e);
        eprintln!("Hint: run 'samtools faidx {}' to create the index", fasta_file);
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

    if list_keys {
        // List all sequence names from the .fai index
        let fai_path = format!("{}.fai", fasta_file);
        let reader = cdbtools::faidx_reader::FaidxReader::open(&fai_path)
            .unwrap_or_else(|e| {
                eprintln!("Error reading .fai index: {}", e);
                process::exit(1);
            });
        for name in reader.names() {
            writeln!(out, "{}", name).ok();
        }
        return;
    }

    // Process region queries
    if !regions.is_empty() {
        // Regions specified on command line
        for region in &regions {
            fetch_region(&mut db, region, &mut out);
        }
    } else {
        // Read region queries from stdin (one per line)
        let stdin = io::stdin();
        for line in stdin.lock().lines() {
            let region = line.unwrap_or_default();
            if region.is_empty() { continue; }
            fetch_region(&mut db, &region, &mut out);
        }
    }
}

fn fetch_region(db: &mut FastaDatabase, region: &str, out: &mut impl Write) {
    // Parse region specification: "name" or "name:start-end"
    let (name, range) = if let Some(colon_pos) = region.find(':') {
        let name = &region[..colon_pos];
        let range_str = &region[colon_pos + 1..];
        let parts: Vec<&str> = range_str.split('-').collect();
        if parts.len() == 2 {
            let start: u64 = parts[0].parse().unwrap_or(0);
            let end: u64 = parts[1].parse().unwrap_or(0);
            (name, Some((start, end)))
        } else {
            (name, None)
        }
    } else {
        (region, None)
    };

    let result = if let Some((start, end)) = range {
        db.fetch_range(name, start, end)
    } else {
        db.fetch(name)
    };

    if let Ok(Some(fasta_text)) = result {
        writeln!(out, "{}", fasta_text).ok();
    } else {
        eprintln!("Warning: sequence '{}' not found in index", name);
    }
}

fn print_usage() {
    eprintln!("Usage: faidx_rust <fasta_file> [regions...] [options]");
    eprintln!();
    eprintln!("Retrieve sequences from a FASTA file using a samtools .fai index.");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  fasta_file       Path to the FASTA file (must have a .fai index)");
    eprintln!("  regions          Sequence regions to retrieve (see formats below)");
    eprintln!();
    eprintln!("Region formats:");
    eprintln!("  <seq_name>                     Retrieve full sequence");
    eprintln!("  <seq_name>:<start>-<end>       Retrieve sub-range (1-based, inclusive)");
    eprintln!();
    eprintln!("Options:");
    eprintln!("  -o <file>       Output file (default: stdout)");
    eprintln!("  -l              List all sequence names in index");
    eprintln!("  -h              Show this help message");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  faidx_rust genome.fa chr1:1000-2000 chr2");
    eprintln!("  echo 'chr1:1000-2000' | faidx_rust genome.fa");
    eprintln!();
    eprintln!("To create the .fai index:");
    eprintln!("  samtools faidx genome.fa");
}
