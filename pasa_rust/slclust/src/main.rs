use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter};
use std::process;

use slclust::Graph;

fn main() {
    let args: Vec<String> = env::args().collect();

    let mut verbose = false;
    let mut want_jaccard = false;
    let mut cutoff: f64 = 0.0;
    let mut input_file: Option<String> = None;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-v" => {
                verbose = true;
            }
            "-vv" => {
                verbose = true;
            }
            "-j" => {
                want_jaccard = true;
                if i + 1 < args.len() {
                    cutoff = args[i + 1].parse().unwrap_or(0.0);
                    i += 1;
                }
            }
            "-h" | "--h" => {
                print_usage();
                process::exit(0);
            }
            _ => {
                input_file = Some(args[i].clone());
            }
        }
        i += 1;
    }

    let stdin = io::stdin();
    let reader: Box<dyn BufRead> = match &input_file {
        Some(path) => {
            let file = File::open(path).unwrap_or_else(|e| {
                eprintln!("Error opening {}: {}", path, e);
                process::exit(1);
            });
            Box::new(BufReader::new(file))
        }
        None => Box::new(stdin.lock()),
    };

    let mut graph = Graph::new();
    let mut pair_count = 0u64;

    for line in reader.lines() {
        let line = line.unwrap_or_default();
        let tokens: Vec<&str> = line.split_whitespace().collect();

        if tokens.len() == 2 && tokens[0] != tokens[1] {
            graph.add_linked_nodes(tokens[0], tokens[1]);
            pair_count += 1;
        }
    }

    if verbose {
        eprintln!("Read {} pairs, {} unique nodes", pair_count, graph.num_nodes());
    }

    let stdout = io::stdout();
    let mut output = BufWriter::new(stdout.lock());

    if want_jaccard {
        let filtered = graph.apply_jaccard_coeff(cutoff);
        let mut filtered = filtered;
        filtered.print_clusters(&mut output).unwrap_or_else(|e| {
            eprintln!("Error writing clusters: {}", e);
            process::exit(1);
        });
    } else {
        graph.print_clusters(&mut output).unwrap_or_else(|e| {
            eprintln!("Error writing clusters: {}", e);
            process::exit(1);
        });
    }
}

fn print_usage() {
    eprintln!("Usage: slclust [options] [input_file]");
    eprintln!();
    eprintln!("Options:");
    eprintln!("  -j <float>  Apply Jaccard coefficient cutoff");
    eprintln!("  -v          Verbose output");
    eprintln!("  -h          Show this help message");
    eprintln!();
    eprintln!("Input: whitespace-separated pairs from stdin or file");
    eprintln!("Output: one line per cluster, space-separated node names");
}
