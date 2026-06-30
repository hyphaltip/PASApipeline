use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process;

use pasa_assembler::alignment_segment::AlignmentSegment;
use pasa_assembler::assembler::CdnaAlignmentAssembler;
use pasa_assembler::cdna_alignment::CdnaAlignment;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: {} inputFile [opts]", args[0]);
        eprintln!("Options:");
        eprintln!("  -F fuzzlength   (bp to discount at alignment termini. default: 20)");
        eprintln!("  -a              illustrate incoming alignments only");
        eprintln!("  -v              verbose");
        process::exit(1);
    }

    let input_file = &args[1];
    let mut fuzzlength: Option<i32> = None;
    let mut illustrate_only = false;

    let mut i = 2;
    while i < args.len() {
        match args[i].as_str() {
            "-F" => {
                if i + 1 < args.len() {
                    let val: i32 = args[i + 1].parse().unwrap_or(20);
                    if val >= 0 && val <= 100 {
                        fuzzlength = Some(val);
                    }
                }
                i += 2;
            }
            "-a" => {
                illustrate_only = true;
                i += 1;
            }
            "-v" => {
                i += 1;
            }
            _ => {
                i += 1;
            }
        }
    }

    let file = File::open(input_file).unwrap_or_else(|e| {
        eprintln!("Error opening {}: {}", input_file, e);
        process::exit(1);
    });

    let reader = BufReader::new(file);
    let mut cdna_list: Vec<CdnaAlignment> = Vec::new();

    println!("//");

    for line in reader.lines() {
        let line = line.unwrap_or_default();
        if !line.contains(',') {
            continue;
        }

        println!("input: {}", line);

        let tokens: Vec<&str> = line.split(',').collect();
        if tokens.len() < 3 {
            continue;
        }

        let acc = tokens[0];
        let orient = tokens[1].chars().next().unwrap_or('+');

        let mut seglist: Vec<AlignmentSegment> = Vec::new();
        for t in &tokens[2..] {
            let coords: Vec<&str> = t.split('-').collect();
            if coords.len() >= 2 {
                let lend: i32 = coords[0].parse().unwrap_or(0);
                let rend: i32 = coords[1].parse().unwrap_or(0);
                seglist.push(AlignmentSegment::new(lend, rend));
            }
        }

        if !seglist.is_empty() {
            let mut align = CdnaAlignment::new(seglist, orient);
            align.set_title(acc);
            cdna_list.push(align);
        }
    }

    let mut assembler = CdnaAlignmentAssembler::new(cdna_list);

    if let Some(fl) = fuzzlength {
        assembler.set_fuzzlength(fl);
    }

    if !illustrate_only {
        assembler.assemble_alignments();
    }

    print!("{}", assembler.to_align_illustration(70));
}
