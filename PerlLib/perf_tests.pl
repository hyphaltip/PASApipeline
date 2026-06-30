#!/usr/bin/env perl

use strict;
use warnings;

use Carp;
use FindBin;
use lib ("$FindBin::Bin/");

use Time::HiRes qw(gettimeofday);
use Fasta_reader;
use Fasta_retriever;
use PSL_parser;
use Gene_obj;
use GFF3_utils;

my $NUM_ITERATIONS = 1000;
my $LARGE_FASTA_SIZE = 10000;

sub run_all_tests {
    print "=" x 70, "\n";
    print "PASA Performance Test Suite\n";
    print "=" x 70, "\n\n";
    
    my @results;
    
    push @results, run_fasta_reader_benchmark();
    push @results, run_gene_obj_benchmark();
    push @results, run_psl_parser_benchmark();
    push @results, run_regex_benchmark();
    push @results, run_string_concat_benchmark();
    
    print_summary(@results);
}

sub run_fasta_reader_benchmark {
    print "Testing Fasta_reader performance...\n";
    
    my $test_fasta = _create_test_fasta();
    
    my $t0 = gettimeofday();
    for (1..$NUM_ITERATIONS) {
        my $reader = Fasta_reader->new($test_fasta);
        while (my $entry = $reader->next()) {
            my $seq = $entry->get_sequence();
        }
    }
    my $elapsed = gettimeofday() - $t0;
    
    unlink $test_fasta if -e $test_fasta;
    
    my $ops_per_sec = sprintf("%.0f", ($NUM_ITERATIONS * $LARGE_FASTA_SIZE) / $elapsed);
    print "  Fasta_reader: ${elapsed}s for $NUM_ITERATIONS iterations\n";
    print "  Throughput: ~$ops_per_sec sequences/sec\n\n";
    
    return {
        name => 'Fasta_reader',
        time => $elapsed,
        ops => $ops_per_sec,
    };
}

sub run_gene_obj_benchmark {
    print "Testing Gene_obj get_exons() caching...\n";
    
    my $gene = _create_test_gene(50);
    
    my $t0 = gettimeofday();
    for (1..$NUM_ITERATIONS) {
        my @exons = $gene->get_exons();
    }
    my $elapsed = gettimeofday() - $t0;
    
    $t0 = gettimeofday();
    for (1..$NUM_ITERATIONS) {
        my @exons = $gene->get_exons();
        my $seq = "";
        for my $exon (@exons) {
            $seq .= $exon->{end5} . "-";
        }
    }
    my $elapsed_with_seq = gettimeofday() - $t0;
    
    my $cached_speedup = $elapsed_with_seq / $elapsed;
    
    print "  get_exons() x$NUM_ITERATIONS: ${elapsed}s\n";
    print "  With sequence creation: ${elapsed_with_seq}s\n";
    print "  Caching benefit: ${cached_speedup}x faster with caching\n\n";
    
    return {
        name => 'Gene_obj get_exons',
        time => $elapsed,
        ops => sprintf("%.0f", $NUM_ITERATIONS / $elapsed),
    };
}

sub run_psl_parser_benchmark {
    print "Testing PSL_parser toString() optimization...\n";
    
    my $psl_line = _create_test_psl_line(100);
    my $entry = PSL_entry->new($psl_line);
    
    my $t0 = gettimeofday();
    for (1..$NUM_ITERATIONS) {
        my $str = $entry->toString();
    }
    my $elapsed = gettimeofday() - $t0;
    
    $t0 = gettimeofday();
    for (1..$NUM_ITERATIONS) {
        my $pid = $entry->get_per_id();
    }
    my $pid_elapsed = gettimeofday() - $t0;
    
    print "  toString() x$NUM_ITERATIONS: ${elapsed}s\n";
    print "  get_per_id() x$NUM_ITERATIONS: ${pid_elapsed}s\n\n";
    
    return {
        name => 'PSL_parser',
        time => $elapsed,
        ops => sprintf("%.0f", $NUM_ITERATIONS / $elapsed),
    };
}

sub run_regex_benchmark {
    print "Testing precompiled regex vs inline...\n";
    
    my $test_line = "ID=gene0001;Name=test_gene;Note=description;Parent=transcript001;Alias=locus1";
    my $iterations = $NUM_ITERATIONS * 10;
    
    my $t0 = gettimeofday();
    for (1..$iterations) {
        $test_line =~ /ID="?([^;\s"]+)"?;?/;
        $test_line =~ /Name="?([^\;"]+)"?/;
        $test_line =~ /Note="?([^\;"]+)"?/;
    }
    my $inline_elapsed = gettimeofday() - $t0;
    
    my $ID_RE = qr/ID="?([^;\s"]+)"?;?/;
    my $NAME_RE = qr/Name="?([^\;"]+)"?/;
    my $NOTE_RE = qr/Note="?([^\;"]+)"?/;
    
    $t0 = gettimeofday();
    for (1..$iterations) {
        $test_line =~ $ID_RE;
        $test_line =~ $NAME_RE;
        $test_line =~ $NOTE_RE;
    }
    my $precompiled_elapsed = gettimeofday() - $t0;
    
    my $speedup = $inline_elapsed / $precompiled_elapsed;
    
    print "  Inline regex x$iterations: ${inline_elapsed}s\n";
    print "  Precompiled x$iterations: ${precompiled_elapsed}s\n";
    print "  Speedup: ${speedup}x faster with precompiled\n\n";
    
    return {
        name => 'Regex precompile',
        time => $precompiled_elapsed,
        ops => sprintf("%.0f", $iterations / $precompiled_elapsed),
    };
}

sub run_string_concat_benchmark {
    print "Testing string concatenation vs join...\n";
    
    my $num_parts = 1000;
    my $iterations = int($NUM_ITERATIONS / 10);
    
    my $t0 = gettimeofday();
    for (1..$iterations) {
        my $str = "";
        for (1..$num_parts) {
            $str .= "part_$_";
        }
    }
    my $concat_elapsed = gettimeofday() - $t0;
    
    $t0 = gettimeofday();
    for (1..$iterations) {
        my @parts;
        for (1..$num_parts) {
            push @parts, "part_$_";
        }
        my $str = join('', @parts);
    }
    my $join_elapsed = gettimeofday() - $t0;
    
    my $speedup = $concat_elapsed / $join_elapsed;
    
    print "  String concat: ${concat_elapsed}s\n";
    print "  Array join: ${join_elapsed}s\n";
    print "  Speedup: ${speedup}x faster with join\n\n";
    
    return {
        name => 'String concat vs join',
        time => $join_elapsed,
        ops => sprintf("%.0f", $iterations / $join_elapsed),
    };
}

sub _create_test_fasta {
    my $fname = "/tmp/perf_test_$$.fa";
    open(my $fh, '>', $fname) or die;
    
    for my $i (1..$LARGE_FASTA_SIZE) {
        print $fh ">$i\n";
        print $fh ("ACGT" x 250 . "\n");
    }
    close $fh;
    return $fname;
}

sub _create_test_gene {
    my $num_exons = shift || 10;
    my $gene = Gene_obj->new();
    
    $gene->{mRNA_exon_objs} = [];
    for my $i (1..$num_exons) {
        push @{$gene->{mRNA_exon_objs}}, { end5 => $i * 100, end3 => $i * 100 + 90 };
    }
    $gene->{strand} = '+';
    
    return $gene;
}

sub _create_test_psl_line {
    my $num_blocks = shift || 10;
    
    my @block_sizes = ((100) x $num_blocks);
    my @q_starts = map { $_ * 100 } 0..$num_blocks-1;
    my @t_starts = map { $_ * 100 } 0..$num_blocks-1;
    
    return join("\t",
        1000, 0, 0, 0, 0, 0, 0, 0, '+',
        "test_accession", 2000, 0, 2000,
        "chr1", 100000,
        1000, 2000,
        $num_blocks,
        join(',', @block_sizes),
        join(',', @q_starts),
        join(',', @t_starts),
    );
}

sub print_summary {
    my @results = @_;
    
    print "=" x 70, "\n";
    print "Summary\n";
    print "=" x 70, "\n\n";
    
    printf "%-30s %10s %15s\n", "Test", "Time (s)", "Ops/sec";
    print "-" x 55, "\n";
    
    for my $r (@results) {
        printf "%-30s %10.3f %15s\n", 
            $r->{name}, $r->{time}, $r->{ops};
    }
    
    print "\nRecommendations:\n";
    print "  - Use precompiled regex patterns in hot loops\n";
    print "  - Cache computed values (get_exons, get_per_id)\n";
    print "  - Use join() instead of concatenation in loops\n";
    print "  - Reuse filehandles when possible\n";
}

if (@ARGV && $ARGV[0] eq '--help') {
    print "Usage: $0 [--help] [NUM_ITERATIONS]\n";
    print "  NUM_ITERATIONS: Number of iterations for benchmarks (default: $NUM_ITERATIONS)\n";
    exit(0);
}

if (@ARGV) {
    $NUM_ITERATIONS = $ARGV[0];
}

run_all_tests();