#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use File::Basename;
use Cwd;

use Carp;
use File::Temp qw(tempfile);
use Getopt::Long qw(:config no_ignore_case bundling pass_through);


my $usage = <<__EOUSAGE__;

######################################################################
#
#  Required:
#  --genome <string>           target genome to align to
#  --transcripts <string>      cdna sequences to align
#
#  Optional:
#  -N <int>                    number of top hits (default: 1)
#  -I <int>                    max intron length
#  --CPU <int>                 number of threads (default: 2)
#  --no_chimera                do not report chimeric alignmetnts
#  --SAM                       output in SAM format
#
#######################################################################


__EOUSAGE__

    ;


my ($genome, $transcriptDB, $max_intron);
my $CPU = 2;

my $help_flag;

my $number_top_hits = 1;
my $no_chimera_flag = 0;
my $SAM_flag = 0;

&GetOptions( 'h' => \$help_flag,
             'genome=s' => \$genome,
             'transcripts=s' => \$transcriptDB,
             'I=i' => \$max_intron,
             'CPU=i' => \$CPU,
             'N=i' => \$number_top_hits,
             'no_chimera' => \$no_chimera_flag,
             'SAM' => \$SAM_flag,
             
             );


unless ($genome && $transcriptDB) {
    die $usage;
}


my $GMAP_CUSTOM_OPTS = $ENV{GMAP_CUSTOM_OPTS} || "";

main: {
	
	my $genomeName = basename($genome);
	my $genomeDir = $genomeName . ".gmap";

	my $genomeBaseDir = dirname($genome);

	my $cwd = cwd();
	
	unless (-d "$genomeBaseDir/$genomeDir") {
		
        #my $cmd = "gmap_build -D $genomeBaseDir -d $genomeBaseDir/$genomeDir -k 13 $genome >&2";
        #my $cmd = "gmap_build -D $genomeBaseDir -T $genomeBaseDir -d $genomeDir -k 13 $genome >&2";
        my $cmd = "gmap_build -D $genomeBaseDir -d $genomeDir -k 13 $genome >&2";
		&process_cmd($cmd);
	}

	
	## run GMAP

    my $num_gmap_top_hits = $number_top_hits;
    if ((! $no_chimera_flag) && $num_gmap_top_hits == 1) {
        $num_gmap_top_hits = 0; # reports two hits if chimera with this setting.
    }
    
    my $format = ($SAM_flag) ? "samse" : "3";

    my $gmap_prog = (-s $genome > 2**32 ) ? "gmapl" : "gmap";
    
	my $cmd = "$gmap_prog -D $genomeBaseDir -d $genomeDir ${GMAP_CUSTOM_OPTS} -f $format -n $num_gmap_top_hits -x 50 -t $CPU -B 5 ";
	if ($max_intron) {
        $cmd .= " --intronlength=$max_intron ";
    }

    ## gmap output goes to a temp file first, so a gmap crash part-way through
    ## never leaves partial alignments in our STDOUT. On success the file is
    ## copied to STDOUT unchanged (same output as before).
    my ($tmp_fh, $tmp_out) = tempfile("gmap_out.XXXXXX", DIR => $cwd, UNLINK => 1);
    close $tmp_fh;
    print STDERR "CMD: $cmd $transcriptDB > $tmp_out\n";
    my $ret = system("$cmd $transcriptDB > $tmp_out");
    if ($ret == 0) {
        &cat_to_stdout($tmp_out);
        exit(0);
    }
    if ($SAM_flag) {
        die "Error, cmd: $cmd $transcriptDB died with ret ($ret)";
    }

    ## gmap can segfault on individual transcripts (seen with gmap 2021-12-17,
    ## 2023-04-28 and 2025-07-31 on a tandem-repeat contig), which kills the
    ## whole run. Fall back to chunks; split failing chunks down to single
    ## transcripts, skip the ones gmap cannot align, and list them.
    print STDERR "WARNING: gmap failed (ret $ret) on the full transcript set; retrying in chunks to skip transcripts that crash gmap.\n";
    my @records = &read_fasta_records($transcriptDB);
    my @failed;
    my $chunk_size = 500;
    for (my $i = 0; $i <= $#records; $i += $chunk_size) {
        my $j = ($i + $chunk_size - 1 > $#records) ? $#records : $i + $chunk_size - 1;
        &run_gmap_chunk($cmd, [ @records[$i..$j] ], $cwd, \@failed);
    }
    my $failed_file = "$cwd/gmap.failed_transcripts.txt";
    open (my $ofh, ">$failed_file") or die "Error, cannot write $failed_file";
    print $ofh "$_\n" foreach @failed;
    close $ofh;
    print STDERR "WARNING: gmap could not align " . scalar(@failed) . " transcript(s); skipped, listed in $failed_file\n";

	exit(0);
}

####
sub run_gmap_chunk {
    my ($cmd, $records_aref, $cwd, $failed_aref) = @_;
    my @records = @$records_aref;
    return unless @records;
    my ($in_fh, $in_file) = tempfile("gmap_chunk.XXXXXX", DIR => $cwd, SUFFIX => ".fa", UNLINK => 1);
    print $in_fh ">$_->{header}\n$_->{seq}\n" foreach @records;
    close $in_fh;
    my ($out_fh, $out_file) = tempfile("gmap_chunk_out.XXXXXX", DIR => $cwd, UNLINK => 1);
    close $out_fh;
    my $ret = system("$cmd $in_file > $out_file 2>> $cwd/gmap.chunk_retry.log");
    if ($ret == 0) {
        &cat_to_stdout($out_file);
    }
    elsif (scalar(@records) == 1) {
        my ($acc) = split(/\s+/, $records[0]->{header});
        push (@$failed_aref, $acc);
    }
    else {
        my $mid = int(scalar(@records) / 2);
        &run_gmap_chunk($cmd, [ @records[0..$mid-1] ], $cwd, $failed_aref);
        &run_gmap_chunk($cmd, [ @records[$mid..$#records] ], $cwd, $failed_aref);
    }
    unlink ($in_file, $out_file);
    return;
}

####
sub read_fasta_records {
    my ($file) = @_;
    my @records;
    open (my $fh, $file) or die "Error, cannot read $file";
    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ /^>(.*)$/) {
            push (@records, { header => $1, seq => "" });
        }
        elsif (@records) {
            $records[-1]->{seq} .= $line;
        }
    }
    close $fh;
    return @records;
}

####
sub cat_to_stdout {
    my ($file) = @_;
    open (my $fh, $file) or die "Error, cannot read $file";
    while (my $line = <$fh>) {
        print $line;
    }
    close $fh;
    return;
}

####
sub process_cmd {
	my ($cmd) = @_;
	
	print STDERR "CMD: $cmd\n";
	#return;

	my $ret = system($cmd);
	if ($ret) {
		die "Error, cmd: $cmd died with ret ($ret)";
	}

	return;
}



