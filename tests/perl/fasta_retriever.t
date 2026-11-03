#!/usr/bin/env perl
#
# Fasta_retriever: same sequences from plain, .fai-indexed, .gz and .bz2 FASTA.
#
# Before this test, a .gz/.bz2 file returned an empty string for every
# accession (get_seq opened a new stream at the file start and stopped at the
# first header), and without a .fai the compressed bytes were indexed as text.
# Now compressed input is read correctly with a one-time warning, and an empty
# result (e.g. a stale .fai) is an error, not a silent empty sequence.
#
# Run: prove -v tests/perl   (or: make perl-test)

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Test::More;
use File::Temp qw(tempdir);

use Fasta_retriever;

my $dir = tempdir(CLEANUP => 1);

my %SEQ = (
    chr1 => "ACGTACGTAC" x 13,           # 130 bp, multi-line
    chr2 => "GGGCCCAAAT" x 5 . "TTA",    # 53 bp
    chr3 => "AC",                        # short last record
);
my @ORDER = qw(chr1 chr2 chr3);

my $plain = "$dir/genome.fa";
open(my $out, ">", $plain) or die $!;
for my $acc (@ORDER) {
    print $out ">$acc description text\n";
    my $s = $SEQ{$acc};
    print $out substr($s, 0, 60, "") . "\n" while length $s;
}
close $out;

sub check_all {
    my ($label, $file) = @_;
    my $r = Fasta_retriever->new($file);
    for my $acc (reverse @ORDER) {   # out of file order on purpose
        is($r->get_seq($acc), $SEQ{$acc}, "$label: $acc");
    }
    return $r;
}

## 1. plain, no index
check_all("plain, no .fai", $plain);

## 2. plain with a .fai (column 3 = byte offset of the first sequence line)
{
    my $fai_plain = "$dir/indexed.fa";
    system("cp", $plain, $fai_plain) == 0 or die;
    open(my $in, "<", $fai_plain) or die $!;
    open(my $fai, ">", "$fai_plain.fai") or die $!;
    my $pos = 0;
    while (my $line = <$in>) {
        $pos += length $line;
        if ($line =~ /^>(\S+)/) {
            print $fai join("\t", $1, length($SEQ{$1}), $pos, 60, 61), "\n";
        }
    }
    close $in; close $fai;
    check_all("plain, .fai", $fai_plain);

    ## stale .fai: chr2 offset points at the >chr3 header -> empty -> error
    ## (offset 0 cannot be used: get_seq treats 0 as "no position")
    open($in, "<", $fai_plain) or die $!;
    my ($p, $chr3_hdr) = (0, undef);
    while (my $line = <$in>) {
        $chr3_hdr = $p if $line =~ /^>chr3/;
        $p += length $line;
    }
    close $in;
    open($fai, ">", "$fai_plain.fai") or die $!;
    print $fai "chr2\t53\t$chr3_hdr\t60\t61\n";
    close $fai;
    my $r = Fasta_retriever->new($fai_plain);
    eval { $r->get_seq("chr2") };
    like($@, qr/empty sequence retrieved for chr2/, "stale .fai: empty result is an error");
}

## 3. gzip and bzip2
for my $c (["gz", "gzip"], ["bz2", "bzip2"]) {
    my ($ext, $tool) = @$c;
    SKIP: {
        skip "$tool not found", 8 unless system("command -v $tool >/dev/null 2>&1") == 0;
        my $file = "$dir/genome_$ext.fa";
        system("cp", $plain, $file) == 0 or die;
        system($tool, "-f", $file) == 0 or die "$tool failed";
        $file .= ".$ext";

        my $warn = "";
        {
            local *STDERR;
            open(STDERR, ">", \$warn) or die;
            my $r = check_all(".$ext", $file);
            $r->get_seq("chr1");
        }
        my $n = () = $warn =~ /WARNING: Fasta_retriever: .* is compressed/g;
        is($n, 1, ".$ext: compressed warning printed once");

        my $r = Fasta_retriever->new($file);
        eval { $r->get_seq("chrX") };
        like($@, qr/no seek pos for acc: chrX|not found/, ".$ext: unknown accession is an error");

        ## same compressed data under a plain .fa name: detected by magic bytes
        my $noext = "$dir/genome_${ext}_noext.fa";
        system("cp", $file, $noext) == 0 or die;
        {
            local *STDERR;
            my $w = "";
            open(STDERR, ">", \$w) or die;
            check_all(".$ext without extension", $noext);
        }
    }
}

done_testing();
