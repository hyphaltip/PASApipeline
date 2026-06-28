

=head1 NAME

package CdbTools

=cut


=head1 DESCRIPTION

    routines for extracting entries from Fasta file using the CDBtools cdbfasta and cdbyank.

    Updated to support Rust-optimized cdbyank_rust and faidx_rust as alternatives.
    - cdbyank(): Uses cdbyank_rust if available, else falls back to cdbyank (C++)
    - get_seq(): Uses faidx_rust with samtools .fai index (preferred for new code)

=cut

    ;

package main;
our $SEE;


package CdbTools;

use strict;
use warnings;
require Exporter;

our @ISA = qw (Exporter);
our @EXPORT = qw (cdbyank linearize cdbyank_linear get_seq);

## cdbfasta and cdbyank must be in path, otherwise the system will die.

## Detect available tools
sub _which {
    my ($tool) = @_;
    for my $path (split /:/, $ENV{PATH} || '') {
        return "$path/$tool" if -x "$path/$tool";
    }
    return undef;
}

my $CDBYANK_RUST = _which('cdbyank_rust');
my $FAIDX_RUST   = _which('faidx_rust');
my $CDBYANK      = _which('cdbyank');
my $SAMTOOLS     = _which('samtools');

=over 4

=item cdbyank()

B<Description:> Retrieves a fasta sequence entry from a fasta database

B<Parameters:> accession, fastaFilename

B<Returns:> fastaEntry

use the linearize method to extract the fasta entry components

=back

=cut


    ;

sub cdbyank {
    my ($accession, $fastaFile) = @_;
    unless (-s "$fastaFile.cidx") {
        ## regenerate index file:
        my $cmd = "cdbfasta -C $fastaFile";
        my $ret = system $cmd;
        if ($ret) {
            die "Error, couldn't create index file: $cmd, ret($ret)\n";
        }
    }

    ## Prefer Rust cdbyank_rust if available (10-100x faster for large records)
    my $cmd;
    if ($CDBYANK_RUST && -x $CDBYANK_RUST) {
        $cmd = "$CDBYANK_RUST -a '$accession' $fastaFile.cidx";
    } elsif ($CDBYANK && -x $CDBYANK) {
        $cmd = "$CDBYANK -a '$accession' $fastaFile.cidx";
    } else {
        die "Error: neither cdbyank_rust nor cdbyank found in PATH\n";
    }

    if ($SEE) {
        print "CMD: $cmd\n";
    }

    my $fastaEntry = `$cmd`;
    if ($?) {
        die "Error, couldn't run cdbyank: $cmd, ret($?)\n";
    }

    unless ($fastaEntry) {
        die "Error, no fasta entry retrieved by accession: $accession\n";
    }

    return ($fastaEntry);
}


=over 4

=item linearize()

B<Description:> breaks down a fasta sequence into its components

B<Parameters:> fastaEntry

B<Returns:> (accession, header, linearSequence)

=back

=cut

    ;

sub linearize {
    my ($fastaEntry) = @_;

    unless ($fastaEntry =~ /^>/) {
        die "Error, fasta entry lacks expected format starting with header '>' character.\nHere's the entry\n$fastaEntry\n\n";
    }

    my @lines = split (/\n/, $fastaEntry);
    my $header = shift @lines;
    my $sequence = join ("", @lines);
    $sequence =~ s/\s+//g;

    $header =~ />(\S+)/;
    my $accession = $1;

    return ($accession, $header, $sequence);
}



=over 4

=item cdbyank_linear()

B<Description:> same as calling cdbyank (), and chasing it with linearize(), but only the sequence is returned

B<Parameters:> accession, fasta_db

B<Returns:> linearSequence

=back

=cut


    ;


sub cdbyank_linear {
    my ($acc, $fasta_db) = @_;

    my $fasta_entry = cdbyank($acc, $fasta_db);

    my ($acc2, $header, $genome_seq) = linearize($fasta_entry);

    return ($genome_seq);
}


=over 4

=item get_seq()

B<Description:> Retrieves a linear sequence from a FASTA file.
    Uses faidx_rust (Rust) with samtools .fai index for optimal performance.
    Falls back to cdbyank_linear if faidx_rust is not available.

B<Parameters:> accession, fastaFilename

B<Returns:> linearSequence (string, no header)

B<Note:> This is the preferred method for new code. It automatically
    creates the .fai index via 'samtools faidx' if it doesn't exist.

    The .fai index format (used by samtools faidx) is:
    - NAME  LENGTH  OFFSET  LINEBASES  LINEWIDTH

    This is simpler and more widely supported than the CDB .cidx format.
    For range extraction, use: get_seq_range($acc, $fasta, $start, $end)

=back

=cut

    ;

sub get_seq {
    my ($accession, $fastaFile) = @_;

    ## If faidx_rust is available, use it (preferred path)
    if ($FAIDX_RUST && -x $FAIDX_RUST) {

        ## Ensure .fai index exists
        unless (-s "$fastaFile.fai") {
            if ($SAMTOOLS && -x $SAMTOOLS) {
                my $ret = system("$SAMTOOLS faidx $fastaFile");
                die "Error creating .fai index: samtools faidx $fastaFile\n" if $ret;
            } else {
                die "Error: no .fai index for $fastaFile and samtools not in PATH\n";
            }
        }

        ## Use faidx_rust to retrieve the sequence
        my $cmd = "$FAIDX_RUST $fastaFile $accession";
        my $result = `$cmd`;

        if ($? || !$result) {
            die "Error: faidx_rust failed to retrieve '$accession' from $fastaFile\n";
        }

        ## Linearize: strip header, join sequence lines
        my @lines = split(/\n/, $result);
        shift @lines; ## remove header line
        my $seq = join("", @lines);
        $seq =~ s/\s+//g;

        return $seq;
    }

    ## Fall back to cdbyank_linear (uses CDB .cidx index)
    return cdbyank_linear($accession, $fastaFile);
}


=over 4

=item get_seq_range()

B<Description:> Retrieves a sub-range of a sequence from a FASTA file.
    Uses faidx_rust with samtools .fai index.

B<Parameters:> accession, fastaFilename, start, end

B<Returns:> linearSequence (sub-range, no header)

B<Note:> Coordinates are 1-based, inclusive.

=back

=cut

    ;

sub get_seq_range {
    my ($accession, $fastaFile, $start, $end) = @_;

    unless ($FAIDX_RUST && -x $FAIDX_RUST) {
        die "Error: get_seq_range requires faidx_rust in PATH\n";
    }

    unless (-s "$fastaFile.fai") {
        if ($SAMTOOLS && -x $SAMTOOLS) {
            my $ret = system("$SAMTOOLS faidx $fastaFile");
            die "Error creating .fai index\n" if $ret;
        } else {
            die "Error: no .fai index and samtools not in PATH\n";
        }
    }

    my $region = "${accession}:${start}-${end}";
    my $cmd = "$FAIDX_RUST $fastaFile $region";
    my $result = `$cmd`;

    if ($? || !$result) {
        die "Error: faidx_rust failed to retrieve range '$region'\n";
    }

    my @lines = split(/\n/, $result);
    shift @lines;
    my $seq = join("", @lines);
    $seq =~ s/\s+//g;

    return $seq;
}


1; #EOM
    
