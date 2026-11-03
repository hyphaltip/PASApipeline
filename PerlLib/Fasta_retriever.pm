package Fasta_retriever;

use strict;
use warnings;
use Carp;
use threads;
use threads::shared;

my $LOCKVAR :shared;
our $DEBUG = 0;

## Compressed FASTA: there is no random access, so each get_seq() reads the
## stream from the start up to the wanted header (slow for many lookups; use
## an uncompressed FASTA for speed). File offsets, including a bgzip .fai,
## are not used for compressed input.
my %DECOMPRESS_CMD = (
    '.gz'  => 'zcat',
    '.bz2' => 'bzcat',
);
my %WARNED_COMPRESSED;

## Returns the decompress command for a compressed file, else undef.
## Detects by extension, then by magic bytes (gzip/bgzip 1f 8b, bzip2 'BZh'),
## so a compressed file without a .gz/.bz2 name is also read correctly.
my %COMPRESSED_CACHE;
sub _is_compressed {
    my $filename = shift;
    return $COMPRESSED_CACHE{$filename} if exists $COMPRESSED_CACHE{$filename};
    my $cmd;
    for my $ext (keys %DECOMPRESS_CMD) {
        if ($filename =~ /\Q$ext\E$/) { $cmd = $DECOMPRESS_CMD{$ext}; last; }
    }
    if (!$cmd && open(my $fh, "<:raw", $filename)) {
        my $magic = "";
        read($fh, $magic, 3);
        close $fh;
        if    ($magic =~ /^\x1f\x8b/) { $cmd = 'zcat'; }
        elsif ($magic eq 'BZh')      { $cmd = 'bzcat'; }
    }
    return $COMPRESSED_CACHE{$filename} = $cmd;
}

sub _open_compressed {
    my $filename = shift;
    my $cmd = _is_compressed($filename) or return (undef, 0);
    open(my $fh, "-|", $cmd, $filename)
        or die "Error, cannot open compressed file: $filename";
    unless ($WARNED_COMPRESSED{$filename}++) {
        print STDERR "WARNING: Fasta_retriever: $filename is compressed; each sequence "
            . "lookup reads the file from the start. Use an uncompressed FASTA for speed.\n";
    }
    return ($fh, 1);
}

sub new {
    my ($packagename) = shift;
    my $filename = shift;
    
    unless ($filename) {
        confess "Error, need filename as param";
    }

    my $self = { filename => $filename,
                 acc_to_pos_index => undef,
                 fh => undef,
    };
    
    my %acc_to_pos_index :shared;

    $self->{acc_to_pos_index} = \%acc_to_pos_index;
        
    bless ($self, $packagename);

    $self->_init();


    return($self);
}


sub _init {
    my $self = shift;
    
    my $filename = $self->{filename};

    if (_is_compressed($filename)) {
        ## index accessions only; get_seq() scans the stream (see top of file)
        my ($fh) = _open_compressed($filename);
        my $n = 0;
        while (<$fh>) {
            if (/^>(\S+)/) {
                $self->{acc_to_pos_index}->{$1} = ++$n;  # order, not a file offset
            }
        }
        close $fh;
        return;
    }

    # use a samtools faidx index if available
    my $index_file = "$filename.fai";
    if (-s $index_file) {
        open(my $fh, $index_file) or die "Error, cannot open file: $index_file";
        while(<$fh>) {
            chomp;
            my @x = split(/\t/);
            my $acc = $x[0];
            my $file_pos = $x[2];
            $self->{acc_to_pos_index}->{$acc} = $file_pos;
        }
        close $fh;
    }
    else {
        print STDERR "-missing faidx file: $index_file, extracting positions directly.\n";
        print STDERR "-Fasta_retriever:: begin initializing for $filename\n";
    
        open (my $fh, $filename) or die $!;
        $self->{fh} = $fh;
        while (<$fh>) {
            if (/>(\S+)/) {
                my $acc = $1;
                my $file_pos = tell($fh);
                $self->{acc_to_pos_index}->{$acc} = $file_pos;
            }
        }
        print STDERR "-Fasta_retriever:: done initializing for $filename\n";
    }
    
    return;
}

sub refresh_fh {
    my $self = shift;
    
    open (my $fh, $self->{filename}) or die "Error, cannot open file : " . $self->{filename};
    $self->{fh} = $fh;
    
    return $fh;
}


sub get_seq {
    my $self = shift;
    my $acc = shift;

    unless (defined $acc) {
        confess "Error, need acc as param";
    }

    {
        lock $LOCKVAR;
    
        my $file_pos = $self->{acc_to_pos_index}->{$acc} or confess "Error, no seek pos for acc: $acc";
        
        my ($fh, $is_compressed) = _open_compressed($self->{filename});
        
        if ($is_compressed) {
            ## no seek on a stream: read up to the wanted header
            my $found = 0;
            while (<$fh>) {
                if (/^>(\S+)/ && $1 eq $acc) { $found = 1; last; }
            }
            confess "Error, accession $acc not found in compressed file $self->{filename}"
                unless $found;
        }
        else {
            $fh = $self->{fh};
            unless ($fh && fileno($fh)) {
                $fh = $self->refresh_fh();
            }
            seek($fh, $file_pos, 0);
        }

        print STDERR "seeking $acc -> $file_pos\n" if $DEBUG;

        my @seq_lines;
        while (<$fh>) {
            if (/^>/) {
                print STDERR "   reached $_, stopping\n" if $DEBUG;
                last;
            }
            push @seq_lines, $_;
        }
        close $fh if $is_compressed;
        print STDERR "-done seeking $acc\n\n" if $DEBUG;

        my $seq = join('', @seq_lines);
        $seq =~ s/\s+//g;
        ## An empty sequence means a wrong offset (e.g. a stale .fai) or an
        ## empty record; never return it silently.
        unless (length $seq) {
            confess "Error, empty sequence retrieved for $acc from $self->{filename} "
                . "(stale or wrong .fai index, or empty FASTA record?)";
        }
        return $seq;
    }
}
    
    
    

1; #EOM
    
