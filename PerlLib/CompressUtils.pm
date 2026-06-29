package CompressUtils;

use strict;
use warnings;
use Carp;
use File::Basename;

## Utility functions for transparent file compression/decompression.
## Supports .gz (gzip) and .bz2 (bzip2) formats.

## Detect if a file is compressed based on extension
sub is_compressed {
    my ($file) = @_;
    return $file =~ /\.(gz|bz2)$/;
}

## Open a file for reading, transparently decompressing .gz/.bz2 files.
## Returns a filehandle.
sub open_read {
    my ($file) = @_;
    my $fh;
    if ($file =~ /\.gz$/) {
        open($fh, "-|", "gzip -dc $file") or confess "Error, cannot gzip -dc $file: $!";
    } elsif ($file =~ /\.bz2$/) {
        open($fh, "-|", "bzip2 -dc $file") or confess "Error, cannot bzip2 -dc $file: $!";
    } else {
        open($fh, "<", $file) or confess "Error, cannot open $file: $!";
    }
    return $fh;
}

## Open a file for writing, transparently compressing to .gz/.bz2.
## The output compression format is determined by the file extension.
## Returns a filehandle.
sub open_write {
    my ($file) = @_;
    my $fh;
    if ($file =~ /\.gz$/) {
        open($fh, "|-", "gzip -c > $file") or confess "Error, cannot write gzip to $file: $!";
    } elsif ($file =~ /\.bz2$/) {
        open($fh, "|-", "bzip2 -c > $file") or confess "Error, cannot write bzip2 to $file: $!";
    } else {
        open($fh, ">", $file) or confess "Error, cannot write to $file: $!";
    }
    return $fh;
}

## Compress a file in-place, adding .gz extension.
## Returns the new filename.
sub gzip_file {
    my ($file) = @_;
    return $file if $file =~ /\.gz$/;
    my $gz_file = "$file.gz";
    system("gzip -c $file > $gz_file") == 0
        or confess "Error compressing $file: $!";
    unlink($file);
    return $gz_file;
}

## Decompress a .gz file in-place, removing .gz extension.
## Returns the new filename.
sub gunzip_file {
    my ($file) = @_;
    return $file unless $file =~ /\.gz$/;
    my $out_file = $file;
    $out_file =~ s/\.gz$//;
    system("gzip -dc $file > $out_file") == 0
        or confess "Error decompressing $file: $!";
    unlink($file);
    return $out_file;
}

1;
