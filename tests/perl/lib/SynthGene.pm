package SynthGene;
## Synthetic multi-exon gene builder for phase tests (no data files needed).
## Segment lengths 5'->3' = (49, 166, 220, 48, 60, 87) = 630 bp = 210 codons.
## trim5 = number of bases removed from the 5' end of the first exon
## (0 = complete ORF with ATG; 1 or 2 = 5'-partial model).
use strict;
use warnings;

our @CDS_LENS = (49, 166, 220, 48, 60, 87);
our $INTRON   = 60;
our $FLANK    = 100;

sub revcomp {
    my ($s) = @_;
    $s = reverse $s;
    $s =~ tr/ACGTacgt/TGCAtgca/;
    return $s;
}

sub make_orf {
    srand(42);
    my @bases = qw(A C G T);
    my @codons;
    while (@codons < 208) {
        my $c = join("", map { $bases[int(rand(4))] } 1 .. 3);
        next if $c =~ /^(TAA|TAG|TGA)$/;
        push @codons, $c;
    }
    return "ATG" . join("", @codons) . "TAA";
}

sub intron { return "GT" . ("C" x ($INTRON - 4)) . "AG"; }

## returns (contig_seq, exons_href end5=>end3, orf_after_trim)
sub build_gene {
    my ($strand, $trim5) = @_;
    $trim5 ||= 0;
    my $orf = make_orf();
    $orf = substr($orf, $trim5);
    my @lens = @CDS_LENS;
    $lens[0] -= $trim5;
    my (@exon_coords);
    my $plus = ("C" x $FLANK);
    my $offset = 0;
    for my $i (0 .. $#lens) {
        my $exon = substr($orf, $offset, $lens[$i]);
        $offset += $lens[$i];
        my $start = length($plus) + 1;
        $plus .= $exon;
        push @exon_coords, [$start, length($plus)];
        $plus .= intron() if $i < $#lens;
    }
    $plus .= ("C" x $FLANK);
    my %exons;
    if ($strand eq '+') {
        $exons{$_->[0]} = $_->[1] for @exon_coords;
        return ($plus, \%exons, $orf);
    }
    my $L = length($plus);
    my $minus_contig = revcomp($plus);
    for my $c (@exon_coords) {
        my ($s, $e) = @$c;
        $exons{$L - $s + 1} = $L - $e + 1;
    }
    return ($minus_contig, \%exons, $orf);
}

## expected GFF3 phases 5'->3' given trim5 (bases to skip in first segment)
sub expected_gff3_phases {
    my ($trim5) = @_;
    $trim5 ||= 0;
    my @lens = @CDS_LENS;
    $lens[0] -= $trim5;
    my $skip0 = (3 - $trim5 % 3) % 3;
    my @ph;
    my $consumed = -$skip0;
    for my $len (@lens) {
        push @ph, (3 - ($consumed % 3)) % 3;
        $consumed += $len;
    }
    return @ph;
}

## CDS phases from to_GFF3_format output, 5'->3'
sub gff3_cds_phases {
    my ($gene_obj) = @_;
    $gene_obj->{com_name} //= "test gene";
    my $strand = $gene_obj->get_orientation();
    my @cds;
    for my $line (split /\n/, $gene_obj->to_GFF3_format()) {
        my @x = split /\t/, $line;
        next unless @x > 7 && $x[2] eq 'CDS';
        push @cds, [$x[3], $x[7]];
    }
    @cds = sort { $a->[0] <=> $b->[0] } @cds;
    @cds = reverse @cds if $strand eq '-';
    return map { $_->[1] } @cds;
}

## internal {phase} values 5'->3'
sub internal_phases {
    my ($gene_obj) = @_;
    my @p;
    for my $exon ($gene_obj->get_exons()) {
        if (my $cds = $exon->get_CDS_obj()) {
            push @p, $cds->{phase};
        }
    }
    return @p;
}

## a GFF3 text for a gene with the given exons and GFF3 phases (5'->3')
sub gff3_text {
    my ($exons_href, $strand, $phases_aref, $id) = @_;
    my @ends = sort { $a <=> $b } keys %$exons_href;
    @ends = reverse @ends if $strand eq '-';
    my ($gl, $gr) = (sort { $a <=> $b } (map { ($_, $exons_href->{$_}) } @ends))[0, -1];
    my $t = "ctg\tsrc\tgene\t$gl\t$gr\t.\t$strand\t.\tID=$id;Name=$id\n";
    $t .= "ctg\tsrc\tmRNA\t$gl\t$gr\t.\t$strand\t.\tID=$id-T1;Parent=$id\n";
    my $i = 0;
    for my $e5 (@ends) {
        my ($l, $r) = sort { $a <=> $b } ($e5, $exons_href->{$e5});
        $t .= "ctg\tsrc\texon\t$l\t$r\t.\t$strand\t.\tID=$id-T1.exon$i;Parent=$id-T1\n";
        $t .= "ctg\tsrc\tCDS\t$l\t$r\t.\t$strand\t$phases_aref->[$i]\tID=$id-T1.cds;Parent=$id-T1\n";
        $i++;
    }
    return $t;
}

1;
