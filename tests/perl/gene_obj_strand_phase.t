#!/usr/bin/env perl
#
# Regression test: Gene_obj exon order, CDS phase and CDS/protein sequence
# for minus-strand multi-exon genes.
#
# Commit f7cede7 ("first tier perl optimization") made two changes to
# PerlLib/Gene_obj.pm:
#   1. get_exons() cached its sorted list in $self->{_sorted_exons}. The
#      strand reversal was applied only when the cache was built, and nothing
#      cleared it. trivial_refinement() calls get_exons() before it sets the
#      strand, so a new minus-strand gene cached ascending (+ strand) order and
#      set_CDS_phases() then counted phases from the 3' end.
#   2. create_CDS_sequence()/create_cDNA_sequence() dropped the ascending sort
#      before joining the exons and reverse-complementing. With a 5'->3' exon
#      list on the minus strand, the CDS and protein came out scrambled.
# Effect seen with funannotate update (Absidia caerulea RSA 138, 2026-11-02):
# 3,535 tbl2asn InternalStop errors, and wrong PASA accept/reject decisions in
# cDNA_annotation_comparer.dbi (protein overlap computed on scrambled proteins).
#
# The gene structure follows FBN33_000004 (6 CDS segments, 630 bp CDS,
# expected phases 5'->3' = 0,2,1,0,0,0). The sequence is synthetic, so the
# test needs no data files.
#
# Run: prove -v tests/perl   (or: make perl-test)

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Test::More;
use Storable qw(nfreeze thaw);

use Gene_obj;
use Exons_to_geneobj;

my @CDS_LENS        = (49, 166, 220, 48, 60, 87);  # 5'->3', sum 630 = 210 codons
my @EXPECTED_PHASES = (0, 2, 1, 0, 0, 0);
my $INTRON          = 60;
my $FLANK           = 100;

sub revcomp {
    my ($s) = @_;
    $s = reverse $s;
    $s =~ tr/ACGTacgt/TGCAtgca/;
    return $s;
}

## Deterministic ORF: ATG + 208 sense codons + TAA. No stop codons inside.
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

## Intron: GT + filler with no ATG + AG. Filler 'C' avoids new ORFs in the cDNA
## (introns are spliced out anyway) and keeps splice sites canonical.
sub intron { return "GT" . ("C" x ($INTRON - 4)) . "AG"; }

## Build a contig carrying the gene on the given strand.
## Returns (contig_seq, exons_href end5=>end3, orf).
sub build_gene {
    my ($strand) = @_;
    my $orf = make_orf();
    my (@pieces, @exon_coords);   # coords on the forward (transcript-sense) string
    my $pos = $FLANK;             # 0-based length so far
    my $plus = ("C" x $FLANK);
    my $offset = 0;
    for my $i (0 .. $#CDS_LENS) {
        my $len = $CDS_LENS[$i];
        my $exon = substr($orf, $offset, $len);
        $offset += $len;
        my $start = length($plus) + 1;
        $plus .= $exon;
        my $end = length($plus);
        push @exon_coords, [$start, $end];
        $plus .= intron() if $i < $#CDS_LENS;
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
        # transcript-sense pos p maps to contig pos L-p+1; end5 > end3 on minus
        $exons{$L - $s + 1} = $L - $e + 1;
    }
    return ($minus_contig, \%exons, $orf);
}

## CDS phases from GFF3 output, in 5'->3' order.
sub gff3_cds_phases {
    my ($gene_obj, $strand) = @_;
    $gene_obj->{com_name} //= "test gene";  # to_GFF3_format warns without one
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

sub internal_stops {
    my ($prot) = @_;
    my $core = $prot;
    $core =~ s/\*$//;
    my $n = () = $core =~ /\*/g;
    return $n;
}

sub check_gene {
    my ($label, $gene_obj, $strand, $seq_ref, $orf) = @_;

    my @end5s = map { $_->{end5} } $gene_obj->get_exons();
    my @want  = sort { $a <=> $b } @end5s;
    @want = reverse @want if $strand eq '-';
    is_deeply(\@end5s, \@want, "$label: get_exons() is 5'->3' order");

    is_deeply([gff3_cds_phases($gene_obj, $strand)], \@EXPECTED_PHASES,
              "$label: GFF3 CDS phases 5'->3' are @EXPECTED_PHASES");

    $gene_obj->create_all_sequence_types($seq_ref);
    is(uc $gene_obj->get_CDS_sequence(), $orf, "$label: CDS sequence equals the ORF");
    my $prot = $gene_obj->get_protein_sequence();
    like($prot, qr/^M/, "$label: protein starts with M");
    is(internal_stops($prot), 0, "$label: protein has no internal stop");

    my $cdna = uc $gene_obj->get_cDNA_sequence();
    is($cdna, $orf, "$label: cDNA sequence equals the ORF (no UTR in this model)");
}

for my $strand ('-', '+') {
    my ($seq, $exons, $orf) = build_gene($strand);

    ## 1. As cDNA_annotation_comparer.dbi builds an update.
    my $g = Exons_to_geneobj::create_gene_obj($exons, \$seq, {});
    ok(defined $g, "strand $strand: gene object created");
    $g->{asmbl_id} = 'ctg'; $g->{TU_feat_name} = 'g1'; $g->{Model_feat_name} = 'g1.t1';
    is($g->get_orientation(), $strand, "strand $strand: orientation");
    check_gene("strand $strand fresh", $g, $strand, \$seq, $orf);

    ## 2. As the stored after_gene_obj: clone -> nfreeze/thaw -> refine
    ##    (read back by dump_valid_annot_updates.dbi and the comparer).
    my $c = $g->clone_gene();
    $c->{asmbl_id} = 'ctg'; $c->{TU_feat_name} = 'g1'; $c->{Model_feat_name} = 'g1.t1';
    my $t = thaw(nfreeze($c));
    $t->refine_gene_object();
    check_gene("strand $strand clone+thaw", $t, $strand, \$seq, $orf);

    ## 3. Phases recomputed on the thawed object must not change.
    $t->set_CDS_phases(\$seq);
    is_deeply([gff3_cds_phases($t, $strand)], \@EXPECTED_PHASES,
              "strand $strand clone+thaw: set_CDS_phases() keeps correct phases");
}

## 4. get_exons() must follow later changes to the exon list and the strand.
{
    my ($seq, $exons) = build_gene('-');
    my $g = Exons_to_geneobj::create_gene_obj($exons, \$seq, {});
    my $n_before = scalar($g->get_exons());
    my @exon_objs = $g->get_exons();
    $g->{mRNA_exon_objs} = [ @exon_objs[0 .. 2] ];   # replace the list directly
    is(scalar($g->get_exons()), 3, "get_exons() follows a replaced exon list (was $n_before)");
}

done_testing();
