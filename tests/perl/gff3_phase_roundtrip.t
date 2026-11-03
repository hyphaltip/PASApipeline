#!/usr/bin/env perl
## Regression test: CDS phase convention between GFF3_utils (reader),
## Gene_obj::set_CDS_phases (computer) and Gene_obj::to_GFF3_format (writer).
##
## Gene_obj {phase} = codon position of the first CDS base (0,1,2).
## GFF3 column 8   = bases to skip before the first full codon = (3 - {phase}) % 3.
##
## Section 1: parse a GFF3 with phases 1 and 2 (+ and - strand) and write it
##            back; the phases must be unchanged.
## Section 2: a 5'-partial model whose first CDS base is codon position 2 or 3;
##            set_CDS_phases + to_GFF3_format must give GFF3 phase 1 / 2.
##
## Run: prove -v tests/perl   (or: make perl-test)
use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib", "$FindBin::Bin/lib");
use Test::More;
use File::Temp qw(tempfile);
use Gene_obj;
use GFF3_utils;
use SynthGene;

diag("Gene_obj.pm: $INC{'Gene_obj.pm'}");

## Section 1: GFF3 round trip
for my $strand ('+', '-') {
    my ($seq, $exons) = SynthGene::build_gene($strand, 0);
    my @in = SynthGene::expected_gff3_phases(0);        # (0 2 1 0 0 0): has a 1 and a 2
    my ($fh, $fn) = tempfile(SUFFIX => ".gff3", UNLINK => 1);
    print $fh SynthGene::gff3_text($exons, $strand, \@in, "g");
    close $fh;
    my %genes;
    GFF3_utils::index_GFF3_gene_objs($fn, \%genes);
    my ($g) = values %genes;
    is($g->get_orientation(), $strand, "strand $strand: parsed orientation");
    my @out = SynthGene::gff3_cds_phases($g);
    is_deeply(\@out, \@in, "strand $strand: GFF3 parse -> to_GFF3_format keeps CDS phases (@in)");
    my @int = SynthGene::internal_phases($g);
    is_deeply(\@int, [ map { (3 - $_) % 3 } @in ],
              "strand $strand: parsed {phase} is codon position (1<->2 of GFF3 phase)");
}

## Section 2: 5'-partial model through set_CDS_phases
for my $strand ('+', '-') {
    for my $trim (1, 2) {
        my ($seq, $exons, $orf) = SynthGene::build_gene($strand, $trim);
        my $g = new Gene_obj();
        $g->populate_gene_obj({%$exons}, {%$exons});
        $g->{asmbl_id} = 'ctg'; $g->{TU_feat_name} = 'g'; $g->{Model_feat_name} = 'g-T1';
        $g->set_CDS_phases(\$seq);
        is(uc $g->get_CDS_sequence(), $orf, "strand $strand trim $trim: CDS sequence is the trimmed ORF");
        ok($g->is_5prime_partial(), "strand $strand trim $trim: flagged 5' partial");
        my @exp = SynthGene::expected_gff3_phases($trim);
        my @out = SynthGene::gff3_cds_phases($g);
        is_deeply(\@out, \@exp, "strand $strand trim $trim: set_CDS_phases -> to_GFF3_format gives (@exp)");
        ## the first segment's phase must also be in the same convention as the later ones
        my @int = SynthGene::internal_phases($g);
        is($int[0], (3 - $exp[0]) % 3, "strand $strand trim $trim: first {phase} is codon position, like later segments");
    }
}

done_testing();
