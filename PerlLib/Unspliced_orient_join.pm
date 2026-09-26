package Unspliced_orient_join;

## Opt-in helper for Launch_PASA_pipeline.pl --UNSPLICED_JOIN_SPLICED.
##
## Unspliced (single-segment) alignments get spliced_orient '?', and PASA
## clusters and subclusters by spliced orientation, so a '?' alignment that
## sits inside a spliced gene never joins it and becomes a separate
## single-exon locus. This module decides, for each '?' alignment, whether
## it can be given the orientation of the spliced alignments that cover it:
## it is assigned only when every spliced alignment covering at least
## $min_pct percent of its span has the same orientation. If it is covered by
## both orientations, or by none, it stays '?'.

use strict;
use warnings;

sub assign_unspliced_orientations {
    my ($alignments_aref, $min_pct) = @_;
    ## each alignment: { id, scaffold, orient ('+', '-' or '?'), lend, rend }
    my %by_scaffold;
    foreach my $aln (@$alignments_aref) {
        push (@{$by_scaffold{$aln->{scaffold}}}, $aln);
    }
    my %assigned;
    foreach my $scaffold (keys %by_scaffold) {
        my @spliced = sort {$a->{lend} <=> $b->{lend}}
                      grep { $_->{orient} eq '+' || $_->{orient} eq '-' } @{$by_scaffold{$scaffold}};
        my @unspliced = grep { $_->{orient} eq '?' } @{$by_scaffold{$scaffold}};
        next unless (@spliced && @unspliced);
        my @starts = map { $_->{lend} } @spliced;
        my $max_span = 0;
        foreach my $s (@spliced) {
            my $span = $s->{rend} - $s->{lend};
            $max_span = $span if $span > $max_span;
        }
        foreach my $u (@unspliced) {
            my $u_len = $u->{rend} - $u->{lend} + 1;
            next unless $u_len > 0;
            my %orients;
            for (my $i = &_lower_bound(\@starts, $u->{lend} - $max_span);
                 $i <= $#spliced && $spliced[$i]->{lend} <= $u->{rend}; $i++) {
                my $s = $spliced[$i];
                next if $s->{rend} < $u->{lend};
                my $overlap = &_min($s->{rend}, $u->{rend}) - &_max($s->{lend}, $u->{lend}) + 1;
                if ($overlap / $u_len * 100 >= $min_pct) {
                    $orients{$s->{orient}} = 1;
                }
            }
            my @o = keys %orients;
            $assigned{$u->{id}} = $o[0] if (scalar(@o) == 1);
        }
    }
    return \%assigned;
}

## Orientation test for pairing two alignments/assemblies: equal
## orientations pair; a '?' pairs with '+' or '-' only if it was assigned
## that orientation.
sub orientations_compatible {
    my ($orient_a, $id_a, $orient_b, $id_b, $assigned_href) = @_;
    return 1 if $orient_a eq $orient_b;
    my $eff_a = ($orient_a eq '?' && exists $assigned_href->{$id_a}) ? $assigned_href->{$id_a} : $orient_a;
    my $eff_b = ($orient_b eq '?' && exists $assigned_href->{$id_b}) ? $assigned_href->{$id_b} : $orient_b;
    return ($eff_a eq $eff_b) ? 1 : 0;
}

sub _lower_bound {
    my ($aref, $value) = @_;
    my ($lo, $hi) = (0, scalar(@$aref));
    while ($lo < $hi) {
        my $mid = int(($lo + $hi) / 2);
        if ($aref->[$mid] < $value) { $lo = $mid + 1; } else { $hi = $mid; }
    }
    return $lo;
}

sub _min { return ($_[0] < $_[1]) ? $_[0] : $_[1]; }
sub _max { return ($_[0] > $_[1]) ? $_[0] : $_[1]; }

1;
