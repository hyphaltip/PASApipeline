#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/..");
use Test::More;
use Unspliced_orient_join;

sub aln { my ($id, $sc, $o, $l, $r) = @_; return { id => $id, scaffold => $sc, orient => $o, lend => $l, rend => $r }; }

## single-exon read inside a plus-strand spliced gene -> '+'
{
    my $a = Unspliced_orient_join::assign_unspliced_orientations(
        [ aln('g1', 'c1', '+', 1000, 3000), aln('u1', 'c1', '?', 1200, 1500) ], 30);
    is($a->{u1}, '+', 'unspliced inside plus gene gets +');
}

## read covered by genes on both strands -> stays '?'
{
    my $a = Unspliced_orient_join::assign_unspliced_orientations(
        [ aln('g1', 'c1', '+', 1000, 3000), aln('g2', 'c1', '-', 1400, 5000),
          aln('u1', 'c1', '?', 1500, 1800) ], 30);
    ok(!exists $a->{u1}, 'unspliced covered by both orientations stays ?');
}

## overlap below threshold does not count
{
    my $a = Unspliced_orient_join::assign_unspliced_orientations(
        [ aln('g1', 'c1', '+', 1000, 2000), aln('u1', 'c1', '?', 1900, 2900) ], 30);
    ok(!exists $a->{u1}, '10% overlap below 30% threshold stays ?');
}

## threshold is inclusive and measured on the unspliced span
{
    my $a = Unspliced_orient_join::assign_unspliced_orientations(
        [ aln('g1', 'c1', '-', 1000, 2000), aln('u1', 'c1', '?', 1701, 2700) ], 30);
    is($a->{u1}, '-', 'exactly 30% of unspliced span covered -> assigned');
}

## a weak overlap with the other strand does not block assignment
{
    my $a = Unspliced_orient_join::assign_unspliced_orientations(
        [ aln('g1', 'c1', '+', 1000, 3000), aln('g2', 'c1', '-', 2950, 6000),
          aln('u1', 'c1', '?', 2000, 2999) ], 30);
    is($a->{u1}, '+', 'minus gene covers only 5% -> ignored');
}

## different scaffolds never interact
{
    my $a = Unspliced_orient_join::assign_unspliced_orientations(
        [ aln('g1', 'c1', '+', 1000, 3000), aln('u1', 'c2', '?', 1200, 1500) ], 30);
    ok(!exists $a->{u1}, 'no assignment across scaffolds');
}

## long spliced gene that starts far upstream is still found (window search)
{
    my @alns = (aln('g1', 'c1', '+', 100, 50000));
    push @alns, aln("s$_", 'c1', '-', 60000 + $_ * 10, 60005 + $_ * 10) for (1 .. 50);
    push @alns, aln('u1', 'c1', '?', 40000, 40500);
    my $a = Unspliced_orient_join::assign_unspliced_orientations(\@alns, 30);
    is($a->{u1}, '+', 'long upstream-starting gene found');
}

## spliced alignments are never reassigned
{
    my $a = Unspliced_orient_join::assign_unspliced_orientations(
        [ aln('g1', 'c1', '+', 1000, 3000), aln('g2', 'c1', '-', 1500, 2500) ], 30);
    is_deeply($a, {}, 'only ? alignments are assigned');
}

## orientations_compatible
{
    my $assigned = { u1 => '+' };
    ok(Unspliced_orient_join::orientations_compatible('+', 'a', '+', 'b', $assigned), '+ with +');
    ok(!Unspliced_orient_join::orientations_compatible('+', 'a', '-', 'b', $assigned), '+ with -');
    ok(Unspliced_orient_join::orientations_compatible('?', 'u1', '+', 'b', $assigned), 'assigned ? with +');
    ok(!Unspliced_orient_join::orientations_compatible('?', 'u1', '-', 'b', $assigned), 'assigned ? with -');
    ok(!Unspliced_orient_join::orientations_compatible('?', 'u2', '+', 'b', $assigned), 'unassigned ? with +');
    ok(Unspliced_orient_join::orientations_compatible('?', 'u2', '?', 'u3', $assigned), '? with ? (unchanged behaviour)');
}

done_testing();
