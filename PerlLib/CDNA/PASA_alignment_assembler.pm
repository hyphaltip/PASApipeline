#!/usr/local/bin/perl

package main;
our $SEE;


package CDNA::PASA_alignment_assembler;

=head1 NAME

CDNA::PASA_alignment_assembler

=cut

=head1 DESCRIPTION

This module is used to assemble compatible cDNA alignments.  The algorithm is as follows:
must describe this here.

=cut

use strict;
use CDNA::CDNA_alignment;
use Data::Dumper;
use Carp;
use File::Temp qw(tempfile);
use Pasa_tmpdir;

## File scoped globals:
my $DELIMETER = "$;,";
our $FUZZLENGTH = 20;


=item new()

=over 4

B<Description:> instantiates a new cDNA assembler obj.

B<Parameters:> none

B<Returns:> $obj_href

$obj_href is the object reference newly instantiated by this new method.

=back

=cut

sub new  {
    my $package_name = shift;
    my $self = {};
    bless ($self, $package_name);
    $self->_init(@_);
    return ($self);
}


sub _init {
    my $self = shift;
    $self->{incoming_alignments} = []; #these are the alignments to be assembled.
    $self->{assemblies} = []; #contains list of all singletons and assemblies.
    $self->{fuzzlength} = $FUZZLENGTH;  #default setting.
    
    my $pasa_bin = `which pasa_rust 2>/dev/null`;
    $pasa_bin =~ s/\s//g;
    unless ($pasa_bin && -x $pasa_bin) {
        $pasa_bin = `which pasa 2>/dev/null`;
        $pasa_bin =~ s/\s//g;
    }
    
    unless (-x $pasa_bin) {
        confess "Error, pasa binary [$pasa_bin] isn't executable or couldn't be found.";
    }
    
    $self->{pasa_bin} = $pasa_bin;

}


=item assemble_alignments()
    
=over 4
    
B<DESCRIPTION:> assembles a series of cDNA aligmnments into one or more cDNA assemblies using a directed acyclic graph.

B<Parameters:> @alignments

@alignments is an array of CDNA::CDNA_alignment objects

B<Returns:> none.

=back

=cut


sub assemble_alignments {
    my $self = shift;
    my @alignments = @_;
    @alignments = sort {$a->{lend}<=>$b->{lend}} @alignments; #keep in order of lend across genomic sequence to provide a layout.
    $self->{incoming_alignments} = [@alignments];
    my %accs;
    my %spliced_orientations; # track so can set later in each assembly based on content.
    my %aligned_orientations;
    my %FL_accs;
    
    foreach my $alignment (@alignments) {
        my $acc = $alignment->get_acc();
        $accs{$acc} = 0;
        my $spliced_orient = $alignment->get_spliced_orientation();
        $spliced_orientations{$acc} = $spliced_orient;
        my $aligned_orient = $alignment->get_orientation() or confess "Error, no orientation set for " . $alignment->toToken();
        $aligned_orientations{$acc} = $aligned_orient;
        $FL_accs{$acc} = $alignment->is_fli();
    }
    $self->{accs_in_assemblies} = \%accs;
    
    my $num_alignments = $#alignments + 1;
    
    $self->force_flexorient('+');
    my @assemblies = $self->pasa_cpp_assemblies('+');
    $self->force_flexorient('-');
    push (@assemblies, $self->pasa_cpp_assemblies('-'));
    
    # sort in order of decreasing score
    @assemblies = reverse sort {$a->{num_contained_aligns}<=>$b->{num_contained_aligns}} @assemblies;  
    
    if ($SEE) {
        print "\n\nScore summary for all assemblies (nr set unchosen):\n\n";
        foreach my $assembly (@assemblies) {
            print "score: " . $assembly->{num_contained_aligns} . ", " . $assembly->toToken . "\n";
        }
        print "\n\n";
    }
    my @report_assemblies;
    my $still_missing = 1;
    foreach my $assembly (@assemblies) {
        print "\nAnalyzing assembly.\n" if $SEE;
        my $contained_aligns_aref = $assembly->{contained_aligns};
        my $have_unseen = 0;
        my $spliced_orient = '?';
        my $is_fli = 0;
        my %aligned_orient_counts;

        foreach my $acc (@$contained_aligns_aref) {
            print "got: $acc\n" if $SEE;
            # check to see if we've encountered this one yet.
            unless ($accs{$acc}) {
                $have_unseen = 1;
            }
            $accs{$acc} = 1;
            my $curr_spliced_orient = $spliced_orientations{$acc};
            print "curr_spliced_orient: $curr_spliced_orient\n" if $SEE;
            if ($curr_spliced_orient ne '?') {
                if ($spliced_orient ne '?' && $spliced_orient ne $curr_spliced_orient) {
                    ## cannot have conflicting spliced orientations in the assembly: corruption.
                    die "Fatal: conflicting spliced orientations in current PASA assembly results (spliced_orient: $spliced_orient, $acc has $curr_spliced_orient).\n";
                }
                $spliced_orient = $curr_spliced_orient; ## retain original spliced orientation.
            }
            if ( (!$is_fli) && $FL_accs{$acc}) {
                $is_fli = 1;
            }
            ## track aligned orientation
            $aligned_orient_counts{ $aligned_orientations{$acc} } ++;

        }
        
        ## set assembly orientation
        $assembly->set_spliced_orientation($spliced_orient);
        if ($spliced_orient eq '?') {
            ## set aligned orientation based on a majority vote
            my @orients = reverse sort { $aligned_orient_counts{$a} <=> $aligned_orient_counts{$b} } keys %aligned_orient_counts;
            my $winning_aligned_orient = shift @orients;
            $assembly->set_orientation($winning_aligned_orient);
        }
        else {
            # got spliced orientation, use it for aligned orientation too.
            $assembly->set_orientation($spliced_orient);
        }
        
        $assembly->set_fli_status($is_fli);
        
        if ($have_unseen) {
            print $assembly->toToken . "\n" if $SEE;
            push (@report_assemblies, $assembly);
        }
        $still_missing = 0;
        foreach my $key (keys %accs) {
            if (! $accs{$key}) {
                $still_missing = 1;
                print "still missing: $key\n" if $SEE;
            }
        }
        if (! $still_missing) {
            last; #got them all.
        }
    }
    
    $self->{assemblies} = \@report_assemblies;
    
    if ($still_missing) {
        die "Didn't obtain assemblies describing all maximal assemblies.\n";
    }
    
    
}


=head2 pasa_cpp_assemblies

Hands the incoming alignments to the external C++/Rust C<pasa> binary and parses
the assemblies back.

B<Temp-file allocation:> this method exchanges data with the child C<pasa>
process through files in a shared temp directory, so those filenames are the
only thing isolating one concurrent caller from another. They are allocated with
C<File::Temp::tempfile>, which creates each file atomically (C<O_EXCL>) and
retries on collision.

Please do not "optimize" this back into a constructed filename. From the 2015
initial import until 2026-07 the names were built as C<time() . "-" . rand()>.
Testing on the funannotate Perl showed that token does in fact differ per
ithread, so it was not observed to collide -- but constructing a name and then
opening it assumes exclusivity, whereas C<O_EXCL> verifies it. That distinction
matters most where a constructed token cannot help at all: a shared or network
C<$TMPDIR>, where pids repeat across hosts. The same reasoning is why this uses
atomic creation rather than a UUID, which would also only lower the probability
of a collision rather than eliminate it.

See docs/THREAD_SAFETY.md, which also records a threading failure that was
initially misattributed to this code and turned out to be a benchmark harness
artifact.

=cut

sub pasa_cpp_assemblies {
    my $self = shift;
    my $forced_orient = shift;

    my $prev_input_sep = $/;
    $/ = "\n";

    my $sequence_ref;
    my $incoming_alignments_aref = $self->{incoming_alignments};

    ## PASA_TMPDIR > TMPDIR > SCRATCH > /tmp > . , each checked for writability.
    ## On a cluster this lands on per-job node-local storage instead of a
    ## node-wide /tmp; see PerlLib/Pasa_tmpdir.pm.
    my $tmpdir = Pasa_tmpdir::get_tmpdir();

    ## These temp files are the ONLY isolation between concurrent assemblers,
    ## and $tmpdir is shared by every thread, every process, and -- if TMPDIR
    ## lives on a network filesystem -- every node.
    ##
    ## The original name was time() . "-" . rand(), which carries no pid and
    ## no thread id. srand() seeds partly from the pid, and the pid is
    ## IDENTICAL across ithreads of one process, so two threads entering here
    ## in the same second could derive the same token and silently clobber
    ## each other's pasa input/output. assemble_clusters.dbi defaults to
    ## -T 2, so this raced by default.
    ##
    ## File::Temp is used rather than a longer/fancier token because it
    ## creates the file atomically with O_EXCL and retries on collision:
    ## uniqueness is enforced by the filesystem instead of assumed from
    ## entropy. That also covers the cross-node case that a pid+tid token
    ## cannot, since pids repeat across hosts.
    my ($in_fh, $pasa_input) = tempfile("pasa.XXXXXXXXXX", DIR => $tmpdir,
                                        SUFFIX => ".$forced_orient.in", UNLINK => 0);
    my ($out_fh, $pasa_output) = tempfile("pasa.XXXXXXXXXX", DIR => $tmpdir,
                                          SUFFIX => ".$forced_orient.out", UNLINK => 0);
    close $out_fh;  ## the pasa binary writes this itself via shell redirection

    my @assemblies;
    foreach my $alignment (@$incoming_alignments_aref) {
        my $acc = $alignment->get_acc();
        ## commas not allowed in acc name:
        if ($acc =~ /,/) { 
            die "ERROR, $acc accession contains comma(s).  This is not allowed.\n";
        }
        my $orient = $alignment->{fixed_orient};
        my $alignText = "$acc,$orient";
        unless (ref $sequence_ref) {
            $sequence_ref = $alignment->get_genomic_seq_ref();
        }
        foreach my $seg ($alignment->get_alignment_segments()) {
            my ($lend, $rend) = $seg->get_coords();
            $alignText .= ",$lend-$rend";
        }
        print $in_fh $alignText . "\n";
    }
    close $in_fh;
    
    if ($SEE) {
        print "PASA_INPUT ($forced_orient):\n====\n";
        system "cat $pasa_input";
        print "====\n";
    }
    my $cmd = $self->{pasa_bin} . " $pasa_input > $pasa_output";
    my $ret = system $cmd;
    if ($ret) {
        ## 'pasa_killer.input' was a fixed name in cwd, so concurrent failing
        ## threads overwrote each other's diagnostic. Keep one per failure.
        my (undef, $killer) = tempfile("pasa_killer.XXXXXXXXXX", DIR => $tmpdir,
                                       SUFFIX => ".input", UNLINK => 0);
        rename($pasa_input, $killer)
            or warn "Error, cannot preserve failing input $pasa_input: $!";
        unlink($pasa_output);
        print STDERR "PASA died on input file.  See $killer\n";
        die "PASA assembler failed; offending input preserved at $killer\n";
    } else {

        # process the output.
        open (my $out_read_fh, '<', $pasa_output) or die "Can't open $pasa_output";
        while (<$out_read_fh>) {
            if (/assembly:\s\(\d+\)\scontains\salignments:\s\[([^\]]+)\]\swith\sstructure\s\[([^\]]+)\]/) {
                print "Extracting assembly output: $_" if $SEE;
                my $acclist = $1;
                my $aligndescript = $2;
                my @x = split (/,/, $aligndescript);
                
                shift @x;
                my $orient = shift @x;
                my @alignSegs;
                my $length = 0;
                foreach my $coordset (@x) {
                    my ($lend, $rend) = sort {$a<=>$b} split (/-/, $coordset);
                    my $seg = new CDNA::Alignment_segment($lend, $rend);
                    $length += ($rend - $lend) + 1;
                    push (@alignSegs, $seg);
                }
                my $assembly = new CDNA::CDNA_alignment($length, \@alignSegs, $sequence_ref);
                
                my @accs = split (/,/, $acclist);
                my $num_accs = $#accs + 1;
                $assembly->{contained_aligns} = [@accs];
                $assembly->{num_contained_aligns} = $num_accs;
                
                $acclist =~ s/,/\//g; #convert list of accessions into a new accession representing a single entry (unity)
                # if we keep the commas, use of this assembly in future PASA runs will break the assembler
                # because of the input file requirements.
                $assembly->set_acc($acclist);
                
                push (@assemblies, $assembly);
            }
            
        }
        close $out_read_fh;
        
        if ($SEE) {
            print "PASA_OUTPUT ($forced_orient):\n####\n";
            system "cat $pasa_output";
            print "####\n";
        }
        unlink ($pasa_input, $pasa_output) unless $SEE;
        
    }
    
    $/ = $prev_input_sep; ## restore
    
    return (@assemblies);
}


sub force_flexorient {
    my $self = shift;
    my $orient = shift;
    my $alignments_aref = $self->{incoming_alignments};
    my $num_alignments = $#{$alignments_aref} + 1;
    ## Fix orientations of fli and multi-segment alignments:
    for (my $i = 0; $i < $num_alignments; $i++) {
        my $alignment = $alignments_aref->[$i];
        my $num_segments = $alignment->get_num_segments();
        my $spliced_orientation = $alignment->get_spliced_orientation();
        
        if ($spliced_orientation =~ /^[+-]$/) { # set specifically
            $alignment->{fixed_orient} = $spliced_orientation;  ## adding tag, using only in this module.
        } else {
            print "Setting $i to flex orient: $orient.\n" if $SEE;
            $alignment->{fixed_orient} = $orient;
        }
    }
    
}


sub unique_entries {
    my @x = @_;
    my %z;
    foreach my $y (@x) {
        $z{$y}=1;
    }
    return (keys %z);
}


=item get_assemblies()
    
=over 4
    
B<Description:> returns all the alignment assemblies resulting from the assembly procedure.

B<Parameters:> none.

B<Returns:> @assemblies

@assemblies is an array of CDNA::CDNA_alignment objects.

use the get_acc() method of the alignment object to retrieve all the accessions of the cDNAs that were merged into the assembly.

=back

=cut


sub get_assemblies {
    my $self = shift;
    return (@{$self->{assemblies}});
}

=item toAlignIllustration()

=over 4

B<Description:> illustrates the individual cDNAs to be assembled along with the final products.

B<Parameters:> $max_line_chars(optional)

$max_line_chars is an integer representing the maximum number of characters in a single line of output to the terminal.  The default is  100.

B<Returns:> $alignment_illustration_text

$alignment_illustration_text is a string containing a paragraph of text which illustrates the alignments and assemblies. An example is below:

 --->    <-->  <----->     <--->    <----------------	(+)gi|1199466

 --->    <-->  <----->     <--->    <------------   (+)gi|1209702
                        
---->    <-->  <----	(+)AV827070

---->    <-->  <---	(+)AV828861

---->    <-->  <---      (+)AV830936

 --->    <-->  <-	(+)H36350

ASSEMBLIES: (1)

---->    <-->  <----->     <--->    <---------------- (+) gi|1199466, gi|1209702, AV827070, AV828861, AV830936, H36350




=back

=cut

    ;    

sub toAlignIllustration () {
    my $self = shift;
    my $max_line_chars = shift;
    $max_line_chars = ($max_line_chars) ? $max_line_chars : 100; #if not specified, 100  chars / line is default.
    
    ## Get minimum coord for relative positioning.
    my @coords;
    my @alignments = @{$self->{incoming_alignments}};
    foreach my $alignment (@alignments) {
        my @c = $alignment->get_coords();
        push (@coords, @c);
    }
    @coords = sort {$a<=>$b} @coords;
    print "coords: @coords\n" if $::SEE;
    my $min_coord = shift @coords;
    my $max_coord = pop @coords;
    my $rel_max = $max_coord - $min_coord;
    my $alignment_text = "";
    ## print each alignment followed by assemblies:
    my $num_alignments = $#alignments + 1;
    $alignment_text .= "Individual Alignments: ($num_alignments)\n";
    my $i = 0;
    foreach my $alignment (@alignments) {
        $alignment_text .= (sprintf ("%3d ", $i)) . $alignment->toAlignIllustration($min_coord, $rel_max, $max_line_chars) . "\n";
        $i++;
    }
    
    my @assemblies = @{$self->{assemblies}};
    my $num_assemblies = $#assemblies + 1;
    $alignment_text .= "\n\nASSEMBLIES: ($num_assemblies)\n";
    foreach my $assembly (@assemblies) {
        $alignment_text .= "    " . $assembly->toAlignIllustration($min_coord, $rel_max, $max_line_chars) . "\n";
    }
    
    return ($alignment_text);
}


1;
