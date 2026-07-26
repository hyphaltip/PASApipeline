#!/usr/bin/env perl

package main;
our $CLUSTERPATH;


package SingleLinkageClusterer;

## package not to be instantiated.  Just provides a namespace.

## Input: Array containing array-refs of pairs:
##               @_ = ( [1,2], [2,3], [6,7], [7,8], ...)
## Output: Array of all clusters as array-refs.
##              return ([1,2,3] , [6,7,8], ...)

use strict;
use warnings;
use File::Temp qw(tempfile);
use Pasa_tmpdir;

__run_test() unless caller;

## Detect available slclust binary
my $SLCLUST_RUST = _which("slclust_rust");
my $SLCLUST      = _which("slclust");

sub _which {
    my ($tool) = @_;
    for my $path (split /:/, $ENV{PATH} || '') {
        return "$path/$tool" if -x "$path/$tool";
    }
    return undef;
}


sub build_clusters {
    my @pairs = @_;
    
    ## Was: hardcoded /tmp with a "$$.<time>.<rand>" name. Two problems on a
    ## cluster -- /tmp is node-wide, often small, and sometimes tmpfs-backed
    ## (charging the job's memory), while schedulers provide faster per-job
    ## node-local scratch; and the name was constructed rather than claimed.
    ## This function is reached from threads via subcluster_builder.dbi and
    ## runs once per cluster, so temp churn here is heavy.
    my $tmpdir = Pasa_tmpdir::get_tmpdir();

    my ($pair_fh, $pairfile) = tempfile("slclust.XXXXXXXXXX", DIR => $tmpdir,
                                        SUFFIX => ".pairs", UNLINK => 0);

    #must do mapping because cluster program doesn't like word chars, just ints.
    my %map_id_to_feat;
    my %map_feat_to_id;
    my $id = 1;

    foreach my $pair (@pairs) {
        my ($a, $b) = @$pair;
        unless ($map_feat_to_id{$a}) {
            $map_feat_to_id{$a} = $id;
            $map_id_to_feat{$id} = $a;
            $id++;
        }
        unless ($map_feat_to_id{$b}) {
            $map_feat_to_id{$b} = $id;
            $map_id_to_feat{$id} = $b;
            $id++;
        }
        
        print $pair_fh "$map_feat_to_id{$a} $map_feat_to_id{$b}\n";
    }
    close $pair_fh;

    my ($cluster_fh, $clusterfile) = tempfile("slclust.XXXXXXXXXX", DIR => $tmpdir,
                                              SUFFIX => ".clusters", UNLINK => 0);
    close $cluster_fh;  ## slclust writes this itself via shell redirection
    unless (-w $clusterfile) { die "Can't write $clusterfile";}
    
    ## Prefer slclust_rust: 3-9x faster than C++ slclust at every scale tested
    ## (plain clustering and Jaccard filtering alike), deterministic output,
    ## and no stack-depth limit to worry about. Fall back to C++ slclust if
    ## slclust_rust isn't installed. Set SLCLUST_BACKEND=cpp to force C++.
    my $cmd;
    my $force_cpp = ($ENV{SLCLUST_BACKEND} || '') eq 'cpp';
    if ($SLCLUST_RUST && -x $SLCLUST_RUST && !$force_cpp) {
        $cmd = "$SLCLUST_RUST < $pairfile > $clusterfile";
    } elsif ($SLCLUST && -x $SLCLUST) {
        $cmd = "ulimit -s unlimited 2>/dev/null; $SLCLUST < $pairfile > $clusterfile";
    } else {
        ## both temp files already exist by this point; don't strand them
        unlink ($pairfile, $clusterfile);
        die "ERROR: Neither slclust nor slclust_rust found in PATH";
    }
    
    my $ret = system ($cmd);
    if ($ret) {
        unlink ($pairfile, $clusterfile);
        die "ERROR: Couldn't run cluster properly via path: $cmd";
    }

    my @clusters;
    open (my $clusters_fh, '<', $clusterfile)
        or die "Error, cannot read cluster output $clusterfile: $!";

    while (my $line = <$clusters_fh>) {
        my @elements;
        while ($line =~ /(\d+)\s?/g) {
            push (@elements, $map_id_to_feat{$1});
        }
        if (@elements) {
            push (@clusters, [@elements]);
        }
    }
    
    close $clusters_fh;
    
    ## clean up
    unlink ($pairfile, $clusterfile);
    
    return (@clusters);
}


############
## Testing
###########

sub __run_test {
    
    my @pairs = ( [1,2], [2,3], [4,5] );

    my @clusters = &SingleLinkageClusterer::build_clusters(@pairs);

    use Data::Dumper;
    
    print "Input: " . Dumper(\@pairs);
    print "Output: " . Dumper(\@clusters);

    exit(0);
}


1;
