package Pasa_tmpdir;

=head1 NAME

Pasa_tmpdir - resolve a writable scratch directory for PASA temp files

=head1 DESCRIPTION

PASA scripts exchange data with external binaries (C<pasa>, C<slclust>) through
temp files. Several call sites hardcoded C</tmp>, which is the wrong choice on
an HPC cluster:

=over

=item *

C</tmp> is often small, shared between every job on the node, or tmpfs-backed
(so it consumes the job's memory allocation).

=item *

Schedulers and site modules usually provide per-job node-local storage that is
faster and automatically cleaned up. On UCR HPCC the C<workspace/scratch> module
sets both C<TMPDIR> and C<SCRATCH> to C</scratch/$USER/$SLURM_JOB_ID>.

=item *

Writing PASA's per-cluster temp churn to shared storage is markedly slower than
node-local disk.

=back

=head2 get_tmpdir()

Returns the first candidate that exists and is writable, in this order:

=over

=item 1. C<$PASA_TMPDIR> - explicit PASA-specific override, wins over all else.

=item 2. C<$TMPDIR> - the POSIX standard, and what schedulers/site modules set.

=item 3. C<$SCRATCH> - common HPC convention where TMPDIR is not set.

=item 4. C</tmp> - last resort.

=item 5. C<.> - current directory, if even /tmp is unusable (e.g. a read-only
container filesystem).

=back

Each candidate is checked for existence B<and> writability, so a stale or
unwritable C<TMPDIR> inherited from the environment degrades to the next option
instead of causing a confusing failure deep inside an assembler.

The result is cached per interpreter. Under ithreads each thread caches
independently, which is harmless -- the value is derived from the environment,
not from shared mutable state.

=cut

use strict;
use warnings;
use Carp;
use Exporter;

our @ISA       = qw(Exporter);
our @EXPORT_OK = qw(get_tmpdir);

my $CACHED_TMPDIR;

sub get_tmpdir {

    return $CACHED_TMPDIR if defined $CACHED_TMPDIR;

    my @candidates = ($ENV{PASA_TMPDIR}, $ENV{TMPDIR}, $ENV{SCRATCH}, '/tmp', '.');

    foreach my $dir (@candidates) {
        next unless defined $dir && length $dir;
        next unless -d $dir && -w $dir;
        $CACHED_TMPDIR = $dir;
        return $CACHED_TMPDIR;
    }

    confess "Error, no writable temp directory found. Tried: "
          . join(", ", map { defined $_ && length $_ ? $_ : "(unset)" } @candidates)
          . ". Set PASA_TMPDIR to a writable path.\n";
}

1; #EOM
