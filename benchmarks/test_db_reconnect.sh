#!/bin/bash
# test_db_reconnect.sh <training_dir>
#
# Verifies DB_connect::reconnect_to_server semantics against a real MariaDB:
#
#   1. when the connection is alive it is KEPT (same underlying $dbh, and the
#      prepare_cached statement cache survives) -- the fast path that makes the
#      per-cluster call in assemble_clusters.dbi cheap;
#   2. when the connection is genuinely dead it is REBUILT and usable again.
#
# Point 1 is the whole reason for the change, so it is asserted rather than
# assumed: the previous implementation rebuilt unconditionally, which silently
# discarded the statement cache on every one of ~16k per-cluster calls.
#
#   sbatch -p short -c 4 --mem 16gb -t 1:00:00 \
#     benchmarks/test_db_reconnect.sh <training_dir>
set -euo pipefail

TRAINING_ROOT="${1:?usage: test_db_reconnect.sh <training_dir>}"

: "${REPO_ROOT:=${SLURM_SUBMIT_DIR:-$PWD}}"
FUNANNOTATE_ENV="${FUNANNOTATE_ENV:-/opt/linux/rocky/8.x/x86_64/pkgs/funannotate/1.9/.pixi/envs/default}"
export PATH="$REPO_ROOT/bin:$FUNANNOTATE_ENV/bin:$PATH"
export PASAHOME="$REPO_ROOT"

module load singularity 2>/dev/null || true
# shellcheck source=/dev/null
source "$REPO_ROOT/benchmarks/mariadb_lib.sh"

SPECIES="$(basename "$TRAINING_ROOT")"
WORKDIR="${WORKDIR:-/scratch/$USER/${SLURM_JOB_ID:-$$}/reconn_${SPECIES}}"
DB_NAME=$(sed -n 's/^DATABASE=//p' "$TRAINING_ROOT/training/pasa/alignAssembly.txt" | head -1)

stage_datadir "$TRAINING_ROOT" "$WORKDIR"
trap 'mariadb_stop' EXIT
mariadb_start "$WORKDIR" "$DB_NAME"

PASACONF="$WORKDIR/mysql_db/conf/pasa-reconn.config.txt"
cat > "$PASACONF" <<CONF
MYSQLSERVER=${MYHOSTNAME}:${PORT}
MYSQL_RW_USER=${DB_USER}
MYSQL_RW_PASSWORD=${DB_PASS}
MYSQL_RO_USER=${DB_USER}
MYSQL_RO_PASSWORD=${DB_PASS}
PASA_ADMIN_EMAIL=test@localhost
PASA_ADMIN_DB=PASA2_admin
USE_PASA_DB_SETUP_HOOK=false
BASE_PASA_URL=http://localhost/cgi-bin/
HOOK_PERL_LIBS=__PASAHOME__/SAMPLE_HOOKS
HOOK_GENE_STRUCTURE_UPDATER=GFF3::GFF3_annot_updater::get_updater_obj
CONF
export PASACONF

echo "=== Testing reconnect_to_server semantics ==="
perl -I "$REPO_ROOT/PerlLib" -e '
    use strict; use warnings;
    use DBI; use Pasa_conf; use DB_connect;
    use Time::HiRes qw(time);

    my $db = shift;
    my $p = DB_connect::connect_to_db(
        Pasa_conf::getParam("MYSQLSERVER"), $db,
        Pasa_conf::getParam("MYSQL_RW_USER"), Pasa_conf::getParam("MYSQL_RW_PASSWORD"));

    my $fail = 0;

    ## 1. alive -> same handle kept
    my $before = "$p->{dbh}";
    $p = DB_connect::reconnect_to_server($p);
    my $after = "$p->{dbh}";
    if ($before eq $after) { print "PASS: live connection retained (same dbh)\n"; }
    else { print "FAIL: live connection was rebuilt ($before -> $after)\n"; $fail++; }

    ## 2. prepare_cached cache survives the call
    my $sth1 = $p->{dbh}->prepare_cached("select count(*) from clusters where annotdb_asmbl_id = ?");
    $p = DB_connect::reconnect_to_server($p);
    my $sth2 = $p->{dbh}->prepare_cached("select count(*) from clusters where annotdb_asmbl_id = ?");
    if ("$sth1" eq "$sth2") { print "PASS: prepare_cached statement survives reconnect call\n"; }
    else { print "FAIL: statement cache discarded ($sth1 -> $sth2)\n"; $fail++; }

    ## 3. cost of the hot path: ~16k calls in assemble_clusters.
    ##    Compare against the previous implementation, reproduced here, which
    ##    rebuilt the connection unconditionally on every call.
    my $N = 2000;

    my $t0 = time();
    $p = DB_connect::reconnect_to_server($p) for (1 .. $N);
    my $new_el = time() - $t0;

    my $legacy = sub {
        my ($dbproc) = @_;
        my $new_dbh = DB_connect::connect_to_db(
            $dbproc->{__server}, $dbproc->{__db}, $dbproc->{__username}, $dbproc->{__password});
        $dbproc->{dbh} = $new_dbh->{dbh};
        return $dbproc;
    };
    my $legacy_N = 200;   ## fewer: each one is a full connect + auth
    my $q = DB_connect::connect_to_db(
        Pasa_conf::getParam("MYSQLSERVER"), $db,
        Pasa_conf::getParam("MYSQL_RW_USER"), Pasa_conf::getParam("MYSQL_RW_PASSWORD"));
    my $t1 = time();
    $q = $legacy->($q) for (1 .. $legacy_N);
    my $old_el = time() - $t1;
    $q->disconnect;

    printf("INFO: new    %6.3f ms/call  => %6.1fs over ~16k clusters\n", 1000*$new_el/$N, 16000*$new_el/$N);
    printf("INFO: legacy %6.3f ms/call  => %6.1fs over ~16k clusters\n", 1000*$old_el/$legacy_N, 16000*$old_el/$legacy_N);
    printf("INFO: connection overhead reduced %.0fx (excludes the separate cost of\n"
         . "      re-preparing every cached statement, which legacy also forced)\n",
           ($old_el/$legacy_N) / ($new_el/$N));

    ## 4. genuinely dead connection -> rebuilt and usable
    $p->{dbh}->disconnect;
    $p = DB_connect::reconnect_to_server($p);
    my $ok = eval {
        my @r = DB_connect::do_sql_2D($p, "select count(*) from clusters");
        print "PASS: dead connection rebuilt, query returned $r[0][0] clusters\n"; 1;
    };
    unless ($ok) { print "FAIL: could not use rebuilt connection: $@\n"; $fail++; }

    $p->disconnect;
    print $fail ? "\nRESULT: $fail check(s) FAILED\n" : "\nRESULT: all checks passed\n";
    exit($fail ? 1 : 0);
' "$DB_NAME"
