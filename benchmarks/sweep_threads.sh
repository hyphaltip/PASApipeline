#!/bin/bash
# sweep_threads.sh <training_dir> <workdir> <threads_csv> [scripts_csv]
#
# Thread-scaling sweep for the PASA steps that dominate real runs.
#
# Rationale (see docs/PERFORMANCE_OPTIMIZATION.md): profiling three real
# Aspergillus flavus funannotate-train runs showed assemble_clusters.dbi at
# 48-66% of PASA wall time and classify_alt_splice_isoforms.dbi at 12-18%,
# and every one of those runs was invoked with -T 2. So before optimizing
# code we need the actual scaling curve for these steps.
#
# This does NOT re-run the whole pipeline. It reuses a completed run's
# populated MariaDB datadir, restarts it in a Singularity instance against a
# scratch COPY (the shared original is never touched), and re-runs just the
# step(s) under test at each thread count.
#
# assemble_clusters.dbi is read-only against the DB -- it only writes
# assemblies/*.assemblies -- so it is re-runnable by deleting that directory.
# classify_alt_splice_isoforms.dbi purges its own tables on entry
# (init_alt_splice_tables), so it is re-runnable in place.
#
# Records per (script, threads): wall seconds, peak RSS, and a checksum of
# the sorted output so a thread count that changes results is caught.
#
# Usage:
#   ./benchmarks/sweep_threads.sh \
#       /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/genome_annotation_training/Aspergillus_flavus_AF36 \
#       /scratch/$USER/sweep_AF36 \
#       1,2,4,7,8,16
set -euo pipefail

TRAINING_ROOT="${1:?usage: sweep_threads.sh <training_dir> <workdir> <threads_csv> [scripts_csv]}"
WORKDIR="${2:?need workdir}"
THREADS_CSV="${3:-1,2,4,8,16}"
SCRIPTS_CSV="${4:-assemble_clusters,classify_alt_splice_isoforms}"

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="$TRAINING_ROOT/training"
MARIADB_SIF="/bigdata/stajichlab/shared/lib/mariadb/mariadb.sif"

[ -d "$SRC/mysql_db/db" ] || { echo "ERROR: no populated datadir at $SRC/mysql_db/db" >&2; exit 1; }
[ -s "$SRC/genome.fasta" ] || { echo "ERROR: no genome at $SRC/genome.fasta" >&2; exit 1; }

DB_NAME=$(sed -n 's/^DATABASE=//p' "$SRC/pasa/alignAssembly.txt" | head -1)
[ -n "$DB_NAME" ] || { echo "ERROR: could not parse DATABASE= from $SRC/pasa/alignAssembly.txt" >&2; exit 1; }

# The PASA scripts are `#!/usr/bin/env perl`, and the system perl has neither
# DBI nor DBD::mysql. Use the same funannotate env the production runs used,
# so the sweep measures the same interpreter/module stack as the real runs.
FUNANNOTATE_ENV="${FUNANNOTATE_ENV:-/opt/linux/rocky/8.x/x86_64/pkgs/funannotate/1.9/.pixi/envs/default}"
[ -x "$FUNANNOTATE_ENV/bin/perl" ] || { echo "ERROR: no perl at $FUNANNOTATE_ENV/bin/perl" >&2; exit 1; }

export PASAHOME="$REPO_ROOT"
export PATH="$REPO_ROOT/bin:$FUNANNOTATE_ENV/bin:$PATH"

echo "=== Staging scratch copy of $DB_NAME into $WORKDIR ==="
rm -rf "$WORKDIR"
mkdir -p "$WORKDIR/mysql_db"
# -L so a symlinked genome is materialized; datadir copy is the bulk (~800M)
cp -rL "$SRC/mysql_db/db"   "$WORKDIR/mysql_db/db"
cp -rL "$SRC/mysql_db/conf" "$WORKDIR/mysql_db/conf"
cp -L  "$SRC/genome.fasta"  "$WORKDIR/genome.fasta"

# The genome indexes are NOT optional for a faithful benchmark. Without
# genome.fasta.fai, Fasta_retriever::_init takes its fallback branch: it scans
# the FASTA itself and leaves an open filehandle in $self->{fh}, which is then
# cloned into every worker thread. With the .fai present it builds only the
# position index and each thread opens its own handle on first get_seq().
#
# Omitting the .fai therefore puts threaded runs on a code path production
# never uses, and it produced corrupted sequence reads that surfaced as
# "spliced orient in db (+) differs from calculated spliced orient (-)" from
# Ath1_cdnas::batch_create_alignment_objs_by_id -- a harness artifact that
# looked exactly like a genuine threading bug.
for idx in .fai .cidx; do
    if [ -s "$SRC/genome.fasta$idx" ]; then
        cp -L "$SRC/genome.fasta$idx" "$WORKDIR/genome.fasta$idx"
    else
        echo "WARNING: $SRC/genome.fasta$idx missing; threaded runs may not match production" >&2
    fi
done
[ -s "$WORKDIR/genome.fasta.fai" ] || {
    echo "ERROR: genome.fasta.fai is required for a faithful benchmark; aborting." >&2
    exit 1
}

chmod -R u+w "$WORKDIR"

cd "$WORKDIR"

module load singularity 2>/dev/null || true

INSTANCE_NAME="pasa_sweep_$$"
PORT=$(shuf -i 14000-14999 -n1)
MYHOSTNAME=$(hostname -s)

# Reuse the credentials the original run already granted inside this datadir.
# We cannot mint a new user: the production datadir's root account has a
# password we do not hold (unlike run_bench_mysql.sh, which mysql_install_db's
# a fresh datadir with passwordless root).
SRC_PASACONF=$(ls "$WORKDIR"/mysql_db/conf/pasa-local*.config.txt 2>/dev/null | head -1)
[ -n "$SRC_PASACONF" ] || { echo "ERROR: no preserved pasa-local*.config.txt in $SRC/mysql_db/conf" >&2; exit 1; }
DB_USER=$(sed -n 's/^MYSQL_RW_USER=//p'     "$SRC_PASACONF" | head -1)
DB_PASS=$(sed -n 's/^MYSQL_RW_PASSWORD=//p' "$SRC_PASACONF" | head -1)
[ -n "$DB_USER" ] && [ -n "$DB_PASS" ] || { echo "ERROR: could not parse credentials from $SRC_PASACONF" >&2; exit 1; }
echo "Reusing preserved credentials for user '$DB_USER'"

sed -i "s/^port = .*/port = $PORT/"          "$WORKDIR/mysql_db/conf/my.cnf"
sed -i "s/^user = .*/user = $USER/"          "$WORKDIR/mysql_db/conf/my.cnf"
sed -i "s/^bind-address.*/bind-address = 0.0.0.0/" "$WORKDIR/mysql_db/conf/my.cnf" || true

stop_mysqldb() { singularity instance stop "$INSTANCE_NAME" 2>/dev/null || true; }
trap "stop_mysqldb; exit 130" SIGHUP SIGINT SIGTERM
trap "stop_mysqldb" EXIT

export SINGULARITY_BINDPATH="$WORKDIR/mysql_db/db:/var/lib/mysql,$WORKDIR/mysql_db/conf/my.cnf:/etc/mysql/my.cnf"

wait_for_mysql() {
  for i in $(seq 1 60); do
    if singularity exec instance://"$INSTANCE_NAME" \
         mysqladmin --socket=/var/run/mysqld/mysqld.sock ping >/dev/null 2>&1; then
      echo "MariaDB up after ${i}s"; return 0
    fi
    sleep 1
  done
  echo "ERROR: MariaDB did not start" >&2; return 1
}

# The datadir's grants are scoped to the node the original run happened to
# land on ('jstajich'@'r22'), and its root account carries a password we do
# not hold -- so neither socket-as-user nor root can widen them. Bring the
# server up once with --skip-grant-tables (no auth at all), re-enable the
# grant system in-session with FLUSH PRIVILEGES, add a '%'-scoped grant, then
# restart normally so the sweep runs against a normally-authenticating server.
# NB: passing --skip-grant-tables as a mysqld_safe argument through
# `singularity instance start` does not reach mysqld (the server still
# enforced privileges and rejected FLUSH for lack of RELOAD). Setting it in
# my.cnf under [mysqld] is honored unambiguously.
echo "=== Pass 1: starting MariaDB with skip-grant-tables to add a host-independent grant ==="
sed -i "/^\[mysqld\]/a skip-grant-tables" "$WORKDIR/mysql_db/conf/my.cnf"
singularity instance start --writable-tmpfs "$MARIADB_SIF" "$INSTANCE_NAME" /usr/bin/mysqld_safe
wait_for_mysql

singularity exec instance://"$INSTANCE_NAME" \
  mysql --socket=/var/run/mysqld/mysqld.sock -u root <<SQL
FLUSH PRIVILEGES;
CREATE USER IF NOT EXISTS '$DB_USER'@'%' IDENTIFIED BY '$DB_PASS';
ALTER USER '$DB_USER'@'%' IDENTIFIED BY '$DB_PASS';
GRANT ALL PRIVILEGES ON \`$DB_NAME\`.* TO '$DB_USER'@'%';
FLUSH PRIVILEGES;
SQL

stop_mysqldb
sleep 5

echo "=== Pass 2: restarting MariaDB normally on port $PORT ==="
sed -i "/^skip-grant-tables$/d" "$WORKDIR/mysql_db/conf/my.cnf"
singularity instance start --writable-tmpfs "$MARIADB_SIF" "$INSTANCE_NAME" /usr/bin/mysqld_safe
wait_for_mysql

PASACONF="$WORKDIR/mysql_db/conf/pasa-sweep.config.txt"
cat > "$PASACONF" <<CONF
MYSQLSERVER=${MYHOSTNAME}:${PORT}
MYSQL_RW_USER=${DB_USER}
MYSQL_RW_PASSWORD=${DB_PASS}
MYSQL_RO_USER=${DB_USER}
MYSQL_RO_PASSWORD=${DB_PASS}
PASA_ADMIN_EMAIL=sweep@localhost
PASA_ADMIN_DB=PASA2_admin
USE_PASA_DB_SETUP_HOOK=false
BASE_PASA_URL=http://localhost/cgi-bin/
HOOK_PERL_LIBS=__PASAHOME__/SAMPLE_HOOKS
HOOK_GENE_STRUCTURE_UPDATER=GFF3::GFF3_annot_updater::get_updater_obj
CONF
export PASACONF

# Fail fast: prove PASA can actually reach the DB over TCP the same way the
# pipeline scripts will, rather than discovering it hours into the sweep.
echo "=== Verifying PASA DB connectivity (${MYHOSTNAME}:${PORT}) ==="
perl -I "$REPO_ROOT/PerlLib" -e '
    use DBI; use Pasa_conf; use DB_connect;
    my $db = shift;
    my $p = DB_connect::connect_to_db(
        Pasa_conf::getParam("MYSQLSERVER"), $db,
        Pasa_conf::getParam("MYSQL_RW_USER"), Pasa_conf::getParam("MYSQL_RW_PASSWORD"));
    my @r = DB_connect::do_sql_2D($p, "select count(*) from clusters");
    print "OK: $r[0][0] clusters visible\n";
    $p->disconnect;
' "$DB_NAME" || { echo "ERROR: PASA cannot connect to the staged DB; aborting sweep." >&2; exit 1; }

RESULTS="$WORKDIR/sweep_results.tsv"
printf "script\tthreads\twall_s\tmax_rss_kb\texit\toutput_checksum\n" > "$RESULTS"
mkdir -p "$WORKDIR/logs"

run_one() {
  local script="$1" threads="$2"
  local tag="${script}.T${threads}"
  local tfile="$WORKDIR/logs/${tag}.time"
  local log="$WORKDIR/logs/${tag}.log"
  local rc=0 checksum="-"

  echo ""
  echo "=== $script  -T $threads ==="

  case "$script" in
    assemble_clusters)
      rm -rf "$WORKDIR/assemblies"
      /usr/bin/time -v -o "$tfile" \
        "$REPO_ROOT/scripts/assemble_clusters.dbi" \
          -G "$WORKDIR/genome.fasta" -M "$DB_NAME" -T "$threads" \
          > "$log" 2>&1 || rc=$?
      # sort so per-contig completion order does not perturb the checksum
      # A failed run may leave no .described files at all; with `set -o
      # pipefail` a failing `cat` here would abort the entire sweep rather
      # than just recording a bad point for this thread count.
      if compgen -G "$WORKDIR/assemblies/*.assemblies.described" >/dev/null 2>&1; then
        checksum=$(cat "$WORKDIR"/assemblies/*.assemblies.described 2>/dev/null \
                   | sort | md5sum | cut -d' ' -f1) || checksum="checksum-failed"
      else
        checksum="no-output"
      fi
      ;;
    classify_alt_splice_isoforms)
      /usr/bin/time -v -o "$tfile" \
        "$REPO_ROOT/scripts/classify_alt_splice_isoforms.dbi" \
          -G "$WORKDIR/genome.fasta" -M "$DB_NAME" -T "$threads" \
          > "$log" 2>&1 || rc=$?
      # Query over TCP, not the socket. The grant we add is '<user>'@'%', but a
      # socket connection authenticates as '<user>'@'localhost' and is
      # rejected -- and under `set -o pipefail` that rejection propagates out
      # of this assignment and kills the whole sweep. That is exactly how a
      # previous run died immediately after a SUCCESSFUL classify -T 1, losing
      # every remaining thread count.
      checksum=$(mysql -h "$MYHOSTNAME" -P "$PORT" -u "$DB_USER" -p"$DB_PASS" -N -B \
                   -e "select type,count(*) from ${DB_NAME}.splice_variation group by type order by type" \
                   2>/dev/null | md5sum | cut -d' ' -f1) || checksum="query-failed"
      ;;
    subcluster_builder)
      # no -T; serial reference point only, run once
      /usr/bin/time -v -o "$tfile" \
        "$REPO_ROOT/scripts/subcluster_builder.dbi" \
          -G "$WORKDIR/genome.fasta" -M "$DB_NAME" -m 50 \
          > "$log" 2>&1 || rc=$?
      checksum=$(md5sum < "$log" | cut -d' ' -f1)
      ;;
    *)
      echo "ERROR: unknown script '$script'" >&2; return 1 ;;
  esac

  local wall rss
  wall=$(awk -F': ' '/Elapsed \(wall clock\)/{print $2}' "$tfile")
  rss=$(awk -F': ' '/Maximum resident set size/{print $2}' "$tfile")
  # h:mm:ss or m:ss -> seconds
  wall=$(awk -v t="$wall" 'BEGIN{n=split(t,a,":"); s=0; for(i=1;i<=n;i++) s=s*60+a[i]; printf "%.1f", s}')

  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$script" "$threads" "$wall" "$rss" "$rc" "$checksum" >> "$RESULTS"
  echo "--> wall=${wall}s rss=${rss}kb exit=$rc checksum=$checksum"
}

IFS=',' read -ra SCRIPTS <<< "$SCRIPTS_CSV"
IFS=',' read -ra THREADS <<< "$THREADS_CSV"

for script in "${SCRIPTS[@]}"; do
  if [ "$script" = "subcluster_builder" ]; then
    run_one "$script" 1
    continue
  fi
  for t in "${THREADS[@]}"; do
    run_one "$script" "$t"
  done
done

echo ""
echo "=== Sweep complete: $RESULTS ==="
column -t "$RESULTS"
python3 "$REPO_ROOT/benchmarks/sweep_summarize.py" "$RESULTS" || true
