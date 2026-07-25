#!/bin/bash
# make_sqlite_copy.sh <training_dir> <out_sqlite_path>
#
# One-off: bring a completed run's MariaDB datadir up on node-local scratch,
# copy every table into a fresh SQLite database, and drop the .sqlite
# somewhere persistent. That file is then the input to the backend/storage
# sweep, so the MariaDB-vs-SQLite comparison runs on identical data rather
# than on two separately-produced runs.
#
#   sbatch -p short -c 8 --mem 32gb -t 2:00:00 \
#     benchmarks/make_sqlite_copy.sh <training_dir> <out.sqlite>
set -euo pipefail

TRAINING_ROOT="${1:?usage: make_sqlite_copy.sh <training_dir> <out_sqlite_path>}"
OUT_SQLITE="${2:?need output .sqlite path}"

: "${REPO_ROOT:=${SLURM_SUBMIT_DIR:-$PWD}}"
[ -f "$REPO_ROOT/benchmarks/mariadb_lib.sh" ] || {
    echo "ERROR: REPO_ROOT=$REPO_ROOT has no benchmarks/mariadb_lib.sh" >&2; exit 1; }

FUNANNOTATE_ENV="${FUNANNOTATE_ENV:-/opt/linux/rocky/8.x/x86_64/pkgs/funannotate/1.9/.pixi/envs/default}"
export PATH="$REPO_ROOT/bin:$FUNANNOTATE_ENV/bin:$PATH"
export PASAHOME="$REPO_ROOT"

module load singularity 2>/dev/null || true
# shellcheck source=/dev/null
source "$REPO_ROOT/benchmarks/mariadb_lib.sh"

SPECIES="$(basename "$TRAINING_ROOT")"
WORKDIR="${WORKDIR:-/scratch/$USER/${SLURM_JOB_ID:-$$}/mk_sqlite_${SPECIES}}"
DB_NAME=$(sed -n 's/^DATABASE=//p' "$TRAINING_ROOT/training/pasa/alignAssembly.txt" | head -1)
[ -n "$DB_NAME" ] || { echo "ERROR: no DATABASE= in alignAssembly.txt" >&2; exit 1; }

stage_datadir "$TRAINING_ROOT" "$WORKDIR"
trap 'mariadb_stop' EXIT
mariadb_start "$WORKDIR" "$DB_NAME"

TMP_SQLITE="$WORKDIR/$(basename "$OUT_SQLITE")"
rm -f "$TMP_SQLITE"

echo "=== Converting $DB_NAME -> $TMP_SQLITE ==="
"$REPO_ROOT/benchmarks/mysql_to_sqlite.pl" \
    --server "${MYHOSTNAME}:${PORT}" \
    --db "$DB_NAME" \
    --user "$DB_USER" \
    --password "$DB_PASS" \
    --sqlite "$TMP_SQLITE"

# Publish as .unverified FIRST. The conversion is the expensive part and the
# workdir is node-local, so a hiccup in the checks below must not discard it.
mkdir -p "$(dirname "$OUT_SQLITE")"
cp "$TMP_SQLITE" "$OUT_SQLITE.unverified"

echo "=== Verifying row counts match ==="
# NB: no bare command substitution here -- under `set -e` a non-zero mysql or
# sqlite3 exit in an assignment kills the script with no output at all, which
# is exactly how an earlier run failed silently after a good conversion.
FAIL=0
for t in clusters align_link alignment cdna_info subcluster_link asmbl_link; do
  # Query over TCP, not the socket: the grant we added is scoped to
  # '<user>'@'%', but a socket connection authenticates as '<user>'@'localhost'
  # and is rejected. TCP is also what the conversion itself used.
  if ! m=$(mysql -h "$MYHOSTNAME" -P "$PORT" -u "$DB_USER" -p"$DB_PASS" -N -B \
             -e "select count(*) from \`$DB_NAME\`.\`$t\`" 2>&1); then
    printf "  %-20s mysql query failed: %s\n" "$t" "$m"; FAIL=1; continue
  fi
  if ! s=$(sqlite3 "$TMP_SQLITE" "select count(*) from $t" 2>&1); then
    printf "  %-20s sqlite query failed: %s\n" "$t" "$s"; FAIL=1; continue
  fi
  m="${m//[[:space:]]/}"; s="${s//[[:space:]]/}"
  if [ "$m" = "$s" ]; then
    printf "  %-20s %10s  OK\n" "$t" "$m"
  else
    printf "  %-20s mysql=%s sqlite=%s  MISMATCH\n" "$t" "$m" "$s"; FAIL=1
  fi
done
[ "$FAIL" -eq 0 ] || {
    echo "ERROR: verification failed; left unverified copy at $OUT_SQLITE.unverified" >&2
    exit 1
}

mv "$OUT_SQLITE.unverified" "$OUT_SQLITE"
ls -lh "$OUT_SQLITE"
echo "=== Done: $OUT_SQLITE ==="
