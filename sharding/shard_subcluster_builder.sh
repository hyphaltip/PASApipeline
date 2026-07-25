#!/bin/bash
# shard_subcluster_builder.sh <alignAssembly.config> <genome.fasta> <PASA_LOG_DIR> [--jobs N]
#
# Runs scripts/subcluster_builder.dbi as one process per contig/scaffold via
# GNU parallel, instead of the script's own serial per-contig loop, then merges
# the results back into exactly the output file + checkpoint file that
# Launch_PASA_pipeline.pl expects for this step -- so it must be run BEFORE
# Launch_PASA_pipeline.pl / funannotate train reaches the subcluster_builder.dbi
# stage. Once this driver succeeds, that stage's checkpoint already exists, and
# Pipeliner.pm (see PerlLib/Pipeliner.pm) skips straight past it to
# subcluster_loader.dbi.
#
# Why process-level sharding instead of in-process threading: subcluster_builder.dbi
# forks pasa/slclust subprocesses per cluster (via SingleLinkageClusterer.pm).
# In-process threading of this script was tried and benchmarked as a regression
# at every thread count (issue #7, closed won't-fix) -- forking from a
# multi-threaded Perl process is more expensive than from a single-threaded
# one. Sharding across separate single-threaded processes sidesteps that cost
# entirely (issue #6).
#
# Why GNU parallel instead of a Slurm array: this is single-node fan-out (the
# production run that motivated this had a single-node allocation with 256
# cores and used only 2 -- no need for multi-node scale-out), and `parallel
# --joblog` gives an exact per-contig exit code + signal, which is simpler and
# more reliable than reconciling `sbatch --wait`'s array-level exit status
# against `sacct` per task index.
#
# Safe to shard: subcluster_loader.dbi (which consumes this step's output)
# tracks a single "current cluster" scalar that resets on every
# "Processing cluster:" header, and every contig's cluster_ids are disjoint --
# so per-contig output chunks may be concatenated in ANY ORDER across contigs,
# as long as each contig's own chunk stays internally intact.
#
# Usage:
#   sharding/shard_subcluster_builder.sh <alignAssembly.config> <genome.fasta> <PASA_LOG_DIR> [--jobs N]
#
# <alignAssembly.config> is the same config file passed to Launch_PASA_pipeline.pl
# via -c (must contain a DATABASE= line). <PASA_LOG_DIR> is the same directory
# Launch_PASA_pipeline.pl uses (normally "pasa_run.log.dir" under the PASA run
# directory). Run this script from that same PASA run directory.
#
# Requires GNU parallel on PATH.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
UTILDIR="$REPO_ROOT/scripts"

usage() {
    echo "Usage: $0 <alignAssembly.config> <genome.fasta> <PASA_LOG_DIR> [--jobs N]" >&2
    exit 1
}

if [ $# -lt 3 ]; then
    usage
fi

CONFIG_FILE="$1"
GENOME_DB="$2"
PASA_LOG_DIR="$3"
shift 3

JOBS=""
while [ $# -gt 0 ]; do
    case "$1" in
        --jobs)
            JOBS="${2:-}"
            [ -n "$JOBS" ] || usage
            shift 2
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage
            ;;
    esac
done

if ! command -v parallel >/dev/null 2>&1; then
    echo "Error: GNU parallel not found on PATH. Load it (e.g. 'module load parallel') before running this driver." >&2
    exit 1
fi

if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: config file not found: $CONFIG_FILE" >&2
    exit 1
fi
if [ ! -f "$GENOME_DB" ]; then
    echo "Error: genome fasta not found: $GENOME_DB" >&2
    exit 1
fi

DATABASE="$(grep -E '^DATABASE=' "$CONFIG_FILE" | head -1 | cut -d= -f2- | tr -d '[:space:]')"
if [ -z "$DATABASE" ]; then
    echo "Error: could not find a DATABASE= line in $CONFIG_FILE" >&2
    exit 1
fi

# Mirrors DB_connect::configure_db_driver + Launch_PASA_pipeline.pl:257,315 --
# same rule the main pipeline uses to pick the checkpoint dir name.
case "$DATABASE" in
    */*) DBI_DRIVER="SQLite" ;;
    *)   DBI_DRIVER="mysql" ;;
esac

CHKPTS_DIR="__pasa_$(basename "$DATABASE")_${DBI_DRIVER}_chkpts"
CHKPT_FILE="$CHKPTS_DIR/alignment_assembly_subclustering.ok"
OUT_FILE="$PASA_LOG_DIR/alignment_assembly_subclustering.out"

if [ -e "$CHKPT_FILE" ]; then
    echo "Checkpoint $CHKPT_FILE already exists -- this step is already done; the main pipeline will skip it. Nothing to do." >&2
    exit 0
fi

mkdir -p "$CHKPTS_DIR" "$PASA_LOG_DIR"

SHARD_DIR="$(mktemp -d "$(pwd)/shard_subcluster_builder.XXXXXX")"
echo "=== Shard workspace: $SHARD_DIR ==="

fail() {
    echo "FAILURE: $1" >&2
    echo "Shard workspace left at $SHARD_DIR for inspection. Checkpoint NOT written -- the main pipeline will still (re-)run this step normally if invoked, or fix the issue and rerun this driver." >&2
    exit 1
}

echo "=== Enumerating contigs for database '$DATABASE' ==="
"$UTILDIR/list_asmbl_ids.dbi" -M "$DATABASE" > "$SHARD_DIR/contigs.list" \
    || fail "list_asmbl_ids.dbi failed to enumerate contigs"

NUM_CONTIGS=$(wc -l < "$SHARD_DIR/contigs.list" | tr -d '[:space:]')
if [ "$NUM_CONTIGS" -eq 0 ]; then
    fail "no contigs found for database $DATABASE"
fi
echo "Found $NUM_CONTIGS contigs."

if [ -z "$JOBS" ]; then
    JOBS=$(nproc 2>/dev/null || echo 2)
    if [ "$JOBS" -gt "$NUM_CONTIGS" ]; then
        JOBS="$NUM_CONTIGS"
    fi
fi
echo "=== Dispatching $NUM_CONTIGS contigs across $JOBS parallel jobs ==="

START_TIME=$(date +%s)
set +e
parallel --jobs "$JOBS" --joblog "$SHARD_DIR/joblog.tsv" --halt now,fail=1 \
    "'$UTILDIR/subcluster_builder.dbi' -G '$GENOME_DB' -M '$DATABASE' -R {} > '$SHARD_DIR'/{}.out 2> '$SHARD_DIR'/{}.err" \
    :::: "$SHARD_DIR/contigs.list"
PARALLEL_STATUS=$?
set -e
END_TIME=$(date +%s)
echo "=== parallel exited with status $PARALLEL_STATUS after $((END_TIME - START_TIME))s ==="

# --- Verification: nothing below this line may write $OUT_FILE/$CHKPT_FILE
# --- until every check below has passed.

if [ "$PARALLEL_STATUS" -ne 0 ]; then
    fail "GNU parallel reported at least one failed job (exit $PARALLEL_STATUS); see $SHARD_DIR/joblog.tsv"
fi

# joblog columns (tab-separated, header row): Seq Host Starttime JobRuntime
# Send Receive Exitval Signal Command
if ! tail -n +2 "$SHARD_DIR/joblog.tsv" | awk -F'\t' '
    { exitval=$7; signal=$8; if (exitval != 0 || signal != 0) { print > "/dev/stderr"; bad=1 } }
    END { if (bad) exit 1 }
'; then
    fail "one or more shard jobs had nonzero Exitval/Signal in $SHARD_DIR/joblog.tsv"
fi

JOBLOG_COUNT=$(tail -n +2 "$SHARD_DIR/joblog.tsv" | wc -l | tr -d '[:space:]')
if [ "$JOBLOG_COUNT" -ne "$NUM_CONTIGS" ]; then
    fail "$SHARD_DIR/joblog.tsv has $JOBLOG_COUNT entries, expected $NUM_CONTIGS -- some contig was never dispatched"
fi

echo "=== Checking every contig produced shard output ==="
while IFS= read -r CONTIG; do
    [ -f "$SHARD_DIR/$CONTIG.out" ] || fail "missing shard output file for contig $CONTIG"
done < "$SHARD_DIR/contigs.list"

echo "=== Structural completeness check (per-contig cluster counts) ==="
while IFS= read -r CONTIG; do
    EMITTED=$(grep -c '^// Processing cluster:' "$SHARD_DIR/$CONTIG.out" || true)
    EXPECTED=$("$UTILDIR/count_clusters_for_asmbl_id.dbi" -M "$DATABASE" -R "$CONTIG")
    if [ "$EMITTED" -ne "$EXPECTED" ]; then
        fail "contig $CONTIG emitted $EMITTED 'Processing cluster:' blocks but $EXPECTED clusters exist in the DB -- likely truncated/partial shard output"
    fi
done < "$SHARD_DIR/contigs.list"

echo "=== All checks passed. Merging $NUM_CONTIGS shard outputs into $OUT_FILE ==="
: > "$OUT_FILE"
while IFS= read -r CONTIG; do
    cat "$SHARD_DIR/$CONTIG.out" >> "$OUT_FILE"
done < "$SHARD_DIR/contigs.list"

touch "$CHKPT_FILE"

echo "=== Done: $NUM_CONTIGS contigs merged, checkpoint written: $CHKPT_FILE ==="
echo "=== Shard workspace retained at $SHARD_DIR (safe to remove once the main pipeline has completed subcluster_loader.dbi) ==="

exit 0
