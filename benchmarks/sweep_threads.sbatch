#!/bin/bash -l
#SBATCH -p epyc
#SBATCH -c 64
#SBATCH --mem 200gb
#SBATCH --time 3-00:00:00
#SBATCH -J pasa_sweep
#SBATCH -o logs/sweep_%x_%j.log
#
# Thread-scaling sweep launcher.
#
#   sbatch benchmarks/sweep_threads.sbatch \
#          <training_dir> <threads_csv> <scripts_csv> [results_dir]
#
# NB: these are positional arguments, NOT --export vars. `sbatch --export`
# splits its value on commas, so --export=ALL,THREADS=1,2,4 silently sets
# THREADS=1 and then tries to set variables named "2" and "4" -- which
# quietly truncated an earlier sweep to a single thread count.
#
# WORKDIR defaults to this node's local /scratch (much faster for MariaDB
# I/O than shared storage) but is node-local and wiped at job end, so results
# are copied to results_dir before exiting.
set -euo pipefail

module load singularity 2>/dev/null || true

# NB: sbatch copies this file to /var/spool/slurmd, so BASH_SOURCE does NOT
# point into the repo. Fall back to the submit directory, which is the repo
# root when submitted as `sbatch benchmarks/sweep_threads.sbatch`.
: "${REPO_ROOT:=${SLURM_SUBMIT_DIR:-$PWD}}"
[ -x "$REPO_ROOT/benchmarks/sweep_threads.sh" ] || {
    echo "ERROR: REPO_ROOT=$REPO_ROOT does not contain benchmarks/sweep_threads.sh." >&2
    echo "Submit from the repo root, or pass REPO_ROOT=/path/to/PASApipeline." >&2
    exit 1
}

TRAINING_DIR="${1:?usage: sweep_threads.sbatch <training_dir> <threads_csv> <scripts_csv> [results_dir]}"
THREADS="${2:-1,2,4,8,16}"
SCRIPTS="${3:-assemble_clusters,classify_alt_splice_isoforms,subcluster_builder}"

SPECIES="$(basename "$TRAINING_DIR")"
: "${RESULTS_DIR:=${4:-$REPO_ROOT/bench_results/sweeps/${SPECIES}_${SLURM_JOB_ID:-$$}}}"
: "${WORKDIR:=/scratch/$USER/${SLURM_JOB_ID:-$$}/sweep_${SPECIES}}"

echo "host=$(hostname -s) cpus=${SLURM_CPUS_ON_NODE:-?} job=${SLURM_JOB_ID:-none}"
echo "TRAINING_DIR=$TRAINING_DIR"
echo "WORKDIR=$WORKDIR"
echo "RESULTS_DIR=$RESULTS_DIR"
echo "THREADS=$THREADS"
echo "SCRIPTS=$SCRIPTS"

mkdir -p "$RESULTS_DIR"

# WORKDIR is node-local and disappears with the job; make sure whatever the
# sweep produced gets copied out even if it died partway through.
save_results() {
    echo "=== Copying results to $RESULTS_DIR ==="
    cp "$WORKDIR/sweep_results.tsv" "$RESULTS_DIR/" 2>/dev/null || true
    cp -r "$WORKDIR/logs" "$RESULTS_DIR/" 2>/dev/null || true
}
trap save_results EXIT

"$REPO_ROOT/benchmarks/sweep_threads.sh" \
    "$TRAINING_DIR" "$WORKDIR" "$THREADS" "$SCRIPTS"
