#!/bin/bash
# submit_sweeps.sh -- launch the thread-scaling sweeps described in
# docs/PERFORMANCE_OPTIMIZATION.md § Thread scaling.
#
# Two genomes chosen to separate the two things that limit per-contig
# threading in assemble_clusters.dbi / classify_alt_splice_isoforms.dbi:
#
#   AF36    7 scaffolds  -- contig count is the hard ceiling. Alignments are
#                           spread 48826/46752/37762/34904/33416/31946/18050,
#                           so the largest contig is 19.4% of the work and no
#                           amount of threading can beat ~5.15x. Sweeping past
#                           T=7 shows the plateau (and any oversubscription
#                           cost) directly.
#   21_Fla  437 scaffolds -- enough contigs that thread count, not contig
#                           count, is the binding constraint. This is where a
#                           real scaling knee (and the ithread clone cost of
#                           437 spawns) should show up.
#
# E1376 (103 scaffolds) is deliberately skipped for now: its assemble_clusters
# step is 16629s at -T 2, so a T=1 point alone is ~9h.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TRAIN_ROOT="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/genome_annotation_training"
cd "$REPO_ROOT"
mkdir -p logs

STEPS="assemble_clusters,classify_alt_splice_isoforms,subcluster_builder"

submit() {
  local species="$1" threads="$2" cores="$3" mem="$4" walltime="$5"
  echo "=== $species  threads=$threads  (-c $cores --mem $mem -t $walltime)"
  # Thread/script lists go as POSITIONAL args: sbatch --export splits its
  # value on commas, so passing them there truncates the list to its first
  # element (see the note in sweep_threads.sbatch).
  sbatch -p epyc -c "$cores" --mem "$mem" --time "$walltime" \
         -J "sweep_${species}" \
         --export=ALL,REPO_ROOT="$REPO_ROOT" \
         benchmarks/sweep_threads.sbatch \
         "$TRAIN_ROOT/$species" "$threads" "$STEPS" \
         "$REPO_ROOT/bench_results/sweeps/$species"
}

# 7 scaffolds: sweep across and past the contig-count ceiling
submit Aspergillus_flavus_AF36   1,2,4,7,8,16      32 128gb 1-00:00:00

# 437 scaffolds: sweep far enough to find the knee
submit Aspergillus_flavus_21_Fla 1,2,4,8,16,32,64  64 200gb 1-12:00:00

echo
echo "Watch with:  squeue -u $USER -n sweep_Aspergillus_flavus_AF36,sweep_Aspergillus_flavus_21_Fla"
echo "Results in:  $REPO_ROOT/bench_results/sweeps/<species>/sweep_results.tsv"
