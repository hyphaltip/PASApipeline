#!/bin/bash
# run_bench.sh <run_dir> <cpu> [extra __run_sample_pipeline.pl args...]
#
# Runs the PASA sample_data align-assembly pipeline (the stage funannotate
# train actually exercises: alignment import, validate_alignments_in_db.dbi,
# clustering/assign_clusters_by_*.dbi, subcluster_builder.dbi,
# alignment_assembly_to_gene_models.dbi, TRANSDECODER via
# pasa_asmbls_to_training_set-equivalent), timing each sub-command via
# Pipeliner's verbose "* [datestamp] Running CMD: ..." log lines.
#
# Usage:
#   ./run_bench.sh run_sqlite_baseline 4
#   ./run_bench.sh run_sqlite_fixed 4
set -euo pipefail

# Reuse the funannotate 1.9 pixi env for blat/gmap/minimap2/samtools/TransDecoder
# and a perl with DBD::SQLite already installed (this checkout's own
# Launch_PASA_pipeline.pl resolves its sibling scripts/ relative to itself via
# FindBin, so only external tool binaries need to come from the env PATH).
FUNANNOTATE_ENV="/opt/linux/rocky/8.x/x86_64/pkgs/funannotate/1.9/.pixi/envs/default"
export PATH="$FUNANNOTATE_ENV/bin:$PATH"

RUN_DIR="$1"
CPU="${2:-4}"
shift 2 || true

# RUN_DIR must live directly at the repo root (sibling of PerlLib/,
# Launch_PASA_pipeline.pl, scripts/) -- __run_sample_pipeline.pl resolves
# PerlLib via "$FindBin::Bin/../PerlLib" and invokes "../Launch_PASA_pipeline.pl"
# as a literal relative path, both of which assume exactly one directory level
# below the repo root (matching sample_data/'s own location).
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"

# Launch_PASA_pipeline.pl (via -r/-C resume+recreate) keeps its OWN internal
# checkpoint state beyond just the sqlite file and the outer __chkpts_* dir
# (e.g. pasa_run.log.dir/, per-step .ok markers) -- selectively deleting only
# the obvious artifacts left stale internal checkpoints that skipped table
# creation after a DB recreate, causing "no such table: clusters" on rerun.
# Always regenerate RUN_DIR fresh from the checked-in sample_data/ fixture
# instead of trying to enumerate every generated file.
rm -rf "$REPO_ROOT/$RUN_DIR"
cp -r "$REPO_ROOT/sample_data" "$REPO_ROOT/$RUN_DIR"
cd "$REPO_ROOT/$RUN_DIR"
sed -i "s#^DATABASE=.*#DATABASE=$REPO_ROOT/$RUN_DIR/bench_pasa.sqlite#" \
  sqlite.confs/alignAssembly.config sqlite.confs/annotCompare.config

echo "=== Running PASA align-assembly benchmark in $RUN_DIR (CPU=$CPU) ==="
# --stringent_alignment_overlap matches production funannotate train's own
# Launch_PASA_pipeline.pl invocation (funannotate/train.py passes 30.0) --
# without it, Launch_PASA_pipeline.pl takes the reassign_clusters_via_valid_align_coords.dbi
# fallback path instead of assign_clusters_by_stringent_alignment_overlap.dbi,
# which is one of the scripts these fixes target.
/usr/bin/time -v ./__run_sample_pipeline.pl \
  --align_assembly_config sqlite.confs/alignAssembly.config \
  --annot_compare_config sqlite.confs/annotCompare.config \
  --TRANSDECODER --ALIGNERS minimap2 --CPU "$CPU" --just_align_assembly \
  --stringent_alignment_overlap 30 "$@" \
  > pasa_run.stdout.log 2> pasa_run.log

echo "=== Done. Parsing per-command timing from pasa_run.log ==="
python3 "$REPO_ROOT/benchmarks/parse_bench_log.py" pasa_run.log
