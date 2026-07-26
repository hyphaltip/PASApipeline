# Chromosome/Contig-Level Sharding for subcluster_builder.dbi

See also: [PERFORMANCE_OPTIMIZATION.md](PERFORMANCE_OPTIMIZATION.md) (in-process
threading of the other per-`asmbl_id` scripts).

## Problem

`scripts/subcluster_builder.dbi` loops serially over every contig/scaffold
(`annotdb_asmbl_id`) in the genome. In a profiled production funannotate/PASA
training run (Aspergillus flavus, ~20 scaffolds, MySQL backend), this step cost
~28 minutes of a ~3h40m PASA run — one of the largest single line items,
alongside `assemble_clusters.dbi` and `classify_alt_splice_isoforms.dbi`
(both already parallelized with `-T`/`Thread_helper`; see
PERFORMANCE_OPTIMIZATION.md § Parallelization).

## Why not in-process threading

Unlike its sibling per-`asmbl_id` scripts, `subcluster_builder.dbi` calls
`SingleLinkageClusterer::build_clusters`, which forks external `slclust`/
`slclust_rust` subprocesses per cluster (`PerlLib/SingleLinkageClusterer.pm`).
In-process threading of this script was implemented and benchmarked (GitHub
issue #7), then closed as won't-fix: it was **slower than the serial baseline
at every thread count tested**, because forking from a large multi-threaded
Perl process is more expensive than forking from a small single-threaded one —
a structural cost that gets worse, not better, at production scale.

## Approach: external process-level sharding

`sharding/shard_subcluster_builder.sh` runs one `subcluster_builder.dbi`
process per contig (each single-threaded, so no fork-cost penalty), fanning
out with GNU `parallel`, then merges the per-contig outputs back into exactly
the output file and checkpoint file `Launch_PASA_pipeline.pl` expects for this
stage — so it must be run **before** `Launch_PASA_pipeline.pl` / `funannotate
train` reaches the `subcluster_builder.dbi` step. `PerlLib/Pipeliner.pm` skips
any step whose checkpoint file already exists, so once the driver succeeds, the
main pipeline proceeds straight to `subcluster_loader.dbi` without ever
invoking `subcluster_builder.dbi` itself. No changes to `Launch_PASA_pipeline.pl`
were needed.

Why GNU `parallel` rather than a Slurm array: this is single-node fan-out — the
production run that motivated this had a single-node allocation with 256 cores
and used only 2 (the real waste was CPU allocation, not a need for multi-node
scale-out). `parallel --joblog` also gives an exact per-contig exit code and
signal, which is simpler and more reliable to verify than reconciling
`sbatch --wait`'s array-level exit status against `sacct` per task index.

### `-R <asmbl_id>` flag

`subcluster_builder.dbi` now accepts `-R <asmbl_id>` to restrict processing to
a single contig, mirroring `cDNA_annotation_comparer.dbi`'s
`--RESTRICT_SINGLE_CONTIG`. This is what each shard process passes.

### Concatenation-order guarantee

`subcluster_loader.dbi` (which consumes `subcluster_builder.dbi`'s output)
tracks a single "current cluster" scalar that resets on every
`Processing cluster:` header, and every contig's `cluster_id`s are disjoint.
So per-contig output chunks may be concatenated in **any order across
contigs** — the driver does not need to preserve the original single-process
contig ordering — as long as each contig's own chunk stays internally intact
(never interleaved with another contig's lines, which the driver guarantees by
writing each shard to its own file).

## Usage

```
sharding/shard_subcluster_builder.sh <alignAssembly.config> <genome.fasta> <PASA_LOG_DIR> [--jobs N]
```

Run from the same PASA run directory `Launch_PASA_pipeline.pl` would use, with
the same `alignAssembly.config` (must contain `DATABASE=`) and the same
`PASA_LOG_DIR` (normally `pasa_run.log.dir`). Requires GNU `parallel` on
`PATH`. `--jobs` defaults to `nproc`, capped at the contig count.

Then run `funannotate train` / `Launch_PASA_pipeline.pl` as usual — it will
detect the checkpoint and skip straight past this step.

## Failure handling

Verification happens **before** the output file or checkpoint is written, so a
failed or partial run never gets treated as done:

1. `parallel`'s own exit status is checked first (nonzero if any job failed).
2. Every row of `joblog.tsv` must show `Exitval == 0` and `Signal == 0` — this
   also catches OOM-killed shards (`Signal` reflects the kill signal), and is
   more precise than treating non-empty shard STDERR as a failure signal
   (`subcluster_builder.dbi`'s normal per-contig progress trace already writes
   to STDERR even on success).
3. Every contig must have a corresponding shard output file.
4. Structural completeness check: the number of `Processing cluster:` blocks a
   contig's shard emitted must match an independent
   `select count(distinct cluster_id) from clusters where annotdb_asmbl_id = ?`
   query for that contig (`scripts/count_clusters_for_asmbl_id.dbi`) — catches
   a shard that exited 0 but was truncated before all clusters were written.

If any check fails, the driver exits non-zero, leaves the shard workspace
(including `joblog.tsv`) in place for inspection, and does **not** touch the
checkpoint — the main pipeline will still correctly run this step itself (or
the driver can be fixed and rerun) if invoked afterward.

## New/changed files

- `scripts/subcluster_builder.dbi` — added `-R <asmbl_id>`.
- `scripts/list_asmbl_ids.dbi` (new) — enumerates contigs for the driver.
- `scripts/count_clusters_for_asmbl_id.dbi` (new) — structural completeness
  check helper.
- `PerlLib/Ath1_cdnas.pm` — added `get_all_annotdb_asmbl_ids($dbproc)`.
- `sharding/shard_subcluster_builder.sh` (new) — the driver described above.
