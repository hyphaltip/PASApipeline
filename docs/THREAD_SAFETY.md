# Thread Safety in PASA Perl Scripts

See also: [PERFORMANCE_OPTIMIZATION.md](PERFORMANCE_OPTIMIZATION.md),
[CHROMOSOME_SHARDING.md](CHROMOSOME_SHARDING.md).

## Status

As of 2026-07 there is **no confirmed threading defect** in
`assemble_clusters.dbi`. An earlier revision of this document claimed one; that
claim was wrong and is retracted below, along with the benchmarking mistake that
produced it. Two hardening changes were kept, and the harness that produced the
false signal was fixed.

## Retracted claim: the temp-file token race

An earlier revision asserted that `PerlLib/CDNA/PASA_alignment_assembler.pm`
raced across threads because its temp files were named:

```perl
srand();
my $uniq_token = time() . "-" . rand();
```

and that `srand()` reseeds from the pid — identical across ithreads — so two
threads in the same second would derive the same token and clobber each other.

**That mechanism was never verified, and testing does not support it.** On the
Perl in the funannotate environment, ithreads receive independent RNG state, and
an explicit `srand()` does not collapse it:

```
no srand():   6 distinct tokens of 6 threads  -> all distinct
with srand(): 8 distinct tokens of 8 threads  -> all distinct
```

The claim was written from reading the code, without running the test and
without reading the actual error message — which was sitting in the benchmark's
own per-run log the whole time.

## What the failures actually were

Sweeping `assemble_clusters.dbi` across thread counts showed `exit=25` and a
different output checksum at every `-T >= 2`, with `-T 1` clean. The real error:

```
Error, spliced orient in db (+) differs from calculated spliced orient (-)
  at PerlLib/Ath1_cdnas.pm line 744
  Ath1_cdnas::batch_create_alignment_objs_by_id(...)
  called at scripts/assemble_clusters.dbi line 206
```

and immediately above it in the same log:

```
-missing faidx file: .../genome.fasta.fai, extracting positions directly.
```

**The benchmark harness had copied `genome.fasta` without `genome.fasta.fai`.**
That is not cosmetic. `Fasta_retriever::_init` branches on it:

* **With `.fai`** (what production runs have): only the position index is built.
  `$self->{fh}` stays `undef`, so each worker thread opens its own filehandle on
  first `get_seq()`.
* **Without `.fai`**: `_init` scans the FASTA itself and leaves an **open
  filehandle in `$self->{fh}`** — created *before* threads are spawned, and
  therefore cloned into every worker.

So every threaded measurement exercised a `Fasta_retriever` code path that
production never uses under threading. The resulting bad sequence reads made the
computed splice orientation disagree with the database value, which is precisely
what `batch_create_alignment_objs_by_id` is designed to `confess` on.

The failures are therefore best explained as an **artifact of the harness**, not
a PASA defect. `benchmarks/sweep_threads.sh` now copies `.fai`/`.cidx` and
aborts if the `.fai` is absent.

## What the corrected sweep showed

With `.fai` staged correctly, every thread count **exits 0** — confirming the
failures above were the harness. But two real findings replaced them.

### `assemble_clusters.dbi` output is nondeterministic even single-threaded

Running `-T 1` twice, and `-T 2` twice, on identical input:

| `-T` | run | wall | exit | output checksum |
|------|-----|--------|------|-----------------|
| 1 | A | 401.3s | 0 | `1bbda95c…` |
| 1 | B | 415.4s | 0 | `88ccca08…` |
| 2 | A | 337.1s | 0 | `1a388ae5…` |
| 2 | B | 337.1s | 0 | `68ffc111…` |

**Two serial runs of the same input produce different output.** So this is not
a threading defect at all — threading merely makes it more visible. The
checksum is taken over order-normalized output, so this is a genuine content
difference, and it means `assemble_clusters.dbi` is not reproducible run to run
regardless of thread count.

Prime suspect is `scripts/assemble_clusters.dbi:193`:

```sql
select al.align_acc, al.align_id, al.lend, al.score from align_link al
where al.cluster_id = ? and al.validate = 1
```

There is **no `ORDER BY`**. SQL guarantees nothing about row order without one,
and InnoDB's actual order can shift with buffer-pool state and concurrency.
`assemble_alignments(@alignments)` consumes them in the order returned, so a
different row order can produce different assemblies. Not yet confirmed by
patching the query -- treat as the leading hypothesis, not a diagnosis.

### `assemble_clusters.dbi` barely benefits from threading

| `-T` | wall | speedup | efficiency | peak RSS |
|------|--------|---------|------------|----------|
| 1 | 647.8s | 1.00x | 100% | 145.2 MB |
| 2 | 543.1s | 1.19x | 60% | 195.9 MB |
| 4 | 509.0s | 1.27x | 32% | 259.2 MB |
| 7 | 545.7s | 1.19x | 17% | 361.7 MB |
| 8 | 543.7s | 1.19x | 15% | 346.7 MB |
| 16 | 553.3s | 1.17x | 7% | 346.1 MB |

Peak is 1.27x at `-T 4`, degrading beyond that. The contig work distribution on
this genome allows roughly 5x, so something is serializing the threads. Leading
candidate is `assemble_clusters.dbi:191`, which calls
`DB_connect::reconnect_to_server` **once per cluster** (~16k times): on MySQL
that opens a brand-new connection and discards the `prepare_cached` statement
cache every iteration.

## Hardening changes that were kept

These are worth keeping on their own merits. Neither is a fix for anything
observed above.

### 1. Atomic temp-file allocation in `pasa_cpp_assemblies`

The temp files exchanged with the external `pasa` binary are now allocated with
`File::Temp::tempfile`, which creates each atomically with `O_EXCL` and retries
on collision, rather than constructing a name and opening it. Bareword
`TMPIN`/`TMPOUT` handles were replaced with lexicals.

Rationale, independent of any observed bug: constructing a name and then opening
it is a *guess* that the name is unclaimed; `O_EXCL` is a *check-and-claim*.
Collisions become impossible rather than improbable, at no cost and with no new
dependency (`File::Temp` is core). It also covers a case no constructed token
can: a shared or network `$TMPDIR`, where pids repeat across hosts.

A UUID would be the weaker option here — it lowers collision probability but
still assumes rather than verifies, and no UUID module is installed in the
funannotate environment. UUIDs are the right tool for *identity* (tracing a work
unit across logs), not for *exclusivity*.

### 2. Unique failure-diagnostic filename

On a failed `pasa` invocation the code did `system "mv $pasa_input
pasa_killer.input"` — a fixed name in the current directory, so concurrent
failing threads destroyed each other's diagnostic. Each failure now gets its own
file.

## Lessons for benchmarking parallel code

1. **Read the error before theorising about it.** The mechanism above was
   invented from code reading while the real message sat unread in a log.
2. **Stage the full input set.** A missing sidecar index silently moved the code
   onto a different branch. If a benchmark copies inputs, it must copy
   *everything the real run had* — and assert on what it cannot do without.
3. **Checksum order-normalized output alongside timings.** This is what revealed
   the runs were failing rather than merely slow;
   `benchmarks/sweep_threads.sh` records it per thread count and
   `benchmarks/sweep_summarize.py` flags divergence from the serial reference.
4. **Make instrumentation failures non-fatal.** Under `set -o pipefail`, a
   rejected checksum query aborted an entire sweep immediately after a
   *successful* 24-minute run, discarding every remaining data point. Verifiers
   should degrade to a recorded failure, not kill the experiment.
5. **`-T` defaults to 2, not 1**, in `assemble_clusters.dbi`,
   `classify_alt_splice_isoforms.dbi`, `alignment_assembly_to_gene_models.dbi`,
   `find_alternate_internal_exons.dbi`,
   `assign_clusters_by_stringent_alignment_overlap.dbi`, and
   `extract_transcript_alignment_clusters.dbi`. These scripts are threaded in
   practice even when no thread option was passed.

## Open items in threaded scripts

Found while auditing; none yet fixed.

| Location | Issue |
|---|---|
| `SingleLinkageClusterer.pm:37-64` | Hardcodes `/tmp`, ignoring `TMPDIR`/`SCRATCH`; constructs its temp name rather than claiming it; bareword `PAIRLIST`/`CLUSTERS` handles. Reached from threads via `subcluster_builder.dbi`. See `PerlLib/Pasa_tmpdir.pm`. |
| `classify_alt_splice_isoforms.dbi:172` | Concat loop reopens the destination with `>` *inside* the loop, so `alt_splicing_analysis.results.out` retains only the last subcluster. |
| `classify_alt_splice_isoforms.dbi:159` | One `threads->create` per subcluster (thousands), each cloning the interpreter. Measured ithread clone cost scales with pre-spawn state: 2.3 ms empty, 131 ms at 100k hash keys, 725 ms at 500k. |
| `assemble_clusters.dbi:191` | `DB_connect::reconnect_to_server` per cluster (~16k times). On MySQL this opens a fresh connection and discards the `prepare_cached` statement cache each iteration. |
