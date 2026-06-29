# PASA Pipeline Performance Optimization Report

## Executive Summary

This document outlines a comprehensive review of the PASApipeline codebase with focus on:
1. Performance bottlenecks and optimization opportunities
2. File size reduction through compression
3. Algorithm improvements suitable for Rust implementation
4. Implemented optimizations and performance testing infrastructure

## Rust Migration - Implemented Tiers

### Tier 1: cdna_alignment_assembler with Interval Tree (IMPLEMENTED)

**Location**: `pasa_rust/pasa-assembler/`

The original C++ `determine_compatibilities_and_encapsulations()` method uses an
O(n²) all-vs-all comparison. For large numbers of alignments (thousands), this
becomes the primary bottleneck.

**Optimization**: Replaced the O(n²) pairwise comparison with an interval tree
that enables O(n log n + k) overlap queries, where k is the number of actually
overlapping pairs. For genomic alignments where most don't overlap, this reduces
the number of `canMerge` calls from O(n²) to O(k).

```rust
// Interval tree built from sorted alignment coordinates
struct IntervalTree {
    sorted_starts: Vec<(i32, usize)>,  // Sorted by lend
}

// For each alignment, query the interval tree to find overlapping alignments
// Only perform the full canMerge check on overlapping pairs
fn determine_compatibilities_and_encapsulations(&mut self) {
    let interval_tree = IntervalTree::new(&self.alignments);
    for i in 0..self.num_alignments {
        let candidates = interval_tree.query_overlaps(
            &self.alignments, self.alignments[i].get_coords()
        );
        for &j in &candidates {
            if j <= i { continue; }
            if self.can_merge(i, j) { /* ... */ }
        }
    }
}
```

### Tier 2: Lobject Containment Tracking with Bitset (IMPLEMENTED)

**Location**: `pasa_rust/pasa-assembler/src/lobject.rs`

The original C++ `Lobject` uses `vector<bool>` for containment tracking, which
is a packed bit representation but with poor cache locality and no SIMD
exploitation.

**Optimization**: Replaced `vector<bool>` with a custom bitset backed by
`Vec<u64>`. The key algorithmic improvement is in `num_unique_contained`:
instead of iterating bit-by-bit in O(n), we use bitwise XOR + hardware popcount
to compute the result in O(n/64) with 64 bits processed per word.

```rust
pub fn num_unique_contained(&self, other: &Lobject) -> i32 {
    let mut num = 0i32;
    for i in 0..min_words {
        // self_word & !other_word → bits in self but not in other
        // .count_ones() → hardware popcount instruction
        let diff = self.contained_cdna_indices[i] & !other.contained_cdna_indices[i];
        num += diff.count_ones() as i32;
    }
    num
}
```

### Tier 3: cdbyank_rust + faidx_rust (IMPLEMENTED)

**Location**: `pasa_rust/cdbtools/`

Two Rust CLI tools replace the C++ `cdbyank` with modern, efficient I/O:

1. **`cdbyank_rust`** — Drop-in replacement for C++ `cdbyank`
   - Reads CDB `.cidx` index format (backward compatible with `cdbfasta`)
   - In-memory hash lookup using djb2 algorithm (matches C++ hash)
   - Single `read` syscall per record (vs byte-by-byte in C++)
   - Supports: `-a` (fetch by accession), `-l` (list keys), `-n` (count)

2. **`faidx_rust`** — Preferred alternative using samtools `.fai` index
   - Reads samtools `.fai` plain-text index format
   - Supports full sequence retrieval and sub-range extraction
   - `faidx_rust genome.fa chr1` → full sequence
   - `faidx_rust genome.fa chr1:1000-2000` → sub-range

**Perl Integration** (`PerlLib/CdbTools.pm`):
- `cdbyank()` auto-detects `cdbyank_rust` in PATH, falls back to C++ `cdbyank`
- `get_seq()` — new API using `faidx_rust` with samtools `.fai` index
- `get_seq_range()` — sub-range extraction via `faidx_rust`
- Auto-creates `.fai` index via `samtools faidx` if missing

**Verified**: Output identical to C++ `cdbyank` for all test cases.

### Tier 4: slclust with HashSet-Based Duplicate Detection (IMPLEMENTED)

**Location**: `pasa_rust/slclust/`

The original C++ `Graph::addLinkedNode()` uses a linear scan through the
adjacency list to check for duplicates, making it O(degree) per edge addition.

**Optimization**: Replaced the `vector<Graphnode*>` adjacency list with
`HashSet<usize>` for O(1) average-case duplicate detection. Also replaced the
recursive DFS (which required `ulimit -s unlimited` for large clusters) with an
iterative DFS using an explicit stack.

```rust
pub struct Graph {
    nodes: Vec<GraphNode>,
    node_lookup: HashMap<String, usize>,  // O(1) lookup vs O(log n)
}

struct GraphNode {
    name: String,
    neighbors: HashSet<usize>,  // O(1) dedup vs O(degree) linear scan
    marked: bool,
}

// Iterative DFS (eliminates stack overflow risk)
fn print_clusters(&mut self) {
    for start in 0..self.nodes.len() {
        if self.nodes[start].marked { continue; }
        let mut stack = vec![start];
        while let Some(idx) = stack.pop() {
            if self.nodes[idx].marked { continue; }
            self.nodes[idx].marked = true;
            for &neighbor in &self.nodes[idx].neighbors {
                if !self.nodes[neighbor].marked {
                    stack.push(neighbor);
                }
            }
        }
    }
}
```

**Known limitation**: slclust_rust is slower than C++ at all tested scales.
The HashMap/HashSet overhead outweighs the O(1) dedup benefit for small-to-medium
degree graphs. Pipeline currently falls back to C++ slclust. The primary benefit
of slclust_rust is iterative DFS eliminating `ulimit -s unlimited` — a reliability
improvement, not a speed improvement.

## Performance Testing

A performance test suite has been created at `PerlLib/perf_tests.pl`.

Run with:
```bash
perl PerlLib/perf_tests.pl [ITERATIONS]
```

Rust component tests:
```bash
cd pasa_rust && cargo test --release
```

## Implemented Optimizations

### 1. Gene_obj Caching (HIGH IMPACT) — DONE

**File**: `PerlLib/Gene_obj.pm`

**Issue**: `get_exons()` was sorting exons on every call, even when the gene object never changed.

**Fix**: Cache the sorted exons on first access:
```perl
sub get_exons {
    my ($self) = shift;
    if ($self->{mRNA_exon_objs} && @{$self->{mRNA_exon_objs}}) {
        unless ($self->{_sorted_exons}) {
            my @exons = @{$self->{mRNA_exon_objs}};
            @exons = sort {$a->{end5}<=>$b->{end5}} @exons;
            if ($self->{strand} eq '-') {
                @exons = reverse @exons;
            }
            $self->{_sorted_exons} = \@exons;
        }
        return @{$self->{_sorted_exons}};
    } else {
        return ();
    }
}
```

**Impact**: ~2.77x faster for repeated get_exons() calls.

Also fixed redundant sorting in:
- `create_cDNA_sequence()` 
- `create_CDS_sequence()`

### 2. Fasta_retriever File Handle Reuse (HIGH IMPACT) — DONE

**File**: `PerlLib/Fasta_retriever.pm`

**Issue**: File was reopened on every `get_seq()` call.

**Fix**: Reuse existing handle unless closed:
```perl
my $fh = $self->{fh};
unless ($fh && fileno($fh)) {
    $fh = $self->refresh_fh();
}
```

Also fixed O(n²) string concatenation with array join:
```perl
my @seq_lines;
while (<$fh>) {
    last if /^>/;
    push @seq_lines, $_;
}
my $seq = join('', @seq_lines);
```

Added compression support (.gz, .bz2):
```perl
sub _open_compressed {
    my $filename = shift;
    if ($filename =~ /\.gz$/) {
        open(my $fh, '-|', "zcat $filename |") or die...;
        return ($fh, 1);
    }
    # ... similar for .bz2
}
```

### 3. GFF3_utils Regex Precompilation (MEDIUM IMPACT) — DONE

**File**: `PerlLib/GFF3_utils.pm`

**Issue**: Regex patterns compiled on every line in hot loop.

**Fix**: Precompile at package level:
```perl
my $ID_RE = qr/ID="?([^;\s"]+)"?;?/;
my $NAME_RE = qr/Name="?([^\;"]+)"?/;
my $NOTE_RE = qr/Note="?([^\;"]+)"?/;
my $ALIAS_RE = qr/Alias=([^;]+)/;
my $PARENT_RE = qr/Parent="?([^;\s"]+)"?;?/;
my $FEAT_TYPE_RE = qr/^(gene|mRNA|transcript|CDS|exon)$/;
```

### 4. PSL_parser Caching and O(n²) Fix (MEDIUM IMPACT) — DONE

**File**: `PerlLib/PSL_parser.pm`

**Issue**: 
- `get_per_id()` recalculated on every call
- `toString()` used `shift` in loop (O(n²))

**Fix**: Cache percentage identity:
```perl
sub get_per_id {
    my $self = shift;
    return $self->{_per_id} //= do {
        my $matches = $self->get_match_count();
        my $mismatches = $self->get_mismatch_count();
        sprintf("%.2f", $matches / ($matches + $mismatches) * 100);
    };
}
```

Fix O(n²) shift loop with index iteration:
```perl
for (my $i = 0; $i < @genome_coords; $i++) {
    push @align_parts, 
        $genome_coords[$i]->[0] . "...";
}
$ret_text .= "\n" . join("....", @align_parts) . "\n";
```

### 5. Database Batch Operations — DONE

**Files**: `PerlLib/DB_connect.pm`, `scripts/populate_alignments_via_btab.dbi`, `scripts/import_spliced_alignments.dbi`

**Issue**: Each accession triggered a separate query in loop (N+1 query pattern).

**Fix**: 
- `DB_connect.pm` `RunMod`: switched from `do()` to `prepare_cached()` for repeated statement handle reuse
- `populate_alignments_via_btab.dbi`: added transaction support (`AutoCommit=0`, batch commits every 1000 records), used `prepare_cached` for alignment inserts
- `import_spliced_alignments.dbi`: used `prepare_cached` for `store_cluster_links` batch operations

### 6. Pipeliner Backtick Removal — DONE

**File**: `PerlLib/Pipeliner.pm`

**Issue**: Using backticks for simple file operations (spawning subprocesses unnecessarily).

**Fix**:
```perl
# Replace: my $errmsg = `cat $tmp_stderr`;
open my $err_fh, '<', $tmp_stderr;
my $errmsg = do { local $/; <$err_fh> };

# Replace: `touch $checkpoint_file`;
open my $cp_fh, '>', $checkpoint_file;
close $cp_fh;
```

### 7. Database Indexes — DONE

**Files**: `schema/cdna_alignment_mysqlschema`, `schema/cdna_alignment_sqliteschema`

**Issue**: Missing composite indexes for common query patterns.

**Fix**: Added composite indexes to both MySQL and SQLite schemas:
- `alignment(align_id, lend, rend)` — for range queries joining on align_id
- `align_link(prog, cluster_id)` — for filtering by both prog and cluster_id

### 8. Compression for Intermediate Files — DONE

**Files**: `PerlLib/CompressUtils.pm` (new), `PerlLib/Pipeliner.pm`, `Launch_PASA_pipeline.pl`

**Issue**: Pipeline creates large uncompressed intermediate files (GFF3, BED, GTF).

**Fix**:
- New `CompressUtils.pm` module with transparent gzip/bzip2 file compression utilities
- Pipeliner gains `compress_intermediates` option and per-Command `set_compress_files`/`get_compress_files` methods
- After a command completes, specified intermediate files are gzipped in-place
- `Launch_PASA_pipeline.pl` gains `--compress-intermediates` CLI flag
- `compress_files` metadata added to GFF3/BED/GTF output commands for valid/failed alignments and pasa_assemblies

## High-Priority Optimizations Remaining

### C++ Components

#### 1. cdna_alignment_assembler.cpp - O(n²) Compatibility Matrix — DONE

**Location**: `pasa_cpp/cdna_alignment_assembler.cpp:660-695`

**Issue**: Nested loop for compatibility checking is O(n²)

**Fix**: Added early termination based on sorted `lend` positions. Since alignments are sorted by `lend`, once `alignments[j].lend > alignments[i].rend`, no further alignments can overlap `i`. This reduces the inner loop from O(n) to O(k) where k is the number of potentially overlapping alignments.

#### 2. map<int,bool> for accountedFor — DONE

**Location**: `pasa_cpp/cdna_alignment_assembler.cpp:122`

**Issue**: Using `map` where vector suffices - O(log n) instead of O(1)

**Fix**: Replaced `map<int,bool>` with `vector<bool>` in:
- `accountedFor` in `assembleAlignments()`
- `tracker` in `forwardTrace()` and `backTrace()`
- `uniqueMap` in `unique_entries()`
- `get_max_missing_Lobj()` signature updated to `vector<bool>&`

#### 3. cdbyank Byte-by-Byte I/O

**Location**: `pasa-plugins/cdbtools/cdbyank.cpp:214-245`

**Issue**: Reading large records byte-by-byte for output

**Fix**: Use buffered I/O with `fread`/`fwrite` or memory-mapped files

**Rust Implementation Value**: MEDIUM-HIGH
- Parallel I/O with async
- Zero-copy parsing with `bytes` crate

#### 4. slclust O(n) Duplicate Detection

**Location**: `pasa-plugins/slclust/graphnode.cpp:16-29`

**Issue**: Linear search for duplicate detection in `addLinkedNode`

**Fix**: Use `std::unordered_set<std::string>`

**Rust Implementation Value**: MEDIUM
- Idiomatic Rust `HashSet` provides O(1) average case

### Database Optimizations

#### 1. Lock Contention in DB_connect — DONE

**File**: `PerlLib/DB_connect.pm`

**Issue**: Global lock (`$LOCKVAR` via `threads::shared`) serialized all database operations across all threads, even when using independent connections.

**Fix**: Removed the global lock entirely from `do_sql_2D` and `RunMod`. Each thread creates its own connection via `connect_to_db`, so each has its own `$dbh` and `prepare_cached` cache. No cross-thread sharing occurs, making the lock unnecessary. Also removed the `use threads; use threads::shared;` imports that were only needed for the lock.

#### 2. Database Schema Compression

**Issue**: SQLite/MySQL databases can grow large without compression

**Fix Options**:
- Enable page-level compression for MySQL
- Use SQLite with compression extension
- Implement streaming compression for exports

### Memory Optimization

#### 1. Gene_obj Memory Leaks

**Issue**: Large sequence strings stored in gene objects

**Fix**: Clear sequences after use or use lazy loading

## Parallelization Opportunities

### 1. Alignment Processing
- Batch alignments across multiple cores
- Use Rust's `rayon` for data parallelism
- Perl: Thread::Semaphore for controlled parallelism

### 2. GFF3 Parsing
- Process chromosomes in parallel
- Merge results at the end

### 3. Database Bulk Operations
- Batch inserts with prepared statements
- Parallel writer threads

## Benchmarking Infrastructure

Created `PerlLib/perf_tests.pl` with tests for:
- Fasta_reader throughput
- Gene_obj caching effectiveness
- PSL_parser performance
- Regex compilation vs inline
- String operations

Run with: `perl PerlLib/perf_tests.pl [ITERATIONS]`

## Benchmark Results

Benchmarks run with `run_benchmarks.sh`. Test data: 500 sequences × 5000 bp,
500 alignments, 4994 clustering pairs.

### Tool-Level Benchmarks (Rust vs C++)

| Tool | Rust Time | C++ Time | Speedup | Notes |
|------|-----------|----------|---------|-------|
| cdbyank (500 lookups) | 1.70s | 2.45s | 1.44x | In-memory djb2 hash, bulk I/O |
| faidx_rust (500 lookups) | 2.29s | 2.45s | 1.07x | .fai index, marginal gain for small seqs |
| pasa (500 alignments) | 30.7ms | 41.4ms | 1.35x | Interval tree O(n log n + k) |
| slclust (4994 pairs) | 28.9ms | 5.5ms | 0.19x | HashSet overhead dominates at small scale |

### PASA Assembly Scaling (Interval Tree vs O(n²))

| Alignments | C++ Time | Rust Time | Speedup |
|------------|----------|-----------|---------|
| 100 | 13.9ms | 10.5ms | 1.33x |
| 500 | 47.9ms | 32.3ms | 1.49x |
| 1000 | 236.6ms | 158.0ms | 1.50x |
| 2000 | 5916.6ms | 3014.0ms | 1.96x |

Speedup increases with alignment count — interval tree's O(n log n + k)
complexity becomes more beneficial as n grows. At n=2000, Rust is ~2x faster.

### slclust Scaling (HashSet vs vector)

| Pairs | C++ Time | Rust Time | Speedup |
|-------|----------|-----------|---------|
| 998 | 6.51ms | 8.73ms | 0.75x |
| 4998 | 5.45ms | 19.41ms | 0.28x |
| 9997 | 5.34ms | 31.20ms | 0.17x |
| 19999 | 5.58ms | 64.55ms | 0.09x |

**Known limitation**: slclust_rust is slower than C++ at all tested scales.
The HashMap/HashSet overhead outweighs the O(1) dedup benefit for small-to-medium
degree graphs. The C++ version's `vector` with linear dedup has better cache locality.

**Pipeline fallback**: `SingleLinkageClusterer.pm` now prefers C++ `slclust` (with
`ulimit -s unlimited`) over `slclust_rust`, since the C++ version is faster at all
tested scales. The Rust version is retained as a fallback for environments where
`ulimit` is unavailable or restricted.

**Primary benefit of slclust_rust**: Iterative DFS eliminates the need for
`ulimit -s unlimited` and prevents stack overflow on large clusters — a reliability
improvement, not a speed improvement.

### Perl-Level Benchmarks

```
Fasta_reader:          ~109,000 sequences/sec
Gene_obj get_exons:    ~526,000 calls/sec (cached)
PSL_parser toString:   ~4,187 ops/sec
```

Caching provides approximately 2.77x speedup for frequently accessed gene data.

---

*Generated: 2026-06-28*
*PASApipeline Version: 2.5.3*
