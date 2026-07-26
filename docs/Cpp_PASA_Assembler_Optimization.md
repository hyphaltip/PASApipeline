# C++ pasa_assembler Optimization Plan

## Goal

Apply algorithmic and micro-optimizations to the C++ `pasa_assembler` (in `pasa_cpp/`) to improve throughput, reduce memory footprint, and eliminate known bugs — without changing the toolchain (staying in C++).

The Rust implementation at `pasa_rust/pasa-assembler/` validates several of these approaches; the C++ code remains the canonical reference for the PASA assembly algorithm.

## Hard constraint: bit-identical output

The assembler must emit exactly the same assemblies, in the same order, as the pre-optimization code. Several loops break score ties by "first candidate at the maximum wins" (`curr_total_score > top_score`), so **traversal order is part of the algorithm's semantics**. Any change that alters the order in which candidate indices are visited is a behavior change, not a refactor. See the "Order-sensitive sites" section below.

## Which assembler runs

The C++ assembler is the default. `make` builds and installs `bin/pasa` and no longer builds the optional Rust binaries; pass `make WITH_RUST=1` (or run `make rust`) to get `pasa_rust`, `slclust_rust`, `cdbyank_rust` and `faidx_rust` as well.

`PASA_alignment_assembler.pm` previously probed for `pasa_rust` first and silently preferred it whenever it was on `$PATH`. It now resolves `pasa` first and falls back to `pasa_rust` only if no C++ binary exists. Set `PASA_ASSEMBLER=pasa_rust` in the environment to force the Rust path — useful for A/B comparison.

Note that `scripts/install.sh` and the conda recipe still build and install the Rust binaries unconditionally; `conda-recipe/meta.yaml` runs `pasa_rust --help` as a build test. Those were left alone. With the runtime preference flipped, a conda install still gets the C++ assembler by default.

## Optimization Summary

| # | Optimization | Area | Impact | Status |
|---|-------------|------|--------|--------|
| 1 | Sparse adjacency (sorted `vector<int>` rows) | `determine_compatibilities_and_encapsulations` | HIGH — ~25% less peak RSS | DONE |
| 2 | Bitset with hardware POPCNT | Lobject containment tracking | MEDIUM | DONE |
| 3 | `twoDarray`/`free2Darray` | Memory management | — | DELETED (see below) |
| 4 | Const-reference sort comparators | alignment and segment sorts | MEDIUM | DONE |
| 5 | Move semantics + enum segment type | Assembly construction | **HIGH — ~4× on its own** | DONE |
| 6 | Allocation-free splice-coordinate sets | `mergeAlignments` | HIGH — ~2× on merge-heavy input | DONE |
| 7 | OpenMP parallelism | compatibility scan only | NEGLIGIBLE — see below | DONE, no measurable gain |
| 8 | Compiler flags | Build system | MEDIUM (`-msse4.2` is the one that matters) | DONE |

Measured on a synthetic 3000-alignment locus (`-O3`, single socket):

| Build | Wall | Peak RSS |
|---|---|---|
| pre-optimization (`master`) | 21.0 s | 67.6 MB |
| after items 1–4, 6–8 | 9.2 s | 50.7 MB |
| after item 5 as well, `OMP_NUM_THREADS=1` | **2.24 s** | 43.9 MB |
| after item 5 as well, `OMP_NUM_THREADS=4` | 2.24 s | 44.0 MB |

**9.4× faster, 35% less peak RSS, output byte-identical to `master` in all cases.**

---

## Order-sensitive sites

Both DP scans pick the **first** candidate that attains the maximum score, so the direction of traversal decides ties:

- `do_full_Fscan` — original scanned `j = i-1 … 0`, so the **highest** index wins a tie. Rows are stored ascending, so this loop must iterate with `rbegin()/rend()`.
- `do_full_Rscan` — original scanned `j = i+1 … n-1`, so the **lowest** index wins a tie. Ascending iteration is correct here.

An `unordered_set` was used for these rows at one point. That is incorrect: hash iteration order is unspecified and unrelated to index order, so tie-breaking became a function of hash-table layout. Rows are now `vector<int>` sorted ascending, which gives a defined order *and* makes containment lookups a `binary_search`. This is not a stylistic preference — an ascending Fscan produces different assemblies from `master` on realistic staggered input.

`populateLobjects` and `forwardTrace`/`backTrace` are order-insensitive in effect (results feed a bitset), so they need no special handling.

---

## 1. Sparse Adjacency Lists

**Was**: `bool**` ragged arrays (n×n) for `compatibilities` and `encapsulations` — 2×n² bytes, ~200 MB at n=10k, never freed.

**Now**: `vector<vector<int>>`, one ascending row per alignment.

Note on complexity: the alignments are already sorted by `lend` in the constructor, so the candidate scan for row `i` is simply `j = i+1` forward, breaking at the first `lend > rend_i`. That is the O(k)-per-row scan; an explicit interval tree or binary search over a `(lend, index)` array adds nothing, because such an array is just the identity permutation of the already-sorted input. An earlier revision built and binary-searched that array but then scanned it from index 0 each time, which reintroduced Θ(n²) skip work — worse than the loop it replaced.

The realistic claim is **memory**, not asymptotics: the dense matrix is gone. Note that `vector<int>` rows are only a win while the compatibility graph is sparse; each entry costs 4 bytes against the old 1 byte per pair, so a fully dense locus would use more, not less. Measured peak RSS on the benchmark above dropped 25%.

**Files**: `cdna_alignment_assembler.h`, `cdna_alignment_assembler.cpp`

---

## 2. Bitset + Hardware POPCNT

**Was**: `vector<bool>` for `contained_cdna_indices`, with `num_unique_contained` iterating all n elements bit by bit.

**Now**: `vector<uint64_t>`; `num_unique_contained` computes `self[w] & ~other[w]` word-wise and sums with `__builtin_popcountll`, 4 words per iteration for instruction-level parallelism. `forwardTrace`/`backTrace` walk set bits with `__builtin_ctzll` + `x &= x - 1`, which preserves ascending index order.

This optimization only pays off if the compiler actually emits `POPCNT` — see item 8. Note also that a locus with fewer than 64 alignments occupies a single word, so the unrolled path never runs on small input; do not extrapolate its benefit from a small benchmark.

**Files**: `Lobject.h`, `Lobject.cpp`, `cdna_alignment_assembler.cpp`

---

## 3. `twoDarray` / `free2Darray` — deleted

These had empty `free2Darray` bodies (the `delete[]` calls were commented out), leaking O(n²) per call. Item 1 removed the last caller, and a repo-wide grep found no others, so all four functions were **deleted** rather than repaired. An allocation helper whose free is a silent no-op is a trap for the next caller, and there was no next caller to serve.

Both `twoDarray` overloads were also wrong independently of the leak: the row loop was bounded by `y` instead of `x`, so a non-square request would under-allocate rows and walk off the array. Another reason not to keep them around.

**Files**: `common_subs.cpp`, `common_subs.h`

---

## 4. Pass by Const Reference in sort Comparators

`sort_CDNA_alignments` took `CDNA_alignment` by value, deep-copying both alignments (including their segment vectors) on every one of O(n log n) comparisons. Now `const CDNA_alignment&`.

`segmentSortCriteria` in `cdna_alignment.cpp` had the same defect and was missed the first time round — it took two `Alignment_segment`s by value, and it is called from `refineAlignment`, i.e. from *every* `CDNA_alignment` construction, i.e. from every merge step. Now `const Alignment_segment&`, which needed a `const` overload of `Alignment_segment::get_coords`.

It also compared with `<=`, which is not a strict weak ordering: `comp(a, a)` returning true is undefined behavior in `std::sort` and can run the unguarded insertion-sort loop off the end of the range. Now `<`. Output is unchanged on all fixtures.

**Files**: `cdna_alignment.cpp`, `alignment_segment.h`, `alignment_segment.cpp`, `cdna_alignment_assembler.cpp`

---

## 5. Move Semantics and Enum Segment Type

`create_assembly` copied `alignments[alignment_index]` and then round-tripped every merge through a named temporary (`newAssembly`), copy-assigning it into the accumulator — O(k²) segment copies for an assembly of k alignments. Profiling the otherwise-optimized build put roughly half of runtime here, with `std::string` operations alone at ~30%:

```
20.4%  CDNA_alignment_assembler::create_assembly
19.9%  std::string::_M_construct<char*>   (+ ~13% more across string ops)
12.8%  __introsort_loop<Alignment_segment*>
 8.2%  __unguarded_linear_insert<Alignment_segment*>
```

The dominant string cost was **not** the alignment `title`, as first assumed. It was `Alignment_segment::type`, a `std::string` holding one of `"first"`, `"last"`, `"internal"`, `"single"`, `"?"` — one per segment, written on every `refineAlignment` and copied on every segment copy. It is now an `enum segment_type`. The accessors `set_type`/`get_type` had no callers at all; the field is only assigned in `refineAlignment` and compared in `toAlignIllustration`, so the change is contained within `pasa_cpp`.

`CDNA_alignment` and `Alignment_segment` declare no destructor or copy operations, so the compiler already generates move constructors and move assignment for both. Nothing needed writing — the copies came from code that never let a move happen:

- `create_assembly` now move-assigns straight from the returned temporary: `assembly = mergeAlignments(assembly, nextAlignment)`.
- `mergeAlignments` returns `CDNA_alignment(std::move(new_seg_list), orientation)` and reserves the segment vector up front.
- The `CDNA_alignment` constructor takes its segment list by value and now `std::move`s it into the member instead of copying.
- `add_alignment_segment` moves its by-value parameter into the vector.

Together with the enum this took the benchmark from 9.2 s to 2.24 s. The remaining profile is `create_assembly` at 36% — inherent merge work now, not copying.

**Files**: `cdna_alignment.h`, `cdna_alignment.cpp`, `alignment_segment.h`, `alignment_segment.cpp`, `cdna_alignment_assembler.cpp`

---

## 6. Allocation-Free Splice-Coordinate Sets

`mergeAlignments` built two `map<int,bool>` (later `unordered_set<int>`) of splice coordinates on every call. Both are the wrong container: the sets hold one entry per alignment segment — typically fewer than ten — so the hash table's bucket array and per-node allocations dominated, and `mergeAlignments` is called once per merge step of every assembly. Profiling put the hash insert at 17% of runtime with another 12% in `malloc` underneath it.

**Now**: two `vector<int>` scratch buffers held as members and `clear()`ed per call (retaining capacity, so zero allocations after warm-up), with a linear `has_coord` scan. Membership order is irrelevant here, and duplicates are harmless, so no sorting or dedup is needed. This change alone took the benchmark from 19.3 s to 9.2 s.

The lesson generalizes: O(1)-vs-O(log n) is the wrong lens at n < 16. Constant factors and allocator traffic decide it.

**Files**: `cdna_alignment_assembler.h`, `cdna_alignment_assembler.cpp`

---

## 7. OpenMP Parallelism

**Only the compatibility scan is parallelized**, and it buys nothing measurable (2.24 s at both 1 and 4 threads) because that scan is not the bottleneck — `create_assembly` is. Keep it for larger loci, or drop it; do not cite it as a speedup.

Implementation note: each unordered pair is discovered exactly once, by the lower index, so worker threads accumulate results into thread-local buffers that are appended under a single `critical` at the end of the parallel region and scattered into the rows serially. Rows are then sorted, which makes the contents independent of thread count. An earlier revision put a `critical` around every individual insert, which serialized the loop body.

The region carries `if(!DEBUG)`, so `-v` runs it serially. The per-candidate tracing writes to `cout` from inside the loop, and concurrent threads interleave it mid-line into unreadable output. With the guard, `-v` output is identical at any thread count.

**`do_full_Fscan` and `do_full_Rscan` must NOT be parallelized.** An earlier version of this plan claimed all three loops were independent; that is wrong. Fscan reads `Lobjects[j].LscoreF` for `j < i`, and that field is *written* by the same loop when the outer index was `j`. It is a sequential dynamic-programming recurrence. Rscan is the mirror case. Adding `#pragma omp parallel for` there would read stale scores depending on thread scheduling.

**Files**: `cdna_alignment_assembler.cpp`, `Makefile`

---

## 8. Compiler Flags

```makefile
CXX ?= g++
ARCHFLAGS ?= -msse4.2 -mtune=native
CFLAGS  = -O3 $(ARCHFLAGS) -flto -funroll-loops -fopenmp -std=c++11
LDFLAGS = -O3 $(ARCHFLAGS) -flto -fopenmp
```

- `CXX ?=` rather than `CXX =` — the build must respect an externally supplied compiler (clang on macOS, a specific gcc on HPC). Hardcoding `g++` reverts an earlier fix.
- **`-msse4.2` is the flag that actually matters.** Baseline `x86-64` has no `POPCNT`, so without it `__builtin_popcountll` compiles to a `__popcountdi2` libgcc call and item 2 is undermined. Verified on this toolchain: no flag → `__popcountdi2`; `-mtune=native` alone → still `__popcountdi2`; `-msse4.2` → `popcnt` instruction.
- **`-mtune=native` vs `-march=x86-64-v2`**: these are not alternatives. `-march=` raises the instruction-set baseline (and so restricts which machines can run the binary); `-mtune=` only changes the scheduling/cost model and never emits new instructions, so it is always safe to leave in. `-march=x86-64-v2` would be a reasonable baseline — it includes SSE4.2 — but GCC only accepts it from version 11, and the cluster's default is GCC 8.5, where it is a hard error. `-msse4.2 -mtune=native` gets the same codegen benefit and builds everywhere.
- **`-march=native` was removed.** It bakes the build host's ISA into the binary. On a heterogeneous cluster, or for a conda/container artifact built elsewhere, that is a `SIGILL` waiting to happen. Override `ARCHFLAGS=-march=native` explicitly for a single-machine build.
- `-flto` is in both `CFLAGS` and `LDFLAGS`, and `-O3` is repeated on the link line so the LTO backend optimizes at the same level.

---

## Verification

Equivalence must be **byte-identical stdout against a binary built from the pre-optimization source**, not against the Rust implementation (a different codebase, whose agreement would not prove the C++ refactor preserved C++ tie-breaking).

```bash
git archive master pasa_cpp | tar -x -C /tmp/ref && make -C /tmp/ref/pasa_cpp
make -C pasa_cpp
for t in tests/*.txt; do
  /tmp/ref/pasa_cpp/pasa "$t" > ref.out
  for th in 1 2 8; do OMP_NUM_THREADS=$th pasa_cpp/pasa "$t" | cmp - ref.out || echo "DIFF $t @ $th"; done
done
```

The input set must cover:

1. **Score ties** — many identical/near-identical spliced alignments over the same span. This is what catches Fscan/Rscan traversal-order regressions. Organic data may not trigger it reliably, so the fixture has to be deliberately tied.
2. **>64 and >128 alignments in one locus** — anything smaller leaves `num_words == 1` and never exercises the multi-word bitset path or its unrolled loop. `sample_data/` alone is almost certainly too small.
3. **Randomized staggered tiling** — overlapping alignments at varied offsets and segment counts. This is the case that actually exposed the Fscan ordering bug; the tie-only fixture did not.
4. **Thread-count sweep** — every input at `OMP_NUM_THREADS` of 1, 2, and 8, to confirm the parallel compatibility scan is order-stable end to end.

All four categories currently pass byte-identical against `master`.
