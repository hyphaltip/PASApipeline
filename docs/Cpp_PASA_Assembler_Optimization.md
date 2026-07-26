# C++ pasa_assembler Optimization Plan

## Goal

Apply algorithmic and micro-optimizations to the C++ `pasa_assembler` (in `pasa_cpp/`) to improve throughput, reduce memory footprint, and eliminate known bugs — without changing the toolchain (staying in C++).

The Rust implementation at `pasa_rust/pasa-assembler/` validates several of these approaches; the C++ code remains the canonical reference for the PASA assembly algorithm.

## Optimization Summary

| # | Optimization | Area | Impact | Status |
|---|-------------|------|--------|--------|
| 1 | Sparse adjacency lists + interval tree | `determine_compatibilities_and_encapsulations` | HIGH — O(n²) matrix → O(n log n + k) | PENDING |
| 2 | Bitset with hardware POPCNT | Lobject containment tracking | HIGH — 64× faster bit ops, 16× less memory | PENDING |
| 3 | Fix memory leak in `free2Darray` | Memory management | HIGH — leaked O(n²) memory per call | PENDING |
| 4 | Const-reference sort comparator | `sort_CDNA_alignments` | MEDIUM — eliminates deep copy per comparison | PENDING |
| 5 | Move semantics in `create_assembly` | Assembly construction | MEDIUM — eliminates O(k²) segment copies | PENDING |
| 6 | `unordered_set` over `map<int,bool>` | `mergeAlignments` splice tracking | MEDIUM — O(1) vs O(log n) | PENDING |
| 7 | OpenMP parallelism | Fscan / Rscan loops | HIGH — near-linear speedup on multi-core | PENDING |
| 8 | Compiler flags: `-march=native -flto` | Build system | LOW — ~10-20% generic CPU gain | PENDING |

---

## 1. Sparse Adjacency Lists + Interval Tree

**Current (`cdna_alignment_assembler.cpp:656-698`)**:
- `bool**` ragged array (n×n) allocated for both `compatibilities` and `encapsulations`
- Each `bool` occupies 1 byte → 2×n² bytes (~200 MB at n=10k)
- Inner loop uses early-break on sorted `lend`, but still O(n²/k) worst-case
- `free2Darray` is a no-op (commented out) — memory leak

**Approach (ports interval tree from Rust implementation)**:
- Build a sorted array of `(lend, alignment_index)` pairs
- For each alignment `i`, binary-search (`std::upper_bound`) for the first candidate with `lend > rend_i`
- Only call `canMerge` on the O(k) truly overlapping candidates
- Store results as `vector<unordered_set<int>>` (sparse) — one set per row, only populated entries consume memory

**Files**: `cdna_alignment_assembler.h`, `cdna_alignment_assembler.cpp`

---

## 2. Bitset + Hardware POPCNT

**Current (`Lobject.h:21`, `Lobject.cpp:34-44`)**:
- `vector<bool>` for `contained_cdna_indices` — packed bits but bit-by-bit iteration
- `num_unique_contained` iterates all n elements one-by-one O(n)

**Approach (ported from `pasa_rust/pasa-assembler/src/lobject.rs`)**:
- Replace `vector<bool>` with `vector<uint64_t>` — each uint64_t holds 64 bits
- `num_unique_contained`: compute `self[w] & ~other[w]` word-wise, sum with `__builtin_popcountll`
- `forwardTrace`/`backTrace`: iterate set bits via `__builtin_ctzll` (count trailing zeros) + `x &= x - 1` (clear lowest set bit)

**Files**: `Lobject.h`, `Lobject.cpp`, `cdna_alignment_assembler.cpp`

---

## 3. Fix Memory Leak in `free2Darray`

**Current (`common_subs.cpp:101-128`)**:
Both `free2Darray(int**)` and `free2Darray(bool**)` have their entire body commented out.

**Fix**: Properly `delete[]` each row, then `delete[]` the row-pointer array.

**Note**: Once the sparse adjacency lists (#1) and bitset (#2) are in place, the `twoDarray`/`free2Darray` functions are no longer called from the assembler, but the fix is applied for any remaining callers.

**File**: `common_subs.cpp`

---

## 4. Pass by Const Reference in sort Comparator

**Current (`cdna_alignment_assembler.cpp:7-19`)**:
```cpp
bool sort_CDNA_alignments(CDNA_alignment a, CDNA_alignment b)
```
Passes entire alignment objects (with segment vectors) by value — deep copy on each of O(n log n) comparisons.

**Fix**: Use `const CDNA_alignment&`.

**File**: `cdna_alignment_assembler.cpp`

---

## 5. Move Semantics in `create_assembly`

**Current (`cdna_alignment_assembler.cpp:769-789`)**:
- `CDNA_alignment assembly = alignments[alignment_index]` — deep copy
- Each `mergeAlignments` returns a new deep-copied alignment
- Results in O(k²) segment copies for an assembly of k alignments

**Fix**: Add explicit move constructor/assignment to `CDNA_alignment` and `Alignment_segment`. Use `assembly = mergeAlignments(std::move(assembly), nextAlignment)` or merge in-place.

**Files**: `cdna_alignment.h`, `cdna_alignment.cpp`, `cdna_alignment_assembler.cpp`

---

## 6. `unordered_set` over `map<int,bool>`

**Current (`cdna_alignment_assembler.cpp:429-430`)**:
```cpp
map<int,bool> leftsplicecoords;
map<int,bool> rightsplicecoords;
```

**Fix**: Use `unordered_set<int>` — O(1) average insert/lookup vs O(log n), and stores only keys (no bool values).

**File**: `cdna_alignment_assembler.cpp`

---

## 7. OpenMP Parallelism

**Current**: All loops are single-threaded.

**Approach**: Parallelize three hot loops:
1. `determine_compatibilities_and_encapsulations` — outer i loop (each i reads alignments[i], writes to compat/encap rows i and j)
2. `do_full_Fscan` — outer i loop (each i reads Lobjects[j] for j < i, writes only to Lobjects[i])
3. `do_full_Rscan` — outer i loop (each i reads Lobjects[j] for j > i, writes only to Lobjects[i])

All three have independent per-iteration work with no cross-iteration writes to shared data.

**Files**: `cdna_alignment_assembler.cpp`, `Makefile` (add `-fopenmp`)

---

## 8. Compiler Flags

**Current (`pasa_cpp/Makefile:22`)**:
```makefile
CFLAGS = -O3
```

**Update**:
```makefile
CFLAGS = -O3 -march=native -flto -funroll-loops
LDFLAGS += -flto
```

- `-march=native` — enables CPU-specific SIMD (SSE4.2, AVX2) for `__builtin_popcountll`, auto-vectorization
- `-flto` — link-time optimization, enables cross-module inlining
- `-funroll-loops` — unrolls hot short-loops
- PGO (profile-guided optimization) available as a follow-up

---

## Verification

After all changes:
```bash
cd pasa_cpp
make clean && make
```

Compare output on `sample_data/` against the baseline `rust_optimize` branch — identical assembly structure expected.
