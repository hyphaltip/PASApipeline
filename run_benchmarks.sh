#!/usr/bin/env bash
#
# PASA Pipeline Performance Benchmark Suite
#
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BIN_DIR="$SCRIPT_DIR/bin"
TMPDIR=$(mktemp -d /tmp/pasa_bench.XXXXXX)
trap "rm -rf $TMPDIR" EXIT

fmt_time() {
    local t=$1
    if (( $(echo "$t > 1" | bc -l) )); then
        printf "%.2fs" "$t"
    elif (( $(echo "$t > 0.001" | bc -l) )); then
        printf "%.1fms" "$(echo "$t * 1000" | bc -l)"
    else
        printf "%.1fus" "$(echo "$t * 1000000" | bc -l)"
    fi
}

fmt_speedup() {
    local s=$1
    if (( $(echo "$s > 10" | bc -l) )); then
        printf "%.0fx" "$s"
    elif (( $(echo "$s > 1.5" | bc -l) )); then
        printf "%.1fx" "$s"
    else
        printf "%.2fx" "$s"
    fi
}

print_header() {
    echo ""
    echo "=================================================================="
    echo "  $1"
    echo "=================================================================="
}

print_result_row() {
    local label=$1 rust_t=$2 cpp_t=$3
    local speedup=$(echo "$cpp_t / $rust_t" | bc -l)
    printf "  %-40s  Rust: %-10s  C++: %-10s  Speedup: %s\n" \
        "$label" "$(fmt_time $rust_t)" "$(fmt_time $cpp_t)" "$(fmt_speedup $speedup)"
}

# ============================================================
#  Generate Test Data
# ============================================================
print_header "Generating Test Data"

NUM_SEQS=500
SEQ_LEN=5000

echo "  Creating FASTA with $NUM_SEQS sequences ($SEQ_LEN bp each)..."
TEST_FASTA="$TMPDIR/test_genome.fa"
perl -e '
srand(42);
my @b = ("A","C","G","T");
for my $i (1..'$NUM_SEQS') {
    print ">seq", sprintf("%04d", $i), "\n";
    my $s = "";
    for (1..'$SEQ_LEN') { $s .= $b[int(rand(4))]; }
    for (my $j = 0; $j < length($s); $j += 60) {
        print substr($s, $j, 60), "\n";
    }
}
' > "$TEST_FASTA"

echo "  Building CDB index (.cidx)..."
bin/cdbfasta -C "$TEST_FASTA" > /dev/null 2>&1

echo "  Building samtools .fai index..."
samtools faidx "$TEST_FASTA" 2>/dev/null

ACC_LIST="$TMPDIR/acc_list.txt"
perl -e '
srand(123);
my @a;
for my $i (1..'$NUM_SEQS') { push @a, sprintf("seq%04d", $i); }
for (my $i = $#a; $i > 0; $i--) {
    my $j = int(rand($i + 1));
    @a[$i, $j] = @a[$j, $i];
}
print join("\n", @a), "\n";
' > "$ACC_LIST"

echo ""
echo "  Test data ready: $NUM_SEQS sequences, $SEQ_LEN bp each"

# ============================================================
#  Benchmark 1: cdbyank_rust vs C++ cdbyank
# ============================================================
print_header "Benchmark 1: Sequence Retrieval (CDB .cidx index)"

CDBYANK_RUST="$BIN_DIR/cdbyank_rust"
CDBYANK_CPP="$BIN_DIR/cdbyank"

echo "  Retrieving $NUM_SEQS sequences (one call per sequence)..."
echo ""

CPP_T=$( { t0=$(date +%s.%N); while IFS= read -r acc; do "$CDBYANK_CPP" -a "$acc" "$TEST_FASTA.cidx" >/dev/null 2>&1; done < "$ACC_LIST"; t1=$(date +%s.%N); echo "$t1 - $t0" | bc -l; } )

RUST_T=$( { t0=$(date +%s.%N); while IFS= read -r acc; do "$CDBYANK_RUST" -a "$acc" "$TEST_FASTA.cidx" >/dev/null 2>&1; done < "$ACC_LIST"; t1=$(date +%s.%N); echo "$t1 - $t0" | bc -l; } )

print_result_row "cdbyank ($NUM_SEQS lookups)" $RUST_T $CPP_T

# ============================================================
#  Benchmark 2: faidx_rust vs cdbyank
# ============================================================
print_header "Benchmark 2: Sequence Retrieval (samtools .fai index)"

FAIDX_RUST="$BIN_DIR/faidx_rust"

echo "  Retrieving $NUM_SEQS sequences via faidx_rust..."
echo ""

FAIDX_T=$( { t0=$(date +%s.%N); while IFS= read -r acc; do "$FAIDX_RUST" "$TEST_FASTA" "$acc" >/dev/null 2>&1; done < "$ACC_LIST"; t1=$(date +%s.%N); echo "$t1 - $t0" | bc -l; } )

print_result_row "faidx_rust vs C++ cdbyank" $FAIDX_T $CPP_T
print_result_row "faidx_rust vs cdbyank_rust" $FAIDX_T $RUST_T

# ============================================================
#  Benchmark 3: pasa_rust vs C++ pasa
# ============================================================
print_header "Benchmark 3: PASA Alignment Assembly"

PASA_RUST="$BIN_DIR/pasa_rust"
PASA_CPP="$BIN_DIR/pasa"

PASA_INPUT="$TMPDIR/pasa_input.txt"
perl -e '
srand(99);
my @lines;
for my $i (1..500) {
    my $acc = sprintf("transcript%04d", $i);
    my $orient = (rand() > 0.5) ? "+" : "-";
    my $num_segs = int(rand(5)) + 1;
    my $pos = int(rand(100000)) + 1;
    my @segs;
    for my $s (1..$num_segs) {
        my $seg_len = int(rand(400)) + 100;
        push @segs, "$pos-" . ($pos + $seg_len);
        $pos += $seg_len + int(rand(450)) + 50;
    }
    push @lines, "$acc,$orient," . join(",", @segs);
}
print join("\n", @lines), "\n";
' > "$PASA_INPUT"

NUM_ALIGNMENTS=$(wc -l < "$PASA_INPUT")
echo "  Input: $NUM_ALIGNMENTS alignments"
echo ""

CPP_T=$( { t0=$(date +%s.%N); "$PASA_CPP" "$PASA_INPUT" >/dev/null 2>&1; t1=$(date +%s.%N); echo "$t1 - $t0" | bc -l; } )
RUST_T=$( { t0=$(date +%s.%N); "$PASA_RUST" "$PASA_INPUT" >/dev/null 2>&1; t1=$(date +%s.%N); echo "$t1 - $t0" | bc -l; } )

print_result_row "pasa ($NUM_ALIGNMENTS alignments)" $RUST_T $CPP_T

# ============================================================
#  Benchmark 4: slclust_rust vs C++ slclust
# ============================================================
print_header "Benchmark 4: Single-Linkage Clustering (slclust)"

SLCLUST_RUST="$BIN_DIR/slclust_rust"
SLCLUST_CPP="$BIN_DIR/slclust"

SLCLUST_INPUT="$TMPDIR/slclust_input.txt"
perl -e '
srand(77);
my $num_nodes = 2000;
for my $i (1..5000) {
    my $a = int(rand($num_nodes)) + 1;
    my $b = int(rand($num_nodes)) + 1;
    print "$a $b\n" if ($a != $b);
}
' > "$SLCLUST_INPUT"

NUM_PAIRS=$(wc -l < "$SLCLUST_INPUT")
echo "  Input: $NUM_PAIRS pairs, up to 2000 nodes"
echo ""

CPP_T=$( { t0=$(date +%s.%N); ulimit -s unlimited 2>/dev/null; "$SLCLUST_CPP" < "$SLCLUST_INPUT" >/dev/null 2>&1; t1=$(date +%s.%N); echo "$t1 - $t0" | bc -l; } )
RUST_T=$( { t0=$(date +%s.%N); "$SLCLUST_RUST" < "$SLCLUST_INPUT" >/dev/null 2>&1; t1=$(date +%s.%N); echo "$t1 - $t0" | bc -l; } )

print_result_row "slclust ($NUM_PAIRS pairs)" $RUST_T $CPP_T

# ============================================================
#  Summary
# ============================================================
print_header "Summary"

echo "  Test data: $NUM_SEQS seqs × ${SEQ_LEN}bp, $NUM_ALIGNMENTS alignments, $NUM_PAIRS pairs"
echo ""
echo "  Rust Tier optimizations:"
echo "    Tier 1: Interval tree O(n log n + k) vs O(n²)  [pasa-assembler]"
echo "    Tier 2: Vec<u64> bitset + hardware popcount      [lobject]"
echo "    Tier 3: Bulk I/O, in-memory djb2 hash            [cdbyank_rust, faidx_rust]"
echo "    Tier 4: HashSet O(1) dedup, iterative DFS        [slclust_rust]"
echo ""
echo "  Perl-level optimizations:"
echo "    Gene_obj::get_exons() caching (~2.77x)"
echo "    Fasta_retriever file handle reuse + compression"
echo "    GFF3_utils precompiled regex patterns"
echo "    PSL_parser O(n²) shift→index + get_per_id caching"
echo ""
echo "  Pipeline integration:"
echo "    18 scripts: cdbyank_linear() → get_seq() with faidx_rust"
echo "    PASA_alignment_assembler.pm: prefers pasa_rust"
echo "    SingleLinkageClusterer.pm: prefers slclust_rust"
