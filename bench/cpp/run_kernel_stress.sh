#!/usr/bin/env bash
# run_kernel_stress.sh - C++ kernel exactness batteries, compiled from the include/ headers without any Rust.
#
# It runs
#   * dbg_stress, the dbruijn counter, graph and unitigs against brute-force truth
#       and the independent bidirected oracle of bench/cpp/dbg_oracle.hpp
#   * minimizer_stress, an exactness battery, and reader_stress, the corner cases of the
#       FASTA, FASTQ and gzip readers, where a truncated gzip must throw
#   * the mphf, spss, delta, dyn_spss and seed_index stresses, and wfa_affine_stress, optik's
#       own gap-affine WFA against an independent Gotoh DP
#   * kernels_stress, which checks the constant-time reverse complement and the bit-parallel
#       edit distance in scalar and AVX2 lanes against a naive DP, the k-order-BWT FM index
#       against the brute-force window set, and SPSS distinct counts and strand-independent ids
#       on both the build and load paths against brute force
#
# Usage: bench/cpp/run_kernel_stress.sh [g++]
#
# It needs a C++17 compiler with OpenMP, by default 'g++' on PATH.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# Nix's compiler wrapper drops -march=native unless NIX_ENFORCE_NO_NATIVE is 0, and these builds are meant for this CPU.
# The variable has no effect outside Nix.
export NIX_ENFORCE_NO_NATIVE=0
ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CXX="${1:-g++}"
if ! command -v "$CXX" >/dev/null 2>&1; then
    echo "  skipping C++ kernel stress: compiler '$CXX' not found (pass the g++ path as the first argument)"
    exit 0
fi
TMP="${TMPDIR:-/tmp}/optik_kernel_stress_$$"
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT

FAIL=0
run_ok() { # name status
    if [ "$2" -eq 0 ]; then echo "  => $1: PASS";
    else echo "  => $1: FAIL"; FAIL=1; fi
}
# run_test NAME LOG CMD... runs CMD with stdout and stderr going to LOG, without letting `set -e` abort the script, then reports PASS or FAIL and the log's last line.
run_test() {
    local name="$1" log="$2"; shift 2
    local st=0
    "$@" > "$log" 2>&1 || st=$?
    run_ok "$name" "$st"
    tail -1 "$log"
}
INC=(-I "$ROOT/include" -I "$ROOT/src" -I "$SCRIPT_DIR")
KC="$ROOT/src/kmer_count.cpp"

echo "=============================================="
echo "  optik C++ kernel exactness stress"
echo "=============================================="

# dBG unitig batteries against the independent bidirected oracle (dbg_oracle.hpp).
# dbg_stress covers counter records and masks, mask queries, the unitig multiset, partition, coverage and thread-count invariance, on dense K=4..7 inputs and genomes or reads at K=15..31.
for t in dbg_stress; do
    if [ -f "$SCRIPT_DIR/$t.cpp" ]; then
        echo "[build] $t"
        "$CXX" -std=c++17 -O2 -march=native -fopenmp "${INC[@]}" \
            "$SCRIPT_DIR/$t.cpp" "$KC" -o "$TMP/$t" -lz
        run_test "$t" "$TMP/$t.log" "$TMP/$t"
    fi
done

# Index exactness batteries
for t in minimizer_stress reader_stress; do
    if [ -f "$SCRIPT_DIR/$t.cpp" ]; then
        echo "[build] $t"
        "$CXX" -std=c++17 -O2 -march=native -fopenmp "${INC[@]}" \
            "$SCRIPT_DIR/$t.cpp" -o "$TMP/$t" -lz
        run_test "$t" "$TMP/$t.log" "$TMP/$t"
    fi
done

# SPSS core, MPHF, spss-update deltas, seed index and optik's own gap-affine WFA
for t in mphf_stress spss_stress delta_stress dyn_spss_stress wfa_affine_stress seed_index_stress kernels_stress; do
    if [ -f "$SCRIPT_DIR/$t.cpp" ]; then
        echo "[build] $t"
        "$CXX" -std=c++17 -O2 -march=native -fopenmp -pthread "${INC[@]}" \
            "$SCRIPT_DIR/$t.cpp" "$KC" -o "$TMP/$t" -lz
        run_test "$t" "$TMP/$t.log" "$TMP/$t"
    fi
done

echo "=============================================="
if [ "$FAIL" -eq 0 ]; then echo "  ALL KERNEL STRESS PASSED";
else echo "  SOME KERNEL STRESS FAILED"; fi
echo "=============================================="
exit "$FAIL"