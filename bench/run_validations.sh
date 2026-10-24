#!/usr/bin/env bash
# run_validations.sh - end-to-end correctness validation of optik
#
# Runs the official exact-equality suites for query, mapping and the de Bruijn graph, the C++ kernel exactness stress tests, the WFA2-lib waveform oracle bench, the long-read suite when its data is present, and the further checks listed per step below.
# Any failure makes the exit status non-zero.
#
# Usage
#     bench/run_validations.sh [dataset] [k]
#
# The defaults are ecoli_k31 and k = 31.
#
# Environment
#     OPTIK_BIN        path to the optik binary  [target/release/optik]
#     MINIMAP2_BIN     path to minimap2          [auto-detect]
#     BIFROST_BIN      path to Bifrost           [auto-detect]
#     CUTTLEFISH_BIN   path to Cuttlefish        [auto-detect]
#     SSHASH_BIN       path to SSHash            [auto-detect]
#     KMC_BIN          path to KMC               [auto-detect]
#     SAMTOOLS_BIN     path to samtools          [auto-detect]
#     NIHSS_REF        stitched-unitig genome for sshash --check
#     NUM_READS        reads for mapping         [10000]
#     THREADS          thread count              [usable CPUs, nproc]
#     SALMONELLA_REF   Salmonella reference such as NC_003197.1, which enables step 7, the long-read suite of 50 kb and more
#     REQUIRE_ALL      1 makes any skipped step or comparison (missing tool or data) a failure.
#                      With the default 0, skips are listed in the summary, which then never claims a full pass.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# Nix's compiler wrapper drops -march=native unless NIX_ENFORCE_NO_NATIVE is 0, and these builds are meant for this CPU.
# The variable has no effect outside Nix.
export NIX_ENFORCE_NO_NATIVE=0

# Defaults
DATASET="${1:-ecoli_k31}"
K="${2:-31}"
OPTIK_BIN="${OPTIK_BIN:-$PROJECT_ROOT/target/release/optik}"
MINIMAP2_BIN="${MINIMAP2_BIN:-}"
BIFROST_BIN="${BIFROST_BIN:-}"
CUTTLEFISH_BIN="${CUTTLEFISH_BIN:-}"
SSHASH_BIN="${SSHASH_BIN:-}"
KMC_BIN="${KMC_BIN:-}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-}"
NUM_READS="${NUM_READS:-10000}"
THREADS="${THREADS:-$(nproc)}"

# Tools are found in temp/ or in parent directories.
TOOLS_DIR="${TOOL_BENCH_TOOLS_DIR:-$PROJECT_ROOT/temp}"
if [ -d "$TOOLS_DIR" ]; then
    [ -z "$MINIMAP2_BIN" ] && MINIMAP2_BIN="$TOOLS_DIR/minimap2/minimap2"
    [ -z "$BIFROST_BIN" ] && BIFROST_BIN="$TOOLS_DIR/bifrost/build/src/Bifrost"
    [ -z "$CUTTLEFISH_BIN" ] && CUTTLEFISH_BIN="$TOOLS_DIR/cuttlefish/target/release/cuttlefish"
    [ -z "$SSHASH_BIN" ] && SSHASH_BIN="$TOOLS_DIR/sshash/build/sshash"
    [ -z "$KMC_BIN" ] && KMC_BIN="$TOOLS_DIR/KMC/bin/kmc"
    [ -z "$SAMTOOLS_BIN" ] && SAMTOOLS_BIN="${SAMTOOLS_BIN:-$(command -v samtools || true)}"
fi

# Wavefront mapper, MashMap with WFA, for mapping speed and sensitivity parity
WFMASH_BIN="${WFMASH_BIN:-}"
[ -z "$WFMASH_BIN" ] && [ -x "$TOOLS_DIR/wfmash/build/bin/wfmash" ] && WFMASH_BIN="$TOOLS_DIR/wfmash/build/bin/wfmash"
[ -z "$WFMASH_BIN" ] && command -v wfmash >/dev/null 2>&1 && WFMASH_BIN="$(command -v wfmash)"

# The older layout keeps tools in the parent directory.
[ -z "$MINIMAP2_BIN" ] && [ -x "$PROJECT_ROOT/../minimap2/minimap2" ] && MINIMAP2_BIN="$PROJECT_ROOT/../minimap2/minimap2"
[ -z "$BIFROST_BIN" ] && [ -x "$PROJECT_ROOT/../bifrost/build/src/Bifrost" ] && BIFROST_BIN="$PROJECT_ROOT/../bifrost/build/src/Bifrost"
[ -z "$CUTTLEFISH_BIN" ] && [ -x "$PROJECT_ROOT/../cuttlefish/target/release/cuttlefish" ] && CUTTLEFISH_BIN="$PROJECT_ROOT/../cuttlefish/target/release/cuttlefish"
[ -z "$SSHASH_BIN" ] && [ -x "$PROJECT_ROOT/../sshash/build/sshash" ] && SSHASH_BIN="$PROJECT_ROOT/../sshash/build/sshash"
[ -z "$KMC_BIN" ] && [ -x "$PROJECT_ROOT/../KMC/bin/kmc" ] && KMC_BIN="$PROJECT_ROOT/../KMC/bin/kmc"
[ -z "$SAMTOOLS_BIN" ] && [ -x "$PROJECT_ROOT/../samtools/samtools" ] && SAMTOOLS_BIN="$PROJECT_ROOT/../samtools/samtools"

# Dataset details from the bench configuration
get_dataset_ref() {
    python3 -c "
import sys; sys.path.insert(0, '$SCRIPT_DIR')
from config import DATASETS
ds = DATASETS.get('$DATASET')
if ds:
    print(ds.get('unitigs', ''))
else:
    print('')
" 2>/dev/null
}

get_dataset_queries() {
    python3 -c "
import sys; sys.path.insert(0, '$SCRIPT_DIR')
from config import DATASETS
ds = DATASETS.get('$DATASET')
if ds:
    print(ds.get('queries', ''))
else:
    print('')
" 2>/dev/null
}

get_dataset_dbg_input() {
    python3 -c "
import sys; sys.path.insert(0, '$SCRIPT_DIR')
from config import DATASETS
ds = DATASETS.get('$DATASET')
if ds:
    print(ds.get('dbg_input', ds.get('queries', '')))
else:
    print('')
" 2>/dev/null
}

get_dataset_k() {
    python3 -c "
import sys; sys.path.insert(0, '$SCRIPT_DIR')
from config import DATASETS
ds = DATASETS.get('$DATASET')
if ds:
    print(ds.get('k', ''))
else:
    print('')
" 2>/dev/null
}

# The dataset's k applies unless k was given.
if [ "$K" = "31" ]; then
    DS_K=$(get_dataset_k)
    [ -n "$DS_K" ] && K="$DS_K"
fi

REF=$(get_dataset_ref)
QUERIES=$(get_dataset_queries)
# The dataset's files come with the SSHash checkout, so a missing one stops the suite here with the command that provides it.
for f in "$REF" "$QUERIES"; do
    if [ -n "$f" ] && [ ! -e "$f" ]; then
        echo "error: $DATASET input $f is missing (bench/build_tools.sh sshash fetches it)" >&2
        exit 1
    fi
done

# Look for a Salmonella reference when SALMONELLA_REF is unset.
if [ -z "${SALMONELLA_REF:-}" ]; then
    for candidate in \
        "$TOOLS_DIR/salmonella_LT2.fa" \
        "$PROJECT_ROOT/temp/salmonella_LT2.fa" \
        "$TOOLS_DIR/salmonella_enterica/NC_003197.1.fa" \
        "$PROJECT_ROOT/../salmonella/NC_003197.1.fa"; do
        if [ -f "$candidate" ]; then
            SALMONELLA_REF="$candidate"
            break
        fi
    done
fi

echo "=============================================="
echo "  optik correctness validation suite"
echo "=============================================="
echo "  dataset:   $DATASET"
echo "  reference: $REF"
echo "  k:         $K"
echo "  optik:     $OPTIK_BIN"
echo "  minimap2:  ${MINIMAP2_BIN:-not found}"
echo "  bifrost:   ${BIFROST_BIN:-not found}"
echo "  cuttlefish: ${CUTTLEFISH_BIN:-not found}"
echo "  sshash:    ${SSHASH_BIN:-not found}"
echo "  kmc:       ${KMC_BIN:-not found}"
echo "  samtools:  ${SAMTOOLS_BIN:-not found}"
echo "  wfmash:    ${WFMASH_BIN:-not found}"
echo "  salmonella: ${SALMONELLA_REF:-not set}"
echo "  reads:     $NUM_READS"
echo "  threads:   $THREADS"
echo "=============================================="

FAIL=0
# Wall time of each step, printed at the end to show where a slow run spent its time.
STEP_NAME=""; STEP_T0=$SECONDS; STEP_TIMES=()
step_mark() {
    [ -n "$STEP_NAME" ] && STEP_TIMES+=("$(printf '%6ds  %s' $((SECONDS - STEP_T0)) "$STEP_NAME")")
    STEP_NAME="$1"; STEP_T0=$SECONDS
}
SKIPPED=()
# skip WHAT WHY records a step or comparison that did not run.
skip() {
    SKIPPED+=("$1 ($2)")
    echo "  SKIPPED: $1 -- $2"
}
# Bookkeeping of required tools for the comparisons the per-step scripts would otherwise drop without a word
[ -x "$SSHASH_BIN" ] || skip "SSHash comparisons (steps 1, 4)" "sshash not found"
[ -x "$MINIMAP2_BIN" ] || skip "minimap2 comparisons (steps 2, 4, 7, 9)" "minimap2 not found"
{ [ -n "$WFMASH_BIN" ] && [ -x "$WFMASH_BIN" ]; } || skip "wfmash comparisons (steps 2, 4, 7B)" "wfmash not found"
[ -x "$BIFROST_BIN" ] || skip "Bifrost comparisons (steps 3, 4)" "Bifrost not found"
[ -x "$CUTTLEFISH_BIN" ] || skip "Cuttlefish comparisons (steps 3, 4)" "Cuttlefish not found"
[ -x "$SAMTOOLS_BIN" ] || skip "samtools round-trips (step 4)" "samtools not found"

# Build optik if needed.
if [ ! -x "$OPTIK_BIN" ]; then
    echo ""
    echo "[setup] Building optik..."
    cargo build --release --manifest-path "$PROJECT_ROOT/Cargo.toml"
fi

# Decompress the reference if it is gzipped.
REF_UNCOMPRESSED="$REF"
CLEANUP_FILES=()
if [[ "$REF" == *.gz ]]; then
    REF_UNCOMPRESSED="/tmp/ref_uncomp_$$.fa"
    if [ ! -f "$REF_UNCOMPRESSED" ]; then
        gunzip -c "$REF" > "$REF_UNCOMPRESSED"
    fi
    CLEANUP_FILES+=("$REF_UNCOMPRESSED")
fi
trap "rm -f ${CLEANUP_FILES[*]:-}" EXIT

# Step 1.
# Query validation
echo ""
echo "----------------------------------------------"
echo "  [1/9] Query validation"
step_mark "1 Query validation"
echo "----------------------------------------------"

QUERY_ARGS=(-r "$REF_UNCOMPRESSED" -k "$K" --optik "$OPTIK_BIN" -s 42)
[ -x "$SSHASH_BIN" ] && QUERY_ARGS+=(--sshash "$SSHASH_BIN")

if python3 "$SCRIPT_DIR/validate_query.py" "${QUERY_ARGS[@]}"; then
    echo "  => QUERY VALIDATION: PASS"
else
    echo "  => QUERY VALIDATION: FAIL"
    FAIL=1
fi

# Step 2.
# Mapping validation
echo ""
echo "----------------------------------------------"
echo "  [2/9] Mapping validation"
step_mark "2 Mapping validation"
echo "----------------------------------------------"

MAP_ARGS=(-r "$REF_UNCOMPRESSED" -k "$K" --optik "$OPTIK_BIN"
          -n "$NUM_READS" -s 42)
[ -x "$MINIMAP2_BIN" ] && MAP_ARGS+=(--minimap2 "$MINIMAP2_BIN")
[ -n "$WFMASH_BIN" ] && [ -x "$WFMASH_BIN" ] && MAP_ARGS+=(--wfmash "$WFMASH_BIN")
# When wgatools is available (temp/wgatools/target/release/wgatools or a system binary), the independent PAF round-trip checks are included, validate for position and CIGAR consistency and pafcov for reference depth.
WGATOOLS_BIN="${WGATOOLS_BIN:-}"
[ -z "$WGATOOLS_BIN" ] && [ -x "$PROJECT_ROOT/temp/wgatools/target/release/wgatools" ] && \
    WGATOOLS_BIN="$PROJECT_ROOT/temp/wgatools/target/release/wgatools"
[ -z "$WGATOOLS_BIN" ] && command -v wgatools >/dev/null 2>&1 && WGATOOLS_BIN="$(command -v wgatools)"
[ -n "$WGATOOLS_BIN" ] && [ -x "$WGATOOLS_BIN" ] && MAP_ARGS+=(--wgatools "$WGATOOLS_BIN")

if python3 "$SCRIPT_DIR/validate_mapping.py" "${MAP_ARGS[@]}"; then
    echo "  => MAPPING VALIDATION: PASS"
else
    echo "  => MAPPING VALIDATION: FAIL"
    FAIL=1
fi
# MD and cs tags are checked with samtools calmd, by replay against reference and SEQ, and with minimap2.
TAGS_REF=$(get_dataset_dbg_input)
if [ -n "$TAGS_REF" ] && [ -n "$SAMTOOLS_BIN" ] && [ -x "$SAMTOOLS_BIN" ]; then
    step_mark "2 MD/cs tags"
    TAGS_DIR="/tmp/optik_tags_$$"
    TAGS_IN="$TAGS_REF"
    if [[ "$TAGS_REF" == *.gz ]]; then TAGS_IN="$TAGS_DIR.ref.fa"; gunzip -c "$TAGS_REF" > "$TAGS_IN"; fi
    if python3 "$SCRIPT_DIR/validate_tags.py" --ref "$TAGS_IN" --workdir "$TAGS_DIR" \
            --samtools "$SAMTOOLS_BIN" --threads "$THREADS" > "$TAGS_DIR.log" 2>&1; then
        echo "  => MD/CS TAG VALIDATION: PASS"
    else
        echo "  => MD/CS TAG VALIDATION: FAIL"
        grep -E "\[FAIL" "$TAGS_DIR.log" | cut -c1-200 | head -5
        FAIL=1
    fi
    rm -rf "$TAGS_DIR" "$TAGS_DIR.ref.fa"
else
    skip "MD/cs tag validation" "needs samtools and a dataset genome"
fi

# Step 3.
# De Bruijn graph validation
echo ""
echo "----------------------------------------------"
echo "  [3/9] De Bruijn graph validation"
step_mark "3 De Bruijn graph validation"
echo "----------------------------------------------"

# -c 1 runs on the dataset's dbg_input, the genome its unitigs came from, as the benchmark suite does.
# -c 2 runs on reads simulated from that genome below.
DBG_INPUT=$(get_dataset_dbg_input)
# Decompress if needed.
if [[ "$DBG_INPUT" == *.gz ]]; then
    DBG_INPUT_UNCOMPRESSED="/tmp/dbg_input_$$.fa"
    if [ ! -f "$DBG_INPUT_UNCOMPRESSED" ]; then
        gunzip -c "$DBG_INPUT" > "$DBG_INPUT_UNCOMPRESSED"
    fi
    CLEANUP_FILES+=("$DBG_INPUT_UNCOMPRESSED")
    trap "rm -f ${CLEANUP_FILES[*]}" EXIT
else
    DBG_INPUT_UNCOMPRESSED="$DBG_INPUT"
fi

DBG_ARGS=(-i "$DBG_INPUT_UNCOMPRESSED" -k "$K"
          --optik "$OPTIK_BIN" --min-count 1 -t "$THREADS")
[ -x "$BIFROST_BIN" ] && DBG_ARGS+=(--bifrost "$BIFROST_BIN")
[ -x "$CUTTLEFISH_BIN" ] && DBG_ARGS+=(--cuttlefish "$CUTTLEFISH_BIN")

if [ "$K" -gt 31 ]; then
    skip "DBG validation" "k=$K exceeds dbruijn's k <= 31"
else
    step_mark "3 dbg -c 1"
    if python3 "$SCRIPT_DIR/validate_dbg.py" "${DBG_ARGS[@]}"; then
        echo "  => DBG VALIDATION (-c 1): PASS"
    else
        echo "  => DBG VALIDATION (-c 1): FAIL"
        FAIL=1
    fi
    # -c 2 is the abundance cutoff users run on real reads with errors, and the solid set, partition, links, tags and node table must stay exact there.
    # The input is 20x of 150 bp reads with 1% substitutions and 0.1% indels drawn from the first 1 Mbp of the same genome, so -c 2 keeps the genome's k-mers and drops most error k-mers (a genome at -c 2 would keep only its repeats).
    # The slice bounds the work of the pure-Python oracle, which took about 8 min on the whole 4.9 Mbp genome at 20x.
    # Coverage, error model and every check are the same.
    step_mark "3 dbg -c 2 (incl. read simulation)"
    DBG_SLICE="/tmp/dbg_slice_$$.fa"
    DBG_READS="/tmp/dbg_reads_$$.fa"
    CLEANUP_FILES+=("$DBG_SLICE" "$DBG_READS")
    trap "rm -f ${CLEANUP_FILES[*]}" EXIT
    python3 - "$DBG_INPUT_UNCOMPRESSED" "$DBG_SLICE" <<'PYEOF'
import sys
src, dst, want = sys.argv[1], sys.argv[2], 1_000_000
out, n = [], 0
with open(src) as fh:
    for line in fh:
        if line.startswith(">"):
            if out:
                break
            continue
        line = line.strip()
        out.append(line[:want - n]); n += len(out[-1])
        if n >= want:
            break
with open(dst, "w") as fh:
    fh.write(">slice\n" + "".join(out) + "\n")
PYEOF
    python3 "$SCRIPT_DIR/gen_reads.py" -r "$DBG_SLICE" -c 20 -l 150 \
        -e 0.01 --ins 0.001 --del 0.001 -s 42 -o "$DBG_READS" >/dev/null
    DBG_ARGS_C2=("${DBG_ARGS[@]}")
    for i in "${!DBG_ARGS_C2[@]}"; do
        [ "${DBG_ARGS_C2[$i]}" = "--min-count" ] && DBG_ARGS_C2[$((i + 1))]=2
        [ "${DBG_ARGS_C2[$i]}" = "-i" ] && DBG_ARGS_C2[$((i + 1))]="$DBG_READS"
    done
    if python3 "$SCRIPT_DIR/validate_dbg.py" "${DBG_ARGS_C2[@]}"; then
        echo "  => DBG VALIDATION (-c 2): PASS"
    else
        echo "  => DBG VALIDATION (-c 2): FAIL"
        FAIL=1
    fi
fi

# Step 4.
# Official exact-equality validation
echo ""
echo "----------------------------------------------"
echo "  [4/9] Official exact-equality validation"
step_mark "4 Official exact-equality validation"
echo "----------------------------------------------"

# optik's k-mer code is limited to k <= 32 (MAX_K).
# optik cannot index datasets with a larger k, so for them the strict official suite is skipped with a warning rather than failed.
OFFICIAL_PASS=0
if [ "$K" -gt 32 ]; then
    skip "official suite" "k=$K exceeds optik MAX_K=32"
else
    OFF_WORKDIR="${OFF_WORKDIR:-/tmp/optik_official_$$}"
    rm -rf "$OFF_WORKDIR"
    mkdir -p "$OFF_WORKDIR"

    # The query suite checks exact k-mer membership against SSHash and the KMC oracle on real reads.
    if [ -z "$KMC_BIN" ] || [ ! -x "$KMC_BIN" ]; then
        skip "official query/dbg suites" "KMC not found (set KMC_BIN or run build_tools.sh)"
    else
    OFF_ARGS=(-r "$REF_UNCOMPRESSED" -k "$K" --queries "$QUERIES" -t "$THREADS"
              --optik "$OPTIK_BIN" --sshash "$SSHASH_BIN" --kmc "$KMC_BIN"
              --sshash-m 15 -s 42 --workdir "$OFF_WORKDIR/query")
    step_mark "4 official query"
    if python3 "$SCRIPT_DIR/validate_official.py" query "${OFF_ARGS[@]}" \
            >"$OFF_WORKDIR/query.json" 2>"$OFF_WORKDIR/query.err"; then
        echo "  => OFFICIAL QUERY VALIDATION: PASS"
    else
        echo "  => OFFICIAL QUERY VALIDATION: FAIL"
        grep "FAILED:" "$OFF_WORKDIR/query.err" | cut -c1-200 | head -5
        OFFICIAL_PASS=1
    fi
    fi

    # The map suite checks exact positions on error-free reads, NM and CIGAR, sensitivity scored against the truth beside minimap2 on the same 10%-error reads, and samtools round trips.
    # Each mapper runs with its own documented defaults (optik map, minimap2 -ax sr or map-ont), and no seed parameters are forced.
    MAP2_ARGS=(-r "$REF_UNCOMPRESSED" -k "$K" -n "$NUM_READS" -t "$THREADS"
               --read-len 150 --err-rate 0.10
               --optik "$OPTIK_BIN" -s 42 --workdir "$OFF_WORKDIR/map")
    [ -x "$MINIMAP2_BIN" ] && MAP2_ARGS+=(--minimap2 "$MINIMAP2_BIN")
    [ -n "$WFMASH_BIN" ] && [ -x "$WFMASH_BIN" ] && MAP2_ARGS+=(--wfmash "$WFMASH_BIN")
    [ -x "$SAMTOOLS_BIN" ] && MAP2_ARGS+=(--samtools "$SAMTOOLS_BIN")
    step_mark "4 official map"
    if python3 "$SCRIPT_DIR/validate_official.py" map "${MAP2_ARGS[@]}" \
            >"$OFF_WORKDIR/map.json" 2>"$OFF_WORKDIR/map.err"; then
        echo "  => OFFICIAL MAPPING VALIDATION: PASS"
    else
        echo "  => OFFICIAL MAPPING VALIDATION: FAIL"
        grep "FAILED:" "$OFF_WORKDIR/map.err" | cut -c1-200 | head -5
        OFFICIAL_PASS=1
    fi

    # The dbg suite requires the unitig set to equal the KMC oracle with no branch crossings.
    if [ -x "$KMC_BIN" ]; then
        DBG2_ARGS=(-i "$QUERIES" -k "$K" --min-count 1 -t "$THREADS"
                  --optik "$OPTIK_BIN" --kmc "$KMC_BIN" -s 42
                  --workdir "$OFF_WORKDIR/dbg")
        [ -x "$BIFROST_BIN" ] && DBG2_ARGS+=(--bifrost "$BIFROST_BIN")
        [ -x "$CUTTLEFISH_BIN" ] && DBG2_ARGS+=(--cuttlefish "$CUTTLEFISH_BIN")
        step_mark "4 official dbg"
        if python3 "$SCRIPT_DIR/validate_official.py" dbg "${DBG2_ARGS[@]}" \
                >"$OFF_WORKDIR/dbg.json" 2>"$OFF_WORKDIR/dbg.err"; then
            echo "  => OFFICIAL DBG VALIDATION: PASS"
        else
            echo "  => OFFICIAL DBG VALIDATION: FAIL"
            grep "FAILED:" "$OFF_WORKDIR/dbg.err" | cut -c1-200 | head -5
            OFFICIAL_PASS=1
        fi
    fi
fi
[ "$OFFICIAL_PASS" -ne 0 ] && FAIL=1

# Step 5.
# C++ kernel exactness
echo ""
echo "----------------------------------------------"
echo "  [5/9] C++ kernel exactness stress"
step_mark "5 C++ kernel exactness stress"
echo "----------------------------------------------"

# Exactness batteries for the C++ kernels (dBG oracle, membership indexes, MPHF, SPSS, dynamic deltas, seed index, optik's own WFA) plus a ThreadSanitizer run.
# They compile the include/ headers directly and need no optik binary.
if [ -x "$SCRIPT_DIR/cpp/run_kernel_stress.sh" ]; then
    if bash "$SCRIPT_DIR/cpp/run_kernel_stress.sh"; then
        echo "  => C++ KERNEL STRESS: PASS"
    else
        echo "  => C++ KERNEL STRESS: FAIL"
        FAIL=1
    fi
else
    skip "C++ kernel stress" "cpp/run_kernel_stress.sh not found"
fi

# Step 6.
# Waveform and WFA2 exactness
echo ""
echo "----------------------------------------------"
echo "  [6/9] Waveform/WFA2 exactness bench"
step_mark "6 Waveform/WFA2 exactness bench"
echo "----------------------------------------------"
echo "  (independently verifies optik's own WFA against WFA2-lib;"
echo "   skipped unless temp/WFA2-lib is built)"

# The oracle bench runs optik's WFA in edit mode, as `map` calls it, both uncapped and capped.
# Its results must equal WFA2-lib ends-free, a naive DP and the constructed edit distance over the 150 bp / 10% regime, read lengths from 30 bp to 4 kb and indel-rich consensus strata.
# The battery is deterministic, and it needs WFA2-lib (build_tools.sh wfa2) and g++ to compile the oracle.
if [ -f "$PROJECT_ROOT/temp/WFA2-lib/lib/libwfa.a" ] && \
   [ -f "$SCRIPT_DIR/waveform_bench.py" ]; then
    if python3 "$SCRIPT_DIR/waveform_bench.py" \
            --workdir "$PROJECT_ROOT/temp/waveform"; then
        echo "  => WAVEFORM ORACLE: PASS"
    else
        echo "  => WAVEFORM ORACLE: FAIL"
        FAIL=1
    fi
else
    skip "waveform/WFA2 oracle (step 6)" "WFA2-lib not built (run bench/build_tools.sh wfa2)"
fi
# Gap-affine score parity of optik's WFA with WFA2-lib's classic WFA and its BiWFA mode (bidirectional, O(s) memory) on every pair, with speed printed as well.
if [ -f "$PROJECT_ROOT/temp/WFA2-lib/lib/libwfa.a" ]; then
    WFAB="$(mktemp /tmp/optik_wfab_XXXXXX)"
    if g++ -std=c++17 -O3 -march=native -I "$PROJECT_ROOT/include" -I "$PROJECT_ROOT/temp/WFA2-lib" \
            "$SCRIPT_DIR/cpp/wfa_affine_bench.cpp" "$PROJECT_ROOT/temp/WFA2-lib/bindings/cpp/WFAligner.cpp" \
            "$PROJECT_ROOT/temp/WFA2-lib/lib/libwfa.a" -fopenmp -o "$WFAB" && "$WFAB"; then
        echo "  => WFA2 / BiWFA SCORE PARITY: PASS"
    else
        echo "  => WFA2 / BiWFA SCORE PARITY: FAIL"
        FAIL=1
    fi
    rm -f "$WFAB"
fi

# Step 7.
# Salmonella long-read suite on reads of 50 kb and more
echo ""
echo "----------------------------------------------"
echo "  [7/9] Salmonella >=50 kb long-read suite"
step_mark "7 Salmonella >=50 kb long-read suite"
echo "----------------------------------------------"

# Runs only when SALMONELLA_REF points at a Salmonella reference (such as NC_003197.1, 4.86 Mb).
# Reads of 500 bp or more take optik's chaining path, which writes one primary (with possibly soft-clipped flanks) plus supplementary records carrying SA tags.
#
#   Suite A maps 30 error-free 50 kb reads from both strands and checks them with
#           validate_mapping.py.  The primary must be anchored at its origin
#           (pos - 1 == gt + leading_S) on the truth strand with NM 0, parity
#           with minimap2 -ax map-ont is scored against the truth, and the
#           wgatools PAF round trip runs as well.
#   Suite B maps 30 50 kb reads at 0.1% error with optik and wfmash, both on
#           their defaults, and scores them against the simulation truth.  optik
#           must place every read and at least as many as wfmash.
if [ -n "${SALMONELLA_REF:-}" ] && [ -f "$SALMONELLA_REF" ]; then
    if [ ! -s "$SALMONELLA_REF.fai" ]; then
        if command -v samtools >/dev/null 2>&1; then
            samtools faidx "$SALMONELLA_REF"
        else
            echo "  warning: no .fai for $SALMONELLA_REF and no samtools; wfmash baseline may fail"
        fi
    fi


    SAL_A_FAIL=0
    SAL_A_ARGS=(-r "$SALMONELLA_REF" -k "$K" --optik "$OPTIK_BIN"
                -n 30 --read-len 50000 -s 42)
    [ -x "$MINIMAP2_BIN" ] && SAL_A_ARGS+=(--minimap2 "$MINIMAP2_BIN")
    [ -n "$WFMASH_BIN" ] && [ -x "$WFMASH_BIN" ] && SAL_A_ARGS+=(--wfmash "$WFMASH_BIN")
    [ -n "$WGATOOLS_BIN" ] && [ -x "$WGATOOLS_BIN" ] && SAL_A_ARGS+=(--wgatools "$WGATOOLS_BIN")
    if python3 "$SCRIPT_DIR/validate_mapping.py" "${SAL_A_ARGS[@]}"; then
        echo "  => LONG-READ SUITE A (error-free 50 kb): PASS"
    else
        echo "  => LONG-READ SUITE A (error-free 50 kb): FAIL"
        SAL_A_FAIL=1
    fi

    SAL_B_FAIL=0
    if [ -n "$WFMASH_BIN" ] && [ -x "$WFMASH_BIN" ] && \
       [ -n "$WGATOOLS_BIN" ] && [ -x "$WGATOOLS_BIN" ]; then
        SAL_B_ARGS=(--ref "$SALMONELLA_REF" --optik "$OPTIK_BIN"
                    --wfmash "$WFMASH_BIN" --wgatools "$WGATOOLS_BIN"
                    -n 30 -l 50000 -e 0.001)
        if python3 "$SCRIPT_DIR/salmonella_long_reads.py" "${SAL_B_ARGS[@]}"; then
            echo "  => LONG-READ SUITE B (0.1% errors, wfmash baseline): PASS"
        else
            echo "  => LONG-READ SUITE B (0.1% errors, wfmash baseline): FAIL"
            SAL_B_FAIL=1
        fi
    else
        skip "long-read suite B" "wfmash and/or wgatools not found"
    fi

    [ "$SAL_A_FAIL" -ne 0 ] && FAIL=1
    [ "$SAL_B_FAIL" -ne 0 ] && FAIL=1
else
    skip "long-read suites (step 7)" "no salmonella reference (run bench/build_tools.sh salmonella, or set SALMONELLA_REF)"
fi

# Step 8.
# Determinism and the huge-page A/B
echo ""
echo "----------------------------------------------"
echo "  [8/9] Mapper determinism (threads, THP)"
step_mark "8 Mapper determinism (threads, THP)"
echo "----------------------------------------------"
# The SAM body must be byte-identical at 1 and N threads, and with huge pages disabled (OPTIK_THP=0) and enabled, since output may never depend on scheduling or memory layout.
# The reads are short and long, from both strands, with errors.
H8_WD="$(mktemp -d /tmp/optik_det_XXXXXX)"
H8_REF="$H8_WD/ref.fa"
python3 - "$REF" "$H8_REF" <<'PY'
import gzip, sys
src, dst = sys.argv[1], sys.argv[2]
open_fn = gzip.open if src.endswith(".gz") else open
with open_fn(src, "rt") as f, open(dst, "w") as o:
    for line in f:
        o.write(line)
PY
H8_FAIL=0
python3 "$SCRIPT_DIR/gen_reads.py" -r "$H8_REF" -n 4000 -l 150 -s 42 -e 0.02 \
    --ins 0.002 --del 0.002 -o "$H8_WD/short.fa" >/dev/null || H8_FAIL=1
python3 "$SCRIPT_DIR/gen_reads.py" -r "$H8_REF" -n 100 -l 8000 -s 43 -e 0.05 \
    --ins 0.01 --del 0.01 -o "$H8_WD/long.fa" >/dev/null || H8_FAIL=1
cat "$H8_WD/short.fa" "$H8_WD/long.fa" > "$H8_WD/reads.fa"
for cfg in "1:1" "$THREADS:1" "$THREADS:0"; do
    IFS=':' read -r t thp <<< "$cfg"
    if ! OPTIK_THP="$thp" "$OPTIK_BIN" map -r "$H8_REF" -q "$H8_WD/reads.fa" -t "$t" \
            -o "$H8_WD/t${t}_thp${thp}.sam" 2>"$H8_WD/t${t}_thp${thp}.err"; then
        echo "  => optik map failed (t=$t THP=$thp): $(tail -1 "$H8_WD/t${t}_thp${thp}.err")"
        H8_FAIL=1
    fi
    grep -v '^@' "$H8_WD/t${t}_thp${thp}.sam" > "$H8_WD/t${t}_thp${thp}.body" 2>/dev/null || true
done
for other in "t${THREADS}_thp1" "t${THREADS}_thp0"; do
    if ! cmp -s "$H8_WD/t1_thp1.body" "$H8_WD/$other.body"; then
        echo "  => SAM body differs: t1_thp1 vs $other"
        H8_FAIL=1
    fi
done
N8=$(wc -l < "$H8_WD/t1_thp1.body" 2>/dev/null || echo 0)
# Identical output is evidence only when the reads map, because an empty or all-unmapped SAM is also identical at every thread count.
M8=$(awk -F'\t' '!and($2,4)' "$H8_WD/t1_thp1.body" 2>/dev/null | wc -l)
NREADS8=$(grep -c '^>' "$H8_WD/reads.fa" 2>/dev/null || echo 0)
[ "$N8" -gt 0 ] && [ "$M8" -ge $((NREADS8 * 9 / 10)) ] || { echo "  => only $M8 mapped records for $NREADS8 reads"; H8_FAIL=1; }
if [ "$H8_FAIL" -eq 0 ]; then
    echo "  => DETERMINISM: PASS ($N8 records, $M8 mapped, identical at t=1, t=$THREADS, THP on/off)"
else
    echo "  => DETERMINISM: FAIL"
    FAIL=1
fi
rm -rf "$H8_WD"

# Step 9.
# Official map suite on a second genome, Salmonella, when present
echo ""
echo "----------------------------------------------"
echo "  [9/9] Official map suite on a second genome"
step_mark "9 Official map suite on a second genome"
echo "----------------------------------------------"
# The truth-scored comparison of optik with minimap2 from step 4, repeated on another genome (2% error by default) so that parity is not tuned to one reference.
if [ -n "${SALMONELLA_REF:-}" ] && [ -f "$SALMONELLA_REF" ] && [ -x "$MINIMAP2_BIN" ]; then
    SAL_WD="$(mktemp -d /tmp/optik_sal_map_XXXXXX)"
    SAL_MAP_ARGS=(-r "$SALMONELLA_REF" -k "$K" -n "$NUM_READS" -t "$THREADS"
                  --read-len 150 --err-rate 0.02 --optik "$OPTIK_BIN" -s 7
                  --minimap2 "$MINIMAP2_BIN" --workdir "$SAL_WD")
    [ -n "$WFMASH_BIN" ] && [ -x "$WFMASH_BIN" ] && SAL_MAP_ARGS+=(--wfmash "$WFMASH_BIN")
    [ -x "$SAMTOOLS_BIN" ] && SAL_MAP_ARGS+=(--samtools "$SAMTOOLS_BIN")
    step_mark "9 official map (salmonella)"
    if python3 "$SCRIPT_DIR/validate_official.py" map "${SAL_MAP_ARGS[@]}" \
            >"$SAL_WD/map.json" 2>"$SAL_WD/map.err"; then
        echo "  => SALMONELLA OFFICIAL MAP: PASS"
    else
        echo "  => SALMONELLA OFFICIAL MAP: FAIL ($(grep FAILED "$SAL_WD/map.err" | head -3))"
        FAIL=1
    fi
else
    skip "second-genome map suite (step 9)" "needs the salmonella reference (bench/build_tools.sh salmonella) and minimap2"
fi

# Step 10.
# Edge cases and negative controls
echo ""
echo "----------------------------------------------"
echo "  [10] Adversarial edge cases and negative controls"
step_mark "10 Adversarial edge cases and negative controls"
echo "----------------------------------------------"
# edge_cases.py covers lowercase, N, IUPAC, CRLF, FASTQ, short and empty records, even k, k = 32, cycles and sequence ends for query, dbg and map against independent oracles.
if python3 "$SCRIPT_DIR/edge_cases.py" --optik "$OPTIK_BIN"; then
    echo "  => EDGE CASES: PASS"
else
    echo "  => EDGE CASES: FAIL"
    FAIL=1
fi
# fuzz_small.py generates random adversarial genomes (reverse-complement palindromes, repeats, homopolymers, N runs, lowercase, short records) at random k.
# It checks membership counts, --ids and -c builds against a Python oracle, dbruijn segments against the solid k-mer set, and thread invariance, over three seeds of 80 rounds each.
FUZZ_OK=1
for seed in 1 2 3; do
    python3 "$SCRIPT_DIR/fuzz_small.py" --optik "$OPTIK_BIN" --rounds 80 --seed "$seed" || FUZZ_OK=0
done
if [ "$FUZZ_OK" -eq 1 ]; then
    echo "  => FUZZ (small adversarial genomes): PASS"
else
    echo "  => FUZZ (small adversarial genomes): FAIL"
    FAIL=1
fi
# mapq_calibration.py --gate compares wrong placements with the number each MAPQ bin implies, on a genome with segmental duplications, tandem arrays, a high-copy family and a satellite block.
# Any over-confident optik bin fails it.
MQC_MM2="$PROJECT_ROOT/temp/minimap2/minimap2"
if [ -x "$MQC_MM2" ]; then
    # Read error rates of 1%, 4% and 8%.
    # The cap and the per-edit ratio were fitted across all three, and 8% failed the gate before the clip-length events were added.
    MQC_OK=1
    for MQC_E in 0.01 0.04 0.08; do
        python3 "$SCRIPT_DIR/mapq_calibration.py" --gate --reads 300000 --error "$MQC_E" \
            --optik "$OPTIK_BIN" --minimap2 "$MQC_MM2" --out "${TMPDIR:-/tmp}/optik_mapq_calibration" --report '' \
            > "/tmp/optik_mqc_$$.log" 2>&1 || { MQC_OK=0; echo "  error $MQC_E:"; grep -E "FAIL|OVERCONF" "/tmp/optik_mqc_$$.log" | grep -v minimap2 | head -4; }
    done
    rm -f "/tmp/optik_mqc_$$.log"
    if [ "$MQC_OK" -eq 1 ]; then
        echo "  => MAPQ CALIBRATION (repeat-rich genome, 1/4/8% error): PASS"
    else
        echo "  => MAPQ CALIBRATION (repeat-rich genome, 1/4/8% error): FAIL"
        FAIL=1
    fi
    # Paired-end mapping (-q R1 -q R2) on the same repeat-rich genome.
    # Per-mate accuracy must reach the single-end level with no over-confident MAPQ band, and proper-pair flags, TLEN and mate fields must be consistent.
    # samtools fixmate and markdup run, and with PICARD_BIN so does Picard ValidateSamFile, next to minimap2 -ax sr R1 R2.
    PB_REF="${TMPDIR:-/tmp}/optik_mapq_calibration/repeats_s17.fa"
    if [ -f "$PB_REF" ] && python3 "$SCRIPT_DIR/paired_bench.py" --gate --ref "$PB_REF" --pairs 50000 \
            --optik "$OPTIK_BIN" --minimap2 "$MQC_MM2" --threads "$THREADS" \
            --out "${TMPDIR:-/tmp}/optik_paired_bench" --report '' > "/tmp/optik_pb_$$.log" 2>&1; then
        echo "  => PAIRED-END (repeat-rich genome, vs minimap2): PASS"
    else
        echo "  => PAIRED-END (repeat-rich genome, vs minimap2): FAIL"
        grep -E "\[FAIL|correct" "/tmp/optik_pb_$$.log" | head -6
        FAIL=1
    fi
    rm -f "/tmp/optik_pb_$$.log"
else
    skip "MAPQ calibration" "minimap2 not built (bench/build_tools.sh minimap2)"
fi
# mapq_calibration.py --multi tests the MAPQ model on fresh reads of several genomes (the seed-303 synthetic genome, Salmonella, and E. coli and hg38 chr21 when present).
# Reads of 100, 150 and 250 bp at 1, 4 and 8% error and HiFi- and ONT-like long reads must show no MAPQ bin with significantly more wrong reads than its MAPQs imply, and no long-read record may carry a false indel of 50 bp or more.
MQM_REFS=""
# E. coli is downloaded into temp/genomes/ on first use, and left out when the download fails.
ECOLI_FA="$(python3 "$SCRIPT_DIR/genomes.py" ecoli 2>/dev/null || true)"
for r in "$ECOLI_FA" "$PROJECT_ROOT/temp/salmonella_LT2.fa" \
         "$PROJECT_ROOT/temp/scaling/chr21.fa"; do
    [ -f "$r" ] && MQM_REFS="${MQM_REFS:+$MQM_REFS,}$r"
done
if python3 "$SCRIPT_DIR/mapq_calibration.py" --multi --gate --refs "$MQM_REFS" --reads 20000 \
        --optik "$OPTIK_BIN" --threads "$THREADS" --out "${TMPDIR:-/tmp}/optik_mapq_multi" --report '' \
        > "/tmp/optik_mqm_$$.log" 2>&1; then
    echo "  => MAPQ CALIBRATION (fresh reads: genomes x read lengths x error, short + long reads): PASS"
else
    echo "  => MAPQ CALIBRATION (fresh reads: genomes x read lengths x error, short + long reads): FAIL"
    grep -E "\[FAIL|\[SKIP" "/tmp/optik_mqm_$$.log" | head -8
    FAIL=1
fi
rm -f "/tmp/optik_mqm_$$.log"
# boundary_cases.py drives every CLI parameter to and past its limits and feeds degenerate inputs.
# None may cause a crash, a panic or a silent default.
if python3 "$SCRIPT_DIR/boundary_cases.py" --optik "$OPTIK_BIN"; then
    echo "  => BOUNDARY CASES: PASS"
else
    echo "  => BOUNDARY CASES: FAIL"
    FAIL=1
fi
# cli_smoke.py runs every command as a user would, with FASTA, FASTQ and gzip input, brute-force query counts, --ids, spss-update, GFA validation, samtools (quickcheck, flagstat, sort, index, calmd) on map's SAM, and the error paths.
if python3 "$SCRIPT_DIR/cli_smoke.py" --optik "$OPTIK_BIN"; then
    echo "  => CLI SMOKE (user workflows): PASS"
else
    echo "  => CLI SMOKE (user workflows): FAIL"
    FAIL=1
fi
# thread_check.py treats -t as a budget for everything a command keeps running, i.e. its workers, its readers decoding .gz in their own thread, and OpenMP.
# At -t 1, 2 and 4 and at the default of half the logical CPUs, no subcommand may have more than t threads runnable at once or a CPU time above t times its wall time.
TC_DIR="/tmp/optik_threads_$$"
if python3 "$SCRIPT_DIR/thread_check.py" --optik "$OPTIK_BIN" --work "$TC_DIR" > "$TC_DIR.log" 2>&1; then
    echo "  => THREAD BUDGET (-t): PASS"
else
    echo "  => THREAD BUDGET (-t): FAIL"
    grep -E "FAIL" "$TC_DIR.log" | cut -c1-200 | head -5
    FAIL=1
fi
rm -rf "$TC_DIR" "$TC_DIR.log"
# spill_check.py forces the k-mer counting buffers to disk and requires byte-identical dbruijn, build-fm and build-spss --compact output at 1 to 16 threads.
# It also checks the --max-mem, graph-over-limit and spill-directory errors.
SPILL_DIR="/tmp/optik_spill_$$"
if python3 "$SCRIPT_DIR/spill_check.py" --optik "$OPTIK_BIN" --workdir "$SPILL_DIR" > "$SPILL_DIR.log" 2>&1; then
    echo "  => DISK SPILL (--max-mem): PASS"
else
    echo "  => DISK SPILL (--max-mem): FAIL"
    grep -E "\[FAIL" "$SPILL_DIR.log" | cut -c1-200 | head -5
    FAIL=1
fi
rm -rf "$SPILL_DIR" "$SPILL_DIR.log"
# test_harness.py runs known-answer tests of the validation and benchmark helpers.
if python3 -W ignore "$SCRIPT_DIR/test_harness.py" >/dev/null 2>"/tmp/optik_harness_$$.log"; then
    echo "  => HARNESS SELF-TESTS: PASS"
else
    echo "  => HARNESS SELF-TESTS: FAIL ($(tail -3 /tmp/optik_harness_$$.log))"
    FAIL=1
fi
rm -f "/tmp/optik_harness_$$.log"
# negative_controls.py requires the real binary to pass every validator, and each of about 30 corruptions of its output or of a baseline tool to make the checks responsible for that aspect fail, so no check can pass vacuously.
if python3 "$SCRIPT_DIR/negative_controls.py" "$REF" "$QUERIES"; then
    echo "  => NEGATIVE CONTROLS: PASS"
else
    echo "  => NEGATIVE CONTROLS: FAIL"
    FAIL=1
fi

# Step 11.
# Optional pipeline parity, which needs bcftools and gfatools
echo ""
echo "----------------------------------------------"
echo "  [11] Pipeline parity (outputs consumed as pipelines consume the reference tools')"
step_mark "11 Pipeline parity"
echo "----------------------------------------------"
BCFTOOLS_BIN="${BCFTOOLS_BIN:-$(command -v bcftools 2>/dev/null || true)}"
[ -z "$BCFTOOLS_BIN" ] && [ -x "$PROJECT_ROOT/temp/bcftools/bcftools" ] && BCFTOOLS_BIN="$PROJECT_ROOT/temp/bcftools/bcftools"
GFATOOLS_BIN="${GFATOOLS_BIN:-$PROJECT_ROOT/temp/gfatools/gfatools}"
PP_GENOME=$(get_dataset_dbg_input)
if [ -z "$BCFTOOLS_BIN" ] || [ ! -x "$BCFTOOLS_BIN" ]; then
    skip "pipeline parity" "bcftools not found (put it on PATH, set BCFTOOLS_BIN, or: bench/build_tools.sh bcftools)"
elif [ ! -x "$GFATOOLS_BIN" ]; then
    skip "pipeline parity" "gfatools not found (bench/build_tools.sh gfatools)"
elif [ -z "$PP_GENOME" ]; then
    skip "pipeline parity" "dataset has no genome (dbg_input)"
else
    PP_DIR="/tmp/optik_pipeline_parity_$$"
    if python3 "$SCRIPT_DIR/pipeline_parity.py" --genome "$PP_GENOME" --workdir "$PP_DIR" \
            --bcftools "$BCFTOOLS_BIN" --gfatools "$GFATOOLS_BIN" --threads "$THREADS" \
            > "$PP_DIR.log" 2>&1; then
        echo "  => PIPELINE PARITY: PASS"
    else
        echo "  => PIPELINE PARITY: FAIL"
        grep -E "\[FAIL" "$PP_DIR.log" | cut -c1-200 | head -8
        FAIL=1
    fi
    rm -rf "$PP_DIR"
fi

# Summary
echo ""
echo "=============================================="
if [ "${#SKIPPED[@]}" -gt 0 ]; then
    echo "  NOT RUN (${#SKIPPED[@]}):"
    for x in "${SKIPPED[@]}"; do echo "    - $x"; done
    [ "${REQUIRE_ALL:-0}" = "1" ] && FAIL=1
fi
step_mark ""
echo ""
echo "  step times:"
for t in "${STEP_TIMES[@]}"; do echo "  $t"; done
echo "  total: $((SECONDS))s"
if [ "$FAIL" -ne 0 ]; then
    echo "  SOME VALIDATIONS FAILED"
elif [ "${#SKIPPED[@]}" -gt 0 ]; then
    echo "  ALL RUN VALIDATIONS PASSED -- ${#SKIPPED[@]} SKIPPED (see above; REQUIRE_ALL=1 fails on skips)"
else
    echo "  ALL VALIDATIONS PASSED"
fi
echo "=============================================="

exit "$FAIL"
