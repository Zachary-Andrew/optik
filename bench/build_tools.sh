#!/usr/bin/env bash
#
# build_tools.sh - fetch and build the external tools the benches compare against
#
# Every external tool the benchmark suite uses is cloned and compiled into temp/, the sibling of bench/ (/path/to/bio/optik/temp/).
#
# Usage
#   ./build_tools.sh              # build all tools
#   ./build_tools.sh sshash       # build only sshash
#   ./build_tools.sh -j4          # pass extra flags to make or cmake --build
#
# tool_bench then runs with
#   TOOL_BENCH_TOOLS_DIR=../temp ./tool_bench.py --dataset ecoli_k31 --all
#
# Without TOOL_BENCH_TOOLS_DIR, tool_bench.py looks for each tool in a sibling directory of optik (../sshash/, ../bifrost/ and so on), the older layout.
#
set -euo pipefail

# Script directory and tool output directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS_DIR="$(cd "$SCRIPT_DIR/.." && pwd)/temp"
# Nix's compiler wrapper drops -march=native unless NIX_ENFORCE_NO_NATIVE is 0, and the tools are compiled for this CPU like optik.
# The variable has no effect outside Nix.
export NIX_ENFORCE_NO_NATIVE=0

# Parallelism comes from the environment or the CPU count.
JOBS="${JOBS:-$(nproc)}"
EXTRA_MAKE_FLAGS=()
SELECTED_TOOLS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        -j*)
            JOBS="${1#-j}"
            shift
            ;;
        *)
            SELECTED_TOOLS+=("$1")
            shift
            ;;
    esac
done

# With no tools named, everything is built.
ALL_TOOLS=(sshash fmsi kmercamel minimap2 bifrost cuttlefish sbwt kmc samtools bwa bcftools picard wfmash wgatools wfa2 gfatools salmonella)
# Targets built only when named, because of their large downloads (build_tools.sh scaling | hg38)
EXTRA_TARGETS=(scaling hg38 zlib-ng rtg)
if [[ ${#SELECTED_TOOLS[@]} -eq 0 ]]; then
    SELECTED_TOOLS=("${ALL_TOOLS[@]}")
fi

info()  { printf "\033[1;34m[build]\033[0m %s\n" "$*"; }
ok()    { printf "\033[1;32m[  ok]\033[0m %s\n" "$*"; }
fail()  { printf "\033[1;31m[fail]\033[0m %s\n" "$*" >&2; exit 1; }

# Stop unless a required dependency is on PATH.
check_dep() {
    command -v "$1" >/dev/null 2>&1 || fail "required dependency not found: $1"
}

# Clone a git repository, or update it when it already exists.
clone_or_update() {
    local url="$1" dest="$2"
    if [[ -d "$dest/.git" && "${NO_UPDATE:-0}" == 1 ]]; then
        info "NO_UPDATE=1: building $dest at $(git -C "$dest" rev-parse --short HEAD)"
    elif [[ -d "$dest/.git" ]]; then
        info "updating $dest"
        git -C "$dest" fetch origin
        git -C "$dest" reset --hard origin/HEAD
        git -C "$dest" submodule update -q --init --recursive
    else
        info "cloning $url -> $dest"
        git clone --recursive "$url" "$dest"
    fi
}

# Run one build function, warning and moving on when it fails, so that one broken tool does not keep the others from being built.
try_build() {
    local name="$1"
    shift
    if "$@"; then
        return 0
    else
        printf "\033[1;33m[skip]\033[0m %s failed to build -- will be unavailable\n" "$name"
        return 0
    fi
}

for dep in git cmake make gcc g++; do
    check_dep "$dep"
done

mkdir -p "$TOOLS_DIR"
info "tools will be built into $TOOLS_DIR"

# CMake flags shared by every tool, which silence warnings-as-errors that older code bases trip with newer compilers.
# Individual tools may add flags of their own.
#
# -Wno-template-body.
# Several tools here (sbwt, fmsi and anything else that pulls in sdsl-lite as a submodule) vendor sdsl-lite (github.com/simongog/sdsl-lite), which has a long-standing upstream bug in louds_tree::swap().
# That function still names the members m_select1/m_select0, which the rest of the class renamed to m_bv_select1/m_bv_select0.
# GCC 13 added the -Wtemplate-body diagnostic, which checks member references inside template bodies more strictly and rightly reports this as a hard error, whereas older GCC deferred the check and never caught it.
# sdsl-lite is a third-party submodule cloned fresh on every build, so there is no local source to patch, and the diagnostic is silenced to restore the older behaviour.
# Compiled semantics do not change, since -Wtemplate-body is only a diagnostic that did not exist before GCC 13.
#
# -DCMAKE_POLICY_DEFAULT_CMP0167=OLD.
# CMake 3.30 deprecated its bundled FindBoost.cmake, the module-mode lookup that finds Boost by searching for its headers and libraries.
# Without that module, find_package(Boost) falls back to config mode, which needs Boost to ship its own BoostConfig.cmake files.
# nixpkgs builds Boost with Boost's own b2 system rather than CMake and never generates those files, so config mode cannot succeed there however Boost is installed.
# The OLD policy brings module mode back, which works with a plain Boost install found through CMAKE_PREFIX_PATH.
#
# Every tool is compiled for the host CPU, as optik is (build.rs uses -march=native), since comparing a native optik with generic builds of the other tools would flatter optik.
# NATIVE=0 builds generic binaries instead.
if [[ "${NATIVE:-1}" == 1 ]]; then OPT_FLAGS="-O3 -march=native"; else OPT_FLAGS="-O3"; fi
_CMAKE_FLAGS=(
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_POLICY_VERSION_MINIMUM=3.5
    -DCMAKE_POLICY_DEFAULT_CMP0167=OLD
    -DCMAKE_CXX_FLAGS="-Wno-error -Wno-unused-but-set-variable -Wno-unknown-attributes -Wno-template-body $OPT_FLAGS"
    -DCMAKE_C_FLAGS="$OPT_FLAGS"
)

# The flags are also exported as CXXFLAGS, so Makefile-based tools that build sdsl-lite through their own cmake or make call (fmsi, whose Makefile builds src/sdsl-lite as a submodule) get the same suppression without going through _CMAKE_FLAGS.
export CXXFLAGS="${CXXFLAGS:-} -Wno-template-body"

# sshash  (https://github.com/jermp/sshash)
#   CMake project that produces build/sshash
build_sshash() {
    local src="$TOOLS_DIR/sshash"
    clone_or_update "https://github.com/jermp/sshash.git" "$src"
    # SSHash sizes its control codewords without the skew-index heavy-load form, so at short m the codeword is truncated and the index is corrupt, as `build --check` reports with "wrong id assigned", which patches/sshash_control_bits.py fixes.
    # Indexes that built correctly before still build identically.
    if [ -f "$SCRIPT_DIR/patches/sshash_control_bits.py" ]; then
        if python3 "$SCRIPT_DIR/patches/sshash_control_bits.py" "$src"; then
            ok "sshash control-codeword width patch applied"
        else
            info "warning: sshash control-codeword patch could not be applied (upstream drift?)"
        fi
    fi
    rm -rf "$src/build"
    mkdir -p "$src/build"
    local native=OFF; [[ "${NATIVE:-1}" == 1 ]] && native=ON
    cmake -S "$src" -B "$src/build" "${_CMAKE_FLAGS[@]}" -DSSHASH_USE_ARCH_NATIVE=$native
    cmake --build "$src/build" -j "$JOBS"
    [[ -x "$src/build/sshash" ]] || fail "sshash build produced no binary"
    ok "sshash -> $src/build/sshash"
}

# bifrost  (https://github.com/pmelsted/bifrost)
#   CMake project that produces build/Bifrost
build_bifrost() {
    local src="$TOOLS_DIR/bifrost"
    clone_or_update "https://github.com/pmelsted/bifrost.git" "$src"
    # Upstream bug, the copy constructor reads o.sz_link instead of o.unitig_cs_link.
    sed -i 's/o\.sz_link\[i\]/o.unitig_cs_link[i]/g' "$src/src/DataStorage.tcc"
    rm -rf "$src/build"
    mkdir -p "$src/build"
    cmake -S "$src" -B "$src/build" "${_CMAKE_FLAGS[@]}"
    cmake --build "$src/build" -j "$JOBS"
    [[ -x "$src/build/src/Bifrost" ]] || fail "bifrost build produced no binary"
    ok "bifrost -> $src/build/src/Bifrost"
}

# cuttlefish  (https://github.com/COMBINE-lab/cuttlefish)
#   Cargo project that produces target/release/cuttlefish
build_cuttlefish() {
    local src="$TOOLS_DIR/cuttlefish"
    clone_or_update "https://github.com/COMBINE-lab/cuttlefish.git" "$src"
    command -v cargo >/dev/null 2>&1 || fail "cuttlefish requires Rust (cargo). Install from https://rustup.rs"
    if [[ "${NATIVE:-1}" == 1 ]]; then
        RUSTFLAGS="-C target-cpu=native" cargo build --release --manifest-path "$src/Cargo.toml"
    else
        cargo build --release --manifest-path "$src/Cargo.toml"
    fi
    [[ -x "$src/target/release/cuttlefish" ]] || fail "cuttlefish build produced no binary"
    ok "cuttlefish -> $src/target/release/cuttlefish"
}

# fmsi  (https://github.com/OndrejSladky/fmsi)
#   make builds the fmsi binary in the repository root
build_fmsi() {
    local src="$TOOLS_DIR/fmsi"
    clone_or_update "https://github.com/OndrejSladky/fmsi.git" "$src"
    rm -rf "$src/src/sdsl-lite/build"
    mkdir -p "$src/src/sdsl-lite/build"
    rm -f "$src/fmsi" "$src/src/fmsi"
    make -C "$src" -j "$JOBS" CXXFLAGS="-g -Wall -Wno-unused-function -std=c++17 $OPT_FLAGS"
    [[ -x "$src/fmsi" ]] || fail "fmsi build produced no binary"
    ok "fmsi -> $src/fmsi"
}

# kmercamel  (https://github.com/OndrejSladky/kmercamel)
#   Built with make into the repository root, kmercamel needs
#   libglpk-dev (apt) or glpk (brew)
build_kmercamel() {
    local src="$TOOLS_DIR/kmercamel"
    # kmercamel needs libglpk-dev, so a missing header gets a message saying which package to install.
    if ! printf '#include <glpk.h>\nint main(){return 0;}\n' \
         | g++ -x c++ - -o /dev/null -lglpk 2>/dev/null; then
        fail "kmercamel requires libglpk-dev.  Install with: apt install libglpk-dev  OR  brew install glpk"
    fi
    clone_or_update "https://github.com/OndrejSladky/kmercamel.git" "$src"
    make -C "$src" clean >/dev/null 2>&1 || true
    make -C "$src" -j "$JOBS" CXXFLAGS="-g -Wall -Wno-unused-function -std=c++17 $OPT_FLAGS"
    [[ -x "$src/kmercamel" ]] || fail "kmercamel build produced no binary"
    ok "kmercamel -> $src/kmercamel"
}


# minimap2  (https://github.com/lh3/minimap2)
#   Its Makefile leaves the minimap2 binary in the repository root
build_minimap2() {
    local src="$TOOLS_DIR/minimap2"
    clone_or_update "https://github.com/lh3/minimap2.git" "$src"
    make -C "$src" clean >/dev/null 2>&1 || true
    make -C "$src" -j "$JOBS" CFLAGS="-g -Wall $OPT_FLAGS -Wc++-compat"
    [[ -x "$src/minimap2" ]] || fail "minimap2 build produced no binary"
    ok "minimap2 -> $src/minimap2"
}

# gfatools  (https://github.com/lh3/gfatools)
#   GFA parser and statistics, which bench/pipeline_parity.py uses to check optik's GFA
#   the way graph tools read it
build_gfatools() {
    local src="$TOOLS_DIR/gfatools"
    clone_or_update "https://github.com/lh3/gfatools.git" "$src"
    make -C "$src" -j "$JOBS"
    [[ -x "$src/gfatools" ]] || fail "gfatools build produced no binary"
    ok "gfatools -> $src/gfatools"
}

# sbwt  (https://github.com/algbio/SBWT)
#   CMake project that produces build/bin/sbwt
build_sbwt() {
    local src="$TOOLS_DIR/SBWT"
    clone_or_update "https://github.com/algbio/SBWT.git" "$src"
    # Upstream bug, cxxopts.hpp uses uint8_t without including <cstdint>.
    # Older libstdc++ pulled the header in through <regex> or <string>, which newer GCC no longer guarantees, leaving the type undefined.
    head -1 "$src/include/sbwt/cxxopts.hpp" | grep -q '<cstdint>' || sed -i '1i #include <cstdint>' "$src/include/sbwt/cxxopts.hpp"
    rm -rf "$src/build"
    mkdir -p "$src/build"
    # No -march=native here.
    # SBWT's vendored KMC builds per-ISA radix-sort objects (-msse2 up to -mavx2) for runtime dispatch, and a global -march=native breaks that link.
    # SBWT itself gates only its Elias-Fano paths on BMI2.
    local sbwt_flags="-Wno-error -Wno-unused-but-set-variable -Wno-unknown-attributes -Wno-template-body -O3"
    [[ "${NATIVE:-1}" == 1 ]] && sbwt_flags+=" -mbmi2 -mpopcnt"
    cmake -S "$src" -B "$src/build" -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
        -DCMAKE_POLICY_DEFAULT_CMP0167=OLD -DCMAKE_CXX_FLAGS="$sbwt_flags" -DMAX_KMER_LENGTH=32
    cmake --build "$src/build" -j "$JOBS"
    [[ -x "$src/build/bin/sbwt" ]] || fail "sbwt build produced no binary"
    ok "sbwt -> $src/build/bin/sbwt"
}


# kmc  (https://github.com/refresh-bio/KMC)
#   make project that produces bin/kmc, bin/kmc_dump and bin/kmc_tools, whose set
#   operations bench/compare_tools_published_sets.py uses to compare unitig k-mer sets.
#   Only these targets are built, linked dynamically to avoid static libc problems on
#   NixOS and similar systems.
build_kmc() {
    local src="$TOOLS_DIR/KMC"
    clone_or_update "https://github.com/refresh-bio/KMC" "$src"
    make -C "$src" -j "$JOBS" kmc kmc_dump kmc_tools STATIC_CFLAGS= STATIC_LFLAGS=
    [[ -x "$src/bin/kmc" ]] || fail "kmc build produced no binary"
    [[ -x "$src/bin/kmc_dump" ]] || fail "kmc build produced no kmc_dump"
    [[ -x "$src/bin/kmc_tools" ]] || fail "kmc build produced no kmc_tools"
    ok "kmc -> $src/bin/kmc"
}

# samtools  (https://github.com/samtools/samtools)
#   A system samtools is preferred and the source build is a fallback, because tools
#   further down may not find a locally built htslib at run time (samtools needs
#   LD_LIBRARY_PATH for it).
build_samtools() {
    if command -v samtools >/dev/null 2>&1; then
        ok "samtools: using system binary $(command -v samtools)"
        return 0
    fi
    local src="$TOOLS_DIR/samtools"
    clone_or_update "https://github.com/samtools/samtools.git" "$src"
    clone_or_update "https://github.com/samtools/htslib.git" "$src/htslib"
    # htslib is built first, and samtools picks up the submodule at run time when present.
    ( cd "$src" && git submodule update --init 2>/dev/null || true )
    ( cd "$src" && autoreconf -i 2>/dev/null && ./configure && make -j "$JOBS" ) \
        || make -C "$src" -j "$JOBS"
    [[ -x "$src/samtools" ]] || fail "samtools build produced no binary"
    # samtools links htslib through -L../htslib, which is exported for later commands.
    export LD_LIBRARY_PATH="$src/htslib:${LD_LIBRARY_PATH:-}"
    ok "samtools -> $src/samtools"
}

# bwa  (https://github.com/lh3/bwa)
#   BWA-MEM, the short-read baseline.  A system bwa is used when present, otherwise the
#   latest source is built, which needs zlib.
build_bwa() {
    if command -v bwa >/dev/null 2>&1; then
        ok "bwa: using system binary $(command -v bwa)"
        return 0
    fi
    local src="$TOOLS_DIR/bwa"
    clone_or_update "https://github.com/lh3/bwa.git" "$src"
    make -C "$src" -j "$JOBS"
    [[ -x "$src/bwa" ]] || fail "bwa build produced no binary"
    ok "bwa -> $src/bwa"
}

# bcftools  (https://github.com/samtools/bcftools)
#   Variant calling in the pipeline-parity validation.  A system bcftools is used when
#   present, otherwise it is built.
build_bcftools() {
    if command -v bcftools >/dev/null 2>&1; then
        ok "bcftools: using system binary $(command -v bcftools)"
        return 0
    fi
    local src="$TOOLS_DIR/bcftools" tag
    clone_or_update "https://github.com/samtools/bcftools.git" "$src"
    clone_or_update "https://github.com/samtools/htslib.git" "$src/htslib"
    # The latest release of each, by numeric tag, since the default branch is development.
    tag=$(git -C "$src" tag | grep -E '^[0-9]+\.[0-9]+(\.[0-9]+)?$' | sort -V | tail -1)
    git -C "$src" checkout -q "$tag"
    tag=$(git -C "$src/htslib" tag | grep -E '^[0-9]+\.[0-9]+(\.[0-9]+)?$' | sort -V | tail -1)
    git -C "$src/htslib" checkout -q "$tag"
    ( cd "$src/htslib" && git submodule update --init -q 2>/dev/null || true )
    # Configured without network I/O (libcurl, GCS, S3), because the validation reads local files and the Makefile's default configuration expects libcurl headers.
    ( cd "$src/htslib" && rm -f config.h && autoreconf -i >/dev/null 2>&1 \
        && ./configure --disable-libcurl --disable-gcs --disable-s3 >/dev/null ) || fail "htslib configure failed"
    # htslib is cloned inside the bcftools tree, not beside it as the Makefile's ../htslib default expects.
    make -C "$src" -j "$JOBS" HTSDIR=htslib
    [[ -x "$src/bcftools" ]] || fail "bcftools build produced no binary"
    ok "bcftools -> $src/bcftools"
}

# picard  (https://github.com/broadinstitute/picard)
#   ValidateSamFile, the SAM spec check of every mapper's output.  A system `picard` is
#   used when present, otherwise the latest release jar, which needs java on PATH.
build_picard() {
    if command -v picard >/dev/null 2>&1; then
        ok "picard: using system binary $(command -v picard)"
        return 0
    fi
    command -v java >/dev/null 2>&1 || fail "picard needs java (a JRE 17+) on PATH"
    local dir="$TOOLS_DIR/picard" url
    mkdir -p "$dir"
    url=$(curl -fsSL https://api.github.com/repos/broadinstitute/picard/releases/latest \
          | python3 -c "import json,sys; print([a['browser_download_url'] for a in json.load(sys.stdin)['assets'] if a['name'] == 'picard.jar'][0])")
    curl -fsSL --retry 3 -o "$dir/picard.jar.part" "$url" && mv "$dir/picard.jar.part" "$dir/picard.jar"
    # Picard prints its version but exits 1 for this invocation, so check the output rather than the status.
    local out
    out=$(java -jar "$dir/picard.jar" ValidateSamFile --version 2>&1) || true
    [[ "$out" == *Version:* ]] || fail "picard.jar does not run: $out"
    ok "picard -> $dir/picard.jar ($url)"
}

# wfmash  (https://github.com/waveygang/wfmash)
#   Long-read pairwise aligner, the successor of the pangenome MashMap pipeline.  It
#   vendors minimap2, htslib and WFA2-lib as submodules and needs Boost (system,
#   iostreams, program_options), zlib and bzip2 at configure time.  It is a benchmark
#   baseline only and never linked into the optik crate.
#
#   Where the dev headers are not on the default include path, as on NixOS, pass
#   CMake cache variables through WFMASH_CMAKE_ARGS, as in
#     WFMASH_CMAKE_ARGS="-DBZIP2_ROOT=/nix/store/...-bzip2-1.0.8-dev \
#                        -DBoost_ROOT=/nix/store/...-boost-1.89.0" \
#     ./build_tools.sh wfmash
#   and export CPATH and LIBRARY_PATH pointing at the dependencies' -dev and lib store
#   paths, so that both CMake and the vendored htslib's ./configure (which cannot see
#   CMAKE_*_PATH) find the headers and libraries.
build_wfmash() {
    check_dep cmake
    local src="$TOOLS_DIR/wfmash"
    clone_or_update "https://github.com/waveygang/wfmash.git" "$src"
    ( cd "$src" && git submodule update --init --recursive )
    # WFA2-lib dereferences misaligned uint64_t pointers in its block-compare extend kernels, which is undefined behaviour.
    # GCC 16 and later vectorize those loops assuming 8-byte alignment, and wfmash then segfaults in BiWFA, which patches/wfa2_unaligned_loads.py fixes.
    if [ -f "$SCRIPT_DIR/patches/wfa2_unaligned_loads.py" ]; then
        if python3 "$SCRIPT_DIR/patches/wfa2_unaligned_loads.py" "$src/deps/WFA2-lib"; then
            ok "wfmash WFA2-lib unaligned-load patch applied"
        else
            info "warning: wfmash WFA2-lib unaligned-load patch could not be applied (upstream drift?)"
        fi
    fi
    local extra=()
    if [[ -n "${WFMASH_CMAKE_ARGS:-}" ]]; then
        read -r -a extra <<< "${WFMASH_CMAKE_ARGS}"
    fi
    # Without a system htslib, since wfmash has no htslib submodule, wfmash is told to download and build htslib itself.
    # This is skipped when VENDOR_HTSLIB is already set or a real htslib is visible to the compiler.
    if [[ " ${extra[*]} " != *"VENDOR_HTSLIB"* ]] \
       && ! echo '#include <htslib/faidx.h>' | cc -E - >/dev/null 2>&1; then
        info "system htslib not found; vendoring htslib (-DVENDOR_HTSLIB=ON)"
        extra+=("-DVENDOR_HTSLIB=ON")
    fi
    ( cd "$src" && cmake -S . -B build -DCMAKE_BUILD_TYPE=Release "${extra[@]}" \
        && cmake --build build -j "$JOBS" )
    [[ -x "$src/build/bin/wfmash" ]] || fail "wfmash build produced no binary (set WFMASH_CMAKE_ARGS for library roots)"
    ok "wfmash -> $src/build/bin/wfmash"
}

# salmonella  NCBI RefSeq NC_003197.1, the S.
# Typhimurium LT2 chromosome
#   Second-genome reference for run_validations.sh step 7, the long-read suite of 50 kb
#   and more with the wfmash baseline, and step 9, the official map suite, both skipped
#   without it.  It is fetched from NCBI E-utilities, checked by accession and its exact
#   length of 4,857,432 bp, and written atomically.
fetch_salmonella() {
    check_dep curl
    local out="$TOOLS_DIR/salmonella_LT2.fa"
    if [[ -s "$out" ]] && head -c 13 "$out" | grep -q '^>NC_003197.1'; then
        ok "salmonella -> $out (present)"
    else
        local tmp="$out.part.$$"
        curl -fsSL --retry 3 -o "$tmp" \
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_003197.1&rettype=fasta&retmode=text" \
            || { rm -f "$tmp"; fail "salmonella download failed"; }
        local len
        len=$(grep -v '^>' "$tmp" | tr -d '\n\r' | wc -c)
        if ! head -c 13 "$tmp" | grep -q '^>NC_003197.1' || [[ "$len" -ne 4857432 ]]; then
            rm -f "$tmp"; fail "salmonella download is not NC_003197.1 (length $len)"
        fi
        mv "$tmp" "$out"
        ok "salmonella -> $out (NC_003197.1, $len bp)"
    fi
    # wfmash needs the .fai, and samtools may not be on PATH yet (it may be built into temp/).
    local st
    st=$(command -v samtools || echo "$TOOLS_DIR/samtools/samtools")
    [[ -x "$st" ]] && "$st" faidx "$out" && ok "salmonella .fai written"
}

# scaling  UCSC hg38 chr21 and chr1 from the GRCh38 primary assembly
#   Inputs of bench/scaling_bench.py, the large-genome comparison.  Each chromosome is
#   checked against its GRCh38 length, N bases included.
# hg38  (UCSC hg38, primary chromosomes chr1-22, X, Y and M, 3.1 Gbp)
#   The whole-genome graph test (the README's human genome with 2.5e9 31-mers, built
#   compact, its unitigs checked against Cuttlefish's with bench/cpp/unitig_digest.cpp).
#   About 1 GB to download.
fetch_hg38() {
    check_dep curl
    check_dep gunzip
    mkdir -p "$TOOLS_DIR/scaling"
    local out="$TOOLS_DIR/scaling/hg38_primary.fa"
    if [[ -s "$out" ]] && [[ "$(grep -c '^>' "$out")" -eq 25 ]]; then
        ok "hg38 -> $out (present)"
        return
    fi
    local tmp="$out.part.$$"
    curl -fsSL --retry 3 "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz" \
        | gunzip -c | awk '/^>/{keep = ($1 ~ /^>chr([0-9]+|X|Y|M)$/)} keep' > "$tmp" \
        || { rm -f "$tmp"; fail "hg38 download failed"; }
    if [[ "$(grep -c '^>' "$tmp")" -ne 25 ]]; then
        rm -f "$tmp"; fail "hg38 download does not hold the 25 primary chromosomes"
    fi
    mv "$tmp" "$out"
    ok "hg38 -> $out"
}

fetch_scaling() {
    check_dep curl
    check_dep gunzip
    mkdir -p "$TOOLS_DIR/scaling"
    local c want
    for c in chr21:46709983 chr1:248956422; do
        want=${c#*:}; c=${c%%:*}
        local out="$TOOLS_DIR/scaling/$c.fa"
        if [[ -s "$out" ]] && [[ "$(grep -v '^>' "$out" | tr -d '\n\r' | wc -c)" -eq "$want" ]]; then
            ok "scaling $c -> $out (present)"
            continue
        fi
        local tmp="$out.part.$$"
        curl -fsSL --retry 3 "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/$c.fa.gz" \
            | gunzip -c > "$tmp" || { rm -f "$tmp"; fail "scaling $c download failed"; }
        local len
        len=$(grep -v '^>' "$tmp" | tr -d '\n\r' | wc -c)
        if ! head -c 7 "$tmp" | grep -q "^>$c" || [[ "$len" -ne "$want" ]]; then
            rm -f "$tmp"; fail "scaling $c download is not hg38 $c (length $len, want $want)"
        fi
        mv "$tmp" "$out"
        ok "scaling $c -> $out ($len bp)"
    done
}

# wgatools  (https://github.com/wjwei-hust/wgatools)
#   Rust toolkit for whole-genome alignment formats, with PAF validation, pafcov and more,
#   built with cargo.  bench/validate_mapping.py uses it to check the PAF converted from
#   optik's SAM and its reference coverage.  A bench-only baseline.
build_wgatools() {
    check_dep cargo
    local src="$TOOLS_DIR/wgatools"
    clone_or_update "https://github.com/wjwei-handsome/wgatools.git" "$src"
    ( cd "$src" && cargo build --release )
    [[ -x "$src/target/release/wgatools" ]] || fail "wgatools build produced no binary"
    ok "wgatools -> $src/target/release/wgatools"
}

# wfa2  (https://github.com/smarco/WFA2-lib)
#   The Makefile produces lib/libwfa.a, lib/libwfacpp.a and the reference
#   binaries bin/align_benchmark and bin/generate_dataset.  BUILD_TOOLS=1 builds the
#   align_benchmark oracle of the waveform exactness bench, and BUILD_EXAMPLES=0 skips the
#   examples target, whose precompiled keras and pybind dependencies are fragile.
build_wfa2() {
    local src="$TOOLS_DIR/WFA2-lib"
    clone_or_update "https://github.com/smarco/WFA2-lib.git" "$src"
    # Apply the boundary-key fix for edit and indel ends-free alignment (patches/wfa2_edit_endsfree.py).
    # clone_or_update hard-resets to origin/HEAD, so this must come after it.
    # The patcher is idempotent and leaves the tree untouched if upstream has drifted, with only a warning, in which case the align_benchmark oracle still builds and ends-free edit alignment falls back to end-to-end.
    if [ -f "$SCRIPT_DIR/patches/wfa2_edit_endsfree.py" ]; then
        if python3 "$SCRIPT_DIR/patches/wfa2_edit_endsfree.py" \
                "$src/wavefront/wavefront_compute.c"; then
            ok "wfa2 edit-endsfree patch applied"
        else
            info "warning: wfa2 ends-free patch could not be applied (upstream drift?)"
        fi
    fi
    # The same fix for misaligned uint64_t loads as in wfmash's copy of WFA2-lib
    if [ -f "$SCRIPT_DIR/patches/wfa2_unaligned_loads.py" ]; then
        if python3 "$SCRIPT_DIR/patches/wfa2_unaligned_loads.py" "$src"; then
            ok "wfa2 unaligned-load patch applied"
        else
            info "warning: wfa2 unaligned-load patch could not be applied (upstream drift?)"
        fi
    fi
    # Build the core library and the align_benchmark tool, skipping the examples.
    # WFA2-lib's Makefile does not order lib_wfa before the tools that link lib/libwfa.a, so a single -j build races with "No rule to make target '../../lib/libwfa.a'".
    # The library is therefore built in parallel first and the tools serially afterwards.
    # Both go through the default `all` target, which carries -O3 -march=native as a target-specific variable.
    make -C "$src" setup
    make -C "$src" -j "$JOBS" BUILD_TOOLS=0 BUILD_EXAMPLES=0
    make -C "$src" BUILD_TOOLS=1 BUILD_EXAMPLES=0
    [[ -f "$src/lib/libwfa.a" ]] || fail "wfa2 build produced no library"
    [[ -x "$src/bin/align_benchmark" ]] || fail "wfa2 build produced no align_benchmark"
    ok "wfa2 -> $src/lib/libwfa.a + $src/bin/align_benchmark"

    # Build the waveform oracle harness (optik's exact edit distance against WFA2-lib), a convenience for run_validations.sh [6/7] and manual use that compiles optik's include/ headers directly.
    # The bench rebuilds it itself when this step is skipped, so a failure here is informational.
    if command -v g++ >/dev/null 2>&1 && \
       [ -f "$SCRIPT_DIR/../include/wfa_affine.hpp" ]; then
        if g++ -std=c++17 -O2 -march=native -I "$SCRIPT_DIR/../include" \
                -I "$src" "$SCRIPT_DIR/cpp/wfa2_oracle.cpp" \
                "$src/lib/libwfa.a" -lm -fopenmp -o "$src/bin/wfa2_oracle"; then
            ok "wfa2 oracle -> $src/bin/wfa2_oracle"
        else
            info "warning: wfa2_oracle harness build failed (bench will compile it as needed)"
        fi
    fi
}

# zlib-ng  (https://github.com/zlib-ng/zlib-ng)
#   A faster zlib, built at its latest release tag in zlib-compatible mode as a static
#   library under temp/zlib-ng.  optik's build.rs links it when it is there, so gzipped
#   input inflates faster in every reader, and falls back to the system zlib otherwise.
build_zlib_ng() {
    local src="$TOOLS_DIR/zlib-ng-src" dest="$TOOLS_DIR/zlib-ng"
    clone_or_update "https://github.com/zlib-ng/zlib-ng.git" "$src"
    local tag
    tag="$(git -C "$src" tag -l '[0-9]*.[0-9]*.[0-9]*' | sort -V | tail -1)"
    [[ -n "$tag" ]] || fail "zlib-ng: no release tag found"
    git -C "$src" checkout -q "$tag"
    rm -rf "$src/build" "$dest"
    cmake -S "$src" -B "$src/build" -DCMAKE_BUILD_TYPE=Release -DZLIB_COMPAT=ON \
        -DBUILD_SHARED_LIBS=OFF -DZLIB_ENABLE_TESTS=OFF -DZLIBNG_ENABLE_TESTS=OFF \
        -DWITH_GTEST=OFF -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_INSTALL_PREFIX="$dest"
    cmake --build "$src/build" -j "$JOBS"
    cmake --install "$src/build"
    [[ -f "$dest/lib/libz.a" || -f "$dest/lib64/libz.a" ]] || fail "zlib-ng build produced no libz.a"
    ok "zlib-ng $tag -> $dest (rebuild optik with cargo build --release to link it)"
}

# rtg  (RTG Tools, https://github.com/RealTimeGenomics/rtg-tools)
#   vcfeval, the haplotype-aware VCF comparison inside GIAB's hap.py, which
#   bench/real_reads.py variants uses.  The latest release's Linux build bundles its
#   own JRE, so no system Java is needed.  Installed as temp/rtg/rtg.
build_rtg() {
    local dest="$TOOLS_DIR/rtg" tag zip
    tag="$(curl -fsSL https://api.github.com/repos/RealTimeGenomics/rtg-tools/releases/latest \
        | python3 -c 'import json,sys; print(json.load(sys.stdin)["tag_name"])')"
    [[ -n "$tag" ]] || fail "rtg: no release found"
    zip="$TOOLS_DIR/rtg-tools-$tag-linux-x64.zip"
    curl -fsSL -o "$zip" \
        "https://github.com/RealTimeGenomics/rtg-tools/releases/download/$tag/rtg-tools-$tag-linux-x64.zip"
    rm -rf "$dest" "$TOOLS_DIR/rtg-tools-$tag"
    # Python's zipfile, since unzip is not on every machine.
    # It drops the permission bits, which are restored from the archive so the launcher and the bundled JRE stay executable.
    python3 - "$zip" "$TOOLS_DIR" <<'PY'
import os, sys, zipfile
z = zipfile.ZipFile(sys.argv[1])
for info in z.infolist():
    path = z.extract(info, sys.argv[2])
    mode = (info.external_attr >> 16) & 0o777
    if mode:
        os.chmod(path, mode)
PY
    mv "$TOOLS_DIR/rtg-tools-$tag" "$dest"
    rm -f "$zip"
    # Answer the first-run usage-logging question once, so later runs are non-interactive.
    printf 'RTG_TALKBACK=false\nRTG_USAGE=false\n' >> "$dest/rtg.cfg"
    "$dest/rtg" version >/dev/null || fail "rtg build produced no working rtg"
    ok "rtg $tag -> $dest/rtg"
}

for tool in "${SELECTED_TOOLS[@]}"; do
    case "$tool" in
        sshash)      build_sshash ;;
        bifrost)     build_bifrost ;;
        cuttlefish)  build_cuttlefish ;;
        fmsi)        build_fmsi ;;
        kmercamel)   build_kmercamel ;;
        minimap2)    build_minimap2 ;;
        sbwt)        build_sbwt ;;
        kmc)         build_kmc ;;
        samtools)    build_samtools ;;
        bwa)         build_bwa ;;
        bcftools)    build_bcftools ;;
        picard)      build_picard ;;
        wfmash)      build_wfmash ;;
        wgatools)    build_wgatools ;;
        wfa2)        build_wfa2 ;;
        gfatools)    build_gfatools ;;
        salmonella)  fetch_salmonella ;;
        scaling)     fetch_scaling ;;
        hg38)        fetch_hg38 ;;
        zlib-ng)     build_zlib_ng ;;
        rtg)         build_rtg ;;
        *)           fail "unknown tool: $tool (valid: ${ALL_TOOLS[*]} ${EXTRA_TARGETS[*]})" ;;
    esac
done

echo ""
info "all selected tools built successfully in $TOOLS_DIR"
info "run benchmarks with:  TOOL_BENCH_TOOLS_DIR=$TOOLS_DIR ./tool_bench.py --dataset ecoli_k31 --all"
