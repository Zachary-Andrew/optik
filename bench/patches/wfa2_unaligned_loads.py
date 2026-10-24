#!/usr/bin/env python3
"""wfa2_unaligned_loads.py - idempotent patch that makes WFA2-lib's 8-byte block compares and stores use well-defined unaligned memory access.

Why this patch exists
WFA2-lib compares sequences 8 characters at a time by casting `char*` positions to `uint64_t*` and dereferencing them.

    uint64_t* pattern_blocks = (uint64_t*)(pattern + v);
    uint64_t cmp = *pattern_blocks ^ *text_blocks;   // v, h: any offset

`pattern + v` has byte alignment, so this dereferences a misaligned `uint64_t*`, which is undefined behaviour in C.
Up to GCC 15 the code still compiled to plain unaligned `mov`s.
GCC 16 (-Ofast -funroll-all-loops, as wfmash builds it) vectorizes the `while (cmp == 0)` loop with early break and peels it for 32-byte alignment on the assumption that a `uint64_t*` is 8-byte aligned.
On a pointer 4 bytes off that peel yields an aligned `vmovdqa` on a misaligned address, so wfmash segfaults (SIGSEGV, exit -11) in wavefront_extend_matches_packed_end2end_avx2 during BiWFA.
In optik's validation that showed up as long-read suite B "wfmash_rc: -11" and step 2 "read_opts: exit -11".

The fix routes every such access through memcpy-based macros (WF_LOAD_U64 / WF_STORE_U64).
They compile to the same single unaligned `mov` but carry no alignment assumption, so the code is correct under any optimizer.
They are macros, not static functions, because the patched kernels are non-static `inline` functions (FORCE_INLINE), which C99 forbids from referencing a static function (GCC reports "'f' is static but used in inline function ... which is not static"), whereas memcpy itself has external linkage.
Files patched by the first version of this script (static wf_load_u64/wf_store_u64 helpers) are upgraded in place.
The patched sites are
  wavefront/wavefront_extend_kernels.c      block-compare extend kernel
  wavefront/wavefront_extend_kernels_avx.c  same kernel (AVX2 build)
  wavefront/wavefront_pcigar.c              pcigar extend loads + "MMMMMMMM" stores
  wavefront/wavefront_backtrace.c           "MMMMMMMM" stores

build_tools.sh applies it in build_wfa2 and build_wfmash, whose WFA2-lib is the deps/WFA2-lib submodule, after the source is reset, so it always runs on a pristine tree.
It is idempotent, so a file that already contains the marker is left alone.
If upstream drift breaks an anchor, it exits 1 and leaves every file untouched.

Usage:  wfa2_unaligned_loads.py <path/to/WFA2-lib> Exit:   0 when patched or already patched, 1 on an anchor mismatch with no file modified
"""

import os
import sys
from pathlib import Path

MARKER = "WF_LOAD_U64"

HELPERS = """
/* optik patch (bench/patches/wfa2_unaligned_loads.py): well-defined
 * unaligned 8-byte access; compiles to a single unaligned mov.  Macros (GNU
 * statement expressions), not static functions: the kernels using them are
 * non-static inline functions. */
#include <string.h>
#define WF_LOAD_U64(p) ({ uint64_t wf_ld_; memcpy(&wf_ld_,(p),sizeof(wf_ld_)); wf_ld_; })
#define WF_STORE_U64(p,x) do { const uint64_t wf_st_ = (x); memcpy((p),&wf_st_,sizeof(wf_st_)); } while (0)
"""

# The first version of this patch used static helper functions, which are upgraded in place.
OLD_HELPERS = """
/* optik patch (bench/patches/wfa2_unaligned_loads.py): well-defined
 * unaligned 8-byte access; compiles to a single unaligned mov. */
#include <string.h>
static inline uint64_t wf_load_u64(const void* const p) {
  uint64_t x; memcpy(&x,p,sizeof(x)); return x;
}
static inline void wf_store_u64(void* const p, const uint64_t x) {
  memcpy(p,&x,sizeof(x));
}
"""

# file -> list of old, new, expected occurrence count
EDITS = {
    "wavefront/wavefront_extend_kernels.c": [
        (
            "uint64_t* pattern_blocks = (uint64_t*)(",
            "const char* pattern_blocks = (const char*)(",
            1,
        ),
        (
            "uint64_t* text_blocks = (uint64_t*)(",
            "const char* text_blocks = (const char*)(",
            1,
        ),
        (
            "*pattern_blocks ^ *text_blocks",
            "WF_LOAD_U64(pattern_blocks) ^ WF_LOAD_U64(text_blocks)",
            2,
        ),
        ("++pattern_blocks;", "pattern_blocks += 8;", 1),
        ("++text_blocks;", "text_blocks += 8;", 1),
    ],
    "wavefront/wavefront_extend_kernels_avx.c": [
        (
            "uint64_t* pattern_blocks = (uint64_t*)(",
            "const char* pattern_blocks = (const char*)(",
            1,
        ),
        (
            "uint64_t* text_blocks = (uint64_t*)(",
            "const char* text_blocks = (const char*)(",
            1,
        ),
        (
            "*pattern_blocks ^ *text_blocks",
            "WF_LOAD_U64(pattern_blocks) ^ WF_LOAD_U64(text_blocks)",
            2,
        ),
        ("++pattern_blocks;", "pattern_blocks += 8;", 1),
        ("++text_blocks;", "text_blocks += 8;", 1),
    ],
    "wavefront/wavefront_pcigar.c": [
        (
            "#define CIGAR_8MATCHES_UINT64 *((uint64_t*)matches_lut)",
            "#define CIGAR_8MATCHES_UINT64 WF_LOAD_U64(matches_lut)",
            1,
        ),
        (
            "uint64_t* pattern_blocks = (uint64_t*)(pattern+v);",
            "const char* pattern_blocks = pattern+v;",
            1,
        ),
        (
            "uint64_t* text_blocks = (uint64_t*)(text+h);",
            "const char* text_blocks = text+h;",
            1,
        ),
        (
            "uint64_t pattern_block = *pattern_blocks;",
            "uint64_t pattern_block = WF_LOAD_U64(pattern_blocks);",
            1,
        ),
        (
            "uint64_t text_block = *text_blocks;",
            "uint64_t text_block = WF_LOAD_U64(text_blocks);",
            1,
        ),
        (
            "    pattern_block = *pattern_blocks;",
            "    pattern_block = WF_LOAD_U64(pattern_blocks);",
            1,
        ),
        (
            "    text_block = *text_blocks;",
            "    text_block = WF_LOAD_U64(text_blocks);",
            1,
        ),
        ("++pattern_blocks;", "pattern_blocks += 8;", 1),
        ("++text_blocks;", "text_blocks += 8;", 1),
        (
            "*((uint64_t*)cigar_buffer) = CIGAR_8MATCHES_UINT64;",
            "WF_STORE_U64(cigar_buffer,CIGAR_8MATCHES_UINT64);",
            2,
        ),
    ],
    "wavefront/wavefront_backtrace.c": [
        (
            "*((uint64_t*)(operations+1)) = matches_lut;",
            "WF_STORE_U64(operations+1,matches_lut);",
            1,
        ),
    ],
}


def insert_helpers(src):
    """Insert the helpers right after the last top-level #include."""
    lines = src.split("\n")
    last = max(i for i, line in enumerate(lines) if line.startswith("#include"))
    return "\n".join(lines[: last + 1] + HELPERS.split("\n") + lines[last + 1 :])


def main():
    """Patch every site under the given WFA2-lib directory, or leave all of them alone."""
    if len(sys.argv) != 2:
        print("usage: " + __doc__.split("Usage:")[1].strip(), file=sys.stderr)
        return 1
    root = sys.argv[1]
    out = {}
    for rel, edits in EDITS.items():
        path = os.path.join(root, rel)
        if not os.path.exists(path):
            print(f"wfa2_unaligned_loads: missing {path}", file=sys.stderr)
            return 1
        src = Path(path).read_text()
        if MARKER in src:
            continue  # already patched
        if OLD_HELPERS in src:  # patched by the first version, so upgrade it
            out[path] = (
                src.replace(OLD_HELPERS, HELPERS)
                .replace("wf_load_u64(", "WF_LOAD_U64(")
                .replace("wf_store_u64(", "WF_STORE_U64(")
            )
            continue
        for old, new, n in edits:
            got = src.count(old)
            if got != n:
                print(
                    f"wfa2_unaligned_loads: {rel}: anchor {old!r} found {got}x, expected {n}x "
                    "(upstream drift; nothing modified)",
                    file=sys.stderr,
                )
                return 1
            src = src.replace(old, new)
        out[path] = insert_helpers(src)
    for path, src in out.items():  # every anchor matched, so write
        Path(path).write_text(src)
    print(
        f"wfa2_unaligned_loads: {'patched ' + str(len(out)) + ' file(s)' if out else 'already patched'}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
