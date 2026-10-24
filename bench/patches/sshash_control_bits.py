#!/usr/bin/env python3
"""sshash_control_bits.py - idempotent patch that widens SSHash's control codewords to fit the heavy-load codeword of the skew index.

Why this patch exists
Each minimizer's control codeword (src/builder/build_sparse_and_skew_index.cpp) is one of three forms.

    singleton  offset << 1                                 | 0
    mid-load   ((list_id << min_l) | (bucket_size - 2)) << 2 | 01
    heavy-load ((heavy_offset << 3) | partition_id)    << 2 | 11

The width is computed as

    max(num_bits_per_offset + 1, 2 + min_l + bits_for_list_id)

which covers the first two forms only.
heavy_offset is the start of the bucket among all skew-index positions, so the third form needs 2 + 3 + ceil(log2(skew-index positions + 1)) bits.
When most positions sit in buckets larger than 2^min_l, as with a short m against a large set, that is more than the width and the codeword is truncated.
The `assert(code < ...)` at that line fires in a build without NDEBUG, while a Release build stores the truncated value.
`build --check` then reports "wrong id assigned" and queries can read out of bounds (on chr1 eulertigs at k = 15 every m <= 11 fails, queries segfault at m = 11, and at m = 7 the width is 29 bits where the heavy codeword needs more).
Threads and the -g budget make no difference.

The fix adds the third form to the max.
Decoding is unchanged, since the query side reads the width from the stored compact vector and takes the low 2, then 3 bits, then the rest (include/sparse_and_skew_index.hpp), so a wider codeword decodes the same way.
Where the old width already sufficed the max is the same, so every index that built correctly before is built identically.

build_tools.sh applies it in build_sshash after the source is reset.
It is idempotent, so a file that already contains the marker is left alone.
If upstream drift breaks the anchor, it exits 1 and leaves the file untouched.

Usage:  sshash_control_bits.py <path/to/sshash> Exit:   0 when patched or already patched, 1 on an anchor mismatch with the file unmodified
"""

import os
import sys
from pathlib import Path

MARKER = "optik patch (bench/patches/sshash_control_bits.py)"
FILE = os.path.join("src", "builder", "build_sparse_and_skew_index.cpp")

OLD = """    const uint64_t num_bits_for_control =
        std::max(num_bits_per_offset + 1, 2 + constants::min_l + bits_for_list_id);
"""
NEW = """    // optik patch (bench/patches/sshash_control_bits.py): the heavy-load
    // codeword ((offset << 3) | partition_id) << 2 | 11 must fit as well
    const uint64_t bits_for_heavy =
        2 + 3 + uint64_t(std::ceil(std::log2(num_minimizer_positions_of_buckets_in_skew_index + 1)));
    const uint64_t num_bits_for_control =
        std::max({num_bits_per_offset + 1, 2 + constants::min_l + bits_for_list_id, bits_for_heavy});
"""


def main():
    """Patch SSHash's builder in the source tree given as the only argument."""
    if len(sys.argv) != 2:
        sys.exit("usage: " + __doc__.split("Usage:")[1].strip())
    path = os.path.join(sys.argv[1], FILE)
    src = Path(path).read_text()
    if MARKER in src:
        print(f"{path}: already patched")
        return 0
    if src.count(OLD) != 1:
        print(
            f"{path}: anchor not found (upstream drift?); not modified", file=sys.stderr
        )
        return 1
    with open(path, "w") as f:
        f.write(src.replace(OLD, NEW))
    print(f"{path}: patched")
    return 0


if __name__ == "__main__":
    sys.exit(main())
