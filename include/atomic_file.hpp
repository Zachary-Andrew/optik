#pragma once

/**
 * @file
 * @brief Write-then-rename file output.
 *
 * Output goes to `<path>.tmp.<pid>` and is renamed over `<path>` only after a clean flush and close.
 * A full disk often surfaces only at that last flush, which is why success is judged after close().
 * If anything fails, the temp file is removed and `<path>` is left untouched.
 * As a result, in-place updates (-i X -o X) never truncate the input while it is still being read.
 */

#include <cstdio>
#include <fstream>
#include <string>

#include <unistd.h>

namespace tpoptoa {

/// Runs @p write on a stream for the temp file, then renames it over @p path.
/// A false return means nothing at @p path changed.
template <typename Writer>
bool atomic_write(const char* path, Writer&& write) {
    const std::string tmp = std::string(path) + ".tmp." + std::to_string(::getpid());
    bool written = false;

    {
        std::ofstream out(tmp, std::ios::binary);

        if (out) {
            write(out);
            out.flush();
            written = static_cast<bool>(out);
            out.close();
            written = written && !out.fail();
        }
    }
    if (written && std::rename(tmp.c_str(), path) == 0) return true;

    std::remove(tmp.c_str());
    return false;
}

}
