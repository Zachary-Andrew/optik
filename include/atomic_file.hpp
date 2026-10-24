#pragma once

/*
 * atomic_file.hpp - write-then-rename file output.
 *
 * The content goes to "<path>.tmp.<pid>", which is flushed and closed before it is renamed over <path>.
 * A full disk often shows up only at that final flush, so success is decided after close.
 * On any failure the temporary is removed and <path> is left as it was, so there is never a truncated output, and an in-place update (-i X -o X) cannot truncate the input it is still reading.
 */

#include <cstdio>
#include <fstream>
#include <string>

#include <unistd.h>

namespace tpoptoa {

/// atomic_write() - call @write on a stream for "<path>.tmp.<pid>" and rename the result over @path.
///
/// Returns true only if the stream stayed good through flush and close and the rename succeeded.
template <typename Writer>
bool atomic_write(const char* path, Writer&& write) {
    const std::string tmp = std::string(path) + ".tmp." + std::to_string(::getpid());
    bool ok = false;
    {
        std::ofstream os(tmp, std::ios::binary);
        if (os) {
            write(os);
            os.flush();
            ok = static_cast<bool>(os);
            os.close();
            ok = ok && !os.fail();
        }
    }
    if (ok && std::rename(tmp.c_str(), path) == 0) return true;
    std::remove(tmp.c_str());
    return false;
}

}
