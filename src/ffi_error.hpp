#pragma once

/*
 * ffi_error.hpp - error reporting through the C ABI's caller buffers.
 *
 * Functions that can fail take (char* errbuf, uint64_t errbuf_len) and write a NUL-terminated message there.
 * The Rust side reads it with CStr::from_ptr.
 * A null or empty buffer is ignored.
 */

#include <cstdint>
#include <cstdio>

namespace tpoptoa::ffi {

/// set_error() - copy @msg into @errbuf, truncated to fit.
inline void set_error(char* errbuf, uint64_t errbuf_len, const char* msg) {
    if (errbuf && errbuf_len) std::snprintf(errbuf, errbuf_len, "%s", msg);
}

}
