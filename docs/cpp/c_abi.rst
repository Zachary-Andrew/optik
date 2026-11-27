C ABI
=====

The Rust front end reaches the core only through these extern "C" functions.
Each one catches C++ exceptions at the boundary and reports failure through its return value or its output arguments.
Builders and loaders also write a message into a caller-owned buffer, while the mapper's kernels return kStatusFailed or an edit distance past the cap.

ffi_error.hpp
-------------

.. doxygenfile:: ffi_error.hpp

ffi.cpp
-------

.. doxygenfile:: ffi.cpp

ffi_dbg.cpp
-----------

.. doxygenfile:: ffi_dbg.cpp

ffi_fm.cpp
----------

.. doxygenfile:: ffi_fm.cpp

ffi_spss.cpp
------------

.. doxygenfile:: ffi_spss.cpp

ffi_map.cpp
-----------

.. doxygenfile:: ffi_map.cpp
