Read mapping and alignment
==========================

The mapper seeds reads with reference minimizers from SeedIndex, verifies short-read candidates four at a time with bit-parallel edit distance, and aligns with gap-affine wavefronts.
BiWFA brings alignment memory down to O(s) for long gaps.

seed_index.hpp
--------------

.. doxygenfile:: seed_index.hpp

bpm_edit.hpp
------------

.. doxygenfile:: bpm_edit.hpp

wfa_affine.hpp
--------------

.. doxygenfile:: wfa_affine.hpp
