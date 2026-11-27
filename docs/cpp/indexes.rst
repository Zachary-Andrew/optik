K-mer indexes
=============

Unitigs or any other spectrum-preserving string set feed the two static indexes, SpssCore and FmIndex.
SpssCore buckets super-k-mers by minimizer behind a minimal perfect hash, and DynSpss layers concurrent insert and erase sets over it.
FmIndex answers membership by backward search over the separator-joined unitigs.

spss_core.hpp
-------------

.. doxygenfile:: spss_core.hpp

spss_choose_m.hpp
-----------------

.. doxygenfile:: spss_choose_m.hpp

spss_scan.hpp
-------------

.. doxygenfile:: spss_scan.hpp

spss_io.hpp
-----------

.. doxygenfile:: spss_io.hpp

dyn_spss.hpp
------------

.. doxygenfile:: dyn_spss.hpp

delta_oa.hpp
------------

.. doxygenfile:: delta_oa.hpp

phobic_mphf.hpp
---------------

.. doxygenfile:: phobic_mphf.hpp

fm_index.hpp
------------

.. doxygenfile:: fm_index.hpp
