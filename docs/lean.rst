Lean proofs
===========

``lean/`` proves optik's exactness and MAPQ claims in Lean 4, about an executable model of each algorithm.
The binary and the Rust tests must reproduce the test vectors each model writes.

``docs/lean_docs.py`` builds this page from the Lean sources.
Each module's section opens with its ``/-!`` header and lists the declarations that have docstrings, with signatures up to ``:=``.
The proofs themselves are in the ``.lean`` file linked from each section.

Building and checking
---------------------

.. include:: ../README.md
   :parser: myst_parser.sphinx_
   :start-after: <!-- lean-setup -->
   :end-before: <!-- /lean-setup -->

Three controls in ``bench/negative_controls.py`` corrupt the binary's unitigs, GFA links and presence answers, and the vector checks must report each one.

.. include:: ../lean/APPENDIX_B.md
   :parser: myst_parser.sphinx_

Modules
-------

Modules are listed in pipeline order, starting with k-mers and ending with SAM fields.

.. lean-module:: Kmer
.. lean-module:: Hash
.. lean-module:: Rolling
.. lean-module:: Sort
.. lean-module:: Minimizer
.. lean-module:: SeedIndex
.. lean-module:: Mphf
.. lean-module:: Spss
.. lean-module:: Dyn
.. lean-module:: Varint
.. lean-module:: Count
.. lean-module:: Dbg
.. lean-module:: Ef
.. lean-module:: Fm
.. lean-module:: Bpm
.. lean-module:: Cigar
.. lean-module:: Align
.. lean-module:: Mapq
.. lean-module:: Paired
