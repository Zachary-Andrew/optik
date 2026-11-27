Input and k-mer counting
========================

Records arrive through the streaming FASTA/FASTQ reader and become 2-bit packed k-mers.
Pass 1 cuts the input into super-k-mers and files them by minimizer into partition buffers, which spill to disk once they pass the memory budget.
Pass 2 counts each partition in a private cache-sized table.
Every super-k-mer keeps its flanking bases, so edges are counted alongside the k-mers and no separate (k+1)-mer pass is needed.

fasta_reader.hpp
----------------

.. doxygenfile:: fasta_reader.hpp

kmer.hpp
--------

.. doxygenfile:: kmer.hpp

kmer_extractor.hpp
------------------

.. doxygenfile:: kmer_extractor.hpp

minimizer.hpp
-------------

.. doxygenfile:: minimizer.hpp

kmer_count.hpp
--------------

.. doxygenfile:: kmer_count.hpp

kmer_count.cpp
--------------

.. doxygenfile:: kmer_count.cpp

count_planning.hpp
------------------

.. doxygenfile:: count_planning.hpp

part_count_table.hpp
--------------------

.. doxygenfile:: part_count_table.hpp

spill.hpp
---------

.. doxygenfile:: spill.hpp
