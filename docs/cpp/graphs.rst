De Bruijn graphs and unitigs
============================

Counted k-mers and their edge masks end up in one of three graph forms, picked by k and by whether a plain table fits the memory limit.
DbgGraph keeps nodes in one open-addressing table, DbgGraphEF stores sorted keys as Elias-Fano, and partitioned compaction builds unitigs without any global table.
Every form emits the same maximal unitigs in the same order.

dbg_graph.hpp
-------------

.. doxygenfile:: dbg_graph.hpp

dbg_graph_ef.hpp
----------------

.. doxygenfile:: dbg_graph_ef.hpp

dbg_part.hpp
------------

.. doxygenfile:: dbg_part.hpp

dbg_unitigs.hpp
---------------

.. doxygenfile:: dbg_unitigs.hpp
