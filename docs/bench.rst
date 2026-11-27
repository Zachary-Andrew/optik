Bench helpers
=============

The bench scripts in ``bench/`` share these modules, so every validator reads files, times commands and reports failures the same way.
Their pages come straight from the docstrings.

Sequences and formats
---------------------

.. automodule:: dna

.. automodule:: sam_spec

Measurement
-----------

.. automodule:: measure

.. automodule:: memcap

Downloads
---------

.. automodule:: rangefetch

Shared values and checks
------------------------

.. automodule:: units

.. literalinclude:: ../bench/units.py
   :language: python
   :start-at: BYTES_PER_KIB

.. automodule:: diagnostics

.. literalinclude:: ../bench/diagnostics.py
   :language: python
   :start-at: # A check's

.. automodule:: checklog
