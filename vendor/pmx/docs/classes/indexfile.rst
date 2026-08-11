IndexFile
=========

Index files are frequently used in Gromacs to select subsets of a simulation
system. In **pmx** you can read, change and write index files with the
:class:`IndexFile` class.

API Reference
-------------

.. autoclass:: pmx.ndx.IndexGroup

.. autoclass:: pmx.ndx.IndexFile

   .. rubric:: Methods

   .. autosummary::
      add_group
      delete_group
      parse
      write
