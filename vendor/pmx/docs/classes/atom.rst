Atom
====

This class stores all properties (attributes) related to an atom.
When **pmx** reads a PDB file, each line starting with "ATOM" or "HETATM"
is converted to an :class:`Atom` instance. In addition, the atom class also has a number
of methods that can be used to calculate distances, angles and dihedrals.

API Reference
-------------

.. currentmodule:: pmx.atom

.. autoclass:: Atom
   :exclude-members: readPDBString, read_mol2_line, make_long_name, copy, get_symbol, get_order, set_resname, set_chain_id

   .. rubric:: Methods

   .. autosummary::
      dist
      dist2
      translate
      angle
      dihedral
      nm2a
      a2nm
