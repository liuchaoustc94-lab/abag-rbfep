Model
=====

The :class:`Model` class stores the structural data of your system, and as such it contains
all :class:`Chain`, :class:`Molecule`, and :class:`Atom` instances of your system.
It is possibly the most useful class for a user in practice, as the function that generates
hybrid systems to be used in alchemical calculations (:func:`pmx.alchemy.mutate`)
acts on this object.

API Reference
-------------

.. autoclass:: pmx.model.Model
   :exclude-members: com, writeFASTA, writePIR, writePDB, atomlistFromTop, append, assign_moltype,
                     make_chains, make_residues

   .. rubric:: Methods

   .. autosummary::
      read
      write
      renumber_residues
      renumber_atoms
      insert_residue
      insert_chain
      remove_atom
      remove_residue
      remove_chain
      replace_residue
      fetch_residue
      fetch_residues
