Chain
=====

This class stores arrays of residues/molecules. Thus, a :class:`Chain`
instance has :attr:`residues` and :attr:`atoms` attributes.

API Reference
-------------

.. autoclass:: pmx.chain.Chain
   :exclude-members: com, atomlistFromTop, copy, create, cterminus, nterminus,
                     get_bonded, insert_chain, make_residue_tree, residue

   .. rubric:: Methods

   .. autosummary::
      set_chain_id
      get_sequence
      append
      insert_residue
      remove_residue
      replace_residue
      nbuild
      cbuild
      add_nterm_cap
      add_cterm_cap
      fetch_residues
      rename_atoms_to_gmx
