Molecule
========

This class stores arrays of atoms. It is also used for residues within
polymers, such as amino acids within a protein. As such, a :class:`Molecule`
instance has :attr:`id`, :attr:`resname`, and :attr:`atoms` attributes.


API Reference
-------------

.. autoclass:: pmx.molecule.Molecule
   :exclude-members: assign_moltype, has_atom, new_aa, get_real_resname,
                     get_psi, set_psi, set_psi_down, get_phi, set_phi, set_phi_down, get_omega, set_omega, set_omega_down, nchi, get_chi, set_chi,
                     set_conformation, set_orig_resid, set_chain, set_molecule, set_chain_id,
                     copy, is_protein_residue, get_mol2_types, get_bonded, atomlistFromTop, com

   .. rubric:: Methods

   .. autosummary::
      set_resname
      set_resid
      insert_atom
      append
      remove_atom
      fetch
      fetchm
      is_hybrid
      rename_atoms_to_gmx
