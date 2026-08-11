.. _example_mutate:

Create Mutant Structure
-----------------------

Example on how to create a structure file containing hybrid residues for
an alchemical transformation. At the moment, **pmx** supports the mutation of
protein, DNA, and RNA residues.

This is how you can mutate a :download:`protein <peptide.pdb>` residue. ::

    >>> from pmx.model import Model
    >>> from pmx.alchemy import mutate

    >>> # load the PDB file
    >>> m = Model('protein.pdb', rename_atoms=True)
    >>> # perform mutation
    >>> m2 = mutate(m=m, mut_resid=9, mut_resname='R', ff='amber99sb-star-ildn-mut')
    >>> # save mutant PDB file
    >>> m2.write('mutant.pdb')


The ``Model`` can also be modified in place: ::

    >>> mutate(m=m, mut_resid=9, mut_resname='R', ff='amber99sb-star-ildn-mut', inplace=True)
    >>> m.write('mutant.pdb')

Similarly, you can also mutate :download:`DNA <dna.pdb>` with the same :func:`pmx.alchemy.mutate` function. ::

    >>> m = Model('dna.pdb', rename_atoms=True)
    >>> m2 = mutate(m=m, mut_resid=2, mut_resname='T', ff='amber99sb-star-ildn-bsc1-mut')
    >>> m2.write('mutant_dna.pdb')

With the ``verbose`` option you can get additional info printed to screen. ::

    >>> m2 = mutate(m=m, mut_resid=2, mut_resname='T', ff='amber99sb-star-ildn-bsc1-mut', verbose=True)
    log_> Residue to mutate: 2 | DA | A
    log_> Mutation to apply: A->T
    log_> Hybrid residue name: DAT
    log_> Inserted hybrid residue DAT at position 2 (chain A)

Note that, by default, :class:`pmx.model.Model` renumbers all residues from 1.
If you want to keep the original residue IDs, you can do this by using the
``renumber_residues`` argument when loading the ``Model``::

    >>> # load the PDB file
    >>> m = Model('peptide.pdb', rename_atoms=True, renumber_residues=False)

In this case, however, your residue IDs might not be unique anymore. If **pmx**
finds that your residue ID selection is not unique, it will raise an error.
Thus, in this scenario it is safer to also provide the chain ID
to :func:`pmx.alchemy.mutate`, e.g.::

    >>> m2 = mutate(m=m, mut_chain='A', mut_resid=9, mut_resname='R', ff='amber99sb-star-ildn-mut')
