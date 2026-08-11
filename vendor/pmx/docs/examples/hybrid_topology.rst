Generate Hybrid Topology
------------------------

Example on the use of **pmx** to prepare an hybrid topology. This is after running
pdb2gmx on the mutant structure and obtaining a topology file from Gromacs.

The simplest case is the following (:download:`topol.top <topol.top>`), where
there are no itp files of other molecules included in the topology.
::

    >>> from pmx.forcefield import Topology
    >>> from pmx.alchemy import gen_hybrid_top

    >>> # load the topology file
    >>> top = Topology('topol.top')
    >>> # fill B states for hybrid residues present
    >>> pmxtop, _ = gen_hybrid_top(topol=top)
    >>> # write topology to a new file
    >>> pmxtop.write('pmxtop.top')


If the topology file includes different itp files (:download:`topol2.top <topol2.top>`,
:download:`topol2.itp <topol2.itp>`),
these too can be modified and returned by setting the ``recursive`` argument to ``True``. ::

    >>> # load the topology file
    >>> top = Topology('topol2.top')
    >>> # fill B states for hybrid residues present in topol2.itp
    >>> pmxtop, pmxitps = gen_hybrid_top(topol=top, recursive=True)
    >>> # write topology and itps to a new file
    >>> pmxtop.write('pmxtop.top')
    >>> for i, itp in enumerate(pmxitps):
    >>>     itp.write(itp.filename)

The above will overwrite the original itp files. If you want to save the new pmx itp file
without losing the original files, you can do as follows::

    >>> for i, itp in enumerate(pmxitps):
    >>>     itp.write('pmxitp_%s.itp' % i)

You will then need to update the ``pmxtop.top`` accordingly so to include the
new itp file.

If you like getting lots of screen output, you can set the ``verbose`` argument
to ``True``. ::

    >>> pmxtop, pmxitps = gen_hybrid_top(topol=top, recursive=True, verbose=True)
    log_> Reading input top file "topol2.top"
    log_> Making bonds for state B -> 0 bonds with perturbed atoms
    log_> Making angles for state B -> 0 angles with perturbed atoms
    log_> Making dihedrals for state B -> 0 dihedrals with perturbed atoms
    log_> Removed 0 fake dihedrals
    log_> Total charge of state A = 0
    log_> Total charge of state B = 0

    log_> Reading input itp file "topol2.itp""
    log_> Scanning database for D2A
    log_> Scanning database for G2R
    log_> Scanning database for P2A
    log_> Hybrid Residue -> 9 | D2A
    log_> Hybrid Residue -> 10 | G2R
    log_> Hybrid Residue -> 12 | P2A
    log_> Making bonds for state B -> 61 bonds with perturbed atoms
    log_> Making angles for state B -> 127 angles with perturbed atoms
    log_> Making dihedrals for state B -> 190 dihedrals with perturbed atoms
    log_> Removed 19 fake dihedrals
    log_> Total charge of state A = 1
    log_> Total charge of state B = 3
