Gentop
-------

The ``pmx gentop`` script can be used to fill out the B states for a topology
file containing hybrid residues. You can obtain such topology by running
the Gromacs command **pdb2gmx** on a structure file generated with the
:ref:`pmx mutate <script_mutate>` command. ::

    $ pmx gentop -h
    usage: pmx [-h] [-p topol] [-o outfile] [-ff ff] [--split] [--scale_mass]
               [--norecursive]

    This script fills in the B state to a topology file (itp or top) according to
    the hybrid residues present in the file. If you provide a top file with
    include statemets, by default the script will run through the included itp
    files too; this can turned off using the --norecursive flag. You need to use
    this script after having mutated a structure file with pmx mutate, and after
    having passed that mutated structure through pdb2gmx.

    optional arguments:
      -h, --help     show this help message and exit
      -p topol       Input topology file (itp or top). Default is "topol.top"
      -o outfile     Output topology file. Default is "pmxtop.top"
      -ff ff         Force field to use. If -p is a top file, it is not necessary
                     to specify the forcefield, as it will be determined
                     automatically. If -p is an itp file, then -ff is needed, and
                     if it not provided a list of available ff will be shown.
      --split        Write separate topologies for the vdW and charge
                     transformations.
      --scale_mass   Scale the masses of morphing atoms so that dummies have a
                     mass of 1.
      --norecursive  Whether to fill the B states also for all itp files included
                     in the provided topology file. Default is True. This flag
                     sets it to False.
