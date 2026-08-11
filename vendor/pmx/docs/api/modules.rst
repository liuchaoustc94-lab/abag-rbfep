Modules
=======

**pmx** is a essentially a collection of classes and functions that allow the
manipulation of molecular structure and topology files. It has been thought to
be mainly used in conjunction with the Gromacs simulation package. In fact,
pmx makes use of some Gromacs functions to read/write trajectory data and to
allow fast neighborsearching from a python script. With time, **pmx** has also
become a library for alchemical free energy  calculations.
Here below are the main Modules that the user might
want to use to facilitate the setup and analysis of alchemical free energy
calculations in Gromacs.

pmx.alchemy
-----------

.. automodule:: pmx.alchemy
    :members:
    :undoc-members:
    :show-inheritance:

pmx.analysis
------------

.. automodule:: pmx.analysis
    :members:
    :undoc-members:
    :show-inheritance:

pmx.estimators
--------------

.. automodule:: pmx.estimators
    :members:
    :undoc-members:
    :show-inheritance:

pmx.gmx
-------

.. automodule:: pmx.gmx
    :members:
    :undoc-members:
    :show-inheritance:

pmx.builder
-----------

.. automodule:: pmx.builder
    :members:
    :undoc-members:
    :show-inheritance:

pmx.utils
---------

.. automodule:: pmx.utils
    :members:
    :undoc-members:
    :show-inheritance:

pmx.atomselection
-----------------

.. automodule:: pmx.atomselection
    :members:
    :undoc-members:
    :show-inheritance:

pmx.ffparser
------------

.. automodule:: pmx.ffparser
    :members:
    :undoc-members:
    :show-inheritance:

pmx.geometry
------------

.. automodule:: pmx.geometry
    :members:
    :undoc-members:
    :show-inheritance:

pmx.library
-----------

.. automodule:: pmx.library
    :members:
    :undoc-members:
    :show-inheritance:

pmx.mutdb
---------

.. automodule:: pmx.mutdb
    :members:
    :undoc-members:
    :show-inheritance:

pmx.parser
----------

.. automodule:: pmx.parser
    :members:
    :undoc-members:
    :show-inheritance:

pmx.rotamer
-----------

.. automodule:: pmx.rotamer
    :members:
    :undoc-members:
    :show-inheritance:

pmx.xdrfile
-----------

.. automodule:: pmx.xdrfile
    :members:
    :undoc-members:
    :show-inheritance:

Notes
-----
.. [1] Jarzynski C (1997) Equilibrium free-energy
        differences from nonequilibrium measurements:
        A master-equation approach.
        Phys Rev E 56:5018–5035

.. [2] Goette M, Grubmüller H (2009) Accuracy
        and convergence of free energy differences calculated
        from nonequilibrium switching processes.
        J Comput Chem 30(3):447–456

.. [3] Shirts MR, Bair E, Hooker G, Pande VS
        (2003) Equilibrium free energies from non-equilibrium
        measurements using maximum-likelihood methods.
        Phys Rev Lett 91 (14):140601

.. [4] Hahn AM, Then H (2010) Measuring the
        convergence of Monte Carlo free-energy calculations.
        Phys Rev E 81(4):041117

.. [5] Boresch S, Tettinger F, Leitgeb M, Karplus M (2003)
       Absolute Binding Free Energies:  A Quantitative Approach for
       Their Calculation.
       J Phys Chem B 107(35):9535-9551

.. [6] Hummer, G (2001)
       Fast-growth thermodynamic integration: Error and efficiency analysis
       J. Chem. Phys. 114
