pmx: alchemistry in gromacs
===========================
Toolkit for free-energy calculation setup/analysis and biomolecular structure handling.

**Warning:** this is a development version of pmx, it is not stable or reliable yet.
You are welcome to try/test it and provide feedback, but use at your own risk. The current
stable version of pmx can be found here: https://github.com/dseeliger/pmx

**pmx** is a python library that allows users to setup and analyse molecular
dynamics simulations with the `Gromacs <http://gromacs.org>`_ package.
Among its main features are the setup and analysis of alchemical free energy
calculations for protein, nucleic acid, and small molecule mutations.


.. toctree::
   :maxdepth: 1
   :caption: User Documentation

   installation
   examples/index
   scripts/index
   tutorials/index


.. toctree::
   :maxdepth: 1
   :caption: API Reference

   classes/index
   api/modules

.. toctree::
   :maxdepth: 1
   :caption: Other Tools

   Protein webserver <http://pmx.mpibpc.mpg.de/webserver.html>
   DNA webserver <http://pmx.mpibpc.mpg.de/dna_webserver.html>

.. toctree::
   :maxdepth: 1
   :caption: For Developers

   developers/index

Citations
---------
If you use **pmx** in scientific publications, please cite the following papers:

* V. Gapsys, S. Michielssens, D. Seeliger, B.L. de Groot. `pmx: Automated protein structure and topology generation for alchemical perturbations <http://onlinelibrary.wiley.com/doi/10.1002/jcc.23804/abstract>`_. *J. Comp. Chem* **36** (2015), 348-354
* D. Seeliger, B.L. de Groot. `Protein Thermostability Calculations Using Alchemical Free Energy Simulations <https://www.sciencedirect.com/science/article/pii/S000634951000216X>`_. *Biophys. J.* **98** (2010), 2309-2316

::

    @article{Gapsys2015pmx,
        title = {pmx: Automated protein structure and topology
        generation for alchemical perturbations},
        author = {Gapsys, Vytautas and Michielssens, Servaas
        and Seeliger, Daniel and de Groot, Bert L.},
        journal = {Journal of Computational Chemistry},
        volume = {36},
        number = {5},
        pages = {348--354},
        year = {2015},
        doi = {10.1002/jcc.23804}
    }

    @article{Seeliger2010pmx,
        title = {Protein Thermostability Calculations Using
        Alchemical Free Energy Simulations},
        author = {Seeliger, Daniel and de Groot, Bert L.},
        journal = {Biophysical Journal},
        volume = {98},
        number = {10},
        pages = {2309--2316},
        year = {2010},
        doi = {10.1016/j.bpj.2010.01.051}
    }


License
-------
**pmx** is licensed under the GNU Lesser General Public License v3.0 (LGPL v3).
