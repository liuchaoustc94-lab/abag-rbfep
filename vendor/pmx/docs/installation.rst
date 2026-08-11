Installation
============

Clone the source code from GitHub with::

    $ git clone https://github.com/deGrootLab/pmx.git

Then::

    $ cd pmx
    $ pip install .

Software requirements:

* python >= 3.6
* numpy_
* scipy_
* matplotlib_

GMXLIB
------
In order to be able to use the hybrid/alchemical force fields available in
**pmx**, the environment variable ``$GMXLIB`` needs to be set. You can either set
``$GMXLIB`` in each of you shell sessions when you want to use **pmx**, or you
can set the variable directly in your ~/.bashrc file.

The ``pmx gmxlib`` command will tell you where the force field libraries
were installed, and optionally allow you to automatically add ``$GMXLIB``
to your bashrc file. ::

    $ pmx gmxlib

    In order to be able to use the hybrid/alchamical force fields
    available in pmx, the environment variable GMXLIB needs to be set.

    The path to your pmx force field library is the following:
    /path/to/your/pmx/data/mutff

    You can either set GMXLIB in your shell session as follows:
    $ export GMXLIB=/path/to/your/pmx/data/mutff

    Or you can add this directly in your bashrc file.

    Do you wish pmx to set the GMXLIB variable in your ~/.bashrc? [yes|no]
    >>>

Alternatively, you can set the variable only for a certain python process using
the ``set_gmxlib()`` function::

    >>> from pmx.gmx import set_gmxlib
    >>> set_gmxlib()


.. _numpy: http://www.numpy.org/
.. _scipy: https://www.scipy.org/
.. _matplotlib: https://matplotlib.org/
