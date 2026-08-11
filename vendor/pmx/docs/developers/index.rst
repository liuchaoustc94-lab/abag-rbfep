Developers Guide
================

Here are some general information on the current structure of the repository, and
how things like documentation and testing are setup. It is a good point where to start
for who would like to contribute to the code, with some links to useful websites
and resources.

Structure
---------
The repository is structured as follows::

    pmx
    |
     -- src/
    |   |
    |    -- pmx/
    |       |
    |        -- data/
    |       |
    |        -- scripts/
    |       |
    |        -- extensions/
    |
     -- tests/
    |   |
    |    -- data/
    |
     -- docs/
    |
     -- setup.py

The ``src/pmx`` folder contains the actual library. In this folder there are all the
**pmx** modules (e.g. ``model.py``, ``forcefield.py``, ``estimators.py``). See the
following blogs on the advantages of such structure:

|  https://blog.ionelmc.ro/2014/05/25/python-packaging/#the-structure
|  https://hynek.me/articles/testing-packaging/
|

The folder ``src/pmx/data`` contains the files used by **pmx** or Gromacs. For instance,
here are stored the forcefield libraries (in ``data/mutff``), and a number of pickled
files used by **pmx**. The folder ``src/pmx/scripts`` contains all scripts that are meant to
be called from the command line. The file ``src/pmx/scripts/cli.py`` defines the scripts that are
part of the command line interface. The folder ``src/pmx/extensions`` contains
extension written in C, which are compiled during installation.

The folder ``tests`` contains all the tests for the library. The folder ``docs`` contains the documentation.
More on these below.

Style
-----
We do not follow strict style rules, but try to adhere to the PEP8 guide:

|  https://www.python.org/dev/peps/pep-0008/
|

A suggestion that makes life easier with this: use a code editor or IDE that supports
syntax highlighting, code formatting, linting (automatically checks the code for errors),
and allows running the code from within the environment (like in a Jupyter notebook). Examples
are `Atom <https://atom.io/>`_ and `PyCharm <https://www.jetbrains.com/pycharm/>`_. For more,
have a look the following post:

|  https://realpython.com/python-ides-code-editors-guide/
|

As another general suggestion: use functions whenever possible, classes only when necessary.
It makes sense to create classes to, *e.g.*, hold data that can be manipulated in different ways (a PDB file for instance),
or sometimes to group a bunch of related functions. But if you only need to *do* something, then
functions might be all that is needed. Generally it is easier to start writing functions, then,
if you find yourself writing functions inside of functions, it might be time to
write a class instead. On the other hand, if you only have one function in a class, then it might
be better to stick with just a function.

|  https://www.reddit.com/r/Python/comments/1wcc18/what_are_the_advantages_to_using_classes_vs_not/
|

Documentation
-------------
Documentation for the library is stored in the ``docs`` folder and it is build
using ``sphinx``, which converts **reStructuredText** markup language into a range
of output formats, including HTML and LaTeX (for printable PDF versions).

|  http://www.sphinx-doc.org/en/master/
|

The API documentation is built automatically from the source files by reading the
**docstrings** written in in **NumPy** style. Examples can be found here as well as directly in the code:

|  https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html
|

The documentation is published online with GitHub pages. It is
deployed automatically from `Travis CI <https://travis-ci.org/>`_ using ``doctr``. In practice, after testing,
Travic CI builds the documentation and sends it to the branch ``gh-pages`` so to make it
available at this web address.

Info can be found here:

|  https://pages.github.com/
|  https://blog.github.com/2016-08-22-publish-your-project-documentation-with-github-pages/
|  https://docs.travis-ci.com/user/deployment/pages/
|  https://drdoctr.github.io/doctr/
|

To build the documentation locally, in the ``docs/build`` folder, you do the following
from within the ``docs`` folder::

  make html

or ::

  make clean html

Testing
-------
For testing we use ``pytest``:

|  https://docs.pytest.org/en/latest/
|

All tests are in the ``tests`` folder. There you will see there is a Python file
for each module in **pmx**, with the prefix "test\_": *e.g.* the tests for the module
``model.py`` are defined in ``tests/test_model.py``. Similarly, test functions are
specified with the same prefix; thus, in ``tests/test_model.py`` you can
find a function called for example ``test_model_moltype``.

It gets a bit more complicated when having to load/write files, and this is handled
by the ``gf`` function defined in ``tests/conftest.py``. ``Pytest`` also allows to create
test classes (that start with ``Test``, e.g. ``TestClass``), and has additional
features like **fixtures** and **parametrized functions**.
A good intro to testing with ``pytest`` is the following:

|  https://semaphoreci.com/community/tutorials/testing-python-applications-with-pytest
|

To run all tests and get a summary report you just type "pytest" in the root
directory of the repo::

  $ pytest
  ============================= test session starts ==============================
  platform linux2 -- Python 2.7.13, pytest-3.4.1, py-1.5.2, pluggy-0.6.0
  rootdir: /home/maldegh/sw/degrootlab-pmx-develop, inifile:
  plugins: cov-2.5.1
  collected 11 items

  tests/test_chain.py .                                                    [  9%]
  tests/test_estimators.py ...                                             [ 36%]
  tests/test_import.py .                                                   [ 45%]
  tests/test_model.py .....                                                [ 90%]
  tests/test_utils.py .                                                    [100%]

  ========================== 11 passed in 0.14 seconds ===========================

Continuous Integration
----------------------
Every time changes are pushed to the repository, the code is built and tested automatically
by `Travis CI <https://travis-ci.org/>`_.

|  https://travis-ci.org/deGrootLab/pmx
|

So if the new commits break the build or fail the tests
you should get notified. The configuration file for Travis CI is ``.travis.yml``.

The repo is also linked to `Codecov <https://codecov.io/>`_ so to see which parts
of the code are covered by tests.

|  https://codecov.io/gh/deGrootLab/pmx
|

Versioning
----------
We follow the following standard sequence-based scheme::

  MAJOR.MINOR[.MICRO]

where ``MAJOR`` designates a major revision number for the software (e.g. ``2`` or ``3`` for Python).
This is done when adding a lot of features, breaking backward-compatibility,
or drastically changing the API.

``MINOR`` releases involve moderate changes,
like bug fixes and minor improvements. End users should be able to upgrade without worrying that their code
will not work anymore. If there are small API changes, the user should be notified with a deprecation warning.

Sometimes the third level ``MICRO`` can be used, primarily
for bug fixes. Development versions can just be suffixed by a *dev* number,
e.g. *1.2dev0*, *1.2dev1*, *1.2dev2*, *etc*.

More info can be found here:

|  https://the-hitchhikers-guide-to-packaging.readthedocs.io/en/latest/specification.html#sequence-based-scheme
|

Tracking versions is handled by ``versioneer``.

|  https://github.com/warner/python-versioneer
|  https://blog.mozilla.org/warner/2012/01/31/version-string-management-in-python-introducing-python-versioneer/
|

Which makes versioning easy by using Git tags:

|  https://git-scm.com/book/en/v2/Git-Basics-Tagging
|

In practice, if you have modified the code and would like to tag it with a new version::

  $ git tag -a 1.4.3 -m "my version 1.4.3"

To see all versions of the code available::

  $ git tag
  1.0dev0
  1.0
  1.1
  1.2
  1.3
  1.4
  1.4.1
  1.4.2
  1.4.3

To see the tag data along with the commit that was tagged::

  $ git show 1.4.3
  tag 1.4.3
  Tagger: John Doe <john@doe.com>
  Date:   Tue Jan 1 12:01:01 2018 -0700

  my version 1.4.3

  commit ca82a6dff817ec66f44342007202690a93763949
  Author: Mark Smith <mark@smith.com>
  Date:   Mon Dec 17 21:52:11 2017 -0700

      changed the version number
