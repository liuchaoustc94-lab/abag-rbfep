Classes
=======

**pmx** stores structure data in Python classes. A :class:`Model` instance
contains lists of :class:`Chain`, :class:`Molecule`, and :class:`Atom`
instances. A :class:`Chain` instance contains lists of :class:`Molecule` and
:class:`Atom` instances. A :class:`Molecule` instance contains a list of
:class:`Atom` instances only.

Here is a list of the most important classes:

.. toctree::
   :maxdepth: 1

   atom
   molecule
   chain
   model
   trajectory
   topology
   indexfile
