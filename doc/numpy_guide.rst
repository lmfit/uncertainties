.. index: NumPy support

=======================
Uncertainties in arrays
=======================

.. index:: unumpy

The unumpy package
==================

This package contains:

1. utilities that help with the **creation and manipulation** of
NumPy_ arrays and matrices of numbers with uncertainties;

2. **generalizations** of multiple NumPy functions so that they also
work with arrays that contain numbers with uncertainties.

While :ref:`basic operations on arrays <simple_array_use>` that
contain numbers with uncertainties can be performed without it, the
:mod:`unumpy` package is useful for more advanced uses.

Operations on arrays (including their cosine, etc.)  can thus be
performed transparently.

These features can be made available with

>>> from uncertainties import unumpy

.. Here, there is no need to mention unumpy.unlinalg, because it is indeed
   made available through "import unumpy".

Creation and manipulation of arrays and matrices
------------------------------------------------

.. index::
   single: arrays; creation and manipulation
   single: creation; arrays

Arrays
^^^^^^

Arrays of numbers with uncertainties can be built from values and
uncertainties:

>>> arr = unumpy.uarray([1, 2], [0.01, 0.002])
>>> print arr
[1.0+/-0.01 2.0+/-0.002]

NumPy arrays of numbers with uncertainties can also be built directly
through NumPy, thanks to NumPy's support of arrays of arbitrary objects:

>>> arr = numpy.array([ufloat(1, 0.1), ufloat(2, 0.002)])

.. index::
   single: matrices; creation and manipulation
   single: creation; matrices

Matrices
^^^^^^^^

Matrices of numbers with uncertainties are best created in one of
two ways.  The first way is similar to using :func:`uarray`:

>>> mat = unumpy.umatrix([1, 2], [0.01, 0.002])

Matrices can also be built by converting arrays of numbers with
uncertainties into matrices through the :class:`unumpy.matrix` class:

>>> mat = unumpy.matrix(arr)

:class:`unumpy.matrix` objects behave like :class:`numpy.matrix`
objects of numbers with uncertainties, but with better support for
some operations (such as matrix inversion).  For instance, regular
NumPy matrices cannot be inverted, if they contain numbers with
uncertainties (i.e., ``numpy.matrix([[ufloat(…), …]]).I`` does not
work).  This is why the :class:`unumpy.matrix` class is provided: both
the inverse and the pseudo-inverse of a matrix can be calculated in
the usual way: if :data:`mat` is a :class:`unumpy.matrix`,

>>> print mat.I

does calculate the inverse or pseudo-inverse of :data:`mat` with
uncertainties.

.. index::
   pair: nominal value; uniform access (array)
   pair: uncertainty; uniform access (array)
   pair: standard deviation; uniform access (array)

Uncertainties and nominal values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Nominal values and uncertainties in arrays (and matrices) can be
directly accessed (through functions that work on pure float arrays
too):

>>> unumpy.nominal_values(arr)
array([ 1.,  2.])
>>> unumpy.std_devs(mat)
matrix([[ 0.1  ,  0.002]])


.. index:: mathematical operation; on an array of numbers

Mathematical functions
----------------------

This module defines uncertainty-aware mathematical functions that
generalize those from :mod:`uncertainties.umath` so that they work on
NumPy arrays of numbers with uncertainties instead of just scalars:

>>> print unumpy.cos(arr)  # Cosine of each array element

NumPy's function names are used, and not those from the :mod:`math`
module (for instance, :func:`unumpy.arccos` is defined, like in NumPy,
and is not named :func:`acos` like in the :mod:`math` module).

The definition of the mathematical quantities calculated by these
functions is available in the documentation for
:mod:`uncertainties.umath` (which is accessible through :func:`help`
or ``pydoc``).


.. index:: saving to file; array
.. index:: reading from file; array

Storing arrays in text format
=============================

Arrays of numbers with uncertainties can be directly :ref:`pickled
<pickling>`, saved to file and read from a file. Pickling has the
advantage of preserving correlations between errors.

Storing instead arrays in **text format** loses correlations between
errors but has the advantage of being both computer- and
human-readable. This can be done through NumPy's :func:`savetxt` and
:func:`loadtxt`.

Writing the array to file can be done by asking NumPy to use the
*representation* of numbers with uncertainties (instead of the default
float conversion):

>>> numpy.savetxt('arr.txt', arr, fmt='%r')

This produces a file `arr.txt` that contains a text representation of
the array::

  1.0+/-0.01
  2.0+/-0.002

The file can then be read back by instructing NumPy to convert all the
columns with :func:`uncertainties.ufloat_fromstr`. The number
:data:`num_cols` of columns in the input file (1, in our example) must
be determined in advance, because NumPy requires a converter for each
column separately:

>>> converters = dict.fromkeys(range(num_cols), uncertainties.ufloat_fromstr)
>>> arr = numpy.loadtxt('arr.txt', converters=converters, dtype=object)

.. index:: linear algebra; additional functions, ulinalg

Additional array functions: unumpy.ulinalg
==========================================

The :mod:`unumpy.ulinalg` module contains more uncertainty-aware
functions for arrays that contain numbers with uncertainties.

It currently offers generalizations of two functions from
:mod:`numpy.linalg` that work on arrays (or matrices) that contain
numbers with uncertainties, the **matrix inverse and pseudo-inverse**:

>>> unumpy.ulinalg.inv([[ufloat(2, 0.1)]])
array([[0.5+/-0.025]], dtype=object)
>>> unumpy.ulinalg.pinv(mat)
matrix([[0.2+/-0.0012419339757],
        [0.4+/-0.00161789987329]], dtype=object)

.. _NumPy: http://numpy.scipy.org/
