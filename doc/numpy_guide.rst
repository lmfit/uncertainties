.. index: NumPy support

===============================
Uncertainties and numpy arrays
===============================

.. index:: unumpy
.. index:: arrays; simple use, matrices; simple use

.. _simple_array_use:

Arrays of uncertainties Variables
====================================

It is possible to put uncertainties Variable  in NumPy_ arrays and
matrices:

>>> arr = numpy.array([ufloat(1, 0.01), ufloat(2, 0.1)])
>>> 2*arr
[2.0+/-0.02 4.0+/-0.2]
>>> print arr.sum()
3.00+/-0.10

Many common operations on NumPy arrays can be performed transparently
even when these arrays contain numbers with uncertainties.


The unumpy package
==================


While :ref:`basic operations on arrays <simple_array_use>` that
contain numbers with uncertainties can be performed without it, the
:mod:`unumpy` package is useful for more advanced uses.

This package contains:

1. utilities that help with the **creation and manipulation** of
NumPy_ arrays and matrices of numbers with uncertainties;

2. **generalizations** of multiple NumPy functions so that they also
work with arrays that contain numbers with uncertainties.


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
>>> print(arr)
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

>>> print(mat.I)

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

>>> print(unumpy.cos(arr))  # Cosine of each array element

NumPy's function names are used, and not those from the :mod:`math`
module (for instance, :func:`unumpy.arccos` is defined, like in NumPy,
and is not named :func:`acos` like in the :mod:`math` module).

The definition of the mathematical quantities calculated by these
functions is available in the documentation for  :mod:`uncertainties.umath`.

.. index::
   pair: testing and operations (in arrays); NaN

NaN testing and NaN-aware operations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One particular function pertains to NaN testing: ``unumpy.isnan()``. It
returns true for each NaN *nominal value* (and false otherwise).

Since NaN±1 is *not* (the scalar) NaN, functions like
``numpy.nanmean()`` do not skip such values. This is where
``unumpy.isnan()`` is useful, as it can be used for masking out numbers
with a NaN nominal value:

>>> nan = float("nan")
>>> arr = numpy.array([nan, uncertainties.ufloat(nan, 1), uncertainties.ufloat(1, nan), 2])
>>> arr
array([nan, nan+/-1.0, 1.0+/-nan, 2], dtype=object)
>>> arr[~unumpy.isnan(arr)].mean()
1.5+/-nan

or equivalently, by using masked arrays:

>>> masked_arr = numpy.ma.array(arr, mask=unumpy.isnan(arr))
>>> masked_arr.mean()
1.5+/-nan

In this case the uncertainty is NaN as it should be, because one of
the numbers does have an undefined uncertainty, which makes the final
uncertainty undefined (but the average is well defined). In general,
uncertainties are not NaN and one obtains the mean of the non-NaN
values.

.. index:: saving to file; array
.. index:: reading from file; array

Storing arrays in text format
=============================

Arrays of numbers with uncertainties can be directly :ref:`pickled
<pickling>`, saved to file and read from a file. Pickling has the
advantage of preserving correlations between errors.

Storing arrays in **text format** loses correlations between errors but has the
advantage of being both computer- and human-readable. This can be done through
NumPy's :func:`savetxt` and :func:`loadtxt`.

Writing the array to file can be done by asking NumPy to use the
*representation* of numbers with uncertainties (instead of the default float
conversion):

>>> numpy.savetxt('arr.txt', arr, fmt='%r')

This produces a file `arr.txt` that contains a text representation of
the array::

  1.0+/-0.01
  2.0+/-0.002

The file can then be read back by instructing NumPy with :meth:`numpy.loadtxt`,
but for object arrays, this requires a converter function for each column
separately.  We can use func:`uncertainties.ufloat_fromstr`, but
:meth:`numpy.loadtxt` passes bytes to converters, they must first be converted
into a string.  In addition the number of maximum number of columns must be
known.  An example of using all of this to unpack the data saved with
:meth:`numpy.savetxt` would be:

>>> from uncertainties import ufloat_fromstr
>>> max_cols = 1
>>> converters = {col: lambda dat: ufloat_fromstr(dat.decode("utf-8"))
....                              for col in range(max_cols)}
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
