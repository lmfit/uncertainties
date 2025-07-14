.. index: NumPy support

===============================
Uncertainties and numpy arrays
===============================

.. index:: unumpy
.. index:: arrays; simple use

.. _simple_array_use:

Arrays of uncertainties Variables
====================================

It is possible to put :class:`UFloat` objects  in NumPy_ arrays:

>>> import numpy as np
>>> from uncertainties import UFloat
>>> arr = np.array([UFloat(1, 0.01), UFloat(2, 0.1)])
>>> print(2*arr)
[2.0+/-0.02 4.0+/-0.2]
>>> print(str(arr.sum()))
3.00+/-0.10

Many common operations on NumPy arrays can be performed transparently
even when these arrays contain numbers with uncertainties.


The unumpy package
==================

It is possible to perform :ref:`basic operations on arrays <simple_array_use>` that
contain :class:`UFloat` objects using the :mod:`numpy` package together with
:class:`UFloat` objects.
However, the :mod:`uncertainties.unumpy` package provides additional array
functionality including:

1. utilities that help with the **creation and manipulation** of
NumPy_ arrays and matrices of numbers with uncertainties;

2. **generalizations** of multiple NumPy functions so that they also
work with arrays that contain numbers with uncertainties.

Operations on arrays (e.g. their cosine, etc.)  can thus be performed transparently.

These features can be made available with

>>> from uncertainties import unumpy

.. Here, there is no need to mention unumpy.unlinalg, because it is indeed
   made available through "import unumpy".

Creation and manipulation of arrays
-----------------------------------

.. index::
   single: arrays; creation and manipulation
   single: creation; arrays

Arrays
^^^^^^

Arrays of :class:`UFloat objects can be built from sequences of nominal values and
standard deviations:

>>> arr = unumpy.uarray([1, 2], [0.01, 0.002])
>>> print(arr)
[1.0+/-0.01 2.0+/-0.002]

NumPy arrays of :class:`UFloat` objects can also be built directly through NumPy, thanks
to NumPy's support of arrays of arbitrary objects:

>>> arr = np.array([UFloat(1, 0.1), UFloat(2, 0.002)])


.. index::
   pair: nominal value; uniform access (array)
   pair: uncertainty; uniform access (array)
   pair: standard deviation; uniform access (array)

Uncertainties and nominal values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Nominal values and uncertainties in arrays can be directly accessed (through functions
that work on pure float arrays too):

>>> unumpy.nominal_values(arr)
array([1., 2.])


.. index:: mathematical operation; on an array of numbers

Mathematical functions
----------------------

This module defines uncertainty-aware mathematical functions that
generalize those from :mod:`uncertainties.umath` so that they work on
NumPy arrays of numbers with uncertainties instead of just scalars:

>>> print(unumpy.cos(arr))  # Cosine of each array element
[0.5403023058681398+/-0.08414709848078966
 -0.4161468365471424+/-0.0018185948536513636]

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
>>> arr = np.array([nan, UFloat(nan, 1), UFloat(1, nan), 2])
>>> print(arr)
[nan nan+/-1.0 1.0+/-nan 2]
>>> print(arr[~unumpy.isnan(arr)].mean())
1.5+/-nan

or equivalently, by using masked arrays:

>>> masked_arr = np.ma.array(arr, mask=unumpy.isnan(arr))
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

Number with uncertainties can easy be cast to strings and back. This means that arrays
of numbers with uncertainties can also be cast to string representations and back.
There are many ways to convert an array of numbers with uncertainties to a string
representation for storage and then convert it back to a python array of numbers with
uncertainties.
Here is one example set of functions to perform this operation.

>>> import json
>>> from uncertainties import ufloat_fromstr
>>> def serialize_unumpy_array(u_arr):
...     string_u_arr = np.vectorize(repr)(u_arr)
...     return json.dumps(string_u_arr.tolist(), indent=4)
>>>
>>> def deserialize_unumpy_arr(serialized_u_arr):
...     string_u_arr = np.array(json.loads(serialized_u_arr))
...     return np.vectorize(ufloat_fromstr)(string_u_arr)

We can use the first function to serialize an array

>>> u_arr = np.array([
...     [UFloat(1, 0.1), UFloat(2, 0.2)],
...     [UFloat(3, 0.3), UFloat(4, 0.4)],
... ])
>>> print(u_arr)
[[1.0+/-0.1 2.0+/-0.2]
 [3.0+/-0.3 4.0+/-0.4]]
>>> serialized_u_arr = serialize_unumpy_array(u_arr)
>>> print(serialized_u_arr)
[
    [
        "1.0+/-0.1",
        "2.0+/-0.2"
    ],
    [
        "3.0+/-0.3",
        "4.0+/-0.4"
    ]
]

This can then of course be stored in a ``.json`` file using ``json.dump``.
We can then deserialize

>>> u_arr_2 = deserialize_unumpy_arr(serialized_u_arr)
>>> print(u_arr_2)
[[1.0+/-0.1 2.0+/-0.2]
 [3.0+/-0.3 4.0+/-0.4]]

Note that the process of serializing and deserializing the array of numbers with
uncertainties has result in all correlations between numbers within one array, and also
between numbers from the original array and its deserialized copy

>>> print(u_arr[0, 0] - u_arr_2[0, 0])
0.00+/-0.14
>>> print(u_arr[0, 0] == u_arr_2[0, 0])
False

A future release of :mod:`uncertainties` may provide functionality for
serializing/deserializing number with uncertainties in such a way that correlations can
be preserved.

.. index:: linear algebra; additional functions, ulinalg

Storing arrays in text format
=============================

Formatting Arrays
=================

Arrays of :class:`UFloat` objects have ``dtype=object``.
By default, :mod:`numpy` does not know how to format entries in these arrays and falls
back to the object ``__repr__`` or ``__str__`` function.
For :class:`UFloat` objects this may be inconvenient for readability.

>>> arr = np.array([UFloat(1/3, 1/9), UFloat(2/3, 2/9)])
>>> print(arr)
[0.3333333333333333+/-0.1111111111111111
 0.6666666666666666+/-0.2222222222222222]

The following function can be used to more nicely format arrays of :class:`UFloat`
objects.

>>> from functools import partial
>>> def format_uarray(uarr, fmt_spec=""):
...     return np.vectorize(lambda x: format(x, fmt_spec))(uarr)
>>>
>>> print(format_uarray(arr, fmt_spec=".3u"))
['0.333+/-0.111' '0.667+/-0.222']

Additional array functions: unumpy.ulinalg
==========================================

The :mod:`unumpy.ulinalg` module contains more uncertainty-aware
functions for arrays that contain numbers with uncertainties.

It currently offers generalizations of two functions from
:mod:`numpy.linalg` that work on arrays (or matrices) that contain
numbers with uncertainties, the **matrix inverse and pseudo-inverse**:

>>> arr = unumpy.ulinalg.inv(np.array([[UFloat(2, 0.1)]]))
>>> print(arr)
[[0.5+/-0.025]]
>>> arr = np.array([[UFloat(1, 0.1), UFloat(2, 0.002)]])
>>> arr_pinv = unumpy.ulinalg.pinv(arr)
>>> print(format_uarray(arr_pinv))
[['0.200+/-0.012']
 ['0.400+/-0.016']]

.. _NumPy: http://numpy.scipy.org/
