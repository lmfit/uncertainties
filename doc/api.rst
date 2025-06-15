.. index:: technical details

=========================
Advanced Topics
=========================


This page gives more in-depth technical description of the
:mod:`uncertainties` package.

 .. index:: api

.. _api_funcs:


API: Application Programming Interface
==============================================

.. module:: uncertainties

The most common and important functions for creating uncertain
:class:`Variables` are :func:`ufloat` and :func:`ufloat_fromstr`.  In addition,
the :func:`wrap` can be used to support the propagation of uncertainties with a
user-supplied function.

.. autofunction:: ufloat

.. autofunction:: ufloat_fromstr

.. autoclass:: UFloat

.. autofunction:: wrap

Special Technical Topics
============================================================

.. index::
   pair: uncertainty; NaN

NaN uncertainty
----------------------

If linear `error propagation theory`_ cannot be applied, the functions
defined by :mod:`uncertainties` internally use a `not-a-number value
<http://en.wikipedia.org/wiki/Not_a_number>`_ (``nan``) for the
derivative.

As a consequence, it is possible for uncertainties to be ``nan``:

>>> from uncertainties import umath
>>> umath.sqrt(ufloat(0, 1))
0.0+/-nan

This indicates that **the derivative required by linear error
propagation theory is not defined** (a Monte-Carlo calculation of the
resulting random variable is more adapted to this specific case).

.. _math_def_num_uncert:

Mathematical definition of numbers with uncertainties
-----------------------------------------------------

.. index:: number with uncertainty; definition
.. index:: probability distribution

Mathematically, **numbers with uncertainties** are, in this package,
**probability distributions**.  They are *not restricted* to normal
(Gaussian) distributions and can be **any distribution**.  These
probability distributions are reduced to two numbers: a nominal value
and an uncertainty.

Thus, both independent variables (:class:`Variable` objects) and the
result of mathematical operations (:class:`AffineScalarFunc` objects)
contain these two values (respectively in their :attr:`nominal_value`
and :attr:`std_dev` attributes).

.. index:: uncertainty; definition

The **uncertainty** of a number with uncertainty is simply defined in
this package as the **standard deviation** of the underlying probability
distribution.

The numbers with uncertainties manipulated by this package are assumed
to have a probability distribution mostly contained around their
nominal value, in an interval of about the size of their standard
deviation.  This should cover most practical cases.

.. index:: nominal value; definition

A good choice of **nominal value** for a number with uncertainty is thus
the median of its probability distribution, the location of highest
probability, or the average value.

Probability distributions (random variables and calculation results)
are printed as::

  nominal value +/- standard deviation

but this does not imply any property on the nominal value (beyond the
fact that the nominal value is normally inside the region of high
probability density), or that the probability distribution of the
result is symmetrical (this is rarely strictly the case).


.. _differentiation method:

Differentiation method
----------------------

The :mod:`uncertainties` package automatically calculates the
derivatives required by linear error propagation theory.

Almost all the derivatives of the fundamental functions provided by
:mod:`uncertainties` are obtained through analytical formulas (the
few mathematical functions that are instead differentiated through
numerical approximation are listed in ``umath_core.num_deriv_funcs``).

The derivatives of mathematical *expressions* are evaluated through a
fast and precise method: :mod:`uncertainties` transparently implements
`automatic differentiation`_ with reverse accumulation. This method
essentially consists in keeping track of the value of derivatives, and
in automatically applying the `chain rule
<http://en.wikipedia.org/wiki/Chain_rule>`_. Automatic differentiation
is faster than symbolic differentiation and more precise than
numerical differentiation.

The derivatives of any expression can be obtained with
:mod:`uncertainties` in a simple way, as demonstrated in the :ref:`User
Guide <derivatives>`.

.. _variable_tracking:

Tracking of random variables
----------------------------

This package keeps track of all the random variables a quantity
depends on, which allows one to perform transparent calculations that
yield correct uncertainties.  For example:

>>> x = ufloat(2, 0.1)
>>> a = 42
>>> poly = x**2 + a
>>> poly
46.0+/-0.4
>>> poly - x*x
42.0+/-0

Even though ``x*x`` has a non-zero uncertainty, the result has a zero
uncertainty, because it is equal to :data:`a`.

If the variable :data:`a` above is modified, the value of :data:`poly`
is not modified, as is usual in Python:

>>> a = 123
>>> print(poly)  # Still equal to x**2 + 42, not x**2 + 123
46.0+/-0.4

Random variables can, on the other hand, have their uncertainty
updated on the fly, because quantities with uncertainties (like
:data:`poly`) keep track of them:

>>> x.std_dev = 0
>>> print(poly)  # Zero uncertainty, now
46.0+/-0

As usual, Python keeps track of objects as long as they are used.
Thus, redefining the value of :data:`x` does not change the fact that
:data:`poly` depends on the quantity with uncertainty previously stored
in :data:`x`:

>>> x = 10000
>>> print(poly)  # Unchanged
46.0+/-0

These mechanisms make quantities with uncertainties behave mostly like
regular numbers, while providing a fully transparent way of handling
correlations between quantities.


.. index:: number with uncertainty; classes, Variable class
.. index::  AffineScalarFunc class

.. _classes:


Python classes for variables and functions with uncertainty
-----------------------------------------------------------

Numbers with uncertainties are represented through two different
classes:

1. a class for independent random variables (:class:`Variable`, which
   inherits from :class:`UFloat`),

2. a class for functions that depend on independent variables
   (:class:`AffineScalarFunc`, aliased as :class:`UFloat`).

Additional documentation for these classes is available in their Python
docstring.

The factory function :func:`ufloat` creates variables and thus returns
a :class:`Variable` object:

>>> x = ufloat(1, 0.1)
>>> type(x)
<class 'uncertainties.core.Variable'>

:class:`Variable` objects can be used as if they were regular Python
numbers (the summation, etc. of these objects is defined).

Mathematical expressions involving numbers with uncertainties
generally return :class:`AffineScalarFunc` objects, because they
represent mathematical functions and not simple variables; these
objects store all the variables they depend on:

>>> type(umath.sin(x))
<class 'uncertainties.core.AffineScalarFunc'>


.. _automatic differentiation: http://en.wikipedia.org/wiki/Automatic_differentiation

.. _error propagation theory: http://en.wikipedia.org/wiki/Error_propagation

.. _soerp: https://pypi.python.org/pypi/soerp
.. _mcerp: https://pypi.python.org/pypi/mcerp
