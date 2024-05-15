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

.. autoclass:: Variable

.. autofunction:: wrap

Testing whether an object is a number with uncertainty
------------------------------------------------------

The recommended way of testing whether :data:`value` carries an
uncertainty handled by this module is by checking whether
:data:`value` is an instance of :class:`UFloat`, through
``isinstance(value, uncertainties.UFloat)``.




Special Technical Topics
============================================================



.. index:: pickling
.. index:: saving to file; number with uncertainty
.. index:: reading from file; number with uncertainty

.. _pickling:

Pickling
--------

The quantities with uncertainties created by the :mod:`uncertainties`
package can be `pickled <http://docs.python.org/library/pickle.html>`_
(they can be stored in a file, for instance).

If multiple variables are pickled together (including when pickling
:doc:`NumPy arrays <numpy_guide>`), their correlations are preserved:

>>> import pickle
>>> x = ufloat(2, 0.1)
>>> y = 2*x
>>> p = pickle.dumps([x, y])  # Pickling to a string
>>> (x2, y2) = pickle.loads(p)  # Unpickling into new variables
>>> y2 - 2*x2
0.0+/-0

The final result is exactly zero because the unpickled variables :data:`x2`
and :data:`y2` are completely correlated.

However, **unpickling necessarily creates new variables that bear no
relationship with the original variables** (in fact, the pickled
representation can be stored in a file and read from another program
after the program that did the pickling is finished: the unpickled
variables cannot be correlated to variables that can disappear).  Thus

>>> x - x2
0.0+/-0.14142135623730953

which shows that the original variable :data:`x` and the new variable :data:`x2`
are completely uncorrelated.


.. index:: comparison operators; technical details

.. _comparison_operators:


Comparison operators
--------------------

Comparison operations (>, ==, etc.) on numbers with uncertainties have
a **pragmatic semantics**, in this package: numbers with uncertainties
can be used wherever Python numbers are used, most of the time with a
result identical to the one that would be obtained with their nominal
value only.  This allows code that runs with pure numbers to also work
with numbers with uncertainties.

.. index:: boolean value

The **boolean value** (``bool(x)``, ``if x …``) of a number with
uncertainty :data:`x` is defined as the result of ``x != 0``, as usual.

However, since the objects defined in this module represent
probability distributions and not pure numbers, comparison operators
are interpreted in a specific way.

The result of a comparison operation is defined so as to be
essentially consistent with the requirement that uncertainties be
small: the **value of a comparison operation** is True only if the
operation yields True for all *infinitesimal* variations of its random
variables around their nominal values, *except*, possibly, for an
*infinitely small number* of cases.

Example:

>>> x = ufloat(3.14, 0.01)
>>> x == x
True

because a sample from the probability distribution of :data:`x` is always
equal to itself.  However:

>>> y = ufloat(3.14, 0.01)
>>> x == y
False

since :data:`x` and :data:`y` are independent random variables that
*almost* always give a different value (put differently,
:data:`x`-:data:`y` is not equal to 0, as it can take many different
values). Note that this is different
from the result of ``z = 3.14; t = 3.14; print(z == t)``, because
:data:`x` and :data:`y` are *random variables*, not pure numbers.

Similarly,

>>> x = ufloat(3.14, 0.01)
>>> y = ufloat(3.00, 0.01)
>>> x > y
True

because :data:`x` is supposed to have a probability distribution largely
contained in the 3.14±~0.01 interval, while :data:`y` is supposed to be
well in the 3.00±~0.01 one: random samples of :data:`x` and :data:`y` will
most of the time be such that the sample from :data:`x` is larger than the
sample from :data:`y`.  Therefore, it is natural to consider that for all
practical purposes, ``x > y``.

Since comparison operations are subject to the same constraints as
other operations, as required by the :ref:`linear approximation
<linear_method>` method, their result should be essentially *constant*
over the regions of highest probability of their variables (this is
the equivalent of the linearity of a real function, for boolean
values).  Thus, it is not meaningful to compare the following two
independent variables, whose probability distributions overlap:

>>> x = ufloat(3, 0.01)
>>> y = ufloat(3.0001, 0.01)

In fact the function (x, y) → (x > y) is not even continuous over the
region where x and y are concentrated, which violates the assumption
of approximate linearity made in this package on operations involving
numbers with uncertainties.  Comparing such numbers therefore returns
a boolean result whose meaning is undefined.

However, values with largely overlapping probability distributions can
sometimes be compared unambiguously:

>>> x = ufloat(3, 1)
>>> x
3.0+/-1.0
>>> y = x + 0.0002
>>> y
3.0002+/-1.0
>>> y > x
True

In fact, correlations guarantee that :data:`y` is always larger than
:data:`x`: ``y-x`` correctly satisfies the assumption of linearity,
since it is a constant "random" function (with value 0.0002, even
though :data:`y` and :data:`x` are random). Thus, it is indeed true
that :data:`y` > :data:`x`.


.. index:: linear propagation of uncertainties
.. _linear_method:

Linear propagation of uncertainties
-----------------------------------

This package calculates the standard deviation of mathematical
expressions through the linear approximation of `error propagation
theory`_.

The standard deviations and nominal values calculated by this package
are thus meaningful approximations as long as **uncertainties are
"small"**. A more precise version of this constraint is that the final
calculated functions must have **precise linear expansions in the region
where the probability distribution of their variables is the largest**.
Mathematically, this means that the linear terms of the final calculated
functions around the nominal values of their variables should be much
larger than the remaining higher-order terms over the region of
significant probability (because such higher-order contributions are
neglected).

For example, calculating ``x*10`` with :data:`x` = 5±3 gives a
*perfect result* since the calculated function is linear. So does
``umath.atan(umath.tan(x))`` for :data:`x` = 0±1, since only the
*final* function counts (not an intermediate function like
:func:`tan`).

Another example is ``sin(0+/-0.01)``, for which :mod:`uncertainties`
yields a meaningful standard deviation since the sine is quite linear
over 0±0.01.  However, ``cos(0+/-0.01)``, yields an approximate
standard deviation of 0 because it is parabolic around 0 instead of
linear; this might not be precise enough for all applications.

**More precise uncertainty estimates** can be obtained, if necessary,
with the soerp_ and mcerp_ packages. The soerp_ package performs
*second-order* error propagation: this is still quite fast, but the
standard deviation of higher-order functions like f(x) = x\ :sup:`3`
for x = 0±0.1 is calculated as being exactly zero (as with
:mod:`uncertainties`). The mcerp_ package performs Monte-Carlo
calculations, and can in principle yield very precise results, but
calculations are much slower than with approximation schemes.

.. index::
   pair: uncertainty; NaN

NaN uncertainty
----------------------

If linear `error propagation theory`_ cannot be applied, the functions
defined by :mod:`uncertainties` internally use a `not-a-number value
<http://en.wikipedia.org/wiki/Not_a_number>`_ (``nan``) for the
derivative.

As a consequence, it is possible for uncertainties to be ``nan``:

>>> umath.sqrt(ufloat(0, 1))
0.0+/-nan

This indicates that **the derivative required by linear error
propagation theory is not defined** (a Monte-Carlo calculation of the
resulting random variable is more adapted to this specific case).

However, even in this case where the derivative at the nominal value
is infinite, the :mod:`uncertainties` package **correctly handles
perfectly precise numbers**:

>>> umath.sqrt(ufloat(0, 0))
0.0+/-0

is thus the correct result, despite the fact that the derivative of
the square root is not defined in zero.

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
42+/-0

Even though ``x*x`` has a non-zero uncertainty, the result has a zero
uncertainty, because it is equal to :data:`a`.

If the variable :data:`a` above is modified, the value of :data:`poly`
is not modified, as is usual in Python:

>>> a = 123
>>> print(poly)
46.0+/-0.4  # Still equal to x**2 + 42, not x**2 + 123

Random variables can, on the other hand, have their uncertainty
updated on the fly, because quantities with uncertainties (like
:data:`poly`) keep track of them:

>>> x.std_dev = 0
>>> print(poly)
46+/-0  # Zero uncertainty, now

As usual, Python keeps track of objects as long as they are used.
Thus, redefining the value of :data:`x` does not change the fact that
:data:`poly` depends on the quantity with uncertainty previously stored
in :data:`x`:

>>> x = 10000
>>> print(poly)
46+/-0  # Unchanged

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
<class 'uncertainties.Variable'>

:class:`Variable` objects can be used as if they were regular Python
numbers (the summation, etc. of these objects is defined).

Mathematical expressions involving numbers with uncertainties
generally return :class:`AffineScalarFunc` objects, because they
represent mathematical functions and not simple variables; these
objects store all the variables they depend on:

>>> type(umath.sin(x))
<class 'uncertainties.AffineScalarFunc'>


.. _automatic differentiation: http://en.wikipedia.org/wiki/Automatic_differentiation

.. _error propagation theory: http://en.wikipedia.org/wiki/Error_propagation

.. _soerp: https://pypi.python.org/pypi/soerp
.. _mcerp: https://pypi.python.org/pypi/mcerp
