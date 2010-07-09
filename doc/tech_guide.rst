.. index:: technical details

Technical Guide
===============


Mathematical definition of numbers with uncertainties
-----------------------------------------------------

.. index:: number with uncertainty; definition

Mathematically, **numbers with uncertainties** are, in this package,
probability distributions.  These probabilities are reduced to two
numbers: a nominal value and an uncertainty.

Thus, both variables (:class:`Variable` objects) and the result of
mathematical operations (:class:`AffineScalarFunc` objects) contain
these two values (respectively in their :attr:`nominal_value`
attribute and through their :meth:`std_dev` method).

.. index:: uncertainty; definition

The **uncertainty** of a number with uncertainty is simply defined in
this package as the standard deviation of the underlying probability
distribution.

.. index:: probability distribution

The numbers with uncertainties manipulated by this package are assumed
to have a **probability distribution** mostly contained around their
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

.. index:: correlations; technical details

Tracking of random variables
----------------------------

This package keeps track of all the random variables a quantity
depends on, which allows one to perform transparent calculations that
yield correct uncertainties.  For example:

  >>> x = ufloat((2, 0.1))
  >>> a = 42
  >>> square = x**2 + a
  >>> square
  46.0+/-0.4
  >>> square - x*x
  42.0

Even though ``x*x`` has a non-zero uncertainty, the result has a zero
uncertainty, because it is equal ``a``.

However, only the dependence of quantities on random variables created
by this module is tracked.  Thus, if the variable ``a`` above is
modified, the value of ``square`` is not modified, as is usual in
Python:

  >>> a = 123
  >>> print square
  46.0+/-0.4  # Still equal to x**2 + 42, not x**2 + 123

Random variables can, on the other hand, have their uncertainty
updated on the fly, because quantities with uncertainties (like
``square``) keep track of them:

  >>> x.set_std_dev(0)
  >>> print square
  0.04  # Zero uncertainty, now

As usual, Python keeps track of objects as long as they are used.
Thus, redefining the value of ``x`` does not change the fact that
``square`` depends on the quantity with uncertainty previously stored
in ``x``:

  >>> x = 10000
  >>> print square
  0.04  # Unchanged

These mechanisms make quantities with uncertainties behave mostly like
regular numbers, while providing a fully transparent way of handling
correlations between quantities.

.. _linear_method:

Linear error propagation theory
-------------------------------

Linear approximations of functions (around the nominal values) are
used for the calculation of the standard deviation of mathematical
expressions with this package.

The calculated standard deviations and nominal values are thus
meaningful approximations as long as the functions involved have
precise linear expansions in the region where the probability
distribution of their variables is the largest.  It is therefore
important that **uncertainties be "small"**.  Mathematically, this
means that the linear term of functions around the nominal values of
their variables should be much larger than the remaining higher-order
terms over the region of significant probability.

For instance, ``sin(0+/-0.01)`` yields a meaningful standard deviation
since it is quite linear over 0±0.01.  However, ``cos(0+/-0.01)``,
yields an approximate standard deviation of 0 (because around 0, the
cosine is parabolic, not linear), which might not be precise enough
for all applications.

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

The **boolean value** (``bool(x)``, ``if x…``) of a number with
uncertainty ``x`` is defined as the result of ``x != 0``, as usual.

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

  >>> x = ufloat((3.14, 0.01))
  >>> x == x
  True

because a sample from the probability distribution of ``x`` is always
equal to itself.  However:

  >>> y = ufloat((3.14, 0.01))
  >>> x != y
  True

since ``x`` and ``y`` are independent random variables that *almost*
always give a different value.

Similarly,

  >>> x = ufloat((3.14, 0.01))
  >>> y = ufloat((3.00, 0.01))
  >>> x > y
  True

because ``x`` is supposed to have a probability distribution largely
contained in the 3.14±~0.01 interval, while ``y`` is supposed to be
well in the 3.00±~0.01 one: random samples of ``x`` and ``y`` will
most of the time be such that the sample from ``x`` is larger than the
sample from ``y``.  Therefore, it is natural to consider that for all
practical purposes, ``x > y``.

Since comparison operations are subject to the same constraints as
other operations, as required by the :ref:`linear approximation
<linear_method>` method, their result should be essentially *constant*
over the regions of highest probability of their variables (this is
the equivalent of the linearity of a real function, for boolean
values).  Thus, it is not meaningful to compare the following two
independent variables, whose probability distributions overlap:

  >>> x = ufloat((3, 0.01))
  >>> y = ufloat((3.0001, 0.01))

In fact the function (x, y) → (x > y) is not even continuous over the
region where x and y are concentrated, which violates the assumption
made in this package about operations involving numbers with
uncertainties.  Comparing such numbers therefore returns a boolean
result whose meaning is undefined.

However, values with largely overlapping probability distributions can
sometimes be compared unambiguously:

  >>> x = ufloat((3, 1))
  >>> x
  3.0+/-1.0
  >>> y = x + 0.0002
  >>> y
  3.0002+/-1.0
  >>> y > x
  True

In fact, correlations guarantee that ``y`` is always larger than
``x`` (by 0.0002).

.. index:: number with uncertainty; classes, Variable class
.. index::  AffineScalarFunc class

Classes
-------

Numbers with uncertainties are represented through two different
classes:

1. a class for independent random variables (:class:`Variable`),

2. a class for functions that depend on independent variables (:class:`AffineScalarFunc`).

Thus, the factory function :func:`ufloat` creates variables and
returns a :class:`Variable` object:

  >>> x = ufloat((1, 0.1))
  >>> type(x)
  <class 'uncertainties.Variable'>

:class:`Variable` objects can be used as if they were regular Python
numbers (the summation, etc. of these objects is defined).

Mathematical expressions involving numbers with uncertainties
generally return :class:`AffineScalarFunc` objects, because they
represent mathematical functions and not simple variables; these
objects store all the variables they depend from:

  >>> type(umath.sin(x))
  <class 'uncertainties.AffineScalarFunc'>

Note that :class:`Variable` objects are also :class:`AffineScalarFunc`
objects (a variable x is simply considered to be the identity function
x → x): testing whether ``value`` carries an uncertainty handled by
this module can therefore be done with ``insinstance(value,
AffineScalarFunc)``.


