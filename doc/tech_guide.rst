.. index:: technical details

Technical Guide
===============

.. !!!!! The paragraphs below are only notes: they should be adapted.


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
important that uncertainties be "small".  Mathematically, this means
that the linear term of functions around the nominal values of their
variables should be much larger than the remaining higher-order terms
over the region of significant probability.

For instance, ``sin(0+/-0.01)`` yields a meaningful standard deviation
since it is quite linear over 0±0.01.  However, ``cos(0+/-0.01)``,
yields an approximate standard deviation of 0 (because the cosine is
not well approximated by a line around 0), which might not be precise
enough for all applications.

.. index:: comparison operators; technical details

Comparison operators
--------------------

Comparison operations (>, ==, etc.) on numbers with uncertainties have
a **pragmatic semantics**, in this package: numbers with uncertainties
can be used wherever Python numbers are used, most of the time with a
result identical to the one that would be obtained with their nominal
value only.  This allows code that runs with pure numbers to also work
with numbers with uncertainties.

However, since the objects defined in this module represent
probability distributions and not pure numbers, comparison operators
are interpreted in a specific way.

The result of a comparison operation ("==", ">", etc.) is defined so as
to be essentially consistent with the requirement that uncertainties
be small: the value of a comparison operation is True only if the
operation yields True for all infinitesimal variations of its random
variables, except, possibly, for an infinitely small number of cases.

Example:

  >>> x = ufloat((3.14, 0.01))
  >>> x == x
  True

because a sample from the probability distribution of ``x`` is always
equal to itself.  However:

  >>> y = ufloat((3.14, 0.01))
  >>> x == y
  False

since ``x`` and ``y`` are independent random variables that *almost*
never give the same value.

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

Comparison operations are subject to the same constraints as other
operations, as required by the :ref:`linear approximation
<linear_method>` method: their result should be linear (i.e. constant,
for boolean values) over the regions of highest probability of their
variables.  Thus, it is not meaningful to compare the following two
values, whose probability distributions overlap:

  >>> x = ufloat((3, 0.01))
  >>> y = ufloat((3.0001, 0.01))

In fact the function (x, y) → (x > y) is not even continuous (and
linear) over the region where x and y are concentrated, which violates
the assumption made in this package about operations.  Comparing such
numbers therefore returns a boolean result whose meaning is undefined.

The boolean value (``bool(x)``, ``if x…``) of a number with
uncertainty ``x`` is the result of ``x != 0``, as usual.

.. index:: number with uncertainty; classes, Variable class
.. index::  AffineScalarFunc class

Classes
-------

Numbers with uncertainties are represented through two different
classes:
1. A class for independent variables (:class:`Variable`),
2. A class for functions that depend on independent variables
(:class:`AffineScalarFunc`).

Thus, the factory function :func:`ufloat` creates variables and
returns a :class:`Variable` object:

  >>> x = ufloat((1, 0.1))
  >>> type(x)
  <class 'uncertainties.Variable'>

:class:`Variable` objects can be used as if they were regular Python
numbers (the summation, etc. of these objects is defined).

Mathematical expressions involving numbers with uncertainties
generally return :class:`AffineScalarFunc` objects; these objects
store all the variables they depend from:

  >>> type(umath.sin(x))
  <class 'uncertainties.AffineScalarFunc'>

Note that :class:`Variable` objects are also :class:`AffineScalarFunc`
objects (a variable x is simply considered to be the identity function
x → x): testing whether ``my_value`` carries an uncertainty handled by
this module should be done with ``insinstance(my_value,
AffineScalarFunc)``.


