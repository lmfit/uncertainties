.. index:: user guide

User guide
==========

Basic mathematical operations involving numbers with uncertainties
only require a simple import:

  >>> from uncertainties import ufloat

The :func:`ufloat` function creates numbers with uncertainties.
Existing calculation code can use such numbers as input, usually run
with no or little modification, and automatically produce results with
uncertainties.

.. index::
   pair: number with uncertainty; creation

Creating and handling numbers with uncertainties
------------------------------------------------

Numbers with uncertainties can be input either numerically, or through
one of many string representations:

  >>> x = ufloat((0.20, 0.01))  # x = 0.20+/-0.01
  >>> x = ufloat("0.20+/-0.01")
  >>> x = ufloat("0.20(1)")
  >>> x = ufloat("0.20")


Basic math
-----------

Calculations can be performed directly, as with regular real numbers:

  >>> square = x**2
  >>> print square
  0.04+/-0.004

.. index::
   pair: nominal value; of scalar
   pair: uncertainty; of scalar

.. index:: correlations; detailed example

Correlations
------------

Correlations between variables are automatically handled whatever the
number of variables involved, and whatever the complexity of the
calculation.  Thus, each calculation result keeps track of how it is
correlated to random variables.  For example, when ``x`` is the number
with uncertainty defined above,

  >>> square = x**2
  >>> print square
  0.04+/-0.004
  >>> y = x*x + 1
  >>> square - x*x
  0.0
  >>> y - square
  1.0

The results above have a zero uncertainty despite the fact that ``x``
has a non-zero uncertainty: the calculated functions give the same
value for all samples of the random variable ``x``.

Thanks to the tracking of dependencies on random variables,
calculations can therefore be performed in as many steps as necessary,
exactly as with simple floats.  When various quantities are combined
through mathematical operations, the result is calculated by taking
into account all the correlations between the quantities involved.
All of this is done completely transparently.

Access to the uncertainty and to the nominal value
-----------------------------------------------

The nominal value and the uncertainty (standard deviation) on the
calculated square can also be accessed independently:

  >>> print square
  0.04+/-0.004
  >>> print square.nominal_value
  0.04
  >>> print square.std_dev()
  0.004

Details on the classes made available by this package can be found in
the :ref:`Technical Guide <classes>`.

The various independent contributions to an uncertainty can be
directly obtained.  This information is more easily usable when the
variables are tagged:

  >>> u = ufloat((1, 0.1), "u variable")  # Tag
  >>> v = ufloat((10, 0.1), "v variable")
  >>> sum_value = u+2*v
  >>> sum_value
  21.0+/-0.22360679774997899
  >>> for (var, error) in sum_value.error_components().items():
  ...     print "%s: %f" % (var.tag, error)
  ...
  u variable: 0.100000
  v variable: 0.200000

The total uncertainty on the result (``sum_value``) is the quadratic
sum of these independent uncertainties, as it should.

.. index:: mathematical operation; on a scalar, umath

Mathematical operations
-----------------------

Besides being able to apply basic mathematical operations to numbers
with uncertainty, this package provides generalizations of most
functions from the standard :mod:`math` module.  These operations are
found in the :mod:`uncertainties.umath` module::

  >>> from uncertainties.umath import *  # Imports sin(), etc.
  >>> sin(x**2)
  0.039989334186634168+/-0.003996800426643912


.. index:: comparison operators

Comparison operators
---------------------

Comparison operators behave in a natural way::

  >>> print x
  0.2+/-0.01
  >>> y = x + 0.0001
  >>> y
  0.2001+/-0.01
  >>> y > x
  True
  >>> y > 0
  True

One important concept to keep in mind is that :func:`ufloat` creates
a random variable:

  >>> y = ufloat((1, 0.1))
  >>> z = ufloat((1, 0.1))
  >>> print y
  1.0+/-0.1
  >>> z
  1.0+/-0.1
  >>> y == y
  True
  >>> y == z
  False

In physical terms, two rods of the same nominal length and uncertainty
on their length generally are of different sizes: ``y`` is different
from ``z``.

More detailed information on the semantics of comparison operators for
numbers with uncertainties can be found in the :ref:`Technical Guide
<comparison_operators>`.

.. index:: arrays; simple use, matrices; simple use

.. _simple_array_use:

Arrays of numbers with uncertainties
------------------------------------

It is possible to put numbers with uncertainties in NumPy_ arrays and
matrices:

  >>> print 2*numpy.array([ufloat((1, 0.01)), ufloat((2, 0.1))])
  [2.0+/-0.02 4.0+/-0.2]

:doc:`More complex operations on NumPy arrays <numpy_guide>` can be
performed through the dedicated :mod:`uncertainties.unumpy` module.

.. index:: covariance matrix

Covariance matrix
-----------------

The covariance matrix between various variables or calculated
quantities can be simply obtained::

  >>> sum_value = u+2*v
  >>> cov_matrix = uncertainties.covariance_matrix([u, v, sum_value])

has value

::

  [[0.01, 0.0,  0.01],
   [0.0,  0.01, 0.02],
   [0.01, 0.02, 0.05]]

In this matrix, the zero covariances indicate that ``u`` and ``v`` are
independent from each other; the last column shows that ``sum_value``
does depend on these variables.  The uncertainties package keeps track
at all times of all correlations between quantities (variables and
functions):

  >>> sum_value - (u+2*v)
  >>> 0.0

.. index:: correlations; correlated variables

Correlated variables
--------------------

Reciprocally, correlated variables can be created transparently,
provided that the NumPy_ package is available::

  >>> (u2, v2, sum2) = uncertainties.correlated_values([1, 10, 21], cov_matrix)

creates three new variables with the indicated values, and correct
uncertainties and correlations::

  >>> sum_value
  21.0+/-0.22360679774997899
  >>> sum2
  21.0+/-0.22360679774997899
  >>> sum2 - (u2+2*v2)
  0.0+/-3.8337185686225597e-09

The theoretical value of the last expression is exactly zero, like for
``sum - (u+2*v)``, but numerical errors yield a small uncertainty
(3e-9 is indeed very small compared to the uncertainty on ``sum2``:
correlations should in fact cancel the uncertainty on ``sum2``).

The correlation matrix is the desired one::

  >>> uncertainties.covariance_matrix([u2, v2, sum2])

reproduces the desired covariance matrix ``cov_matrix`` (up to
rounding errors).

.. index::
   single: C code; wrapping
   single: Fortran code; wrapping
   single: wrapping (C, Fortran,…) functions

Making functions accept numbers with uncertainties
--------------------------------------------------

This package allows calculations that are performed through non-Python
code (Fortran, C, etc.) to handle numbers with uncertainties instead
of floats.  Similarly, functions that do not have a simple analytical
form can be automatically wrapped so as to also work on float
parameters that contain uncertainties.

It is thus possible to define a function :func:`f` that takes any
number of real numbers, and to automatically generalize it so that it
also works with numbers with uncertainty:

  >>> wrapped_f = uncertainties.wrap(f)

The new function :func:`wrapped_f` can be given numbers with
uncertainties.  It returns the same values as :func:`f`, but with
uncertainties.

Miscellaneous utilities
-----------------------

.. index:: standard deviation; on the fly modification

It is sometimes useful to modify the error on certain parameters so as
to study its impact on a final result.  With this package, the
**uncertainty of a variable can be changed** on the fly:

  >>> sum_value = u+2*v
  >>> sum_value
  21.0+/-0.22360679774997899
  >>> prev_uncert = u.std_dev()
  >>> u.set_std_dev(10)
  >>> sum_value
  21.0+/-10.001999800039989
  >>> u.set_std_dev(prev_uncert)

The relevant concept is that ``sum_value`` does depend on the
variables ``u`` and ``v``: the :mod:`uncertainties` package keeps
track of this fact, as detailed in the :ref:`Technical Guide
<variable_tracking>`, and uncertainties can thus be updated any time.

.. index::
   pair: nominal value; uniform access (scalar)
   pair: uncertainty; uniform access (scalar)
   pair: standard deviation; uniform access (scalar)

When manipulating ensembles of numbers, *some* of which contain
uncertainties while others are simple floats, it can be useful to
access the **nominal value and uncertainty of all numbers in a uniform
manner**.  This is what the :func:`nominal_value` and
:func:`std_dev` functions do:

  >>> print uncertainties.nominal_value(x)
  0.2
  >>> print uncertainties.std_dev(x)
  0.01
  >>> uncertainties.nominal_value(3)
  3
  >>> uncertainties.std_dev(3)
  0.0

Finally, a utility method is provided that directly yields the
**number of standard deviations** between a number and a result with
uncertainty: with ``x`` equal to 0.20±0.01,

  >>> x.position_in_sigmas(0.17)
  -3.0

.. index:: derivatives

.. _derivatives:

Derivatives
-----------

Since the application of :ref:`linear error propagation theory
<linear_method>` involves the calculation of **derivatives**, this
package automatically performs such calculations; users can thus
easily get the derivative of an expression with respect to any of its
variables:

  >>> u = ufloat((1, 0.1))
  >>> v = ufloat((10, 0.1))
  >>> sum_value = u+2*v
  >>> sum_value.derivatives[u]
  1.0
  >>> sum_value.derivatives[v]
  2.0

These values are obtained with a :ref:`fast differentiation algorithm
<differentiation method>`.

Additional information
----------------------

The capabilities of the :mod:`uncertainties` package in terms of array
handling are detailed in :doc:`numpy_guide`.

Details about the theory behind this package are given in the
:doc:`tech_guide`.

.. _NumPy: http://numpy.scipy.org/

