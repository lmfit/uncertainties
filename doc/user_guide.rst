.. index:: user guide
.. _user guide:

==========
User Guide
==========


Basic setup
===========

Basic mathematical operations involving numbers with uncertainties
only require a simple import:

  >>> from uncertainties import ufloat

The :func:`ufloat` function creates numbers with uncertainties. Existing 
calculation code can usually run with no or little modification and 
automatically produce results with uncertainties.

.. The "import uncertainties" is put here because some examples requires
   uncertainties to have been imported (and not only ufloat).

The :mod:`uncertainties` module contains other features, which can be
made accessible through

  >>> import uncertainties

The :mod:`uncertainties` package also contains sub-modules for
:ref:`advanced mathematical functions <advanced math operations>`, and
:doc:`arrays and matrices <numpy_guide>`.

.. index::
   pair: number with uncertainty; creation

Creating and handling numbers with uncertainties
================================================

Numbers with uncertainties can be input either numerically, or through
one of many string representations, so that files containing numbers
with uncertainties can easily be parsed.  Thus, x = 0.20±0.01 can be
expressed in many convenient ways:

  >>> x = ufloat((0.20, 0.01))  # x = 0.20+/-0.01
  >>> x = ufloat("0.20+/-0.01")
  >>> x = ufloat("0.20(1)")
  >>> x = ufloat("20(1)e-2")  # Exponential notation supported
  >>> x = ufloat("0.20")  # Automatic uncertainty of +/-1 on last digit

The available representations can be listed with ``pydoc
uncertainties.ufloat``.  Representations that are invalid raise a
:class:`ValueError` exception.

Basic math
==========

Calculations can be performed directly, as with regular real numbers:

  >>> square = x**2
  >>> print square
  0.04+/-0.004


.. index:: mathematical operation; on a scalar, umath

.. _advanced math operations:

Mathematical operations
=======================

Besides being able to apply basic mathematical operations to numbers
with uncertainty, this package provides generalizations of most of the
functions from the standard :mod:`math` module.  These mathematical
functions are found in the :mod:`uncertainties.umath` module::

  >>> from uncertainties.umath import *  # Imports sin(), etc.
  >>> sin(x**2)
  0.039989334186634168+/-0.003996800426643912

The list of available mathematical functions can be obtained with the
``pydoc uncertainties.umath`` command.

.. index:: arrays; simple use, matrices; simple use

.. _simple_array_use:

Arrays of numbers with uncertainties
====================================

It is possible to put numbers with uncertainties in NumPy_ arrays and
matrices:

  >>> arr = numpy.array([ufloat((1, 0.01)), ufloat((2, 0.1))])
  >>> 2*arr
  [2.0+/-0.02 4.0+/-0.2]
  >>> print arr.sum()
  3.0+/-0.100498756211

Thus, usual operations on NumPy arrays can be performed transparently
even when these arrays contain numbers with uncertainties.

:doc:`More complex operations on NumPy arrays and matrices 
<numpy_guide>` can be
performed through the dedicated :mod:`uncertainties.unumpy` module.

.. index:: correlations; detailed example


Correlations
============

Correlations between variables are automatically handled whatever the
number of variables involved, and whatever the complexity of the
calculation.  Thus, each calculation result keeps track of how it is
correlated to random variables.  For example, when :data:`x` is the number
with uncertainty defined above,

  >>> square = x**2
  >>> print square
  0.04+/-0.004
  >>> square - x*x
  0.0
  >>> y = x*x + 1
  >>> y - square
  1.0

The last two printed results above have a zero uncertainty despite the
fact that :data:`x`, :data:`y` and :data:`square` have a non-zero uncertainty: the
calculated functions give the same value for all samples of the random
variable :data:`x`.

Thanks to the tracking of dependencies on random variables,
calculations can therefore be performed in as many steps as necessary,
exactly as with simple floats.  When various quantities are combined
through mathematical operations, the result is calculated by taking
into account all the correlations between the quantities involved.
All of this is done completely transparently.

.. index::
   pair: nominal value; of scalar
   pair: uncertainty; of scalar

Access to the uncertainty and to the nominal value
==================================================

The nominal value and the uncertainty (standard deviation) can also be
accessed independently:

  >>> print square
  0.04+/-0.004
  >>> print square.nominal_value
  0.04
  >>> print square.std_dev
  0.004

Access to the individual sources of uncertainty
===============================================

The various contributions to an uncertainty can be obtained through
the :func:`error_components` method, which maps the **independent
variables a quantity depends on** to their **contribution to the total
uncertainty**. According to the :ref:`linear error propagation theory
<linear_method>` implemented in :mod:`uncertainties`, the sum of the
squares of these contributions is the squared uncertainty.

The individual contributions to the uncertainty are more easily usable
when the variables are **tagged**:

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

The variance (i.e. squared uncertainty) of the result
(:data:`sum_value`) is the quadratic sum of these independent
uncertainties, as it should be (``0.1**2 + 0.2**2``).

The tags *do not have to be distinct*. For instance, *multiple* random
variables can be tagged as ``"systematic"``, and their contribution to
the total uncertainty of :data:`result` can simply be obtained as:

  >>> syst_error = math.sqrt(sum(  # Error from *all* systematic errors
  ...     error**2
  ...     for (var, error) in result.error_components().items()
  ...     if var.tag == "systematic"))
          
The remaining contribution to the uncertainty is:

  >>> other_error = math.sqrt(result.std_dev()**2 - syst_error**2)

The variance of :data:`result` is in fact simply the quadratic sum of
these two errors, since the variables from
:func:`result.error_components` are independent.

.. index:: comparison operators

Comparison operators
====================

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

One important concept to keep in mind is that :func:`ufloat` creates a
random variable, so that two numbers with the same nominal value and
standard deviation are generally different:

  >>> y = ufloat((1, 0.1))
  >>> z = ufloat((1, 0.1))
  >>> print y
  1.0+/-0.1
  >>> print z
  1.0+/-0.1
  >>> y == y
  True
  >>> y == z
  False

In physical terms, two rods of the same nominal length and uncertainty
on their length are generally of different sizes: :data:`y` is different
from :data:`z`.

More detailed information on the semantics of comparison operators for
numbers with uncertainties can be found in the :ref:`Technical Guide
<comparison_operators>`.


.. index:: covariance matrix

Covariance and correlation matrices
=================

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

In this matrix, the zero covariances indicate that :data:`u` and :data:`v` are
independent from each other; the last column shows that :data:`sum_value`
does depend on these variables.  The :mod:`uncertainties` package
keeps track at all times of all correlations between quantities
(variables and functions):

  >>> sum_value - (u+2*v)
  >>> 0.0

Correlation matrix
------------------

If the NumPy_ package is available, the correlation matrix can be
obtained as well:

  >>> corr_matrix = uncertainties.correlation_matrix([u, v, sum_value])
  >>> corr_matrix
  array([[ 1.        ,  0.        ,  0.4472136 ],
         [ 0.        ,  1.        ,  0.89442719],
         [ 0.4472136 ,  0.89442719,  1.        ]])
  
.. index:: correlations; correlated variables

Correlated variables
====================

Reciprocally, **correlated variables can be created** transparently,
provided that the NumPy_ package is available.

Use of a covariance matrix
--------------------------

Correlated variables can be obtained through the *covariance* matrix::

  >>> (u2, v2, sum2) = uncertainties.correlated_values([1, 10, 21], cov_matrix)

creates three new variables with the listed nominal values, and the given
covariance matrix::

  >>> sum_value
  21.0+/-0.22360679774997899
  >>> sum2
  21.0+/-0.22360679774997899
  >>> sum2 - (u2+2*v2)
  0.0+/-3.8337185686225597e-09

The theoretical value of the last expression is exactly zero, like for
``sum - (u+2*v)``, but numerical errors yield a small uncertainty
(3e-9 is indeed very small compared to the uncertainty on :data:`sum2`:
correlations should in fact cancel the uncertainty on :data:`sum2`).

The covariance matrix is the desired one::

  >>> uncertainties.covariance_matrix([u2, v2, sum2])

reproduces the original covariance matrix :data:`cov_matrix` (up to
rounding errors).

Use of a correlation matrix
---------------------------

Alternatively, correlated values can be defined through a
*correlation* matrix (the correlation matrix is the covariance matrix
normalized with individual standard deviations; it has ones on its
diagonal), along with a list of nominal values and standard deviations::

  >>> (u3, v3, sum3) = uncertainties.correlated_values_norm(
  ...     [(1, 0.1), (10, 0.1), (21, 0.22360679774997899)], corr_matrix)
  >>> print u3
  1.0+/-0.1

The three returned numbers with uncertainties have the correct
uncertainties and correlations (:data:`corr_matrix` can be recovered
through :func:`correlation_matrix`).

.. index::
   single: C code; wrapping
   single: Fortran code; wrapping
   single: wrapping (C, Fortran,…) functions

Making custom functions accept numbers with uncertainties
=========================================================

This package allows **code which is not meant to be used with numbers
with uncertainties to handle them anyway**. This is for instance
useful when calling external functions (which are out of the user's
control), including functions written in C or Fortran.  Similarly,
**functions that do not have a simple analytical form** can be
automatically wrapped so as to also work with arguments that contain
uncertainties.

It is thus possible to have a function :func:`f` that returns a single
float, and to automatically generalize it so that it also works with
numbers with uncertainties:

  >>> wrapped_f = uncertainties.wrap(f)

The new function :func:`wrapped_f` accepts numbers with uncertainties
as arguments wherever a Python float was used for :func:`f`.
:func:`wrapped_f` returns the same values as :func:`f`, but with
uncertainties.

With a simple wrapping call like above, uncertainties in the function
result are automatically calculated numerically. **Analytical
uncertainty calculations can be performed** if derivatives are
provided to :func:`wrap` (for details, see the documentation string of
:func:`wrap` with the ``pydoc`` command, or :func:`help`).

Miscellaneous utilities
=======================

.. index:: standard deviation; on the fly modification

It is sometimes useful to modify the error on certain parameters so as
to study its impact on a final result.  With this package, the
**uncertainty of a variable can be changed** on the fly:

  >>> sum_value = u+2*v
  >>> sum_value
  21.0+/-0.22360679774997899
  >>> prev_uncert = u.std_dev
  >>> u.std_dev = 10
  >>> sum_value
  21.0+/-10.001999800039989
  >>> u.std_dev = prev_uncert

The relevant concept is that :data:`sum_value` does depend on the
variables :data:`u` and :data:`v`: the :mod:`uncertainties` package keeps
track of this fact, as detailed in the :ref:`Technical Guide
<variable_tracking>`, and uncertainties can thus be updated at any time.

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
`standard score <http://en.wikipedia.org/wiki/Standard_score>`_
(number of standard deviations) between a number and a result with
uncertainty: with :data:`x` equal to 0.20±0.01,

  >>> x.std_score(0.17)
  -3.0

.. index:: derivatives

.. _derivatives:

Derivatives
===========

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
======================

The capabilities of the :mod:`uncertainties` package in terms of array
handling are detailed in :doc:`numpy_guide`.

Details about the theory behind this package and implementation 
information are given in the
:doc:`tech_guide`.

.. _NumPy: http://numpy.scipy.org/
