.. index:: user guide
.. _user guide:

==========
User Guide
==========


Basic usage
===========

Basic mathematical operations involving numbers with uncertainties requires
importing the :func:`ufloat` function which creates a :class:`Variable`:
number with both a nominal value and an uncertainty.

     >>> from uncertainties import ufloat
     >>> x = ufloat(2.7, 0.01)   #  a Variable with a value 2.7+/-0.01

The :mod:`uncertainties` module contains sub-modules for :ref:`advanced
mathematical functions <advanced math operations>`, and :doc:`arrays and
matrices <numpy_guide>`, which can be accessed as well.

.. index::
   pair: number with uncertainty; creation

Creating Variables: numbers with uncertainties
================================================

To create a number with uncertainties or *Variable*, use the :func:`ufloat`
function, which takes a *nominal value* (which can be interpreted as the most
likely value, or the mean or central value of the distribution of values), a
*standard error* (the standard deviation or :math:`1-\sigma` uncertainty), and
an optional *tag*:

>>> x = ufloat(2.7, 0.01)  # x = 2.7+/-0.01
>>> y = ufloat(4.5,  1.2, tag='y_variable')  # x = 4..5+/-1.2

.. index::
   pair: nominal value; scalar
   pair: uncertainty; scalar

You can access the nominal value and standard deviation for any Variable with
the `nominal_value` and `std_dev` attributes:

>>> print(x.nominal_value,  x.std_dev)
2.7 0.01


Because these are fairly long to type, for convenience,  `nominal_value` can be
abbreviated as `n` and `std_dev` as `s`:

>>> print(x.n,  x.s)
2.7 0.01

uncertainties Variables can also be created from one of many string
representations.  The following forms will all create Variables with the same
value:

>>> from uncertainties import ufloat_fromstr
>>> x = ufloat(0.2, 0.01)
>>> x = ufloat_fromstr("0.20+/-0.01")
>>> x = ufloat_fromstr("(2+/-0.1)e-01")  # Factored exponent
>>> x = ufloat_fromstr("0.20(1)")  # Short-hand notation
>>> x = ufloat_fromstr("20(1)e-2")  # Exponent notation
>>> x = ufloat_fromstr(u"0.20±0.01")  # Pretty-print form
>>> x = ufloat_fromstr("0.20")  # Automatic uncertainty of +/-1 on last digit

More details on the :func:`ufloat` and :func:`ufloat_from_str` can be found in
:ref:`api_funcs`.

Basic math with uncertain Variables
=========================================

Uncertainties variables created in :func:`ufloat` or :func:`ufloat_fromstr` can
be used in basic mathematical calculations (``+``, ``-``, ``*``, ``/``, ``**``)
as with other Python numbers and variables.

>>> t = ufloat(0.2, 0.01)
>>> double = 2.0*t
>>> print(double)
0.4+/-0.02
>>> square = t**2
>>> print(square)
0.040+/-0.004

When adding two Variables, the uncertainty in the result is the quadrature sum
(square-root of the sum of squares) of the uncertainties of the two Variables:

>>> x = ufloat(20, 4)
>>> y = ufloat(12, 3)
>>> print(x+y)
32.0+/-5.0

We can check that error propagation when adding two independent variables
(using the abbreviation `.s` for the standard error):

>>> from math import sqrt
>>> (x+y).s == sqrt(x.s**2 + y.s**2)
True


Multiplying two Variables will properly propagate those
uncertainties too:

>>> print(x*y)
240.0+/-76.83749084919418
>>> (x*y).s == (x*y).n *  sqrt((x.s/x.n)**2 + (y.s/y.n)**2 )
True

But note that adding a Variable to itself does not add its uncertainties in
quadrature, but are simply scaled:

>>> print(x+x)
40.0+/-8.0
>>> print(3*x + 10)
70.0+/-12.0


It is important to understand that calculations done with Variable know about
the correlation between the Variables.   Variables created with :func:`ufloat`
(and  :func:`ufloat_fromstr`) are completely uncorrelated with each other, but
are known to be completely correlated with themselves.  This means that


>>> x = ufloat(5, 0.5)
>>> y = ufloat(5, 0.5)
>>> x - y
0.0+/-0.7071067811865476
>>> x - x
0.0+/-0

For two *different* Variables, uncorrelated uncertainties will be propagated.
But when doing a calculation with a single Variable, the uncertainties are
correlated, and calculations will reflect that.


.. index:: mathematical operation; on a scalar, umath

.. _advanced math operations:

Mathematical operations with uncertain Variables
=====================================================

Besides being able to apply basic mathematical operations to uncertainties
Variables, this package provides generalized versions of 40 of the the
functions from the standard :mod:`math` *module*.  These mathematical functions
are found in the :mod:`uncertainties.umath` module:

    >>> from uncertainties.umath import sin, exp, sqrt
    >>> x   = ufloat(0.2, 0.01)
    >>> sin(x)
    0.19866933079506122+/-0.009800665778412416
    >>> sin(x*x)
    0.03998933418663417+/-0.003996800426643912
    >>> exp(-x/3.0)
    0.9355069850316178+/-0.003118356616772059
    >>> sqrt(230*x + 3)
    7.0+/-0.16428571428571428


The functions in the :mod:`uncertainties.umath` module include:

    ``acos``, ``acosh``, ``asin``, ``asinh``, ``atan``, ``atan2``, ``atanh``,
    ``ceil``, ``copysign``, ``cos``, ``cosh``, ``degrees``, ``erf``, ``erfc``,
    ``exp``, ``expm1``, ``fabs``, ``factorial``, ``floor``, ``fmod``,
    ``frexp``, ``fsum``, ``gamma``, ``hypot``, ``isinf``, ``isnan``,
    ``ldexp``, ``lgamma``, ``log``, ``log10``, ``log1p``, ``modf``,
    ``pow``, ``radians``, ``sin``, ``sinh``, ``sqrt``, ``tan``, ``tanh``, ``trunc``


Comparison operators
====================

Comparison operators (``==``, ``!=``, ``>``, ``<``, ``>=``, and ``<=``) for Variables with
uncertainties are somewhat complicated, and need special attention.  As we
hinted at above, and will explore in more detail below and in the
:ref:`Technical Guide <comparison_operators>`, this relates to the correlation
between Variables.



Equality and inequality comparisons
------------------------------------

If we compare the equality of two Variables with the same nominal value and
uncertainty, we see

>>> x = ufloat(5, 0.5)
>>> y = ufloat(5, 0.5)
>>> x == x
True
>>> x == y
False

The difference here is that although the two Python objects have the same
nominal value and uncertainty, these are independent, uncorrelated values.  It
is not exactly true that the difference is based on identity, note that

>>> x == (1.0*x)
True
>>> x is (1.0*x)
False

In order for the result of two calculations with uncertainties to be considered
equal, the :mod:`uncertainties` package does not test whether the nominal value
and the uncertainty have the same value.  Instead it checks whether the
difference of the two calculations has a nominal value of 0 *and* an
uncertainty of 0.

>>> (x -x)
0.0+/-0
>>> (x -y)
0.0+/-0.7071067811865476


Comparisons of magnitude
------------------------------------

The concept of comparing the magnitude of values with uncertainties is a bit
complicated.  That is, a Variable with a value of 25 +/- 10 might be greater
than a Variable with a value of 24 +/- 8 most of the time, but *sometimes* it
might be less than it.   The :mod:`uncertainties` package takes the simple
approach of comparing nominal values.  That is

>>> a = ufloat(25, 10)
>>> b = ufloat(24, 8)
>>> a > b
True

Note that combining this comparison and the above discussion of `==` and `!=`
can lead to a result that maybe somewhat surprising:


>>> a = ufloat(25, 10)
>>> b = ufloat(25, 8)
>>> a >= b
False
>>> a > b
False
>>> a == b
False
>>> a.nominal_value >= b.nominal_value
True

That is, since `a` is neither greater than `b` (nominal value only) nor equal to
`b`, it cannot be greater than or equal to `b`.


 .. index::
   pair: testing (scalar); NaN


Handling NaNs and infinities
===============================

NaN values can appear in either the nominal value or uncertainty of a
Variable.  As is always the case, care must be exercised when handling NaN
values.

While :func:`math.isnan` and :func:`numpy.isnan` will raise `TypeError`
exceptions for uncertainties Variables (because an uncertainties Variable is
not a float), the function :func:`umath.isnan` will return whether the nominal
value of a Variable is NaN.  Similarly, :func:`umath.isinf` will return whether
the nominal value of a Variable is infinite.

To check whether the uncertainty is NaN or Inf, use one of :func:`math.isnan`,
:func:`math.isinf`, :func:`nupmy.isnan`, or , :func:`nupmy.isinf` on the
``std_dev`` attribute.


.. index:: correlations; detailed example


Automatic correlations
======================

Correlations between variables are **automatically handled** whatever
the number of variables involved, and whatever the complexity of the
calculation. For example, when :data:`x` is the number with
uncertainty defined above,

>>> square = x**2
>>> print(square)
0.040+/-0.004
>>> square - x*x
0.0+/-0
>>> y = x*x + 1
>>> y - square
1.0+/-0

The last two printed results above have a zero uncertainty despite the
fact that :data:`x`, :data:`y` and :data:`square` have a non-zero uncertainty: the
calculated functions give the same value for all samples of the random
variable :data:`x`.

Thanks to the automatic correlation handling, calculations can be
performed in as many steps as necessary, exactly as with simple
floats.  When various quantities are combined through mathematical
operations, the result is calculated by taking into account all the
correlations between the quantities involved.  All of this is done
transparently.



Access to the individual sources of uncertainty
===============================================

The various contributions to an uncertainty can be obtained through the
:func:`error_components` method, which maps the **independent variables
a quantity depends on** to their **contribution to the total
uncertainty**. According to :ref:`linear error propagation theory
<linear_method>` (which is the method followed by :mod:`uncertainties`),
the sum of the squares of these contributions is the squared
uncertainty.

The individual contributions to the uncertainty are more easily usable
when the variables are **tagged**:

>>> u = ufloat(1, 0.1, "u variable")  # Tag
>>> v = ufloat(10, 0.1, "v variable")
>>> sum_value = u+2*v
>>> sum_value
21.0+/-0.223606797749979
>>> for (var, error) in sum_value.error_components().items():
...     print("{}: {}".format(var.tag, error))
...
u variable: 0.1
v variable: 0.2

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

>>> other_error = math.sqrt(result.std_dev**2 - syst_error**2)

The variance of :data:`result` is in fact simply the quadratic sum of
these two errors, since the variables from
:func:`result.error_components` are independent.

.. index:: comparison operators


.. index:: covariance matrix

Covariance and correlation matrices
===================================

Covariance matrix
-----------------

The covariance matrix between various variables or calculated
quantities can be simply obtained:

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
0.0+/-0

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

Correlated variables can be obtained through the *covariance* matrix:

>>> (u2, v2, sum2) = uncertainties.correlated_values([1, 10, 21], cov_matrix)

creates three new variables with the listed nominal values, and the given
covariance matrix:

>>> sum_value
21.0+/-0.223606797749979
>>> sum2
21.0+/-0.223606797749979
>>> sum2 - (u2+2*v2)
0.0+/-3.83371856862256e-09

The theoretical value of the last expression is exactly zero, like for
``sum - (u+2*v)``, but numerical errors yield a small uncertainty
(3e-9 is indeed very small compared to the uncertainty on :data:`sum2`:
correlations should in fact cancel the uncertainty on :data:`sum2`).

The covariance matrix is the desired one:

>>> uncertainties.covariance_matrix([u2, v2, sum2])

reproduces the original covariance matrix :data:`cov_matrix` (up to
rounding errors).

Use of a correlation matrix
---------------------------

Alternatively, correlated values can be defined through:

- a sequence of nominal values and standard deviations, and
- a *correlation* matrix between each variable of this sequence
  (the correlation matrix is the covariance matrix
  normalized with individual standard deviations; it has ones on its
  diagonal)—in the form of a NumPy array-like object, e.g. a
  list of lists, or a NumPy array.

Example:

>>> (u3, v3, sum3) = uncertainties.correlated_values_norm(
...     [(1, 0.1), (10, 0.1), (21, 0.22360679774997899)], corr_matrix)
>>> print(u3)
1.00+/-0.10

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

It is thus possible to take a function :func:`f` *that returns a
single float*, and to automatically generalize it so that it also
works with numbers with uncertainties:

>>> wrapped_f = uncertainties.wrap(f)

The new function :func:`wrapped_f` (optionally) *accepts a number
with uncertainty* in place of any float *argument* of :func:`f` (note
that floats contained instead *inside* arguments of :func:`f`, like
in a list or a NumPy array, *cannot* be replaced by numbers with
uncertainties).
:func:`wrapped_f` returns the same values as :func:`f`, but with
uncertainties.

With a simple wrapping call like above, uncertainties in the function
result are automatically calculated numerically. **Analytical
uncertainty calculations can be performed** if derivatives are
provided to :func:`wrap`.


Miscellaneous utilities
=======================

.. index:: standard deviation; on the fly modification

It is sometimes useful to modify the error on certain parameters so as
to study its impact on a final result.  With this package, the
**uncertainty of a variable can be changed** on the fly:

>>> sum_value = u+2*v
>>> sum_value
21.0+/-0.223606797749979
>>> prev_uncert = u.std_dev
>>> u.std_dev = 10
>>> sum_value
21.0+/-10.00199980003999
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

>>> print(uncertainties.nominal_value(x))
0.2
>>> print(uncertainties.std_dev(x))
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

>>> u = ufloat(1, 0.1)
>>> v = ufloat(10, 0.1)
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

.. |minus2html| raw:: html

   <sup>-2</sup>
