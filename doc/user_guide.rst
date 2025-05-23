.. index:: user guide
.. _user guide:

==========
User Guide
==========

The :mod:`uncertainties` package is built around the :class:`UFloat` class.
:class:`UFloat` objects model statistical random variables.
A :class:`UFloat` object has a nominal value which models the mean of a random variable
and an uncertainty which can be used to calculate the standard deviation of the random
variable and correlations between one random variable and another.
It is possible to construct new random variables by applying mathematical operations to
and between existing random variables.
Likewise, it is possible to apply mathematical operations such as arithmetic or
trigonometric functions to :class:`UFloat` objects to get new :class:`UFloat` objects.
The nominal values pass through the mathematical operations as if they were regular
:class:`float` objects, but the uncertainties propagate according to the approximate
rules of linear error propagation based on local derivatives of the mathematical
operations.

In addition to the :class:`UFloat` class, the :mod:`uncertainties` package also provides
sub-modules for :ref:`advanced mathematical functions <advanced math operations>`, and
:doc:`arrays and matrix <numpy_guide>` operations.

.. index::
   pair: number with uncertainty;

Creating UFloat objects: numbers with uncertainties
===================================================

To create a number with uncertainties or :class:`UFloat`, use the :func:`ufloat`
function, which takes a *nominal value* (which can be interpreted as the most
likely value, or the mean or central value of the distribution of values), a
*standard error* (the standard deviation or :math:`1-\sigma` uncertainty), and
an optional *tag*:

>>> x = ufloat(2.7, 0.01)  # x = 2.7+/-0.01
>>> y = ufloat(4.5,  1.2, tag='y_variable')  # x = 4..5+/-1.2

.. index::
   pair: nominal value; scalar
   pair: uncertainty; scalar

You can access the nominal value and standard deviation for any :class:`UFloat` object
with the `nominal_value` and `std_dev` attributes:

>>> print(x.nominal_value,  x.std_dev)
2.7 0.01


Because these are fairly long to type, for convenience,  `nominal_value` can be
abbreviated as `n` and `std_dev` as `s`:

>>> print(x.n,  x.s)
2.7 0.01

uncertainties :class:`UFloat` objects can also be created from one of many string
representations.  The following forms will all create :class:`UFloat` objects with
the same values:

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

Basic math with uncertain UFloat objects
========================================

Uncertainties :class:`UFloat` objects created using :func:`ufloat` or
:func:`ufloat_fromstr` can be used in basic mathematical calculations
(``+``, ``-``, ``*``, ``/``, ``**``)
as with other Python numbers.

>>> t = ufloat(0.2, 0.01)
>>> double = 2.0*t
>>> print(double)
0.400+/-0.020
>>> square = t**2
>>> print(square)
0.040+/-0.004

When adding two uncorrelated :class:`UFloat` objects, the standard deviation in the
resulting :class:`UFloat` object is the quadrature sum (square-root of the sum of
squares) of the standard deviations of the two input :class:`UFloat` objects.

>>> x = ufloat(20, 4)
>>> y = ufloat(12, 3)
>>> print(x+y)
32+/-5

Note that calls :func:`ufloat` create :class:`UFloat` objects which have no correlation
to any previously created :class:`UFloat` objects.
We can check the correctness of the error propagation for adding two uncorrelated
:class:`UFloat` objects:

>>> from math import sqrt
>>> (x+y).s == sqrt(x.s**2 + y.s**2)
True

We can also multiply two :class:`UFloat` objects:

>>> print(x*y)
(2.4+/-0.8)e+02
>>> (x*y).s == (x*y).n *  sqrt((x.s/x.n)**2 + (y.s/y.n)**2 )
True

However, consider the behavior when we add a :class:`UFloat` object to itself:

>>> print(x+x)
40+/-8
>>> print(3*x + 10)
70+/-12

In this case the resulting standard deviation is not the quadrature sum of the
standard deviation on the inputs.
This is because the :class:`UFloat` ``x`` is correlated with itself so the error
propagation must be handled accordingly.
This begins to demonstrate that calculations performed using :class:`UFloat` objects are
aware of correlations between :class:`UFloat` objects.
This is demonstrated in these two simple examples:

>>> x = ufloat(5, 0.5)
>>> y = ufloat(5, 0.5)
>>> x - y
0.0+/-0.7071067811865476
>>> x - x
0.0+/-0

.. index:: mathematical operation; on a scalar, umath

.. _advanced math operations:

Mathematical operations with uncertain UFloat objects
=====================================================

Besides being able to apply basic arithmetic operations to uncertainties
:class`UFloat` objects, this package provides generalized versions of 40 of the the
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

Comparison operators (``==``, ``!=``, ``>``, ``<``, ``>=``, and ``<=``) for
:class:`UFloat` objects with uncertainties are somewhat complicated, and need special
attention.  As we hinted at above, and will explore in more detail below and in the
:ref:`Technical Guide <comparison_operators>`, this relates to the correlation
between :class:`UFloat` objects.

Equality and inequality comparisons
------------------------------------

If we compare the equality of two :class:`UFloat` objects with the same nominal value
and standard deviation, we see

>>> x = ufloat(5, 0.5)
>>> y = ufloat(5, 0.5)
>>> x == x
True
>>> x == y
False

The difference here is that although the two :class:`UFloat` objects have the same
nominal value and standard deviation, they model two *uncorrelated* random variables.
This can be demonstrated with the :meth:`UFloat.covariance` method.

>>> print(x.covariance(x))
0.25
>>> print(x.covariance(y))
0

In order for the result of two calculations with uncertainties to be considered
equal, the :mod:`uncertainties` package does not test whether the nominal value
and the standard deviation have the same value.  Instead it checks whether the
difference of the two calculations has a nominal value of 0 *and* a standard deviation
of 0.

>>> (x -x)
0.0+/-0
>>> (x -y)
0.0+/-0.7071067811865476

 .. index::
   pair: testing (scalar); NaN


Handling NaNs and infinities
===============================

NaN values can appear in either the nominal value or uncertainty of a
Variable.  As is always the case, care must be exercised when handling NaN
values.

While :func:`math.isnan` and :func:`numpy.isnan` will raise `TypeError`
exceptions for uncertainties :class:`UFloat` objects (because an uncertainties
:class:`UFloat` object is not a float), the function :func:`umath.isnan` will return
whether the nominal value of a Variable is NaN.  Similarly, :func:`umath.isinf` will
return whether the nominal value of a Variable is infinite.

To check whether the uncertainty is NaN or Inf, use one of :func:`math.isnan`,
:func:`math.isinf`, :func:`nupmy.isnan`, or , :func:`nupmy.isinf` on the
``std_dev`` attribute.


.. index:: correlations; detailed example


Power Function Behavior
=======================

The value of one :class:`UFloat` raised to the power of another can be calculated in two
ways:

>>> from uncertainties import umath
>>>
>>> x = ufloat(4.5, 0.2)
>>> y = ufloat(3.4, 0.4)
>>> print(x**y)
(1.7+/-1.0)e+02
>>> print(umath.pow(x, y))
(1.7+/-1.0)e+02

The function ``x**y`` is defined for all ``x != 0`` and for ``x == 0`` as long as
``y > 0``.
There is not a unique definition for ``0**0``, however python takes the convention for
:class:`float` that ``0**0 == 1``.
If the power operation is performed on an ``(x, y)`` pair for which ``x**y`` is
undefined then an exception will be raised:

>>> x = ufloat(0, 0.2)
>>> y = ufloat(-3.4, 0.4)
>>> print(x**y)
Traceback (most recent call last):
 ...
ZeroDivisionError: 0.0 cannot be raised to a negative power

On the domain where it is defined, ``x**y`` is always real for ``x >= 0``.
For ``x < 0`` it is real for all integer values of ``y``.
If ``x<0`` and ``y`` is not an integer then ``x**y`` has a non-zero imaginary component.
The :mod:`uncertainties` module does not handle complex values:

>>> x = ufloat(-4.5, 0.2)
>>> y = ufloat(-3.4, 0.4)
>>> print(x**y)
Traceback (most recent call last):
 ...
ValueError: The uncertainties module does not handle complex results

The ``x`` derivative is real anywhere ``x**y`` is real except along ``x==0`` for
non-integer ``y``.
At these points the ``x`` derivative would be complex so a NaN value is used:

>>> x = ufloat(0, 0.2)
>>> y=1.5
>>> print((x**y).error_components())
{0.0+/-0.2: nan}

The ``y`` derivative is real anywhere ``x**y`` is real as long as ``x>=0``.
For ``x < 0`` the ``y`` derivative is always complex valued so a NaN value is used:

>>> x = -2
>>> y = ufloat(1, 0.2)
>>> print((x**y).error_components())
{1.0+/-0.2: nan}

Automatic correlations
======================

Correlations between :class:`UFloat` objects are **automatically handled** whatever
the number of :class:`UFloat` objects involved, and whatever the complexity of the
calculation. For example, when :data:`x` is the number with
uncertainty defined above,

>>> x = ufloat(0.2, 0.01)
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



UAtoms: How Uncertainty and Correlations are Tracked
====================================================

The basic, indivisibile, unit of uncertainty in the :mod:`uncertainties` package is the
:class:`UAtom`.
A :class:`UAtom` models a random variable with zero mean and unity standard deviation.
Every :class:`UAtom` is unique and uncorrelated with every other :class:`UAtom`.
The uncertainty of a :class:`UFloat` object is a :class:`UCombo` object which models a
:class:`float`-weighted linear combination of :class:`UAtom` objects.
A :class:`UFloat` object can be thought of as a sum of a fixed nominal value together
with a zero-mean :class:`UCombo` object.

The standard deviation of a single :class:`UFloat` object is calculated by taking the
sum-of-squares of the weights for all the :class:`UAtom` objects that make up the
corresponding :attribute:`uncertainty` attribute for that :class:`UFloat` object.
The correlation between two :class:`UFloat` objects is calculated by taking the sum
of products of the weights of shared :class:`UAtom` objects between the two
:class:`UFloat` :attribute:`uncertainty` attributes.

Every time a new :class:`UFloat` is instantiated via the :func:`ufloat` function a
single new independent :class:`UAtom` is also instantiated (and given the optional tag
passed into :func:`ufloat`) and paired with the new :class:`UFloat`.
When :class:`UFloat` objects are combined together using mathematical operations the
resulting :class:`UFloat` objects inherit dependence on the :class:`UAtom` objects
on which the input :class:`UFloat` objects depend in accordance with
:ref:`linear error propagation theory <linear_method>`.
In this way, the correlation between :class:`UFloat` objects can be tracked.

We can get access to the :class:`UAtom` objects on which a given :class:`UFloat`
depends, and their corresponding weights using the :attribute:`UFloat.error_components`
attribute:


.. testsetup:: uuid

   from uncertainties import ufloat
   from unittest.mock import patch
   import uuid
   import random

   class FakeUUID4:
       def __init__(self):
           self.seed = 0
           self.rnd = random.Random()

       def __call__(self):
           self.rnd.seed(self.seed)
           fake_uuid = uuid.UUID(int=self.rnd.getrandbits(128), version=4)
           self.seed += 1
           return fake_uuid
   fake_uuid4 = FakeUUID4()

   p = patch('uncertainties.ucombo.uuid.uuid4', fake_uuid4)
   p.start()

.. doctest:: uuid

   >>> x = ufloat(1, 0.1)
   >>> y = ufloat(2, 0.3)
   >>> z = x * y
   >>> print(x.error_components)
   {UAtom(e3e70682-c209-4cac-a29f-6fbed82c07cd): 0.1}
   >>> print(y.error_components)
   {UAtom(cd613e30-d8f1-4adf-91b7-584a2265b1f5): 0.3}
   >>> print(z.error_components)
   {UAtom(cd613e30-d8f1-4adf-91b7-584a2265b1f5): 0.3, UAtom(e3e70682-c209-4cac-a29f-6fbed82c07cd): 0.2}

The standard deviation of each :class:`UFloat` is given by the sum of squares of the
weights for all the :class:`UAtom` objects on which that :class:`UFloat` depends

.. doctest:: uuid

   >>> print(x.std_dev)
   0.1
   >>> print(y.std_dev)
   0.3
   >>> print(z.std_dev)
   0.36055512754639896

The :func:`ufloat` function accepts a ``tag`` argument.
If a string is passed in as the ``tag`` then this ``tag`` gets added to the new
:class:`UAtom` object that is instantiated together with the new :class:`UFloat`.
Note that :class:`UFloat` objects do not carry tags, only the underlying :class`UAtom`
objects do.
The tags on :class:`UAtom` objects can be used to keep track of the statistical
relationships in a more human-readable way:

.. doctest:: uuid

   >>> x = ufloat(1, 0.1, tag='x')
   >>> y = ufloat(2, 0.3, tag='y')
   >>> z = x * y
   >>>
   >>> for uatom, weight in z.error_components.items():
   ...     if uatom.tag is not None:
   ...         label = uatom.tag
   ...     else:
   ...         label = uatom.uuid
   ...     print(f"{label}: {weight}")
   y: 0.3
   x: 0.2


.. testcleanup :: uuid

   p.stop()

The tags *do not have to be distinct*. For instance, *multiple* :class:`UFloat` objects
can be tagged as ``"systematic"``, and their contribution to the total uncertainty of
:data:`result` can simply be obtained as:

>>> syst_error = math.sqrt(sum(  # Error from *all* systematic errors
...     error**2
...     for (uatom, error) in result.error_components().items()
...     if uatom.tag == "systematic"))

The remaining contribution to the uncertainty is:

>>> other_error = math.sqrt(result.std_dev**2 - syst_error**2)

The variance of :data:`result` is in fact simply the quadratic sum of these two errors,
since the :class:`UAtom` objects from :func:`result.error_components` are independent.

.. index:: comparison operators


.. index:: covariance matrix

Covariance and correlation matrices
===================================

Covariance matrix
-----------------

The covariance matrix between various variables or calculated
quantities can be simply obtained:

>>> from uncertainties import covariance_matrix
>>> sum_value = u+2*v
>>> cov_matrix = covariance_matrix([u, v, sum_value])

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

>>> from uncertainties import correlation_matrix
>>> corr_matrix = correlation_matrix([u, v, sum_value])
>>> print(corr_matrix)
[[1.         0.         0.4472136 ]
 [0.         1.         0.89442719]
 [0.4472136  0.89442719 1.        ]]

.. index:: correlations; correlated variables

Correlated variables
====================

Reciprocally, **correlated variables can be created** transparently,
provided that the NumPy_ package is available.

Use of a covariance matrix
--------------------------

Correlated variables can be obtained through the *covariance* matrix:

>>> from uncertainties import correlated_values
>>> (u2, v2, sum2) = correlated_values([1, 10, 21], cov_matrix)

creates three new variables with the listed nominal values, and the given
covariance matrix:

>>> print(sum_value)
21.00+/-0.22
>>> print(sum2)
21.00+/-0.22
>>> print(format(sum2 - (u2+2*v2), ".6f"))
0.000000+/-0.000000

The theoretical value of the last expression is exactly zero, like for
``sum - (u+2*v)``, but numerical errors yield a small uncertainty
(3e-9 is indeed very small compared to the uncertainty on :data:`sum2`:
correlations should in fact cancel the uncertainty on :data:`sum2`).

The covariance matrix is the desired one:

>>> import numpy as np
>>> print(np.array_str(np.array(covariance_matrix([u2, v2, sum2])), suppress_small=True))
[[0.01 0.   0.01]
 [0.   0.01 0.02]
 [0.01 0.02 0.05]]

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

>>> from uncertainties import correlated_values_norm
>>> (u3, v3, sum3) = correlated_values_norm(
...     [(1, 0.1), (10, 0.1), (21, 0.22360679774997899)],
...     corr_matrix,
... )
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

>>> from scipy.special import jv
>>> from uncertainties import wrap as u_wrap
>>> x = ufloat(2, 0.01)
>>> jv(0, x)
Traceback (most recent call last):
 ...
TypeError: ufunc 'jv' not supported for the input types, and the inputs could not be safely coerced to any supported types according to the casting rule ''safe''
>>> print(u_wrap(jv)(0, x))
0.224+/-0.006

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

>>> from uncertainties import nominal_value, std_dev
>>> x = ufloat(0.2, 0.01)
>>> print(nominal_value(x))
0.2
>>> print(std_dev(x))
0.01
>>> print(nominal_value(3))
3
>>> print(std_dev(3))
0.0

Finally, a utility method is provided that directly yields the
`standard score <http://en.wikipedia.org/wiki/Standard_score>`_
(number of standard deviations) between a number and a result with
uncertainty:

>>> x = ufloat(0.20, 0.01)
>>> print(x.std_score(0.17))
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
