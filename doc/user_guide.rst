from uncertainties import correlation_matrixfrom uncertainties import covariance_matrix==========
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
:doc:`arrays <numpy_guide>` operations.

.. index::
   pair: number with uncertainty; creation

Creating UFloat Objects
=======================

:class:`UFloat` objects are directly instantiated by passing in a :class:`float`
*nominal value* and positive :class:`float` *standard deviation*.
The nominal value can be interpreted as the most likely, mean, or central value of a
distribution of possible outcomes for the random variables, while the standard deviation
is the standard deviation (the :math:`1-sigma` uncertainty range) for the distribution.
An optional *tag* can be passed which can be used to keep track of different sources of
uncertainty.

>>> from uncertainties import UFloat
>>> x = UFloat(2.7, 0.01)
>>> print(x)
2.700+/-0.010
>>> y = UFloat(4.5,  1.2, tag='y_variable')
>>> print(y)
4.5+/-1.2

.. index::
   pair: nominal value; scalar
   pair: uncertainty; scalar

You can access the nominal value and standard deviation for any :class:`UFloat` object
using the `nominal_value` and `std_dev` attributes:

>>> print(x.nominal_value,  x.std_dev)
2.7 0.01


Because these are fairly long to type, for convenience, `nominal_value` can be
abbreviated as `n` and `std_dev` as `s`:

>>> print(x.n,  x.s)
2.7 0.01

uncertainties :class:`UFloat` objects can also be created from one of many string
representations.  The following forms will all create :class:`UFloat` objects with
the same values:

>>> from uncertainties import ufloat_fromstr
>>> x = UFloat(0.2, 0.01)
>>> x = ufloat_fromstr("0.20+/-0.01")
>>> x = ufloat_fromstr("(2+/-0.1)e-01")  # Factored exponent
>>> x = ufloat_fromstr("0.20(1)")  # Shorthand notation
>>> x = ufloat_fromstr("20(1)e-2")  # Exponent notation
>>> x = ufloat_fromstr(u"0.20±0.01")  # Pretty-print form
>>> x = ufloat_fromstr("0.20")  # Automatic uncertainty of +/-1 on last digit

Historically :class:`UFloat` objects were primary constructed using the :func:`ufloat`
factory method:

>>> from uncertainties import ufloat
>>> x = ufloat(2.7, 0.01)
>>> print(x)
2.700+/-0.010
>>> y = ufloat(4.5,  1.2, tag='y_variable')
>>> print(y)
4.5+/-1.2

However, it is now encouraged to instantiate :class:`UFloat` objects directly using the
class constructor.

More details on the :class:`UFloat` class and :func:`ufloat_from_str` :func:`ufloat`
functions can be found in :ref:`api_funcs`.

Basic Arithmetic with UFloat Objects
====================================

Uncertainties :class:`UFloat` objects can be used in basic mathematical calculations
(``+``, ``-``, ``*``, ``/``, ``**``)
as with other Python numbers.

>>> x = UFloat(0.2, 0.01)
>>> print(2 * x)
0.400+/-0.020
>>> print(x**2)
0.040+/-0.004
>>> y = UFloat(0.1, 0.02)
>>> print(x + y)
0.300+/-0.022
>>> print(x - y)
0.100+/-0.022
>>> print(x * y)
0.020+/-0.004


So we see that we can perform basic mathematical operations between :class:`UFloat` and
:class:`float` objects and also between two :class:`UFloat` objects.
We can also see that :mod:`uncertainties` handles the mathematical propagation of the
uncertainty to the final result.

.. _linear_uncertainty_math:

Linear Uncertainty Propagation
==============================

The :mod:`uncertainties` package uses :class:`UFloat` objects apply the theory of
linear error propagation.
Suppose ``A`` and ``B`` are real random variables which can be expressed as::

   A = A_0 + w_Ax dx + w_Ay dy = A_0 + dA
   B = B_0 + w_By dy + w_Bz dz = B_0 + dB


Here ``A_0`` and ``B_0`` are just real numbers, ``dx``, ``dy``, and ``dz`` are
independent, zero mean, unity variance random variables, and ``w_Ax``, ``w_Ay``,
``w_By`` and ``w_Bz`` are positive real number weights.
Since ``dx``, ``dy`` and ``dz`` are zero mean we can see that ``A_0`` and ``B_0`` are
the means of the random variables ``A`` and ``B`` respectively.

Because ``dx``, ``dy`` and ``dz`` have unity variance, it is easy to calculate the variance
of ``A`` and ``B`` by::

   Var(A) = w_Ax^2 + w_Ay^2
   Var(B) = w_By^2 + w_Bz^2


The theory of linear error propagation allows us to calculate the uncertainty on random
variable ``C = f(A, B)`` of the random variables ``A`` and ``B``::

   C = f(A, B) = f(A_0, B_0) + df/dA dA + df/dB dB
               = f(A_0, B_0) + df/dA w_Ax dx + df/dA w_Ay dy + df/dB w_By dy + df/dB w_Bz dz
               = f(A_0, B_0) + df/dA w_Ax dx + (df/dA w_Ay + df/dB w_By) dy + df/dB w_Bz dz
               = C_0 + dC


From this point we could calculate the variance and standard deviation of ``f(A, B)``,
but we will skip that calculation here.
Here, we will simply observe how, using a simple Taylor expansion, we can track the
dependence of ``C = f(A, B)`` on the random variables ``dx``, ``dy``, and ``dz`` on
which ``A`` and ``B`` depend.
We will also note that ``C = f(A, B)`` has dependence on ``dy`` due to both ``A`` and
``B``.
In other words, ``A`` and ``B`` have non-zero correlation and a proper uncertainty
propagation calculation of ``C`` must take this correlation into account.

Error Components, `UAtom` Objects, and Uncertainty Propagation
==============================================================

We can begin to see how the :mod:`uncertainties` modules performs linear uncertainty
propagation by studying the `UFloat.error_components` property.
A :class:`UFloat` object is like a random variable

.. doctest::
   :hide:

   >>> import random
   >>> random.seed(42)


>>> A = UFloat(10, 0.1, tag="A special tag")
>>> print(A.n)
10.0

``A`` is the random variable and ``A.n == 10.0``, the nominal value, is like the mean
of the random variable ``A``.
The :class:`UFloat` object has an :attr:`error_components` property

>>> print(A.error_components)
{UAtom(1c80317fa3b1799d, tag="A special tag")): 0.1}

We see that the :attr:`error_components` property returns a dict whose keys are
:class:`UAtom` objects and whose values are floats.
A :class:`UAtom` object is like the ``dx``, ``dy``, or ``dz`` random variables above.
It is like an independent random variable with zero mean and unity variance.
Whenver a new :class:`UFloat` object is directly instantiated, a new  :class:`UAtom`
is generated with a unique identifer.
Let's study the single :class:`UAtom` object responsible for the uncertainty in ``A``:

>>> single_uatom = list(A.error_components.keys())[0]
>>> print(type(single_uatom))
<class 'uncertainties.ucombo.UAtom'>
>>> print(type(single_uatom.uuid))
<class 'int'>
>>> print(format(single_uatom.uuid, 'x'))
1c80317fa3b1799d
>>> print(single_uatom.tag)
A special tag

We see that the :class:`UAtom` object has an integer :attr:`uuid` attribute which
appears in hex format in the :class:`UAtom` object's string representations.
There is also an optional ``str`` :attr:`tag` attribute.
We will demonstrate usage of the :attr:`tag` attribute below, but for now, it is
important to know that the :attr:`tag` attribute is not unique between :class:`UAtom`
instances and it in no way replaces the :attr:`uuid` attribute.

We can now reproduce the manipulations in the :ref:`linear_uncertainty_math` section.

>>> dx = UFloat(0, 1)
>>> print(dx.error_components)
{UAtom(bdd640fb06671ad1): 1.0}
>>> dy = UFloat(0, 1)
>>> print(dy.error_components)
{UAtom(3eb13b9046685257): 1.0}
>>> dz = UFloat(0, 1)
>>> print(dz.error_components)
{UAtom(23b8c1e9392456de): 1.0}

Note that we are defining :class:`UFloat` objects with the names ``dx``, ``dy``, and
``dz``, but we should really think of the corresponding :class:`UAtom` objects as the
units of uncertainty.

>>> A0 = 10
>>> A = A0 + 0.1 * dx + 0.2 * dy
>>> print(A.error_components)
{UAtom(3eb13b9046685257): 0.2, UAtom(bdd640fb06671ad1): 0.1}
>>> B0 = 20
>>> B = B0 + 0.3 * dy + 0.4 * dz
>>> print(B.error_components)
{UAtom(23b8c1e9392456de): 0.4, UAtom(3eb13b9046685257): 0.3}

Here we see that ``A`` and ``B`` each contain the appropriate weighting of the
corresponding :class:`UAtom` objects.
Now suppose ``C = f(A, B) = A * B``.
Then

>>> C = A * B
>>> print(C.n)
200.0

Thinking about the error components of ``C``, we expect that ``C`` has dependence on
the ``dx``, ``dy``, and ``dz`` :class:`UAtom` objects.
The dependence of ``C`` on ``dx` only comes through ``A``.
So we expect the weight for ``dx`` on ``C`` to be ``df/dA = B_0 = 20`` times the weight
of ``dx`` on ``A``, 0.1.
So we expect the total weight to be 2.

>>> print(C.error_components)
{UAtom(23b8c1e9392456de): 4.0, UAtom(3eb13b9046685257): 7.0, UAtom(bdd640fb06671ad1): 2.0}

Indeed, this is what we find.
The reader can verify the dependence on ``dz`` using a similar calculation.
The dependence on ``dy`` can also be verified, but, this time it is necessary to take
into account the fact that both ``A`` and ``B`` depend on ``dy``.
The :mod:`uncertainties` packages, can, of course, easily report the total standard
deviation of ``C`` given its error components:

>>> print(C.s)
8.306623862918075

This bookkeeping makes it easy for the :mod:`uncertainties` package to report the
`covariance <https://en.wikipedia.org/wiki/Covariance>`_ and
`correlation <https://en.wikipedia.org/wiki/Correlation>`_
between two :class:`UFloat` objects

>>> print(A.covariance(A))
0.05000000000000001
>>> print(A.covariance(B))
0.06
>>> print(C.covariance(A))
1.6
>>> print(C.covariance(B))
3.7
>>> print(A.correlation(A))
1.0
>>> print(A.correlation(B))
0.5366563145999494
>>> print(C.correlation(A))
0.8614110432930647
>>> print(C.correlation(B))
0.8908553128346921

We plainly see how the :mod:`uncertainties` package is aware of the correlation, or
lack of correlation, between :class:`UFloat` objects by looking at the following simple
example

>>> x = UFloat(5, 0.5)
>>> y = UFloat(5, 0.5)
>>> print(x - y)
0.0+/-0.7
>>> print(x - x)
0.0+/-0

We can calculate the covariance and correlation between ``x`` and ``y``

>>> print(x.covariance(y))
0.0
>>> print(x.correlation(y))
0.0
>>> print(x.covariance(x))
0.25
>>> print(x.correlation(x))
1.0

Here is one more set of examples:

>>> x = UFloat(0.2, 0.01)
>>> square = x**2
>>> print(square)
0.040+/-0.004
>>> print(square - x*x)
0.0+/-0
>>> y = x*x + 1
>>> print(y - square)
1.0+/-0

.. index:: mathematical operation; on a scalar, umath

.. _advanced math operations:

Mathematical operations with UFloat objects
===========================================

Besides being able to apply basic arithmetic operations to uncertainties
:class`UFloat` objects, this package provides generalized versions of 40 of the the
functions from the standard :mod:`math` *module*.  These mathematical functions
are found in the :mod:`uncertainties.umath` module:

.. doctest::
   :hide:

   >>> import random
   >>> random.seed(123)

>>> from uncertainties.umath import sin, exp, sqrt
>>> x = ufloat(0.2, 0.01)
>>> sin(x)
0.19866933079506122+/-0.009800665778412416
>>> sin(x*x)
0.03998933418663417+/-0.003996800426643912
>>> exp(-x/3.0)
0.9355069850316178+/-0.003118356616772059
>>> sqrt(230*x + 3)
7.0+/-0.16428571428571428

We can verify the ``sin(x)`` example follows the linear error propagation formula above.
We know the derivative of ``sin(x)`` is ``cos(x)``. So, if ``x`` only depends on a
single :class:`UAtom` as it does in this example, we expect the corresponding error
contribution to be ``cos(x0)`` times the weight of that :class:`UAtom` for ``x``, 0.01.

>>> from uncertainties.umath import cos
>>> print(x.error_components)
{UAtom(44867db30d67b366): 0.01}
>>> print(0.01 * cos(0.2))
0.009800665778412416
>>> print(sin(x).error_components)
{UAtom(44867db30d67b366): 0.009800665778412416}

We see the expected weighting.

The functions in the :mod:`uncertainties.umath` module include:

    ``acos``, ``acosh``, ``asin``, ``asinh``, ``atan``, ``atan2``, ``atanh``, ``cos``,
    ``cosh``, ``degrees``, ``erf``, ``erfc``, ``exp``, ``expm1``, ``fsum``, ``gamma``,
    ``hypot``, ``isinf``, ``isnan``, ``lgamma``, ``log``, ``log10``, ``log1p``, ``pow``,
    ``radians``, ``sin``, ``sinh``, ``sqrt``, ``tan``, ``tanh``,


Equality Comparison
===================

Two :class:`UFloat` objects are equal if their nominal values are equal as
:class:`float` objects and their :attr:`error_components` dictionaries are equal.
It is not sufficient for the two :class:`UFloat` to have equal :attr:`nominal_value`
and :attr:`std_dev` attributes.

.. doctest::
   :hide:

   >>> import random
   >>> random.seed(1)

>>> x = ufloat(5, 0.5)
>>> print(x == x)
True
>>> y = ufloat(5, 0.5)
>>> print(x == y)
False

We can see that this is because ``x`` and ``y`` depend on independent :class:`UAtom`
objects.

>>> print(x.error_components)
{UAtom(91b7584a2265b1f5): 0.5}
>>> print(y.error_components)
{UAtom(cd613e30d8f16adf): 0.5}

Note that if  a :class:`UFloat` object is ever found to have dependence on a
:class:`UAtom` object with a weight of 0 then that :class:`UAtom` is excluded from the
:attr:`error_components`.

Recall that :class:`UFloat` objects model real random variables.
It is not conventional to define an ordering on random variables.
For example, suppose ``X`` is a normally distributed random variable with mean 1 and
standard deviation 10 and ``Y`` is also a random variable with mean 2 and standard
deviation 10 then, over all possible samples, it will usually be true that ``Y>X``, but
for a large fraction of samples it will be ``X>Y``.
For this reason, the ordering operations ``<, <=, >=, >`` are not defined on the
:class:`UFloat` class.

In previous versions of the :mod:`uncertainties` package these ordering operations were
defined based on comparison of the :class:`UFloat` :attr:`nominal_value` attribute.
Now, if users want to compare :class:`UFloat` objects based on the ordering of the
:attr:`nominal_value` attribute they can do so explicitly

>>> x = UFloat(1, 10)
>>> y = UFloat(2, 10)
>>> print(y > x)
Traceback (most recent call last):
    ...
TypeError: '>' not supported between instances of 'UFloat' and 'UFloat'
>>> print(y.n > x.n)
True

.. index:: covariance matrix

Covariance and correlation matrices
===================================

To lowest order, a size ``N`` set of random variables can be described by a length ``N``
sequence of mean values together with an ``NxN`` matrix capturing the pairwise
covariance or correlation matrix between the random variables.
We've seen above that the :mod:`uncertainties` package supports calculating the
covariance and correlation between two :class:`UFloat` objects.
The :mod:`uncertainties` package also provides utility functions for calculating the
``NxN`` covariance or correlation matrix for a sequence of :class:`UFloat` objects.

Furthermore, given a length ``N`` sequence of nominal values together with a valid
``NxN`` covariance orcorrelation matrix, :mod:`uncertainties` provides functions to
construct a sequence of :class:`UFloat` objects whose statistics match those inputs.

For ``N`` random variables ``X_i`` the elements of the covariance and correlation
matrices are given by:

   Cov_{i, j} = E[(X_i - E[X_i])(X_j - E[X_j])]
   Corr_{i, j} = Cov_{i, j} / sqrt(Cov_{i, i} * Cov_{j, j})

Calculating the Covariance and Correlation Matrices
---------------------------------------------------

We calculate the covariance and correlation matrices for a sequence of :class:`UFloat`
objects:

>>> from uncertainties import covariance_matrix, correlation_matrix
>>> x = ufloat(1, 0.1)
>>> y = ufloat(10, 0.1)
>>> z = x + 2 * y
>>> cov_mat = covariance_matrix([x, y, z])
>>> corr_mat = correlation_matrix(([x, y, z]))

We can view the matrices

>>> import numpy as np
>>> np.set_printoptions(precision=3)
>>> print(cov_mat)
[[0.01 0.   0.01]
 [0.   0.01 0.02]
 [0.01 0.02 0.05]]

The diagonal elements are the variances, the squares of the standard deviations, of the
three :class:`UFloat` objects.
The two 0 off-diagonal elements indicate that ``x`` and ``y`` are uncorrelated.
The non-zero off-diagonal elements show that ``z`` is correlated to ``x`` and ``y``.
The correlation matrix is a rescaled version of the covariance matrix whose entries
range from 0, for totally uncorrelated random variables, to 1, for perfectly correlated
variables:

>>> print(corr_mat)
[[1.    0.    0.447]
 [0.    1.    0.894]
 [0.447 0.894 1.   ]]

Generating `UFloat` Objects from a Covariance or Correlation Matrix
-------------------------------------------------------------------

Above we generated a covariance or correlation matrix from a sequence of :class:`UFloat`
objects.
The :mod:`uncertainties` package expose the reverse functionality.
Given a covariance matrix and a sequence of nominal values, it is
possible to construct a sequence of :class:`UFloat` with nominal values and correlations
matching the covariance matrix passed in.

>>> from uncertainties import correlated_values, correlated_values_norm
>>> x0 = 1
>>> y0 = 10
>>> z0 = x0 + 2 * y0
>>> nominal_values = [x0, y0, z0]

With this we can generate a sequence of :class:`UFloat` objects given a covariance
matrix

>>> x2, y2, z2 = correlated_values(nominal_values, cov_mat)
>>> print(x2)
1.00+/-0.10
>>> print(y2)
10.00+/-0.10
>>> print(z2)
21.00+/-0.22
>>> cov_mat_2 = covariance_matrix([x2, y2, z2])
>>> print(np.all(np.isclose(cov_mat_2, cov_mat)))
True

We can do the same with a correlation matrix, but because the correlation matrix is
normalized, we must independently supply information about the standard deviations of
the resulting :class:`UFloat` objects.
or given a correlation matrix

>>> dx = 0.1
>>> dy = 0.1
>>> dz = 0.22
>>> x3, y3, z3 = correlated_values_norm(((x0, dx), (y0, dy), (z0, dz)), corr_mat)
>>> print(x3)
1.00+/-0.10
>>> print(y3)
10.00+/-0.10
>>> print(z3)
21.00+/-0.22
>>> corr_mat_3 = correlation_matrix([x3, y3, z3])
>>> print(np.all(np.isclose(corr_mat_3, corr_mat)))
True

 .. index::
   pair: testing (scalar); NaN

Handling NaNs and infinities
===============================

NaN values can appear in either the nominal value or uncertainty of a
Variable.  As is always the case, care must be exercised when handling NaN
values.

The standard library :func:`math.isnan` and :func:`numpy.isnan` functions will raise
``TypeError`` exceptions for :mod:`uncertainties` :class:`UFloat` objects since these
functions can only handle :class:`float` input.
The :mod:`uncertainties` package provides the :func:`umath.isnan` function which reports
if the :attr:`nominal_value` attribute of a :class:`UFloat` object is NaN or not.

>>> from uncertainties import umath
>>> print(umath.isnan(UFloat(float("nan"), float("nan"))))
True
>>> print(umath.isnan(UFloat(float("nan"), 0.1)))
True
>>> print(umath.isnan(UFloat(1.0, float("nan"))))
False
>>> print(umath.isnan(UFloat(1.0, 0.1)))
False

The :func:`umath.isinf` function detects if the :attr:`nominal_value` is infinite.

To check if the :attr:`std_dev` attribute of a :class:`UFloat` object is
NaN or infinite, you must explicitly apply the :func:`math.isnan` or :func:`math.isinf`
function to the :attr:`std_dev` attribute of the :class:`UFloat` object.

>>> import math
>>> print(math.isinf(UFloat(1, float("inf")).s))
True
>>> print(math.isnan(UFloat(1, float("nan")).s))
True

.. index:: correlations; detailed example

Power Function Behavior
=======================

.. doctest::
   :hide:

   >>> import random
   >>> random.seed(2)

The value of one :class:`UFloat` raised to the power of another can be calculated in two
ways:

>>> x = UFloat(4.5, 0.2)
>>> y = UFloat(3.4, 0.4)
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
>>> print((x**y).error_components)
{UAtom(da94e3e8ab73738f): nan}

The ``y`` derivative is real anywhere ``x**y`` is real as long as ``x>=0``.
For ``x < 0`` the ``y`` derivative is always complex valued so a NaN value is used:

>>> x = -2
>>> y = ufloat(1, 0.2)
>>> print((x**y).error_components)
{UAtom(4067c3584ee207f8): nan}

.. index::
   single: C code; wrapping
   single: Fortran code; wrapping
   single: wrapping (C, Fortran,…) functions

Wrapping Functions to Support ``UFloat`` Input
==============================================

Users may have their own functions which manipulate :class:`float` input and produce
:class:`float` output.
With the :mod:`uncertainties` package, users can wrap those functions so that they
accept :class:`UFloat` input and track the uncertainty according to the rules of linear
error propagation.
This is realized using the :func:`wrap` function.
Consider calculating a Bessel function

>>> from scipy.special import jv
>>> x = ufloat(2, 0.01)
>>> jv(0, x)
Traceback (most recent call last):
 ...
TypeError: ufunc 'jv' not supported for the input types, and the inputs could not be safely coerced to any supported types according to the casting rule ''safe''

We see that this naive approach fails because the :func:`jv` function does not support
:class:`UFloat` input.
We can remedy this using the :func:`wrap` function

>>> from uncertainties import wrap
>>> wrapped = jv
>>> wrapper = wrap(wrapped)
>>> print(wrapper(0, x))
0.224+/-0.006

The wrapped function must return exactly one :class:`float`.
The resulting wrapper function can accept either a :class:`float` or :class:`UFloat`
input for any parameter for which the wrapped function accepted a :class:`float` input.
The derivatives of the function with respect to its inputs are, by default, calculated
numerically.
However, if the user has alternative functions available to calculate the derivatives
these can be used instead of the default numerical calculation.
The :mod:`uncertainties` package needs to calculate the derivative for any parameter
into which a :class:`UFloat` object is passed as an argument.
Derivatives for positional arguments can be passed in as a tuple into the
``derivatives_args`` parameter of the :func:`wraps` function and derivativs for keyword
arguments can be passed in as a dictionary into the ``derivatives_kwargs`` parameter.
Any parameters for which a user does not supply an alternative derivative calculation
function will use the default numerical derivative calculation.

>>> def wrapped(x, a, b, c):
...     return a * x**2 + b * x + c
>>>
>>> def x_deriv(x, a, b, c):
...     return 2 * a * x + b
>>>
>>> def a_deriv(x, a, b, c):
...     return x**2
>>>
>>> def c_deriv(x, a , b, c):
...     return 0
>>>
>>>
>>> wrapper_analytic = wrap(
...     wrapped,
...     derivatives_args=(x_deriv, a_deriv),
...     derivatives_kwargs={"c": c_deriv},
... )
>>> wrapper_numerical = wrap(wrapped)
>>> x = ufloat(1, 0.1)
>>> a = ufloat(3, 0.2)
>>> b = ufloat(2, 0.5)
>>> c = ufloat(8, 2)
>>>
>>> print(wrapper_analytic(x, a, b, c) == wrapper_numerical(x, a, b, c))
True

Standard Score
==============

A utility method is provided that directly yields the
`standard score <http://en.wikipedia.org/wiki/Standard_score>`_ (number of standard
deviations) between a number and a :class:`UFloat` object

>>> x = ufloat(0.20, 0.01)
>>> print(x.std_score(0.17))
-3.0

This means that 0.17 is -3.0 standard deviations away from 0.20.

Pickling
========

:class:`UFloat` objects retain all correlations before and after pickling

>>> import pickle
>>> x = UFloat(1, 0.1)
>>> y = x**2
>>> print(x.correlation(y))
1.0
>>> x2 = pickle.loads(pickle.dumps(x))
>>> print(x2.correlation(y))
1.0
>>> print(x2 == x)
True

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
