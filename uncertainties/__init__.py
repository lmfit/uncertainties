#!! Whenever the documentation below is updated, setup.py should be
# checked for consistency.

'''
Calculations with full error propagation for quantities with uncertainties.
Derivatives can also be calculated.

Web user guide: https://pythonhosted.org/uncertainties/.

Example of possible calculation: (0.2 +/- 0.01)**2 = 0.04 +/- 0.004.

Correlations between expressions are correctly taken into account (for
instance, with x = 0.2+/-0.01, 2*x-x-x is exactly zero, as is y-x-x
with y = 2*x).


Examples:

  import uncertainties
  from uncertainties import ufloat
  from uncertainties.umath import *  # sin(), etc.

  # Mathematical operations:
  x = ufloat(0.20, 0.01)  # x = 0.20+/-0.01
  x = ufloat_fromstr("0.20+/-0.01")  # Other representation
  x = ufloat_fromstr("0.20(1)")  # Other representation
  # Implicit uncertainty of +/-1 on the last digit:
  x = ufloat_fromstr("0.20")
  print x**2  # Square: prints "0.040+/-0.004"
  print sin(x**2)  # Prints "0.0399...+/-0.00399..."

  print x.std_score(0.17)  # Prints "-3.0": deviation of -3 sigmas

  # Access to the nominal value, and to the uncertainty:
  square = x**2  # Square
  print square  # Prints "0.040+/-0.004"
  print square.nominal_value  # Prints "0.04"
  print square.std_dev  # Prints "0.004..."

  print square.derivatives[x]  # Partial derivative: 0.4 (= 2*0.20)

  # Correlations:
  u = ufloat(1, 0.05, "u variable")  # Tag
  v = ufloat(10, 0.1, "v variable")
  sum_value = u+v

  u.std_dev = 0.1  # Standard deviations can be updated on the fly
  print sum_value - u - v  # Prints "0+/-0" (exact result)

  # List of all sources of error:
  print sum_value  # Prints "11.00+/-0.14"
  for (var, error) in sum_value.error_components().iteritems():
      print "%s: %f" % (var.tag, error)  # Individual error components

  # Covariance matrices:
  cov_matrix = uncertainties.covariance_matrix([u, v, sum_value])
  print cov_matrix  # 3x3 matrix

  # Correlated variables can be constructed from a covariance matrix, if
  # NumPy is available:
  (u2, v2, sum2) = uncertainties.correlated_values([1, 10, 11],
                                                   cov_matrix)
  print u2  # Value and uncertainty of u: correctly recovered (1.00+/-0.10)
  print uncertainties.covariance_matrix([u2, v2, sum2])  # == cov_matrix

- The main function provided by this module is ufloat, which creates
numbers with uncertainties (Variable objects).  Variable objects can
be used as if they were regular Python numbers.  The main attributes
and methods of Variable objects are defined in the documentation of
the Variable class.

- Valid operations on numbers with uncertainties include basic
mathematical functions (addition, etc.).

Most operations from the standard math module (sin, etc.) can be applied
on numbers with uncertainties by using their generalization from the
uncertainties.umath module:

  from uncertainties.umath import sin
  print sin(ufloat_fromstr("1+/-0.01"))  # 0.841+/-0.005
  print sin(1)  # umath.sin() also works on floats, exactly like math.sin()

Logical operations (>, ==, etc.) are also supported.

Basic operations on NumPy arrays or matrices of numbers with
uncertainties can be performed:

  2*numpy.array([ufloat(1, 0.01), ufloat(2, 0.1)])

More complex operations on NumPy arrays can be performed through the
dedicated uncertainties.unumpy sub-module (see its documentation).

Calculations that are performed through non-Python code (Fortran, C,
etc.) can handle numbers with uncertainties instead of floats through
the provided wrap() wrapper:

  import uncertainties

  # wrapped_f is a version of f that can take arguments with
  # uncertainties, even if f only takes floats:
  wrapped_f = uncertainties.wrap(f)

If some derivatives of the wrapped function f are known (analytically,
or numerically), they can be given to wrap()--see the documentation
for wrap().

- Utility functions are also provided: the covariance matrix between
random variables can be calculated with covariance_matrix(), or used
as input for the definition of correlated quantities (correlated_values()
function--defined only if the NumPy module is available).

- Mathematical expressions involving numbers with uncertainties
generally return AffineScalarFunc objects, which also print as a value
with uncertainty.  Their most useful attributes and methods are
described in the documentation for AffineScalarFunc.  Note that
Variable objects are also AffineScalarFunc objects.  UFloat is an
alias for AffineScalarFunc, provided as a convenience: testing whether
a value carries an uncertainty handled by this module should be done
with insinstance(my_value, UFloat).

- Mathematically, numbers with uncertainties are, in this package,
probability distributions.  These probabilities are reduced to two
numbers: a nominal value and an uncertainty.  Thus, both variables
(Variable objects) and the result of mathematical operations
(AffineScalarFunc objects) contain these two values (respectively in
their nominal_value and std_dev attributes).


The uncertainty of a number with uncertainty is simply defined in
this package as the standard deviation of the underlying probability
distribution.

The numbers with uncertainties manipulated by this package are assumed
to have a probability distribution mostly contained around their
nominal value, in an interval of about the size of their standard
deviation.  This should cover most practical cases.  A good choice of
nominal value for a number with uncertainty is thus the median of its
probability distribution, the location of highest probability, or the
average value.

- When manipulating ensembles of numbers, some of which contain
uncertainties, it can be useful to access the nominal value and
uncertainty of all numbers in a uniform manner:

  x = ufloat_fromstr("3+/-0.1")
  print nominal_value(x)  # Prints 3
  print std_dev(x)  # Prints 0.1
  print nominal_value(3)  # Prints 3: nominal_value works on floats
  print std_dev(3)  # Prints 0: std_dev works on floats

- Probability distributions (random variables and calculation results)
are printed as:

  nominal value +/- standard deviation

but this does not imply any property on the nominal value (beyond the
fact that the nominal value is normally inside the region of high
probability density), or that the probability distribution of the
result is symmetrical (this is rarely strictly the case).

- Linear approximations of functions (around the nominal values) are
used for the calculation of the standard deviation of mathematical
expressions with this package.

The calculated standard deviations and nominal values are thus
meaningful approximations as long as the functions involved have
precise linear expansions in the region where the probability
distribution of their variables is the largest.  It is therefore
important that uncertainties be small.  Mathematically, this means
that the linear term of functions around the nominal values of their
variables should be much larger than the remaining higher-order terms
over the region of significant probability.

For instance, sin(0+/-0.01) yields a meaningful standard deviation
since it is quite linear over 0+/-0.01.  However, cos(0+/-0.01) yields
an approximate standard deviation of 0 (because the cosine is not well
approximated by a line around 0), which might not be precise enough
for all applications.

- Comparison operations (>, ==, etc.) on numbers with uncertainties
have a pragmatic semantics, in this package: numbers with
uncertainties can be used wherever Python numbers are used, most of
the time with a result identical to the one that would be obtained
with their nominal value only.  However, since the objects defined in
this module represent probability distributions and not pure numbers,
comparison operator are interpreted in a specific way.

The result of a comparison operation ("==", ">", etc.) is defined so as
to be essentially consistent with the requirement that uncertainties
be small: the value of a comparison operation is True only if the
operation yields True for all infinitesimal variations of its random
variables, except, possibly, for an infinitely small number of cases.

Example:

  "x = 3.14; y = 3.14" is such that x == y

but

  x = ufloat(3.14, 0.01)
  y = ufloat(3.14, 0.01)

is not such that x == y, since x and y are independent random
variables that almost never give the same value.  However, x == x
still holds.

The boolean value (bool(x), "if x...") of a number with uncertainty x
is the result of x != 0.

- The uncertainties package is for Python 2.3 and above.

- This package contains tests.  They can be run either manually or
automatically with the nose unit testing framework (nosetests).

(c) 2009-2016 by Eric O. LEBIGOT (EOL) <eric.lebigot@normalesup.org>.
Please send feature requests, bug reports, or feedback to this address.

Please support future development by donating $10 or more through PayPal!

This software is released under a dual license.  (1) The BSD license.
(2) Any other license, as long as it is obtained from the original
author.'''

from builtins import map
from .core import *
from .core import __all__  # For a correct help(uncertainties)

# Numerical version:
__version_info__ = (3, 1, 7)
__version__ = '.'.join(map(str, __version_info__))

__author__ = 'Eric O. LEBIGOT (EOL) <eric.lebigot@normalesup.org>'
