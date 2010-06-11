Technical Guide
===============

.. !!!!! The paragraphs below are only notes: they should be adapted.


Mathematical definition of numbers with uncertainties
-----------------------------------------------------

.. index:: nominal value; definition, correlations; definition

Mathematically, numbers with uncertainties are defined in this package
as probability distributions that are described here by two numbers: a
*nominal value* and an *uncertainty*.  Thus, both variables (Variable
objects) and the result of mathematical operations (AffineScalarFunc
objects) contain these two values (respectively in their nominal_value
attribute and through their std_dev() method).

Uncertainties are simply defined here as the standard deviation of the
underlying probability distribution.

Nominal values are normally contained well inside the region of
highest probability of their underlying distribution.  The nominal
value of an expression is the expression evaluated at the nominal
values of its variables.

Good choices of nominal values for random variables (ufloat objects)
are thus their median, the location of highest probability, or their
average value: the nominal value of expressions should then be well
inside their region of highest probability.




- Probability distributions (random variables and calculation results)
are printed as:

  nominal value +/- standard deviation

but this does not imply any property on the nominal value (beyond the
fact that the nominal value is normally inside the region of high
probability density), or that the probability distribution of the
result is symmetrical (this is rarely strictly the case).


Linear error propagation theory
-------------------------------


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
(over the region of significant probability).

For instance, sin(0+/-0.01) yields a meaningful standard deviation
since it is quite linear over 0+/-0.01.  However, cos(0+/-0.01), which
is parabolic around 0, yields an approximate standard deviation of 0,
which might not be precise enough for all applications.



Comparison operators
--------------------


- Logical operations (>, ==, etc.) on numbers with uncertainties have
a pragmatic semantics, in this module: numbers with uncertainties can
be used wherever Python numbers are used, most of the time with a
result identical to the one that would be obtained with their nominal
value only.  However, since the objects defined in this module
represent probability distributions and not pure numbers, logical
operations are defined in a specific way.

The result of a logical operation ("==", ">", etc.) is defined so as
to be essentially consistent with the requirement that uncertainties
be small: the value of a logical operation is True only if the
operation yields True for all infinitesimal variations of its random
variables, except, possibly, for an infinitely small number of cases.

Example:

  "x = 3.14; y = 3.14" is such that x == y

but

  x = ufloat((3.14, 0.01))
  y = ufloat((3.14, 0.01))

is not such that x == y, since x and y are independent random
variables that almost never give the same value.  However, x == x
still holds.

The boolean value (bool(x), "if x...") of a number with uncertainty x
is the result of x != 0.


Classes
-------

- The main function provided by this module is ufloat, which creates
numbers with uncertainties (Variable objects).  Variable objects can
be used as if they were regular Python numbers.  The main attributes
and methods of Variable objects are defined in the documentation of
the Variable class.

- Mathematical expressions involving numbers with uncertainties
(Variable objects) generally return AffineScalarFunc objects, which
also print as a value with uncertainty.  Their most useful attributes
and methods are described in the documentation for AffineScalarFunc.
Note that Variable objects are also AffineScalarFunc objects: testing
whether a value carries an uncertainty handled by this module should
be done with insinstance(my_value, AffineScalarFunc).
