#!! Whenever the documentation below is updated, setup.py should be
# checked for consistency.

'''
Calculations with full error propagation for quantities with uncertainties.
Derivatives can also be calculated.

Web user guide: http://packages.python.org/uncertainties/.

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
  x = ufloat_fromstr("0.20")  # Implicit uncertainty of +/-1 on the last digit
  print x**2  # Square: prints "0.04+/-0.004"
  print sin(x**2)  # Prints "0.0399...+/-0.00399..."

  print x.std_score(0.17)  # Prints "-3.0": deviation of -3 sigmas

  # Access to the nominal value, and to the uncertainty:
  square = x**2  # Square
  print square  # Prints "0.04+/-0.004"  
  print square.nominal_value  # Prints "0.04"
  print square.std_dev  # Prints "0.004..."

  print square.derivatives[x]  # Partial derivative: 0.4 (= 2*0.20)

  # Correlations:
  u = ufloat(1, 0.05, "u variable")  # Tag
  v = ufloat(10, 0.1, "v variable")
  sum_value = u+v
  
  u.std_dev = 0.1  # Standard deviations can be updated on the fly
  print sum_value - u - v  # Prints "0.0" (exact result)

  # List of all sources of error:
  print sum_value  # Prints "11+/-0.1414..."
  for (var, error) in sum_value.error_components().iteritems():
      print "%s: %f" % (var.tag, error)  # Individual error components

  # Covariance matrices:
  cov_matrix = uncertainties.covariance_matrix([u, v, sum_value])
  print cov_matrix  # 3x3 matrix

  # Correlated variables can be constructed from a covariance matrix, if
  # NumPy is available:
  (u2, v2, sum2) = uncertainties.correlated_values([1, 10, 11],
                                                   cov_matrix)
  print u2  # Value and uncertainty of u: correctly recovered (1+/-0.1)
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
  print sin(ufloat_fromstr("1+/-0.01"))  # 0.841...+/-0.005...
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
their nominal_value attribute and through their std_dev() method).

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

- The uncertainties package is for Python 2.5 and above.

- This package contains tests.  They can be run either manually or
automatically with the nose unit testing framework (nosetests).

(c) 2009-2013 by Eric O. LEBIGOT (EOL) <eric.lebigot@normalesup.org>.
Please send feature requests, bug reports, or feedback to this address.

Please support future development by donating $5 or more through PayPal!

This software is released under a dual license.  (1) The BSD license.
(2) Any other license, as long as it is obtained from the original
author.'''

# The idea behind this module is to replace the result of mathematical
# operations by a local approximation of the defining function.  For
# example, sin(0.2+/-0.01) becomes the affine function
# (AffineScalarFunc object) whose nominal value is sin(0.2) and
# whose variations are given by sin(0.2+delta) = 0.98...*delta.
# Uncertainties can then be calculated by using this local linear
# approximation of the original function.

from __future__ import division  # Many analytical derivatives depend on this

import sys
import re
import math
from math import sqrt, log  # Optimization: no attribute look-up
import copy
import warnings
import itertools
import collections
import inspect

# Numerical version:
__version_info__ = (2, 3, 6)
__version__ = '.'.join(map(str, __version_info__))

__author__ = 'Eric O. LEBIGOT (EOL) <eric.lebigot@normalesup.org>'

# Attributes that are always exported (some other attributes are
# exported only if the NumPy module is available...):
__all__ = [

    # All sub-modules and packages are not imported by default,
    # in particular because NumPy might be unavailable.

    'ufloat',  # Main function: returns a number with uncertainty

    # Uniform access to nominal values and standard deviations:
    'nominal_value',
    'std_dev',
    
    # Utility functions (more are exported if NumPy is present):
    'covariance_matrix',
    
    # Class for testing whether an object is a number with
    # uncertainty.  Not usually created by users (except through the
    # Variable subclass), but possibly manipulated by external code
    # ['derivatives()' method, etc.].
    'UFloat',

    # Wrapper for allowing non-pure-Python function to handle
    # quantitities with uncertainties:
    'wrap',

    # The documentation for wrap() indicates that numerical
    # derivatives are calculated through partial_derivative().  The
    # user might also want to change the size of the numerical
    # differentiation step.
    'partial_derivative'

    ]

###############################################################################

def set_doc(doc_string):
    """
    Decorator function that sets the docstring to the given text.

    It is useful for functions whose docstring is calculated
    (including string substitutions).
    """
    def set_doc_string(func):
        func.__doc__ = doc_string
        return func
    return set_doc_string

# Some types known to not depend on Variable objects are put in
# CONSTANT_TYPES.  The most common types can be put in front, as this
# may slightly improve the execution speed.
#
#! In Python 2.6+, numbers.Number could be used instead, here:
FLOAT_LIKE_TYPES = (float, int, long)
CONSTANT_TYPES = FLOAT_LIKE_TYPES+(complex,)

###############################################################################
# Utility for issuing deprecation warnings

def deprecation(message):
    '''
    Warns the user with the given message, by issuing a
    DeprecationWarning.
    '''

    # stacklevel = 3 points to the original user call (not to the
    # function from this module that called deprecation()).
    # DeprecationWarning is ignored by default: not used.
    
    warnings.warn('Obsolete: %s Code can be automatically updated with'
                  ' python -m uncertainties.1to2 -w ProgramDirectory.'
                  % message, stacklevel=3)

###############################################################################

def isnan(x):
    '''
    Equivalent to the math.isnan() of Python 2.6+.
    '''
    return x != x
    
###############################################################################

## Definitions that depend on the availability of NumPy:


try:
    import numpy
except ImportError:
    pass
else:

    # NumPy numbers do not depend on Variable objects:
    FLOAT_LIKE_TYPES += (numpy.number,)
    CONSTANT_TYPES += FLOAT_LIKE_TYPES[-1:]
    
    # Entering variables as a block of correlated values.  Only available
    # if NumPy is installed.

    #! It would be possible to dispense with NumPy, but a routine should be
    # written for obtaining the eigenvectors of a symmetric matrix.  See
    # for instance Numerical Recipes: (1) reduction to tri-diagonal
    # [Givens or Householder]; (2) QR / QL decomposition.
    
    def correlated_values(nom_values, covariance_mat, tags=None):
        """
        Returns numbers with uncertainties (AffineScalarFunc objects)
        that correctly reproduce the given covariance matrix, and have
        the given (float) values as their nominal value.

        The correlated_values_norm() function returns the same result,
        but takes a correlation matrix instead of a covariance matrix.
        
        The list of values and the covariance matrix must have the
        same length, and the matrix must be a square (symmetric) one.

        The numbers with uncertainties returned depend on newly
        created, independent variables (Variable objects).

        If 'tags' is not None, it must list the tag of each new
        independent variable.

        nom_values -- sequence with the nominal (real) values of the
        numbers with uncertainties to be returned.

        covariance_mat -- full covariance matrix of the returned
        numbers with uncertainties (not the statistical correlation
        matrix, i.e., not the normalized covariance matrix). For
        example, the first element of this matrix is the variance of
        the first returned number with uncertainty.
        """

        # If no tags were given, we prepare tags for the newly created
        # variables:
        if tags is None:
            tags = (None,) * len(nom_values)

        # The covariance matrix is diagonalized in order to define
        # the independent variables that model the given values:

        (variances, transform) = numpy.linalg.eigh(covariance_mat)

        # Numerical errors might make some variances negative: we set
        # them to zero:
        variances[variances < 0] = 0.
        
        # Creation of new, independent variables:

        # We use the fact that the eigenvectors in 'transform' are
        # special: 'transform' is unitary: its inverse is its transpose:

        variables = tuple(
            # The variables represent "pure" uncertainties:
            Variable(0, sqrt(variance), tag)
            for (variance, tag) in zip(variances, tags))

        # Representation of the initial correlated values:
        values_funcs = tuple(
            AffineScalarFunc(value, dict(zip(variables, coords)))
            for (coords, value) in zip(transform, nom_values))

        return values_funcs

    __all__.append('correlated_values')

    def correlated_values_norm(values_with_std_dev, correlation_mat,
                               tags=None):
        '''
        Returns correlated values like correlated_values(), but takes
        instead as input:

        - nominal (float) values along with their standard deviation, and
        
        - a correlation matrix (i.e. a normalized covariance matrix
          normalized with individual standard deviations).

        values_with_std_dev -- sequence of (nominal value, standard
        deviation) pairs. The returned, correlated values have these
        nominal values and standard deviations.

        correlation_mat -- correlation matrix (i.e. the normalized
        covariance matrix, a matrix with ones on its diagonal).
        '''

        (nominal_values, std_devs) = numpy.transpose(values_with_std_dev)

        return correlated_values(
            nominal_values,
            correlation_mat*std_devs*std_devs[numpy.newaxis].T,
            tags)
        
    __all__.append('correlated_values_norm')
    
###############################################################################

# Mathematical operations with local approximations (affine scalar
# functions)

class NotUpcast(Exception):
    'Raised when an object cannot be converted to a number with uncertainty'

def to_affine_scalar(x):
    """
    Transforms x into a constant affine scalar function
    (AffineScalarFunc), unless it is already an AffineScalarFunc (in
    which case x is returned unchanged).

    Raises an exception unless 'x' belongs to some specific classes of
    objects that are known not to depend on AffineScalarFunc objects
    (which then cannot be considered as constants).
    """

    if isinstance(x, AffineScalarFunc):
        return x

    if isinstance(x, CONSTANT_TYPES):
        # No variable => no derivative:
        return AffineScalarFunc(x, {})

    # Case of lists, etc.
    raise NotUpcast("%s cannot be converted to a number with"
                    " uncertainty" % type(x))

# !! It would be possible to split the partial derivative calculation
# into two functions: one for positional arguments (case of integer
# arg_ref) and one for keyword arguments (case of string
# arg_ref). However, this would either duplicate the code for the
# numerical differentiation, or require a call, which is probably more
# expensive in time than the tests done here.
def partial_derivative(f, arg_ref):
    """
    Returns a function that numerically calculates the partial
    derivative of function f with respect to its argument arg_ref.

    arg_ref -- describes which variable to use for the
    differentiation. If f is called with f(*args, **kwargs) arguments,
    an integer represents the index of an argument in args, and a
    string represents the name of an argument in kwargs.
    """

    # Which set of function parameter contains the variable to be
    # changed? the positional or the optional keyword arguments?
    change_kwargs = isinstance(arg_ref, basestring)
    
    def partial_derivative_of_f(*args, **kwargs):
        """
        Partial derivative, calculated with the (-epsilon, +epsilon)
        method, which is more precise than the (0, +epsilon) method.
        """

        # args_with_var contains the arguments (either args or kwargs)
        # that contain the variable that must be shifted, as a mutable
        # object (because the variable contents will be modified):

        # The values in args need to be modified, for the
        # differentiation: it is converted to a list:
        args_with_var = kwargs if change_kwargs else list(args)
       
        # The step is relative to the parameter being varied, so that
        # shifting it does not suffer from finite precision limitations:
        step = 1e-8*abs(args_with_var[arg_ref])
        if not step:
            # Arbitrary, but "small" with respect to 1, and of the
            # order of the square root of the precision of double
            # precision floats:
            step = 1e-8

        args_with_var[arg_ref] += step

        shifted_f_plus = (
            f(*args, **args_with_var) if change_kwargs
            else f(*args_with_var, **kwargs))

        args_with_var[arg_ref] -= 2*step  # Optimization: only 1 list copy
        
        shifted_f_minus = (
            f(*args, **args_with_var) if change_kwargs
            else f(*args_with_var, **kwargs))

        return (shifted_f_plus - shifted_f_minus)/2/step

    return partial_derivative_of_f

class NumericalDerivatives(object):
    """
    Convenient access to the partial derivatives of a function,
    calculated numerically.
    """
    # This is not a list because the number of arguments of the
    # function is not known in advance, in general.

    def __init__(self, function):
        """
        'function' is the function whose derivatives can be computed.
        """
        self._function = function

    def __getitem__(self, n):
        """
        Returns the n-th numerical derivative of the function.
        """
        return partial_derivative(self._function, n)

class IndexableIter(object):
    '''
    Iterable whose values can also be accessed through indexing.

    The input iterable values are cached.

    Some attributes:

    iterable -- iterable used for returning the elements one by one.
    
    returned_elements -- list with the elements directly accessible.
    through indexing. Additional elements are obtained from self.iterable.

    none_converter -- function that takes an index and returns the
    value to be returned when None is obtained form the iterable
    (instead of None).
    '''

    def __init__(self, iterable, none_converter=lambda index: None):
        '''
        iterable -- iterable whose values will be returned.
        
        none_converter -- function applied to None returned
        values. The value that replaces None is none_converter(index),
        where index is the index of the element.
        '''
        self.iterable = iterable
        self.returned_elements = []
        self.none_converter = none_converter
        
    def __getitem__(self, index):

        returned_elements = self.returned_elements
        
        try:
            
            return returned_elements[index]
        
        except IndexError:  # Element not yet cached
            
            for pos in range(len(returned_elements), index+1):

                # ! Python 2.6+: next(...)
                value = self.iterable.next()

                if value is None:
                    value = self.none_converter(pos)
                    
                returned_elements.append(value)
            
            return returned_elements[index]

    def __str__(self):
        return '<%s: [%s...]>' % (
            self.__class__.__name__,
            ', '.join(map(str, self.returned_elements)))
    
def wrap(f, derivatives_args=[], derivatives_kwargs={}):
    """
    Wraps a function f into a function that also accepts numbers with
    uncertainties (UFloat objects); the wrapped function returns the
    value of f with the correct uncertainty and correlations. The
    wrapped function is intended to be used as a drop-in replacement
    for the original function: they can be called in the exact same
    way, the only difference being that numbers with uncertainties can
    be given to the wrapped function where f accepts float arguments.

    Doing so may be necessary when function f cannot be expressed
    analytically (with uncertainties-compatible operators and
    functions like +, *, umath.sin(), etc.).

    f must return a float-like (i.e. a float, an int, etc., not a
    list, etc.), unless when called with no number with
    uncertainty. This is because the wrapped function generally
    returns numbers with uncertainties: they represent a probability
    distribution of real numbers.

    If the wrapped function is called with no argument that has an
    uncertainty, the value of f is returned.

    Parameters: the derivatives_* parameters can be used for defining
    some of the partial derivatives of f. All the (non-None)
    derivatives must have the same signature as f.

    derivatives_args --
    
        Iterable that, when iterated over, returns either derivatives
        (functions) or None. derivatives_args can in particular be a
        simple sequence (list or tuple) that gives the derivatives of
        the first positional parameters of f.

        Each function must be the partial derivative of f with respect
        to the corresponding positional parameters.  These functions
        take the same arguments as f.

        The positional parameters of a function are usually
        positional-or-keyword parameters like in the call func(a,
        b=None). However, they also include var-positional parameters
        given through the func(a, b, *args) *args syntax. In the last
        example, derivatives_args can be an iterable that returns the
        derivative with respect to a, b and then to each optional
        argument in args.

        A value of None (instead of a function) obtained when
        iterating over derivatives_args is automatically replaced by
        the relevant numerical derivative. This derivative is not used
        if the corresponding argument is not a number with
        uncertainty. A None value can therefore be used for non-scalar
        arguments of f (like string arguments).

        If the derivatives_args iterable yields fewer derivatives than
        needed, wrap() automatically sets the remaining unspecified
        derivatives to None (i.e. to the automatic numerical
        calculation of derivatives).

        An indefinite number of derivatives can be specified by having
        derivatives_args be an infinite iterator; this can for
        instance be used for specifying the derivatives of functions
        with a undefined number of argument (like sum(), whose partial
        derivatives all return 1).

    derivatives_kwargs --

        Dictionary that maps keyword parameters to their derivatives,
        or None (as in derivatives_args).

        Keyword parameters are defined as being those of kwargs when f
        has a signature of the form f(..., **kwargs). In Python 3,
        these keyword parameters also include keyword-only parameters.
        
        Non-mapped keyword parameters are replaced automatically by
        None: the wrapped function will use, if necessary, numerical
        differentiation for these parameters (as with
        derivatives_args).

        Note that this dictionary only maps keyword *parameters* from
        the *signature* of f. The way the wrapped function is called
        is immaterial: for example, if f has signature f(a, b=None),
        then derivatives_kwargs should be the empty dictionary, even
        if the wrapped f can be called a wrapped_f(a=123, b=42).

    Example (for illustration purposes only, as
    uncertainties.umath.sin() runs faster than the examples that
    follow): wrap(math.sin) is a sine function that can be applied to
    numbers with uncertainties.  Its derivative will be calculated
    numerically.  wrap(math.sin, [None]) would have produced the same
    result.  wrap(math.sin, [math.cos]) is the same function, but with
    an analytically defined derivative.
        
    Numerically calculated derivatives are meaningless when the
    function is not differentiable (e.g., math.hypot(x, y) in (x, y) =
    (0, 0), and sqrt(x) in x = 0). The corresponding uncertainties are
    either meaningless (case of hypot) or raise an exception when
    calculated (case of sqrt). In such cases, it is recommended (but
    not mandatory) to supply instead a derivative function that
    returns NaN where the function is not differentiable. This
    function can still numerically calculate the derivative where
    defined, for instance by using the partial_derivative() function.
        
    The correctness of the supplied analytical derivatives an be
    tested by setting them to None instead and comparing the
    analytical and the numerical differentiation results.
    
    Note on efficiency: the wrapped function assumes that f cannot
    accept numbers with uncertainties as arguments. If f actually does
    handle some arguments even when they have an uncertainty, the
    wrapped function ignores this fact, which might lead to a
    performance hit: wrapping a function that actually accepts numbers
    with uncertainty is likely to make it slower.
    """

    derivatives_args_index = IndexableIter(
        # Automatic addition of numerical derivatives in case the
        # supplied derivatives_args is shorter than the number of
        # arguments in *args:
        itertools.chain(derivatives_args, itertools.repeat(None)))
        

    # Derivatives for keyword arguments (includes var-keyword
    # parameters **kwargs, but also var-or-keyword parameters, and
    # keyword-only parameters (Python 3):
    
    derivatives_all_kwargs = dict(
        
        # Python 2.7+: {name: ...}
        (name,
         # Optimization: None keyword-argument derivatives are converted
         # right away to derivatives (instead of doing this every time a
         # None derivative is encountered when calculating derivatives):
         partial_derivative(f, name) if derivative is None
         else derivative)

        for (name, derivative) in derivatives_kwargs.iteritems()
    )

    # When the wrapped function is called with keyword arguments that
    # map to positional-or-keyword parameters, their derivative is
    # looked for in derivatives_all_kwargs.  We define these
    # additional derivatives:

    try:
        argspec = inspect.getargspec(f)
    except TypeError:
        # Some functions do not provide meta-data about their
        # arguments (see PEP 362). One cannot use keyword arguments
        # for positional-or-keyword parameters with them: nothing has
        # to be done:
        pass
    else:
        # With Python 3, there is no need to handle keyword-only
        # arguments (and therefore to use inspect.getfullargspec())
        # because they are already handled by derivatives_kwargs.

        # ! Python 2.6+: argspec.args
        for (index, name) in enumerate(argspec[0]):
            
            # The following test handles the case of
            # positional-or-keyword parameter for which automatic
            # numerical differentiation is used: when the wrapped
            # function is called with a keyword argument for this
            # parameter, the numerical derivative must be calculated
            # with respect to the parameter name. In the other case,
            # where the wrapped function is called with a positional
            # argument, the derivative with respect to its index must
            # be used:
            
            derivative = derivatives_args_index[index]

            derivatives_all_kwargs[name] = (
                partial_derivative(f, name)
                if derivative is None else derivative)

    # Optimization: None derivatives for the positional arguments are
    # converted to the corresponding numerical differentiation
    # function (instead of doing this over and over later every time a
    # None derivative is found):

    none_converter = lambda index: partial_derivative(f, index)
    
    derivatives_args_index.returned_elements = [
        none_converter(index) if derivative is None
        else derivative
        for (index, derivative) in enumerate(
            derivatives_args_index.returned_elements)]

    # Future None values are also automatically converted:
    derivatives_args_index.none_converter = none_converter

    
    ## Wrapped function:

    #! Setting the doc string after "def f_with...()" does not
    # seem to work.  We define it explicitly:
    @set_doc("""\
    Version of %s(...) that returns an affine approximation
    (AffineScalarFunc object), if its result depends on variables
    (Variable objects).  Otherwise, returns a simple constant (when
    applied to constant arguments).
    
    Warning: arguments of the function that are not AffineScalarFunc
    objects must not depend on uncertainties.Variable objects in any
    way.  Otherwise, the dependence of the result in
    uncertainties.Variable objects will be incorrect.
    
    Original documentation:
    %s""" % (f.__name__, f.__doc__))    
    def f_with_affine_output(*args, **kwargs):

        ########################################
        # The involved random variables must first be gathered, so
        # that they can be independently updated.
        
        # The arguments that contain an uncertainty (AffineScalarFunc
        # objects) are gathered, as positions or names; they will be
        # replaced by simple floats.

        pos_w_uncert = [index for (index, value) in enumerate(args)
                        if isinstance(value, AffineScalarFunc)]
        names_w_uncert = set(key for (key, value) in kwargs.iteritems()
                             if isinstance(value, AffineScalarFunc))

        ########################################
        # Value of f() at the nominal value of the arguments with
        # uncertainty:

        # The usual behavior of f() is kept, if no number with
        # uncertainty is provided:
        if (not pos_w_uncert) and (not names_w_uncert):
            return f(*args, **kwargs)
                    
        ### Nominal values of the (scalar) arguments:

        # !! Possible optimization: If pos_w_uncert is empty, there
        # is actually no need to create a mutable version of args and
        # one could do args_values = args.  However, the wrapped
        # function is typically called with numbers with uncertainties
        # as positional arguments (i.e., pos_w_uncert is not emtpy),
        # so this "optimization" is not implemented here.
        
        ## Positional arguments:
        args_values = list(args)  # Now mutable: modified below
        # Arguments with an uncertainty are converted to their nominal
        # value:
        for index in pos_w_uncert:
            args_values[index] = args[index].nominal_value

        ## Keyword arguments:

        # For efficiency reasons, kwargs is not copied. Instead, its
        # values with uncertainty are modified:

        # The original values with uncertainties are needed: they are
        # saved in the following dictionary:

        kwargs_uncert_values = {}

        for name in names_w_uncert:
            value_with_uncert = kwargs[name]
            # Saving for future use:
            kwargs_uncert_values[name] = value_with_uncert
            # The original dictionary is modified (for efficiency reasons):
            kwargs[name] = value_with_uncert.nominal_value
            
        f_nominal_value = f(*args_values, **kwargs)

        # If the value is not a float, then this code cannot provide
        # the result, as it returns a UFloat, which represents a
        # random real variable. This happens for instance when
        # ufloat()*numpy.array() is calculated: the
        # AffineScalarFunc.__mul__ operator, obtained through wrap(),
        # returns a NumPy array, not a float:
        if not isinstance(f_nominal_value, FLOAT_LIKE_TYPES):
            return NotImplemented
        
        ########################################

        # The chain rule will be applied.  In the case of numerical
        # derivatives, this method gives a better-controlled numerical
        # stability than numerically calculating the partial
        # derivatives through '[f(x + dx, y + dy, ...) -
        # f(x,y,...)]/da' where dx, dy,... are calculated by varying
        # 'a' by 'da'.  In fact, this allows the program to control
        # how big the dx, dy, etc. are, which is numerically more
        # precise.
        
        ########################################
            
        # Calculation of the derivatives with respect to the variables
        # of f that have a number with uncertainty.
        
        # Mappings of each relevant argument reference (numerical
        # index in args, or name in kwargs to the value of the
        # corresponding partial derivative of f (only for those
        # arguments that contain a number with uncertainty).
        
        derivatives_num_args = {}

        for pos in pos_w_uncert:
            derivatives_num_args[pos] = derivatives_args_index[pos](
                *args_values, **kwargs)

        derivatives_num_kwargs = {}
        for name in names_w_uncert:

            # Optimization: caching of the automatic numerical
            # derivatives for keyword arguments that are
            # discovered. This gives a speedup when the original
            # function is called repeatedly with the same keyword
            # arguments:
            derivative = derivatives_all_kwargs.setdefault(
                name,
                # Derivative never needed before:
                partial_derivative(f, name))

            derivatives_num_kwargs[name] = derivative(
                *args_values, **kwargs)

        ########################################
        # Calculation of the derivative of f with respect to all the
        # variables (Variable objects) involved.

        # Involved variables (Variable objects):
        variables = set()
        
        for expr in itertools.chain(
            (args[index] for index in pos_w_uncert),  # From args
            kwargs_uncert_values.itervalues()):  # From kwargs
            
            # !! In Python 2.7+: |= expr.derivatives.viewkeys()
            variables |= set(expr.derivatives)
        
        # Initial value for the chain rule (is updated below):
        derivatives_wrt_vars = dict((var, 0.) for var in variables)

        ## The chain rule is used...

        # ... on args:
        for (pos, f_derivative) in derivatives_num_args.iteritems():
            for (var, arg_derivative) in args[pos].derivatives.iteritems():
                derivatives_wrt_vars[var] += f_derivative * arg_derivative
        # ... on kwargs:
        for (name, f_derivative) in derivatives_num_kwargs.iteritems():
            for (var, arg_derivative) in (kwargs_uncert_values[name]
                                          .derivatives.iteritems()):
                derivatives_wrt_vars[var] += f_derivative * arg_derivative

                
        # The function now returns an AffineScalarFunc object:        
        return AffineScalarFunc(f_nominal_value, derivatives_wrt_vars)

    # It is easier to work with f_with_affine_output, which represents
    # a wrapped version of 'f', when it bears the same name as
    # 'f'.
    f_with_affine_output.__name__ = f.__name__
    # !! Note: Setting __name__ is however not fully sufficient: the
    # name f_with_affine_output is stored in
    # f_with_affine_output.__code__.co_name; being able to change it
    # would be useful for instance when f_with_affine_output() is
    # called with unexpected arguments (unexpected keyword argument,
    # etc.). co_name is read-only, though.

    return f_with_affine_output

def _force_aff_func_args(func):
    """
    Takes an operator op(x, y) and wraps it.

    The constructed operator returns func(x, to_affine_scalar(y)) if y
    can be upcast with to_affine_scalar(); otherwise, it returns
    NotImplemented.

    Thus, func() is only called on two AffineScalarFunc objects, if
    its first argument is an AffineScalarFunc.
    """

    def op_on_upcast_args(x, y):
        """
        Returns %s(self, to_affine_scalar(y)) if y can be upcast
        through to_affine_scalar.  Otherwise returns NotImplemented.
        """ % func.__name__
        
        try:
            y_with_uncert = to_affine_scalar(y)
        except NotUpcast:
            # This module does not know how to handle the comparison:
            # (example: y is a NumPy array, in which case the NumPy
            # array will decide that func() should be applied
            # element-wise between x and all the elements of y):
            return NotImplemented
        else:
            return func(x, y_with_uncert)

    return op_on_upcast_args

########################################

# Definition of boolean operators, that assume that self and
# y_with_uncert are AffineScalarFunc.

# The fact that uncertainties must be small is used, here: the 
# comparison functions are supposed to be constant for most values of 
# the random variables.

# Even though uncertainties are supposed to be small, comparisons 
# between 3+/-0.1 and 3.0 are handled correctly (even though x == 3.0 is 
# not a constant function in the 3+/-0.1 interval).  The comparison 
# between x and x is handled too, when x has an uncertainty.  In fact, 
# as explained in the main documentation, it is possible to give a 
# useful meaning to the comparison operators, in these cases.

def _eq_on_aff_funcs(self, y_with_uncert):
    """
    __eq__ operator, assuming that both self and y_with_uncert are
    AffineScalarFunc objects.
    """
    difference = self - y_with_uncert
    # Only an exact zero difference means that self and y are
    # equal numerically:
    return not(difference._nominal_value or difference.std_dev)

def _ne_on_aff_funcs(self, y_with_uncert):
    """
    __ne__ operator, assuming that both self and y_with_uncert are
    AffineScalarFunc objects.
    """

    return not _eq_on_aff_funcs(self, y_with_uncert)

def _gt_on_aff_funcs(self, y_with_uncert):
    """
    __gt__ operator, assuming that both self and y_with_uncert are
    AffineScalarFunc objects.
    """
    return self._nominal_value > y_with_uncert._nominal_value

def _ge_on_aff_funcs(self, y_with_uncert):
    """
    __ge__ operator, assuming that both self and y_with_uncert are
    AffineScalarFunc objects.
    """

    return (_gt_on_aff_funcs(self, y_with_uncert)
            or _eq_on_aff_funcs(self, y_with_uncert))

def _lt_on_aff_funcs(self, y_with_uncert):
    """
    __lt__ operator, assuming that both self and y_with_uncert are
    AffineScalarFunc objects.
    """
    return self._nominal_value < y_with_uncert._nominal_value

def _le_on_aff_funcs(self, y_with_uncert):
    """
    __le__ operator, assuming that both self and y_with_uncert are
    AffineScalarFunc objects.
    """

    return (_lt_on_aff_funcs(self, y_with_uncert)
            or _eq_on_aff_funcs(self, y_with_uncert))

########################################

class CallableStdDev(float):
    '''
    Class for standard deviation results, which used to be
    callable. Provided for compatibility with old code. Issues an
    obsolescence warning upon call.
    '''
    
    def __call__ (self):
        deprecation('the std_dev attribute should not be called'
                    ' anymore: use .std_dev instead of .std_dev().')
        return self
        
class AffineScalarFunc(object):
    """
    Affine functions that support basic mathematical operations
    (addition, etc.).  Such functions can for instance be used for
    representing the local (linear) behavior of any function.

    This class is mostly meant to be used internally.

    This class can also be used to represent constants.

    The variables of affine scalar functions are Variable objects.

    AffineScalarFunc objects include facilities for calculating the
    'error' on the function, from the uncertainties on its variables.

    Main attributes and methods:
    
    - nominal_value, std_dev: value at the origin / nominal value, and
      standard deviation.  The standard deviation can be NaN.

    - error_components(): error_components()[x] is the error due to
      Variable x.

    - derivatives: derivatives[x] is the (value of the) derivative
      with respect to Variable x.  This attribute is a dictionary
      whose keys are the Variable objects on which the function
      depends.
      
      All the Variable objects on which the function depends are in
      'derivatives'.

    - std_score(x): position of number x with respect to the
      nominal value, in units of the standard deviation.
    """

    # To save memory in large arrays:
    __slots__ = ('_nominal_value', 'derivatives')
    
    #! The code could be modify in order to accommodate for non-float
    # nominal values.  This could for instance be done through
    # the operator module: instead of delegating operations to
    # float.__*__ operations, they could be delegated to
    # operator.__*__ functions (while taking care of properly handling
    # reverse operations: __radd__, etc.).

    def __init__(self, nominal_value, derivatives):
        """
        nominal_value -- value of the function at the origin.
        nominal_value must not depend in any way of the Variable
        objects in 'derivatives' (the value at the origin of the
        function being defined is a constant).

        derivatives -- maps each Variable object on which the function
        being defined depends to the value of the derivative with
        respect to that variable, taken at the nominal value of all
        variables.
 
        Warning: the above constraint is not checked, and the user is
        responsible for complying with it.
        """

        # Defines the value at the origin:

        # Only float-like values are handled.  One reason is that it
        # does not make sense for a scalar function to be affine to
        # not yield float values.  Another reason is that it would not
        # make sense to have a complex nominal value, here (it would
        # not be handled correctly at all): converting to float should
        # be possible.

        self._nominal_value = float(nominal_value)
        self.derivatives = derivatives

    # The following prevents the 'nominal_value' attribute from being
    # modified by the user:
    @property
    def nominal_value(self):
        "Nominal value of the random number."
        return self._nominal_value
    
    ############################################################

        
    ### Operators: operators applied to AffineScalarFunc and/or
    ### float-like objects only are supported.  This is why methods
    ### from float are used for implementing these operators.

    # Operators with no reflection:

    ########################################
        
    # __nonzero__() is supposed to return a boolean value (it is used
    # by bool()).  It is for instance used for converting the result
    # of comparison operators to a boolean, in sorted().  If we want
    # to be able to sort AffineScalarFunc objects, __nonzero__ cannot
    # return a AffineScalarFunc object.  Since boolean results (such
    # as the result of bool()) don't have a very meaningful
    # uncertainty unless it is zero, this behavior is fine.
    
    def __nonzero__(self):
        """
        Equivalent to self != 0.
        """
        #! This might not be relevant for AffineScalarFunc objects
        # that contain values in a linear space which does not convert
        # the float 0 into the null vector (see the __eq__ function:
        # __nonzero__ works fine if subtracting the 0 float from a
        # vector of the linear space works as if 0 were the null
        # vector of that space):
        return self != 0.  # Uses the AffineScalarFunc.__ne__ function

    ########################################
    
    ## Logical operators: warning: the resulting value cannot always
    ## be differentiated.

    # The boolean operations are not differentiable everywhere, but
    # almost...

    # (1) I can rely on the assumption that the user only has "small"
    # errors on variables, as this is used in the calculation of the
    # standard deviation (which performs linear approximations):

    # (2) However, this assumption is not relevant for some
    # operations, and does not have to hold, in some cases.  This
    # comes from the fact that logical operations (e.g. __eq__(x,y))
    # are not differentiable for many usual cases.  For instance, it
    # is desirable to have x == x for x = n+/-e, whatever the size of e.
    # Furthermore, n+/-e != n+/-e', if e != e', whatever the size of e or
    # e'.

    # (3) The result of logical operators does not have to be a
    # function with derivatives, as these derivatives are either 0 or
    # don't exist (i.e., the user should probably not rely on
    # derivatives for his code).
 
    # !! In Python 2.7+, it may be possible to use functools.total_ordering.
   
    # __eq__ is used in "if data in [None, ()]", for instance.  It is
    # therefore important to be able to handle this case too, which is
    # taken care of when _force_aff_func_args(_eq_on_aff_funcs)
    # returns NotImplemented.
    __eq__ = _force_aff_func_args(_eq_on_aff_funcs)
    
    __ne__ = _force_aff_func_args(_ne_on_aff_funcs)
    __gt__ = _force_aff_func_args(_gt_on_aff_funcs)

    # __ge__ is not the opposite of __lt__ because these operators do
    # not always yield a boolean (for instance, 0 <= numpy.arange(10)
    # yields an array).
    __ge__ = _force_aff_func_args(_ge_on_aff_funcs)

    __lt__ = _force_aff_func_args(_lt_on_aff_funcs)
    __le__ = _force_aff_func_args(_le_on_aff_funcs)

    ########################################

    # Uncertainties handling:
    
    def error_components(self):
        """
        Individual components of the standard deviation of the affine
        function (in absolute value), returned as a dictionary with
        Variable objects as keys. The returned variables are the
        independent variables that the affine function depends on.

        This method assumes that the derivatives contained in the
        object take scalar values (and are not a tuple, like what
        math.frexp() returns, for instance).
        """
    
        # Calculation of the variance:
        error_components = {}

        for (variable, derivative) in self.derivatives.iteritems():            
            # Individual standard error due to variable:
            error_components[variable] = (
                0.
                # 0 is returned even for a NaN derivative, since an
                # exact number has a 0 uncertainty:
                if variable._std_dev == 0
                else abs(derivative*variable._std_dev))
            
        return error_components
    
    @property
    def std_dev(self):
        """
        Standard deviation of the affine function.

        This method assumes that the function returns scalar results.

        This returned standard deviation depends on the current
        standard deviations [std_dev] of the variables (Variable
        objects) involved.
        """
        #! It would be possible to not allow the user to update the
        #std dev of Variable objects, in which case AffineScalarFunc
        #objects could have a pre-calculated or, better, cached
        #std_dev value (in fact, many intermediate AffineScalarFunc do
        #not need to have their std_dev calculated: only the final
        #AffineScalarFunc returned to the user does).
        return CallableStdDev(sqrt(sum(
            delta**2 for delta in self.error_components().itervalues())))

    def _general_representation(self, to_string):
        """
        Uses the to_string() conversion function on both the nominal
        value and the standard deviation, and returns a string that
        describes them.

        to_string() is typically repr() or str().
        """

        (nominal_value, std_dev) = (self._nominal_value, self.std_dev)

        # String representation:

        # Not putting spaces around "+/-" helps with arrays of
        # Variable, as each value with an uncertainty is a
        # block of signs (otherwise, the standard deviation can be
        # mistaken for another element of the array).

        return ("%s+/-%s" % (to_string(nominal_value), to_string(std_dev))
                if std_dev
                else to_string(nominal_value))

    def __repr__(self):
        return self._general_representation(repr)
                    
    def __str__(self):
        return self._general_representation(str)
    
    def std_score(self, value):
        """
        Returns 'value' - nominal value, in units of the standard
        deviation.

        Raises a ValueError exception if the standard deviation is zero.
        """
        try:
            # The ._nominal_value is a float: there is no integer division,
            # here:
            return (value - self._nominal_value) / self.std_dev
        except ZeroDivisionError:
            raise ValueError("The standard deviation is zero:"
                             " undefined result.")

    def __deepcopy__(self, memo):
        """
        Hook for the standard copy module.

        The returned AffineScalarFunc is a completely fresh copy,
        which is fully independent of any variable defined so far.
        New variables are specially created for the returned
        AffineScalarFunc object.
        """
        return AffineScalarFunc(
            self._nominal_value,
            dict((copy.deepcopy(var), deriv)
                 for (var, deriv) in self.derivatives.iteritems()))

    def __getstate__(self):
        """
        Hook for the pickle module.

        The slot attributes of the parent classes are returned, as
        well as those of the __dict__ attribute of the object (if
        any).
        """

        # In general (case where this class is subclassed), data
        # attribute are stored in two places: possibly in __dict_, and
        # in slots. Data from both locations is returned by this
        # method.
        
        all_attrs = {}

        # Support for subclasses that do not use __slots__ (except
        # through inheritance): instances have a __dict__
        # attribute. The keys in this __dict__ are shadowed by the
        # slot attribute names (reference:
        # http://stackoverflow.com/questions/15139067/attribute-access-in-python-first-slots-then-dict/15139208#15139208).
        # The method below not only preserves this behavior, but also
        # saves the full contents of __dict__. This is robust:
        # unpickling gives back the original __dict__ even if __dict__
        # contains keys that are shadowed by slot names:

        try:
            all_attrs['__dict__'] = self.__dict__
        except AttributeError:
            pass
        
        # All the slot attributes are gathered.

        # Classes that do not define __slots__ have the __slots__ of
        # one of their parents (the first parent with their own
        # __slots__ in MRO). This is why the slot names are first
        # gathered (with repetitions removed, in general), and their
        # values obtained later.
        
        all_slots = set()
        
        for cls in type(self).mro():

            # In the diamond inheritance pattern, some parent classes
            # may not have __slots__:
            slot_names = getattr(cls, '__slots__', ())

            # Slot names can be given in various forms (string,
            # sequence, iterable):
            if isinstance(slot_names, str):
                all_slots.add(slot_names)  # Single name
            else:
                all_slots.update(slot_names)

        # The slot values are stored:
        for name in all_slots:
            try:
                # !! It might happen that '__dict__' is itself a slot
                # name. In this case, its value is saved
                # again. Alternatively, the loop could be done on
                # all_slots - set(('__dict__',)):
                all_attrs[name] = getattr(self, name)
            except AttributeError:
                pass  # Undefined slot attribute
                
        return all_attrs

    def __setstate__(self, data_dict):
        """
        Hook for the pickle module.
        """        
        for (name, value) in data_dict.iteritems():
            setattr(self, name, value)

# Nicer name, for users: isinstance(ufloat(...), UFloat) is
# True. Also: isinstance(..., UFloat) is the test for "is this a
# number with uncertainties from the uncertainties package?".
UFloat = AffineScalarFunc

###############################################################################

# Some operators can have undefined derivatives but still give
# meaningful values when some of their arguments have a zero
# uncertainty. Such operators return NaN when their derivative is
# not finite. This way, if the uncertainty of the associated
# variable is not 0, a NaN uncertainty is produced, which
# indicates an error; if the uncertainty is 0, then the total
# uncertainty can be returned as 0.

# Exception catching is used so as to not slow down regular
# operation too much:

def nan_if_exception(f):
    '''
    Wrapper around f(x, y) that let f return NaN when f raises one of
    a few numerical exceptions.
    '''

    def wrapped_f(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except (ValueError, ZeroDivisionError, OverflowError):
            return float('nan')
        
    return wrapped_f

def get_ops_with_reflection():

    """
    Returns operators with a reflection, along with their derivatives
    (for float operands).
    """
    
    # Operators with a reflection:

    # We do not include divmod().  This operator could be included, by
    # allowing its result (a tuple) to be differentiated, in
    # derivative_value().  However, a similar result can be achieved
    # by the user by calculating separately the division and the
    # result.

    # {operator(x, y): (derivative wrt x, derivative wrt y)}:

    # Note that unknown partial derivatives can be numerically
    # calculated by expressing them as something like
    # "partial_derivative(float.__...__, 1)(x, y)":

    # String expressions are used, so that reversed operators are easy
    # to code, and execute relatively efficiently:
    
    derivatives_list = {
        'add': ("1.", "1."),
        # 'div' is the '/' operator when __future__.division is not in
        # effect.  Since '/' is applied to
        # AffineScalarFunc._nominal_value numbers, it is applied on
        # floats, and is therefore the "usual" mathematical division.
        'div': ("1/y", "-x/y**2"),
        'floordiv': ("0.", "0."),  # Non exact: there is a discontinuity
        # The derivative wrt the 2nd arguments is something like (..., x//y),
        # but it is calculated numerically, for convenience:
        'mod': ("1.", "partial_derivative(float.__mod__, 1)(x, y)"),
        'mul': ("y", "x"),
        # The case x**y is constant one the line x = 0 and in y = 0;
        # the corresponding derivatives must be zero in these
        # cases. If the function is actually not defined (e.g. 0**-3),
        # then an exception will be raised when the nominal value is
        # calculated.  These derivatives are transformed to NaN if an
        # error happens during their calculation:
        'pow': ("0. if y == 0"
                " else y*x**(y-1) if x != 0 or y % 1 == 0"
                " else float('nan')",
                "0. if (x == 0) and (y > 0) else log(x)*x**y"),
        'sub': ("1.", "-1."),
        'truediv': ("1/y", "-x/y**2")
        }

    # Conversion to Python functions:
    ops_with_reflection = {}
    for (op, derivatives) in derivatives_list.iteritems():
        ops_with_reflection[op] = [
            eval("lambda x, y: %s" % expr) for expr in derivatives ]

        ops_with_reflection["r"+op] = [
            eval("lambda y, x: %s" % expr) for expr in reversed(derivatives)]

    # Undefined derivatives are converted to NaN when the function
    # itself can be calculated:
    for op in ['pow']:
        ops_with_reflection[op] = map(nan_if_exception,
                                      ops_with_reflection[op])
        
        ops_with_reflection['r'+op] = map(nan_if_exception,
                                          ops_with_reflection['r'+op])
        
    return ops_with_reflection

# Operators that have a reflection, along with their derivatives:
_ops_with_reflection = get_ops_with_reflection()

# Some effectively modified operators (for the automated tests):
_modified_operators = []
_modified_ops_with_reflection = []

# Custom versions of some operators (instead of extending some float
# __*__ operators to AffineScalarFunc, the operators in _custom_ops
# are used):
if sys.version_info < (3,):

    _custom_ops = {}

else:

    def _no_complex_result(func):
        '''
        Returns a function that does like func, but that raises a
        ValueError if the result is complex.
        '''
        def no_complex_func(*args, **kwargs):
            '''
            Like %s, but raises a ValueError exception if the result
            is complex.
            ''' % func.__name__
            
            value = func(*args, **kwargs)
            if isinstance(value, complex):
                raise ValueError('The uncertainties module does not handle'
                                 ' complex results')
            else:
                return value

        return no_complex_func

    # This module does not handle uncertainties on complex numbers:
    # complex results for the nominal value of some operations cannot
    # be calculated with an uncertainty:
    _custom_ops = {
        'pow': _no_complex_result(float.__pow__),
        'rpow': _no_complex_result(float.__rpow__)
        }

def add_operators_to_AffineScalarFunc():
    """
    Adds many operators (__add__, etc.) to the AffineScalarFunc class.
    """
    
    ########################################

    #! Derivatives are set to return floats.  For one thing,
    # uncertainties generally involve floats, as they are based on
    # small variations of the parameters.  It is also better to
    # protect the user from unexpected integer result that behave
    # badly with the division.

    ## Operators that return a numerical value:

    # Single-argument operators that should be adapted from floats to
    # AffineScalarFunc objects, associated to their derivative:
    simple_numerical_operators_derivatives = {
        'abs': lambda x: 1. if x>=0 else -1.,
        'neg': lambda x: -1.,
        'pos': lambda x: 1.,
        'trunc': lambda x: 0.
        }

    for (op, derivative) in (
        simple_numerical_operators_derivatives.iteritems()):

        attribute_name = "__%s__" % op
        
        # float objects don't exactly have the same attributes between
        # different versions of Python (for instance, __trunc__ was
        # introduced with Python 2.6):
        try:
            setattr(AffineScalarFunc, attribute_name,
                    wrap(getattr(float, attribute_name), [derivative]))
        except AttributeError:
            pass
        else:
            _modified_operators.append(op)
            
    ########################################
    # Final definition of the operators for AffineScalarFunc objects:
            
    # Reversed versions (useful for float*AffineScalarFunc, for instance):
    for (op, derivatives) in _ops_with_reflection.iteritems():
        attribute_name = '__%s__' % op

        # float objects don't exactly have the same attributes between
        # different versions of Python (for instance, __div__ and
        # __rdiv__ were removed, in Python 3):

        # float objects don't exactly have the same attributes between
        # different versions of Python (for instance, __trunc__ was
        # introduced with Python 2.6):
        try:

            func_to_wrap = (getattr(float, attribute_name)
                            if op not in _custom_ops
                            else _custom_ops[op])

        except AttributeError:
            pass
        else:
            setattr(AffineScalarFunc, attribute_name,
                    wrap(func_to_wrap, derivatives))
            _modified_ops_with_reflection.append(op)            

    ########################################
    # Conversions to pure numbers are meaningless.  Note that the
    # behavior of float(1j) is similar.
    for coercion_type in ('complex', 'int', 'long', 'float'):
        def raise_error(self):
            raise TypeError("can't convert an affine function (%s)"
                            ' to %s; use x.nominal_value'
                            # In case AffineScalarFunc is sub-classed:
                            % (self.__class__, coercion_type))

        setattr(AffineScalarFunc, '__%s__' % coercion_type, raise_error)

add_operators_to_AffineScalarFunc()  # Actual addition of class attributes

class Variable(AffineScalarFunc):    
    """
    Representation of a float-like scalar random variable, along with
    its uncertainty.

    Objects are meant to represent variables that are independent from
    each other (correlations are handled through the AffineScalarFunc
    class).
    """

    # To save memory in large arrays:
    __slots__ = ('_std_dev', 'tag')

    def __init__(self, value, std_dev, tag=None):
        """
        The nominal value and the standard deviation of the variable
        are set.  These values must be scalars.

        'tag' is a tag that the user can associate to the variable.  This
        is useful for tracing variables.

        The meaning of the nominal value is described in the main
        module documentation.
        """

        #! The value, std_dev, and tag are assumed by __copy__() not to
        # be copied.  Either this should be guaranteed here, or __copy__
        # should be updated.

        # Only float-like values are handled.  One reason is that the
        # division operator on integers would not produce a
        # differentiable functions: for instance, Variable(3, 0.1)/2
        # has a nominal value of 3/2 = 1, but a "shifted" value
        # of 3.1/2 = 1.55.
        value = float(value)

        # If the variable changes by dx, then the value of the affine
        # function that gives its value changes by 1*dx:

        # ! Memory cycles are created.  However, they are garbage
        # collected, if possible.  Using a weakref.WeakKeyDictionary
        # takes much more memory.  Thus, this implementation chooses
        # more cycles and a smaller memory footprint instead of no
        # cycles and a larger memory footprint.
        super(Variable, self).__init__(value, {self: 1.})

        self.std_dev = std_dev  # Assignment through a Python property
        
        self.tag = tag

    # !! In Python 2.6+, the std_dev property would be more simply
    # implemented with @property(getter), then @std_dev.setter(setter).

    # Standard deviations can be modified (this is a feature).
    # AffineScalarFunc objects that depend on the Variable have their
    # std_dev automatically modified (recalculated with the new
    # std_dev of their Variables):
    def _set_std_dev(self, std_dev):
    
        # We force the error to be float-like.  Since it is considered
        # as a standard deviation, it must be positive:
        assert std_dev >= 0 or isnan(std_dev), (
            "the error must be a positive number, or NaN")

        self._std_dev = CallableStdDev(std_dev)
    
    def _get_std_dev(self):
        return self._std_dev
        
    std_dev = property(_get_std_dev, _set_std_dev)
    
    # Support for legacy method:
    def set_std_dev(self, value):  # Obsolete
        deprecation('instead of set_std_dev(), please use'
                    ' .std_dev = ...')
        self.std_dev = value
        
    # The following method is overridden so that we can represent the tag:
    def _general_representation(self, to_string):
        """
        Uses the to_string() conversion function on both the nominal
        value and standard deviation and returns a string that
        describes the number.

        to_string() is typically repr() or str().
        """
        num_repr  = super(Variable, self)._general_representation(to_string)
        
        # Optional tag: only full representations (to_string == repr)
        # contain the tag, as the tag is required in order to recreate
        # the variable.  Outputting the tag for regular string ("print
        # x") would be too heavy and produce an unusual representation
        # of a number with uncertainty.
        return (num_repr if ((self.tag is None) or (to_string != repr))
                else "< %s = %s >" % (self.tag, num_repr))

    def __hash__(self):
        # All Variable objects are by definition independent
        # variables, so they never compare equal; therefore, their
        # id() are therefore allowed to differ
        # (http://docs.python.org/reference/datamodel.html#object.__hash__):
        return id(self)
            
    def __copy__(self):
        """
        Hook for the standard copy module.
        """
        
        # This copy implicitly takes care of the reference of the
        # variable to itself (in self.derivatives): the new Variable
        # object points to itself, not to the original Variable.

        # Reference: http://www.doughellmann.com/PyMOTW/copy/index.html

        #! The following assumes that the arguments to Variable are
        # *not* copied upon construction, since __copy__ is not supposed
        # to copy "inside" information:
        return Variable(self.nominal_value, self.std_dev, self.tag)

    def __deepcopy__(self, memo):
        """
        Hook for the standard copy module.

        A new variable is created.
        """
        
        # This deep copy implicitly takes care of the reference of the
        # variable to itself (in self.derivatives): the new Variable
        # object points to itself, not to the original Variable.

        # Reference: http://www.doughellmann.com/PyMOTW/copy/index.html
        
        return self.__copy__()

        
###############################################################################

# Utilities

def nominal_value(x):
    """
    Returns the nominal value of x if it is a quantity with
    uncertainty (i.e., an AffineScalarFunc object); otherwise, returns
    x unchanged.

    This utility function is useful for transforming a series of
    numbers, when only some of them generally carry an uncertainty.
    """

    return x.nominal_value if isinstance(x, AffineScalarFunc) else x

def std_dev(x):
    """
    Returns the standard deviation of x if it is a quantity with
    uncertainty (i.e., an AffineScalarFunc object); otherwise, returns
    the float 0.

    This utility function is useful for transforming a series of
    numbers, when only some of them generally carry an uncertainty.
    """

    return x.std_dev if isinstance(x, AffineScalarFunc) else 0.

def covariance_matrix(nums_with_uncert):
    """
    Returns a matrix that contains the covariances between the given
    sequence of numbers with uncertainties (AffineScalarFunc objects).
    The resulting matrix implicitly depends on their ordering in
    'nums_with_uncert'.

    The covariances are floats (never int objects).

    The returned covariance matrix is the exact linear approximation
    result, if the nominal values of the numbers with uncertainties
    and of their variables are their mean.  Otherwise, the returned
    covariance matrix should be close to its linear approximation
    value.

    The returned matrix is a list of lists.
    """
    # See PSI.411 in EOL's notes.

    covariance_matrix = []
    # ! In Python 2.6+, this could be replaced by enumerate(..., 1),
    # along with updating places where i1 is used (i1+1 => i1).
    for (i1, expr1) in enumerate(nums_with_uncert):
        derivatives1 = expr1.derivatives  # Optimization
        vars1 = set(derivatives1)
        coefs_expr1 = []
        for expr2 in nums_with_uncert[:i1+1]:
            derivatives2 = expr2.derivatives  # Optimization
            coef = 0.
            for var in vars1.intersection(derivatives2):
                # var is a variable common to both numbers with
                # uncertainties:
                coef += (derivatives1[var]*derivatives2[var]*var._std_dev**2)
            coefs_expr1.append(coef)
        covariance_matrix.append(coefs_expr1)

    # We symmetrize the matrix:
    for (i, covariance_coefs) in enumerate(covariance_matrix):
        covariance_coefs.extend(covariance_matrix[j][i]
                                for j in range(i+1, len(covariance_matrix)))

    return covariance_matrix

try:
    import numpy
except ImportError:
    pass
else:
    def correlation_matrix(nums_with_uncert):
        '''
        Returns the correlation matrix of the given sequence of
        numbers with uncertainties, as a NumPy array of floats.
        '''

        cov_mat = numpy.array(covariance_matrix(nums_with_uncert))

        std_devs = numpy.sqrt(cov_mat.diagonal())
        
        return cov_mat/std_devs/std_devs[numpy.newaxis].T

    __all__.append('correlation_matrix')
        
###############################################################################
# Parsing of values with uncertainties:

POSITIVE_DECIMAL_UNSIGNED = r'(\d+)(\.\d*)?'

# Regexp for a number with uncertainty (e.g., "-1.234(2)e-6"), where the
# uncertainty is optional (in which case the uncertainty is implicit):
NUMBER_WITH_UNCERT_RE_STR = '''
    ([+-])?  # Sign
    %s  # Main number
    (?:\(%s\))?  # Optional uncertainty
    ([eE][+-]?\d+)?  # Optional exponent
    ''' % (POSITIVE_DECIMAL_UNSIGNED, POSITIVE_DECIMAL_UNSIGNED)

NUMBER_WITH_UNCERT_RE_SEARCH = re.compile(
    "^%s$" % NUMBER_WITH_UNCERT_RE_STR, re.VERBOSE).search

class NotParenForm(ValueError):
    '''
    Raised when a string representing an exact number or a number with
    an uncertainty indicated between parentheses was expected but not
    found.
    '''
    
def parse_error_in_parentheses(representation):
    """
    Returns (value, error) from a string representing a number with
    uncertainty like 12.34(5), 12.34(142), 12.5(3.4) or 12.3(4.2)e3.
    If no parenthesis is given, an uncertainty of one on the last
    digit is assumed.

    Raises ValueError if the string cannot be parsed.    
    """

    match = NUMBER_WITH_UNCERT_RE_SEARCH(representation)

    if match:
        # The 'main' part is the nominal value, with 'int'eger part, and
        # 'dec'imal part.  The 'uncert'ainty is similarly broken into its
        # integer and decimal parts.
        (sign, main_int, main_dec, uncert_int, uncert_dec,
         exponent) = match.groups()
    else:
        raise NotParenForm("Unparsable number representation: '%s'."
                           " Was expecting a string of the form 1.23(4)"
                           " or 1.234" % representation)

    # The value of the number is its nominal value:
    value = float(''.join((sign or '',
                           main_int,
                           main_dec or '.0',
                           exponent or '')))
                  
    if uncert_int is None:
        # No uncertainty was found: an uncertainty of 1 on the last
        # digit is assumed:
        uncert_int = '1'

    # Do we have a fully explicit uncertainty?
    if uncert_dec is not None:
        uncert = float("%s%s" % (uncert_int, uncert_dec or ''))
    else:
        # uncert_int represents an uncertainty on the last digits:

        # The number of digits after the period defines the power of
        # 10 than must be applied to the provided uncertainty:
        num_digits_after_period = (0 if main_dec is None
                                   else len(main_dec)-1)
        uncert = int(uncert_int)/10**num_digits_after_period

    # We apply the exponent to the uncertainty as well:
    uncert *= float("1%s" % (exponent or ''))

    return (value, uncert)


_cannot_parse_ufloat_msg_pat = (
    'Cannot parse %s: see the documentation of ufloat_fromstr() for a'
    ' list of accepted formats')

# The following function is not exposed because it can in effect be
# obtained by doing x = ufloat_fromstr(representation) and reading
# x.nominal_value and x.std_dev:
def _str_to_number_with_uncert(representation):
    """
    Given a string that represents a number with uncertainty, returns the
    nominal value and the uncertainty.

    See the documentation of ufloat_fromstr() for a list of accepted
    formats.

    When no numerical error is given, an uncertainty of 1 on the last
    digit is implied.

    Raises ValueError if the string cannot be parsed.
    """

    try:
        # Simple form 1234.45+/-1.2:
        (value, uncert) = representation.split('+/-')
    except ValueError:
        # Form with parentheses or no uncertainty:
        try:
            parsed_value = parse_error_in_parentheses(representation)
        except NotParenForm:
            raise ValueError(_cannot_parse_ufloat_msg_pat % representation)
    else:
        try:
            parsed_value = (float(value), float(uncert))
        except ValueError:
            raise ValueError(_cannot_parse_ufloat_msg_pat % representation)
        
    return parsed_value

def ufloat_fromstr(representation, tag=None):
    """
    Returns a new random variable (Variable object) from a string.
    
    Strings 'representation' of the form '12.345+/-0.015',
    '12.345(15)', or '12.3' are recognized (see more complete list
    below).  In the last case, an uncertainty of +/-1 is assigned to
    the last digit.
    
    Examples of valid string representations:

        12.3e10+/-5e3
        -1.23(3.4)
        -1.34(5)
        1(6)
        3(4.2)
        -9(2)
        1234567(1.2)
        12.345(15)
        -12.3456(78)e-6
        12.3(0.4)e-5        
        0.29
        31.
        -31.
        31
        -3.1e10
        169.0(7)
        169.1(15)
    """

    #! The special ** syntax is for Python 2.5 and before (Python 2.6+
    # understands tag=tag):
    (nominal_value, std_dev) = _str_to_number_with_uncert(representation)
    return ufloat(nominal_value, std_dev, tag)

def _ufloat_obsolete(representation, tag=None):
    '''
    Legacy version of ufloat(). Will eventually be removed.

    representation -- either a (nominal_value, std_dev) tuple, or a
    string representation of a number with uncertainty, in a format
    recognized by ufloat_fromstr().
    '''
    return (ufloat(representation[0], representation[1], tag)
            if isinstance(representation, tuple)
            else ufloat_fromstr(representation, tag))

# The arguments are named for the new version, instead of bearing
# names that are closer to their obsolete use (e.g., std_dev could be
# instead std_dev_or_tag, since it can be the tag, in the obsolete
# ufloat((3, 0.14), "pi") form). This has the advantage of allowing
# new code to use keyword arguments as in ufloat(nominal_value=3,
# std_dev=0.14), without breaking when the obsolete form is not
# supported anymore.
def ufloat(nominal_value, std_dev=None, tag=None):
    """
    Returns a new random variable (Variable object).
    
    The only non-obsolete use is:

    - ufloat(nominal_value, std_dev),
    - ufloat(nominal_value, std_dev, tag=...).

    Other input parameters are temporarily supported:

    - ufloat((nominal_value, std_dev)),
    - ufloat((nominal_value, std_dev), tag),
    - ufloat(str_representation),
    - ufloat(str_representation, tag).

    Valid string representations str_representation are listed in
    the documentation for ufloat_fromstr().

    'tag' is an optional string tag for the variable.  Variables
    don't have to have distinct tags.  Tags are useful for tracing
    what values (and errors) enter in a given result (through the
    error_components() method).
    """

    try:
        # Standard case:

        #! The special ** syntax is for Python 2.5 and before (Python 2.6+
        # understands tag=tag):
        return Variable(nominal_value, std_dev, **{'tag': tag})
    # Exception types raised by, respectively: tuple, string that
    # cannot be converted through float(), and string that can be
    # converted through float() (case of a number with no uncertainty):
    except (TypeError, ValueError, AssertionError):
        # Obsolete, two-argument call:
        deprecation('either use ufloat(nominal_value, std_dev),'
                    ' ufloat(nominal_value, std_dev, tag), or the'
                    ' ufloat_fromstr() function, for string representations.')
        return _ufloat_obsolete(nominal_value,  # Tuple or string
                                # tag keyword used:
                                tag if tag is not None
                                # 2 positional arguments form:
                                else std_dev)

