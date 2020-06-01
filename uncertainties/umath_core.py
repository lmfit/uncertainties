# !!!!!!!!!!! Add a header to the documentation, that starts with something
# like "uncertainties.UFloat-compatible version of...", for all functions.

"""
Implementation of umath.py, with internals.
"""

# This module exists so as to define __all__, which in turn defines
# which functions are visible to the user in umath.py through from
# umath import * and Python shell completion.

from __future__ import division  # Many analytical derivatives depend on this

# Standard modules
from builtins import map
import math
import sys
import itertools

# Local modules
import uncertainties.core as uncert_core
from uncertainties.core import (to_affine_scalar, AffineScalarFunc,
                                LinearCombination)

###############################################################################

# We wrap the functions from the math module so that they keep track of
# uncertainties by returning a AffineScalarFunc object.

# Some functions from the math module cannot be adapted in a standard
# way so to work with AffineScalarFunc objects (either as their result
# or as their arguments):

# (1) Some functions return a result of a type whose value and
# variations (uncertainties) cannot be represented by AffineScalarFunc
# (e.g., math.frexp, which returns a tuple).  The exception raised
# when not wrapping them with wrap() is more obvious than the
# one obtained when wrapping them (in fact, the wrapped functions
# attempts operations that are not supported, such as calculation a
# subtraction on a result of type tuple).

# (2) Some functions don't take continuous scalar arguments (which can
# be varied during differentiation): math.fsum, math.factorial...
# Such functions can either be:

# - wrapped in a special way.

# - excluded from standard wrapping by adding their name to
# no_std_wrapping

# Math functions that have a standard interface: they take
# one or more float arguments, and return a scalar:
many_scalars_to_scalar_funcs = []

# Some functions require a specific treatment and must therefore be
# excluded from standard wrapping.  Functions
# no_std_wrapping = ['modf', 'frexp', 'ldexp', 'fsum', 'factorial']

# Functions with numerical derivatives:
#
# !! Python2.7+: {..., ...}
num_deriv_funcs = set(['fmod', 'gamma', 'lgamma'])

# Functions are by definition locally constant (on real
# numbers): their value does not depend on the uncertainty (because
# this uncertainty is supposed to lead to a good linear approximation
# of the function in the uncertainty region). The type of their output
# for floats is preserved, as users should not care about deviations
# in their value: their value is locally constant due to the nature of
# the function (0 derivative). This situation is similar to that of
# comparisons (==, >, etc.).
#
# !! Python 2.7+: {..., ...}
locally_cst_funcs = set(['ceil', 'floor', 'isinf', 'isnan', 'trunc'])

# Functions that do not belong in many_scalars_to_scalar_funcs, but
# that have a version that handles uncertainties. These functions are
# also not in numpy (see unumpy/core.py).
non_std_wrapped_funcs = []

# Function that copies the relevant attributes from generalized
# functions from the math module:
# This is a copy&paste job from the functools module, changing
# the default arugment for assigned
def wraps(wrapper,
          wrapped,
          assigned=('__doc__',),
          updated=('__dict__',)):
    """Update a wrapper function to look like the wrapped function.

    wrapper -- function to be updated
    wrapped -- original function
    assigned -- tuple naming the attributes assigned directly
    from the wrapped function to the wrapper function
    updated -- tuple naming the attributes of the wrapper that
    are updated with the corresponding attribute from the wrapped
    function.
    """
    for attr in assigned:
        setattr(wrapper, attr, getattr(wrapped, attr))
    for attr in updated:
        getattr(wrapper, attr).update(getattr(wrapped, attr, {}))
    # Return the wrapper so this can be used as a decorator via partial()
    return wrapper


########################################
# Wrapping of math functions:

# Fixed formulas for the derivatives of some functions from the math
# module (some functions might not be present in all version of
# Python).  Singular points are not taken into account.  The user
# should never give "large" uncertainties: problems could only appear
# if this assumption does not hold.

# Functions not mentioned in _fixed_derivatives have their derivatives
# calculated numerically.

# Functions that have singularities (possibly at infinity) benefit
# from analytical calculations (instead of the default numerical
# calculation) because their derivatives generally change very fast.
# Even slowly varying functions (e.g., abs()) yield more precise
# results when differentiated analytically, because of the loss of
# precision in numerical calculations.

#def log_1arg_der(x):
#    """
#    Derivative of log(x) (1-argument form).
#    """
#    return 1/x

def log_der0(*args):
    """
    Derivative of math.log() with respect to its first argument.

    Works whether 1 or 2 arguments are given.
    """
    if len(args) == 1:
        return 1/args[0]
    else:
        return 1/args[0]/math.log(args[1])  # 2-argument form

    # The following version goes about as fast:

    ## A 'try' is used for the most common case because it is fast when no
    ## exception is raised:
    #try:
    #    return log_1arg_der(*args)  # Argument number check
    #except TypeError:
    #    return 1/args[0]/math.log(args[1])  # 2-argument form

def _deriv_copysign(x,y):
    if x >= 0:
        return math.copysign(1, y)
    else:
        return -math.copysign(1, y)

def _deriv_fabs(x):
    if x >= 0:
        return 1
    else:
        return -1

def _deriv_pow_0(x, y):
    if y == 0:
        return  0.
    elif x != 0 or y % 1 == 0:
        return y*math.pow(x, y-1)
    else:
        return float('nan')

def _deriv_pow_1(x, y):
    if x == 0 and y > 0:
        return 0.
    else:
        return math.log(x) * math.pow(x, y)

erf_coef = 2/math.sqrt(math.pi)  # Optimization for erf()

fixed_derivatives = {
    # In alphabetical order, here:
    'acos': [lambda x: -1/math.sqrt(1-x**2)],
    'acosh': [lambda x: 1/math.sqrt(x**2-1)],
    'asin': [lambda x: 1/math.sqrt(1-x**2)],
    'asinh': [lambda x: 1/math.sqrt(1+x**2)],
    'atan': [lambda x: 1/(1+x**2)],
    'atan2': [lambda y, x: x/(x**2+y**2),  # Correct for x == 0
              lambda y, x: -y/(x**2+y**2)],  # Correct for x == 0
    'atanh': [lambda x: 1/(1-x**2)],
    'copysign': [_deriv_copysign,
                 lambda x, y: 0],
    'cos': [lambda x: -math.sin(x)],
    'cosh': [math.sinh],
    'degrees': [lambda x: math.degrees(1)],
    'erf': [lambda x: math.exp(-x**2)*erf_coef],
    'erfc': [lambda x: -math.exp(-x**2)*erf_coef],
    'exp': [math.exp],
    'expm1': [math.exp],
    'fabs': [_deriv_fabs],
    'hypot': [lambda x, y: x/math.hypot(x, y),
              lambda x, y: y/math.hypot(x, y)],
    'log': [log_der0,
            lambda x, y: -math.log(x, y)/y/math.log(y)],
    'log10': [lambda x: 1/x/math.log(10)],
    'log1p': [lambda x: 1/(1+x)],
    'pow': [_deriv_pow_0, _deriv_pow_1],
    'radians': [lambda x: math.radians(1)],
    'sin': [math.cos],
    'sinh': [math.cosh],
    'sqrt': [lambda x: 0.5/math.sqrt(x)],
    'tan': [lambda x: 1+math.tan(x)**2],
    'tanh': [lambda x: 1-math.tanh(x)**2]
    }

# Many built-in functions in the math module are wrapped with a
# version which is uncertainty aware:

this_module = sys.modules[__name__]

def wrap_locally_cst_func(func):
    '''
    Return a function that returns the same arguments as func, but
    after converting any AffineScalarFunc object to its nominal value.

    This function is useful for wrapping functions that are locally
    constant: the uncertainties should have no role in the result
    (since they are supposed to keep the function linear and hence,
    here, constant).
    '''
    def wrapped_func(*args, **kwargs):
        args_float = map(uncert_core.nominal_value, args)
        # !! In Python 2.7+, dictionary comprehension: {argname:...}
        kwargs_float = dict(
            (arg_name, uncert_core.nominal_value(value))
            for (arg_name, value) in kwargs.items())
        return func(*args_float, **kwargs_float)
    return wrapped_func

# for (name, attr) in vars(math).items():
for name in dir(math):

    if name in fixed_derivatives:  # Priority to functions in fixed_derivatives
        derivatives = fixed_derivatives[name]
    elif name in num_deriv_funcs:
        # Functions whose derivatives are calculated numerically by
        # this module fall here (isinf, fmod,...):
        derivatives = []  # Means: numerical calculation required
    elif name not in locally_cst_funcs:
        continue  # 'name' not wrapped by this module (__doc__, e, etc.)

    func = getattr(math, name)

    if name in locally_cst_funcs:
        wrapped_func = wrap_locally_cst_func(func)
    else:  # Function with analytical or numerical derivatives:
        # Errors during the calculation of the derivatives are converted
        # to a NaN result: it is assumed that a mathematical calculation
        # that cannot be calculated indicates a non-defined derivative
        # (the derivatives in fixed_derivatives must be written this way):
        wrapped_func = uncert_core.wrap(
            func, map(uncert_core.nan_if_exception, derivatives))

    # !! The same effect could be achieved with globals()[...] = ...
    setattr(this_module, name, wraps(wrapped_func, func))

    many_scalars_to_scalar_funcs.append(name)

###############################################################################

########################################
# Special cases: some of the functions from no_std_wrapping:

##########
# The math.factorial function is not converted to an uncertainty-aware
# function, because it does not handle non-integer arguments: it does
# not make sense to give it an argument with a numerical error
# (whereas this would be relevant for the gamma function).

##########

# fsum takes a single argument, which cannot be differentiated.
# However, each of the arguments inside this single list can
# be a variable.  We handle this in a specific way:

# Only for Python 2.6+:

# For drop-in compatibility with the math module:
factorial = math.factorial
non_std_wrapped_funcs.append('factorial')


# We wrap math.fsum

original_func = math.fsum  # For optimization purposes

# The function below exists so that temporary variables do not
# pollute the module namespace:
def wrapped_fsum():
    """
    Return an uncertainty-aware version of math.fsum, which must
    be contained in _original_func.
    """

    # The fsum function is flattened, in order to use the
    # wrap() wrapper:

    flat_fsum = lambda *args: original_func(args)

    flat_fsum_wrap = uncert_core.wrap(
        flat_fsum, itertools.repeat(lambda *args: 1))

    return wraps(lambda arg_list: flat_fsum_wrap(*arg_list),
                 original_func)

# !!!!!!!! Documented?
fsum = wrapped_fsum()
non_std_wrapped_funcs.append('fsum')

##########

# Some functions that either return multiple arguments (modf, frexp)
# or take some non-float arguments (which should not be converted to
# numbers with uncertainty).

# ! The arguments have the same names as in the math module
# documentation, so that the docstrings are consistent with them.

@uncert_core.set_doc(math.modf.__doc__)
def modf(x):
    """
    Version of modf that works for numbers with uncertainty, and also
    for regular numbers.
    """

    # The code below is inspired by uncert_core.wrap().  It is
    # simpler because only 1 argument is given, and there is no
    # delegation to other functions involved (as for __mul__, etc.).

    aff_func = to_affine_scalar(x)  # Uniform treatment of all numbers

    (frac_part, int_part) = math.modf(aff_func.nominal_value)

    if aff_func._linear_part:  # If not a constant
        # The derivative of the fractional part is simply 1: the
        # linear part of modf(x)[0] is the linear part of x:
        return (AffineScalarFunc(frac_part, aff_func._linear_part), int_part)
    else:
        # This function was not called with an AffineScalarFunc
        # argument: there is no need to return numbers with uncertainties:
        return (frac_part, int_part)
many_scalars_to_scalar_funcs.append('modf')

@uncert_core.set_doc(math.ldexp.__doc__)
def ldexp(x, i):
    # Another approach would be to add an additional argument to
    # uncert_core.wrap() so that some arguments are automatically
    # considered as constants.

    aff_func = to_affine_scalar(x)  # y must be an integer, for math.ldexp

    if aff_func._linear_part:
        return AffineScalarFunc(
            math.ldexp(aff_func.nominal_value, i),
            LinearCombination([(2**i, aff_func._linear_part)]))
    else:
        # This function was not called with an AffineScalarFunc
        # argument: there is no need to return numbers with uncertainties:

        # aff_func.nominal_value is not passed instead of x, because
        # we do not have to care about the type of the return value of
        # math.ldexp, this way (aff_func.nominal_value might be the
        # value of x coerced to a difference type [int->float, for
        # instance]):
        return math.ldexp(x, i)
many_scalars_to_scalar_funcs.append('ldexp')

@uncert_core.set_doc(math.frexp.__doc__)
def frexp(x):
    """
    Version of frexp that works for numbers with uncertainty, and also
    for regular numbers.
    """

    # The code below is inspired by uncert_core.wrap().  It is
    # simpler because only 1 argument is given, and there is no
    # delegation to other functions involved (as for __mul__, etc.).

    aff_func = to_affine_scalar(x)

    if aff_func._linear_part:
        (mantissa, exponent) = math.frexp(aff_func.nominal_value)
        return (
            AffineScalarFunc(
                mantissa,
                # With frexp(x) = (m, e), x = m*2**e, so m = x*2**-e
                # and therefore dm/dx = 2**-e (as e in an integer that
                # does not vary when x changes):
                LinearCombination([2**-exponent, aff_func._linear_part])),
            # The exponent is an integer and is supposed to be
            # continuous (errors must be small):
            exponent)
    else:
        # This function was not called with an AffineScalarFunc
        # argument: there is no need to return numbers with uncertainties:
        return math.frexp(x)
non_std_wrapped_funcs.append('frexp')

###############################################################################
# Exported functions:

__all__ = many_scalars_to_scalar_funcs + non_std_wrapped_funcs
