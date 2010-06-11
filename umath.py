'''
Mathematical operations that generalize many operations from the
standard math module so that they also work on numbers with
uncertainties.

Examples:

  from umath import sin
  
  # Manipulation of numbers with uncertainties:
  x = uncertainties.ufloat((3, 0.1))
  print sin(x)  # prints 0.141120008...+/-0.098999...

  # The umath functions also work on regular Python floats:
  print sin(3)  # prints 0.141120008...  This is a Python float.

Importing all the functions from this module into the global namespace
is possible.  This is encouraged when using a Python shell as a
calculator.  Example:

  import uncertainties
  from uncertainties.umath import *  # Imports tan(), etc.
  
  x = uncertainties.ufloat((3, 0.1))
  print tan(x)  # tan() is the uncertainties.umath.tan function

The numbers with uncertainties handled by this module are objects from
the uncertainties module, from either the Variable or the
AffineScalarFunc class.

(c) 2009-2010 by Eric O. LEBIGOT (EOL) <eric.lebigot@normalesup.org>.
Please send feature requests, bug reports, or feedback to this address.

This software is released under a dual license.  (1) The GNU General
Public License version 2.  (2) Any other license, as long as it is
obtained from the original author.'''

from __future__ import division  # Many analytical derivatives depend on this

# Standard modules
import math
import sys
import itertools
import functools
import inspect

# Local modules
import uncertainties

from uncertainties import __author__

###############################################################################

# We wrap the functions from the math module so that they keep track of
# uncertainties by returning a AffineScalarFunc object.

# Some functions from the math module cannot be simply adapted to work
# with AffineScalarFunc objects (either as their result or as their
# arguments):

# (1) Some functions return a result of a type whose value and
# variations (uncertainties) cannot be represented by AffineScalarFunc
# (e.g., math.frexp, which returns a tuple).  The exception raised
# when not wrapping them with wrap() is more obvious than the
# one obtained when wrapping them (in fact, the wrapped functions
# attempts operations that are not supported, such as calculation a
# subtraction on a result of type tuple).

# (2) Some functions don't take scalar arguments (which can be varied
# during differentiation): math.fsum...  Such functions can either be:

# - wrapped in a special way in wrap_math_functions()

# - excluded from wrapping by adding their name to no_std_wrapping

# - wrapped in a general way in wrap_math_functions(); in this case,
# the function does not have to be mentioned in any location in this
# code.  The function should function normally, except possibly when
# given AffineScalarFunc arguments.

no_std_wrapping = ['frexp', 'modf', 'fsum']  # Exclude from standard wrapping

std_wrapped_math_funcs = []  # Math functions wrapped in the standard way

non_std_wrapped_math_funcs = []  # Math functions wrapped in a non-standard way

# Function that copies the relevant attributes from generalized
# functions from the math module:
wraps = functools.partial(functools.update_wrapper,
                          assigned=('__doc__', '__name__'))

########################################
# Special cases:

# fsum takes a single argument, which cannot be differentiated.
# However, each of the arguments inside this single list can
# be a variable.  We handle this in a specific way:

if sys.version_info[:2] >= (2, 6):    

    original_func = math.fsum  # Shortcut

    # We wrap math.fsum

    def wrapped_fsum():
        """
        Return an uncertainty aware version of math.fsum, which must
        be contained in _original_func.
        """

        # The fsum function is flattened, in order to use the
        # wrap() wrapper:

        flat_fsum = lambda *args: original_func(args)

        flat_fsum_wrap = uncertainties.wrap(
            flat_fsum, itertools.repeat(lambda *args: 1))

        return wraps(lambda arg_list: flat_fsum_wrap(*arg_list),
                     original_func)

    fsum = wrapped_fsum()

    # Wrapping already done:
    non_std_wrapped_math_funcs.append('fsum')

########################################
# Wrapping of built-in math functions not in no_std_wrapping:

# Fixed formulas for the derivatives of some functions from the math
# module (some functions might not be present in all version of
# Python).  Singular points are not taken into account.  The user
# should never give "large" uncertainties: problems could only appear
# if this assumption does not hold.

# Functions not mentioned in _fixed_derivatives have their derivatives
# calculated numerically.

# Functions that have singularities (possibly at infinity) benefit
# from analytical calculations (instead of the default numerical
# calculation).  Even slowly varying functions (e.g., abs()) yield
# more precise results when differentiated analytically, because of
# the loss of precision in numerical calculations.

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
    'ceil': [lambda x: 0],
    'copysign': [lambda x, y: (1 if x >= 0 else -1) * math.copysign(1, y),
                 lambda x, y: 0],
    'cos': [lambda x: -math.sin(x)],
    'cosh': [math.sinh],
    'degrees': [lambda x: math.degrees(1)],
    'exp': [math.exp],
    'fabs': [lambda x: 1 if x >= 0 else -1],
    'floor': [lambda x: 0],
    'hypot': [lambda x, y: x/math.hypot(x, y),
              lambda x, y: y/math.hypot(x, y)],
    'ldexp': [lambda x, y: 2**y,
              # math.ldexp only accepts an integer as its second
              # argument:
              None],
    'log': [log_der0,
            lambda x, y: -math.log(x, y)/y/math.log(y)],
    'log10': [lambda x: 1/x/math.log(10)],
    'log1p': [lambda x: 1/(1+x)],
    'pow': [lambda x, y: y*math.pow(x, y-1),
            lambda x, y: math.log(x) * math.pow(x, y)],
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

# We do not want to wrap module attributes such as __doc__, etc.:
for (name, func) in inspect.getmembers(math, inspect.isbuiltin):

    if name in no_std_wrapping:
        continue

    if name in fixed_derivatives:
        derivatives = fixed_derivatives[name]
    else:
        derivatives = None  # Means: numerical calculation required
    setattr(this_module, name,
            wraps(uncertainties.wrap(func, derivatives), func))
    std_wrapped_math_funcs.append(name)

###############################################################################
# Exported functions:

__all__ = std_wrapped_math_funcs + non_std_wrapped_math_funcs

