# coding=utf-8

"""
Main module for the uncertainties package, with internal functions.
"""

# The idea behind this module is to replace the result of mathematical
# operations by a local approximation of the defining function.  For
# example, sin(0.2+/-0.01) becomes the affine function
# (AffineScalarFunc object) whose nominal value is sin(0.2) and
# whose variations are given by sin(0.2+delta) = 0.98...*delta.
# Uncertainties can then be calculated by using this local linear
# approximation of the original function.

from __future__ import division  # Many analytical derivatives depend on this

from builtins import str, next, map, zip, range, object
import math
from math import sqrt, log, isnan, isinf  # Optimization: no attribute look-up
import re
import sys
try:
    from math import isinfinite  # !! Python 3.2+
except ImportError:
    def isinfinite(x):
        return isinf(x) or isnan(x)

import copy
import warnings
import itertools
import inspect
import numbers
import collections

# The following restricts the local function getargspec() to the common
# features of inspect.getargspec() and inspect.getfullargspec():
if sys.version_info < (3,):  # !! Could be removed when moving to Python 3 only
    from inspect import getargspec
else:
    from inspect import getfullargspec as getargspec

from . compat import (deprecation, FLOAT_LIKE_TYPES, CONSTANT_TYPES, CallableStdDev,
)
from . util import (set_doc, covariance_matrix, partial_derivative,
NumericalDerivatives, IndexableIter, basestring,
)
from . formatting import (first_digit, PDG_precision, format_num, signif_dgt_to_limit
)
from . parsing import (str_to_number_with_uncert,)

try:
    import numpy
except ImportError:
    pass
else:
    from . util import correlation_matrix, correlated_values, correlated_values_norm

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

    Raises an exception unless x belongs to some specific classes of
    objects that are known not to depend on AffineScalarFunc objects
    (which then cannot be considered as constants).
    """

    if isinstance(x, AffineScalarFunc):
        return x

    if isinstance(x, CONSTANT_TYPES):
        # No variable => no derivative:
        return AffineScalarFunc(x, LinearCombination({}))

    # Case of lists, etc.
    raise NotUpcast("%s cannot be converted to a number with"
                    " uncertainty" % type(x))


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
    distribution over the real numbers.

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
    defined, for instance by using the
    uncertainties.core.partial_derivative() function.

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

    derivatives_all_kwargs = {}

    for (name, derivative) in derivatives_kwargs.items():

        # Optimization: None keyword-argument derivatives are converted
        # right away to derivatives (instead of doing this every time a
        # None derivative is encountered when calculating derivatives):

        if derivative is None:
            derivatives_all_kwargs[name] = partial_derivative(f, name)
        else:
            derivatives_all_kwargs[name] = derivative

    # When the wrapped function is called with keyword arguments that
    # map to positional-or-keyword parameters, their derivative is
    # looked for in derivatives_all_kwargs.  We define these
    # additional derivatives:

    try:
        argspec = getargspec(f)
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

        for (index, name) in enumerate(argspec.args):

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

            if derivative is None:
                derivatives_all_kwargs[name] = partial_derivative(f, name)
            else:
                derivatives_all_kwargs[name] = derivative

    # Optimization: None derivatives for the positional arguments are
    # converted to the corresponding numerical differentiation
    # function (instead of doing this over and over later every time a
    # None derivative is found):

    none_converter = lambda index: partial_derivative(f, index)

    for (index, derivative) in enumerate(
        derivatives_args_index.returned_elements):
        if derivative is None:
            derivatives_args_index.returned_elements[index] = (
                none_converter(index))

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
        # replaced by their nominal value in order to calculate
        # the necessary derivatives of f.

        pos_w_uncert = [index for (index, value) in enumerate(args)
                        if isinstance(value, AffineScalarFunc)]
        names_w_uncert = [key for (key, value) in kwargs.items()
                          if isinstance(value, AffineScalarFunc)]

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
        # saved in the following dictionary (which only contains
        # values with uncertainty):

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

        # Calculation of the linear part of the function value,
        # defined by (coefficient, argument) pairs, where 'argument'
        # is an AffineScalarFunc (for all AffineScalarFunc found as
        # argument of f):
        linear_part = []

        for pos in pos_w_uncert:
            linear_part.append((
                # Coefficient:
                derivatives_args_index[pos](*args_values, **kwargs),
                # Linear part of the AffineScalarFunc expression:
                args[pos]._linear_part))

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

            linear_part.append((
                # Coefficient:
                derivative(*args_values, **kwargs),
                # Linear part of the AffineScalarFunc expression:
                kwargs_uncert_values[name]._linear_part))

        # The function now returns the necessary linear approximation
        # to the function:
        return AffineScalarFunc(
            f_nominal_value, LinearCombination(linear_part))

    f_with_affine_output = set_doc("""\
    Version of %s(...) that returns an affine approximation
    (AffineScalarFunc object), if its result depends on variables
    (Variable objects).  Otherwise, returns a simple constant (when
    applied to constant arguments).

    Warning: arguments of the function that are not AffineScalarFunc
    objects must not depend on uncertainties.Variable objects in any
    way.  Otherwise, the dependence of the result in
    uncertainties.Variable objects will be incorrect.

    Original documentation:
    %s""" % (f.__name__, f.__doc__))(f_with_affine_output)

    # It is easier to work with f_with_affine_output, which represents
    # a wrapped version of 'f', when it bears the same name as 'f':
    # ! __name__ is read-only, in Python 2.3:
    f_with_affine_output.name = f.__name__

    return f_with_affine_output


def force_aff_func_args(func):
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
        Return %s(self, to_affine_scalar(y)) if y can be upcast
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

def eq_on_aff_funcs(self, y_with_uncert):
    """
    __eq__ operator, assuming that both self and y_with_uncert are
    AffineScalarFunc objects.
    """
    difference = self - y_with_uncert
    # Only an exact zero difference means that self and y are
    # equal numerically:
    return not(difference._nominal_value or difference.std_dev)

def ne_on_aff_funcs(self, y_with_uncert):
    """
    __ne__ operator, assuming that both self and y_with_uncert are
    AffineScalarFunc objects.
    """

    return not eq_on_aff_funcs(self, y_with_uncert)

def gt_on_aff_funcs(self, y_with_uncert):
    """
    __gt__ operator, assuming that both self and y_with_uncert are
    AffineScalarFunc objects.
    """
    return self._nominal_value > y_with_uncert._nominal_value

def ge_on_aff_funcs(self, y_with_uncert):
    """
    __ge__ operator, assuming that both self and y_with_uncert are
    AffineScalarFunc objects.
    """

    return (gt_on_aff_funcs(self, y_with_uncert)
            or eq_on_aff_funcs(self, y_with_uncert))

def lt_on_aff_funcs(self, y_with_uncert):
    """
    __lt__ operator, assuming that both self and y_with_uncert are
    AffineScalarFunc objects.
    """
    return self._nominal_value < y_with_uncert._nominal_value

def le_on_aff_funcs(self, y_with_uncert):
    """
    __le__ operator, assuming that both self and y_with_uncert are
    AffineScalarFunc objects.
    """

    return (lt_on_aff_funcs(self, y_with_uncert)
            or eq_on_aff_funcs(self, y_with_uncert))



class LinearCombination(object):
    """
    Linear combination of Variable differentials.

    The linear_combo attribute can change formally, but its value
    always remains the same. Typically, the linear combination can
    thus be expanded.

    The expanded form of linear_combo is a mapping from Variables to
    the coefficient of their differential.
    """

    # ! Invariant: linear_combo is represented internally exactly as
    # the linear_combo argument to __init__():
    __slots__ = "linear_combo"

    def __init__(self, linear_combo):
        """
        linear_combo can be modified by the object, during its
        lifetime. This allows the object to change its internal
        representation over time (for instance by expanding the linear
        combination and replacing the original expression with the
        expanded one).

        linear_combo -- if linear_combo is a dict, then it represents
        an expanded linear combination and must map Variables to the
        coefficient of their differential. Otherwise, it should be a
        list of (coefficient, LinearCombination) pairs (that
        represents a linear combination expression).
        """

        self.linear_combo = linear_combo

    def __bool__(self):
        """
        Return True only if the linear combination is non-empty, i.e. if
        the linear combination contains any term.
        """
        return bool(self.linear_combo)

    def expanded(self):
        """
        Return True if and only if the linear combination is expanded.
        """
        return isinstance(self.linear_combo, dict)

    def expand(self):
        """
        Expand the linear combination.

        The expansion is a collections.defaultdict(float).

        This should only be called if the linear combination is not
        yet expanded.
        """

        # The derivatives are built progressively by expanding each
        # term of the linear combination until there is no linear
        # combination to be expanded.

        # Final derivatives, constructed progressively:
        derivatives = collections.defaultdict(float)

        while self.linear_combo:  # The list of terms is emptied progressively

            # One of the terms is expanded or, if no expansion is
            # needed, simply added to the existing derivatives.
            #
            # Optimization note: since Python's operations are
            # left-associative, a long sum of Variables can be built
            # such that the last term is essentially a Variable (and
            # not a NestedLinearCombination): popping from the
            # remaining terms allows this term to be quickly put in
            # the final result, which limits the number of terms
            # remaining (and whose size can temporarily grow):
            (main_factor, main_expr) = self.linear_combo.pop()

            # print "MAINS", main_factor, main_expr

            if main_expr.expanded():
                for (var, factor) in main_expr.linear_combo.items():
                    derivatives[var] += main_factor*factor

            else:  # Non-expanded form
                for (factor, expr) in main_expr.linear_combo:
                    # The main_factor is applied to expr:
                    self.linear_combo.append((main_factor*factor, expr))

            # print "DERIV", derivatives

        self.linear_combo = derivatives

    def __getstate__(self):
        # Not false, otherwise __setstate__() will not be called:
        return (self.linear_combo,)

    def __setstate__(self, state):
        (self.linear_combo,) = state

class AffineScalarFunc(object):
    """
    Affine functions that support basic mathematical operations
    (addition, etc.).  Such functions can for instance be used for
    representing the local (linear) behavior of any function.

    This class can also be used to represent constants.

    The variables of affine scalar functions are Variable objects.

    AffineScalarFunc objects include facilities for calculating the
    'error' on the function, from the uncertainties on its variables.

    Main attributes and methods:

    - nominal_value, std_dev: value at the origin / nominal value, and
      standard deviation.  The standard deviation can be NaN or infinity.

    - n, s: abbreviations for nominal_value and std_dev.

    - error_components(): error_components()[x] is the error due to
      Variable x.

    - derivatives: derivatives[x] is the (value of the) derivative
      with respect to Variable x.  This attribute is a Derivatives
      dictionary whose keys are the Variable objects on which the
      function depends. The values are the numerical values of the
      derivatives.

      All the Variable objects on which the function depends are in
      'derivatives'.

    - std_score(x): position of number x with respect to the
      nominal value, in units of the standard deviation.
    """

    # To save memory in large arrays:
    __slots__ = ('_nominal_value', '_linear_part')

    # !! Fix for mean() in NumPy 1.8.0:
    class dtype(object):
        type = staticmethod(lambda value: value)

    #! The code could be modified in order to accommodate for non-float
    # nominal values.  This could for instance be done through
    # the operator module: instead of delegating operations to
    # float.__*__ operations, they could be delegated to
    # operator.__*__ functions (while taking care of properly handling
    # reverse operations: __radd__, etc.).

    def __init__(self, nominal_value, linear_part):
        """
        nominal_value -- value of the function when the linear part is
        zero.

        linear_part -- LinearCombination that describes the linear
        part of the AffineScalarFunc.
        """

        # ! A technical consistency requirement is that the
        # linear_part can be nested inside a NestedLinearCombination
        # (because this is how functions on AffineScalarFunc calculate
        # their result: by constructing nested expressions for them).

        # Defines the value at the origin:

        # Only float-like values are handled.  One reason is that it
        # does not make sense for a scalar function to be affine to
        # not yield float values.  Another reason is that it would not
        # make sense to have a complex nominal value, here (it would
        # not be handled correctly at all): converting to float should
        # be possible.

        self._nominal_value = float(nominal_value)

        # In order to have a linear execution time for long sums, the
        # _linear_part is generally left as is (otherwise, each
        # successive term would expand to a linearly growing sum of
        # terms: efficiently handling such terms [so, without copies]
        # is not obvious, when the algorithm should work for all
        # functions beyond sums).
        self._linear_part = linear_part

    # The following prevents the 'nominal_value' attribute from being
    # modified by the user:
    @property
    def nominal_value(self):
        "Nominal value of the random number."
        return self._nominal_value

    # Abbreviation (for formulas, etc.):
    n = nominal_value

    ############################################################

    # Making derivatives a property gives the user a clean syntax,
    # which is consistent with derivatives becoming a dictionary.
    @property
    def derivatives(self):
        """
        Return a mapping from each Variable object on which the function
        (self) depends to the value of the derivative with respect to
        that variable.

        This mapping should not be modified.

        Derivative values are always floats.

        This mapping is cached, for subsequent calls.
        """

        if not self._linear_part.expanded():
            self._linear_part.expand()
            # Attempts to get the contribution of a variable that the
            # function does not depend on raise a KeyError:
            self._linear_part.linear_combo.default_factory = None

        return self._linear_part.linear_combo

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

    def __bool__(self):
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
    # taken care of when force_aff_func_args(eq_on_aff_funcs)
    # returns NotImplemented.
    __eq__ = force_aff_func_args(eq_on_aff_funcs)

    __ne__ = force_aff_func_args(ne_on_aff_funcs)
    __gt__ = force_aff_func_args(gt_on_aff_funcs)

    # __ge__ is not the opposite of __lt__ because these operators do
    # not always yield a boolean (for instance, 0 <= numpy.arange(10)
    # yields an array).
    __ge__ = force_aff_func_args(ge_on_aff_funcs)

    __lt__ = force_aff_func_args(lt_on_aff_funcs)
    __le__ = force_aff_func_args(le_on_aff_funcs)

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

        for (variable, derivative) in self.derivatives.items():

            # print "TYPE", type(variable), type(derivative)

            # Individual standard error due to variable:

            # 0 is returned even for a NaN derivative (in this case no
            # multiplication by the derivative is performed): an exact
            # variable obviously leads to no uncertainty in the
            # functions that depend on it.
            if variable._std_dev == 0:
                # !!! Shouldn't the errors always be floats, as a
                # convention of this module?
                error_components[variable] = 0
            else:
                error_components[variable] = abs(derivative*variable._std_dev)

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
            delta**2 for delta in self.error_components().values())))

    # Abbreviation (for formulas, etc.):
    s = std_dev

    def __repr__(self):
        # Not putting spaces around "+/-" helps with arrays of
        # Variable, as each value with an uncertainty is a
        # block of signs (otherwise, the standard deviation can be
        # mistaken for another element of the array).

        std_dev = self.std_dev  # Optimization, since std_dev is calculated

        # A zero standard deviation is printed because otherwise,
        # ufloat_fromstr() does not correctly parse back the value
        # ("1.23" is interpreted as "1.23(1)"):

        if std_dev:
            std_dev_str = repr(std_dev)
        else:
            std_dev_str = '0'

        return "%r+/-%s" % (self.nominal_value, std_dev_str)

    def __str__(self):
        # An empty format string and str() usually return the same
        # string
        # (http://docs.python.org/2/library/string.html#format-specification-mini-language):
        return self.format('')

    def __format__(self, format_spec):
        '''
        Formats a number with uncertainty.

        The format specification are the same as for format() for
        floats, as defined for Python 2.6+ (restricted to what the %
        operator accepts, if using an earlier version of Python),
        except that the n presentation type is not supported. In
        particular, the usual precision, alignment, sign flag,
        etc. can be used. The behavior of the various presentation
        types (e, f, g, none, etc.) is similar. Moreover, the format
        is extended: the number of digits of the uncertainty can be
        controlled, as is the way the uncertainty is indicated (with
        +/- or with the short-hand notation 3.14(1), in LaTeX or with
        a simple text string,...).

        Beyond the use of options at the end of the format
        specification, the main difference with floats is that a "u"
        just before the presentation type (f, e, g, none, etc.)
        activates the "uncertainty control" mode (e.g.: ".6u").  This
        mode is also activated when not using any explicit precision
        (e.g.: "g", "10f", "+010,e" format specifications).  If the
        uncertainty does not have a meaningful number of significant
        digits (0 and NaN uncertainties), this mode is automatically
        deactivated.

        The nominal value and the uncertainty always use the same
        precision. This implies trailing zeros, in general, even with
        the g format type (contrary to the float case). However, when
        the number of significant digits of the uncertainty is not
        defined (zero or NaN uncertainty), it has no precision, so
        there is no matching. In this case, the original format
        specification is used for the nominal value (any "u" is
        ignored).

        Any precision (".p", where p is a number) is interpreted (if
        meaningful), in the uncertainty control mode, as indicating
        the number p of significant digits of the displayed
        uncertainty. Example: .1uf will return a string with one
        significant digit in the uncertainty (and no exponent).

        If no precision is given, the rounding rules from the
        Particle Data Group are used, if possible
        (http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf). For
        example, the "f" format specification generally does not use
        the default 6 digits after the decimal point, but applies the
        PDG rules.

        A common exponent is used if an exponent is needed for the
        larger of the nominal value (in absolute value) and the
        standard deviation, unless this would result in a zero
        uncertainty being represented as 0e... or a NaN uncertainty as
        NaNe.... Thanks to this common exponent, the quantity that
        best describes the associated probability distribution has a
        mantissa in the usual 1-10 range. The common exponent is
        factored (as in "(1.2+/-0.1)e-5"). unless the format
        specification contains an explicit width (" 1.2e-5+/- 0.1e-5")
        (this allows numbers to be in a single column, when printing
        numbers over many lines). Specifying a minimum width of 1 is a
        way of forcing any common exponent to not be factored out.

        The fill, align, zero and width parameters of the format
        specification are applied individually to each of the nominal
        value and standard deviation or, if the shorthand notation is
        used, globally.

        The sign parameter of the format specification is only applied
        to the nominal value (since the standard deviation is
        positive).

        In the case of a non-LaTeX output, the returned string can
        normally be parsed back with ufloat_fromstr(). This however
        excludes cases where numbers use the "," thousands separator,
        for example.

        Options can be added, at the end of the format
        specification. Multiple options can be specified:

        - When "P" is present, the pretty-printing mode is activated: "±"
          separates the nominal value from the standard deviation, exponents
          use superscript characters, etc.
        - When "S" is present (like in .1uS), the short-hand notation 1.234(5)
          is used, indicating an uncertainty on the last digits; if the digits
          of the uncertainty straddle the decimal point, it uses a fixed-point
          notation, like in 12.3(4.5).
        - When "L" is present, the output is formatted with LaTeX.
        - "p" ensures that there are parentheses around the …±… part (no
          parentheses are added if some are already present, for instance
          because of an exponent or of a trailing % sign, etc.). This produces
          outputs like (1.0±0.2) or (1.0±0.2)e7, which can be useful for
          removing any ambiguity if physical units are added after the printed
          number.
    
        An uncertainty which is exactly zero is represented as the
        integer 0 (i.e. with no decimal point).

        The "%" format type forces the percent sign to be at the end
        of the returned string (it is not attached to each of the
        nominal value and the standard deviation).

        Some details of the formatting can be customized as described
        in format_num().
        '''

        # Convention on limits "between" digits: 0 = exactly at the
        # decimal point, -1 = after the first decimal, 1 = before the
        # units digit, etc.

        # Convention on digits: 0 is units (10**0), 1 is tens, -1 is
        # tenths, etc.

        # This method does the format specification parsing, and
        # calculates the various parts of the displayed value
        # (mantissas, exponent, position of the last digit). The
        # formatting itself is delegated to format_num().

        ########################################

        # Format specification parsing:

        match = re.match(r'''
            (?P<fill>[^{}]??)(?P<align>[<>=^]?)  # fill cannot be { or }
            (?P<sign>[-+ ]?)
            (?P<zero>0?)
            (?P<width>\d*)
            (?P<comma>,?)
            (?:\.(?P<prec>\d+))?
            (?P<uncert_prec>u?)  # Precision for the uncertainty?
            # The type can be omitted. Options must not go here:
            (?P<type>[eEfFgG%]??)  # n not supported
            (?P<options>[PSLp]*)  # uncertainties-specific flags
            $''',
            format_spec,
            re.VERBOSE)

        # Does the format specification look correct?
        if not match:
            raise ValueError(
                'Format specification %r cannot be used with object of type'
                ' %r. Note that uncertainties-specific flags must be put at'
                ' the end of the format string.'
                # Sub-classes handled:
                % (format_spec, self.__class__.__name__))

        # Effective format presentation type: f, e, g, etc., or None,
        # like in
        # https://docs.python.org/3.4/library/string.html#format-specification-mini-language. Contrary
        # to what is written in the documentation, it is not true that
        # None is "the same as 'g'": "{}".format() and "{:g}" do not
        # give the same result, on 31415000000.0. None is thus kept as
        # is instead of being replaced by "g".
        pres_type = match.group('type') or None

        # Shortcut:
        fmt_prec = match.group('prec')  # Can be None

        ########################################

        # Since the '%' (percentage) format specification can change
        # the value to be displayed, this value must first be
        # calculated. Calculating the standard deviation is also an
        # optimization: the standard deviation is generally
        # calculated: it is calculated only once, here:
        nom_val = self.nominal_value
        std_dev = self.std_dev

        # 'options' is the options that must be given to format_num():
        options = set(match.group('options'))

        ########################################

        # The '%' format is treated internally as a display option: it
        # should not be applied individually to each part:
        if pres_type == '%':
            # Because '%' does 0.0055*100, the value
            # 0.5499999999999999 is obtained, which rounds to 0.5. The
            # original rounded value is 0.006. The same behavior is
            # found in Python 2.7: '{:.1%}'.format(0.0055) is '0.5%'.
            # If a different behavior is needed, a solution to this
            # problem would be to do the rounding before the
            # multiplication.
            std_dev *= 100
            nom_val *= 100
            pres_type = 'f'
            options.add('%')

        # At this point, pres_type is in eEfFgG or None (not %).

        ########################################

        # Non-real values (nominal value or standard deviation) must
        # be handled in a specific way:
        real_values = [value for value in [abs(nom_val), std_dev]
                       if not isinfinite(value)]

        # Calculation of digits_limit, which defines the precision of
        # the nominal value and of the standard deviation (it can be
        # None when it does not matter, like for NaN±NaN):

        # Reference value for the calculation of a possible exponent,
        # if needed:
        if pres_type in (None, 'e', 'E', 'g', 'G'):
            # Reference value for the exponent: the largest value
            # defines what the exponent will be (another convention
            # could have been chosen, like using the exponent of the
            # nominal value, irrespective of the standard deviation):
            try:
                exp_ref_value = max(real_values)
            except ValueError:  # No non-NaN value: NaN±NaN…
                # No meaningful common exponent can be obtained:
                pass
            ## else:
            ##     print "EXP_REF_VAL", exp_ref_value

        # Should the precision be interpreted like for a float, or
        # should the number of significant digits on the uncertainty
        # be controlled?
        if ((
            # Default behavior: number of significant digits on the
            # uncertainty controlled (if useful, i.e. only in
            # situations where the nominal value and the standard
            # error digits are truncated at the same place):
            (not fmt_prec and len(real_values)==2)
             or match.group('uncert_prec'))  # Explicit control
            # The number of significant digits of the uncertainty must
            # be meaningful, otherwise the position of the significant
            # digits of the uncertainty does not have a clear
            # meaning. This gives us the *effective* uncertainty
            # control mode:
            and std_dev
            and not isinfinite(std_dev)):

            # The number of significant digits on the uncertainty is
            # controlled.

            # The limit digits_limit on the digits of nom_val and std_dev
            # to be displayed is calculated. If the exponent notation is
            # used, this limit is generally different from the finally
            # displayed limit (e.g. 314.15+/-0.01 has digits_limit=-2, but
            # will be displayed with an exponent as (3.1415+/-0.0001)e+02,
            # which corresponds to 4 decimals after the decimal point, not
            # 2).

            # Number of significant digits to use:
            if fmt_prec:
                num_signif_d = int(fmt_prec)  # Can only be non-negative
                if not num_signif_d:
                    raise ValueError("The number of significant digits"
                                     " on the uncertainty should be positive")
            else:
                (num_signif_d, std_dev) = PDG_precision(std_dev)

            digits_limit = signif_dgt_to_limit(std_dev, num_signif_d)

        else:

            # No control of the number of significant digits on the
            # uncertainty.

            ## print "PRECISION NOT BASED ON UNCERTAINTY"

            # The precision has the same meaning as for floats (it is
            # not the uncertainty that defines the number of digits).

            # The usual default precision is used (this is useful for
            # 3.141592±NaN with an "f" format specification, for
            # example):
            #
            # prec is the precision for the main parts of the final
            # format (in the sense of float formatting):
            #
            # https://docs.python.org/3.4/library/string.html#format-specification-mini-language
            if fmt_prec:
                prec = int(fmt_prec)
            elif pres_type is None:
                prec = 12
            else:
                prec = 6

            if pres_type in ('f', 'F'):

                digits_limit = -prec

            else:  # Format type in None, eEgG

                # We first calculate the number of significant digits
                # to be displayed (if possible):

                if pres_type in ('e', 'E'):
                    # The precision is the number of significant
                    # digits required - 1 (because there is a single
                    # digit before the decimal point, which is not
                    # included in the definition of the precision with
                    # the e/E format type):
                    num_signif_digits = prec+1

                else:  # Presentation type in None, g, G

                    # Effective format specification precision: the rule
                    # of
                    # http://docs.python.org/2.7/library/string.html#format-specification-mini-language
                    # is used:

                    # The final number of significant digits to be
                    # displayed is not necessarily obvious: trailing
                    # zeros are removed (with the gG presentation
                    # type), so num_signif_digits is the number of
                    # significant digits if trailing zeros were not
                    # removed. This quantity is relevant for the
                    # rounding implied by the exponent test of the g/G
                    # format:

                    # 0 is interpreted like 1 (as with floats with a
                    # gG presentation type):
                    num_signif_digits = prec or 1

                # The number of significant digits is important for
                # example for determining the exponent:

                ## print "NUM_SIGNIF_DIGITS", num_signif_digits

                digits_limit = (
                    signif_dgt_to_limit(exp_ref_value, num_signif_digits)
                    if real_values
                    else None)

                ## print "DIGITS_LIMIT", digits_limit

        #######################################

        # Common exponent notation: should it be used? use_exp is set
        # accordingly. If a common exponent should be used (use_exp is
        # True), 'common_exp' is set to the exponent that should be
        # used.

        if pres_type in ('f', 'F'):
            use_exp = False
        elif pres_type in ('e', 'E'):
            if not real_values:
                use_exp = False
            else:
                use_exp = True
                # !! This calculation might have been already done,
                # for instance when using the .0e format:
                # signif_dgt_to_limit() was called before, which
                # prompted a similar calculation:
                common_exp = first_digit(round(exp_ref_value, -digits_limit))

        else:  # None, g, G

            # The rules from
            # https://docs.python.org/3.4/library/string.html#format-specification-mini-language
            # are applied.

            # Python's native formatting (whose result could be parsed
            # in order to determine whether a common exponent should
            # be used) is not used because there is shared information
            # between the nominal value and the standard error (same
            # last digit, common exponent) and extracting this
            # information from Python would entail parsing its
            # formatted string, which is in principle inefficient
            # (internally, Python performs calculations that yield a
            # string, and the string would be parsed back into
            # separate parts and numbers, which is in principle
            # unnecessary).

            # Should the scientific notation be used? The same rule as
            # for floats is used ("-4 <= exponent of rounded value <
            # p"), on the nominal value.

            if not real_values:
                use_exp = False
            else:
                # Common exponent *if* used:
                common_exp = first_digit(round(exp_ref_value, -digits_limit))

                # print "COMMON EXP TEST VALUE", common_exp
                # print "LIMIT EXP", common_exp-digits_limit+1
                # print "WITH digits_limit", digits_limit

                # The number of significant digits of the reference value
                # rounded at digits_limit is exponent-digits_limit+1:
                if -4 <= common_exp < common_exp-digits_limit+1:
                    use_exp = False
                else:
                    use_exp = True

        ########################################

        # Calculation of signif_limit (position of the significant
        # digits limit in the final fixed point representations; this
        # is either a non-positive number, or None), of
        # nom_val_mantissa ("mantissa" for the nominal value,
        # i.e. value possibly corrected for a factorized exponent),
        # and std_dev_mantissa (similarly for the standard
        # deviation). common_exp is also set to None if no common
        # exponent should be used.

        if use_exp:

            # Not 10.**(-common_exp), for limit values of common_exp:
            factor = 10.**common_exp

            nom_val_mantissa = nom_val/factor
            std_dev_mantissa = std_dev/factor
            # Limit for the last digit of the mantissas:
            signif_limit = digits_limit - common_exp

        else:  # No common exponent

            common_exp = None

            nom_val_mantissa = nom_val
            std_dev_mantissa = std_dev
            signif_limit = digits_limit

        ## print "SIGNIF_LIMIT", signif_limit

        ########################################

        # Format of the main (i.e. with no exponent) parts (the None
        # presentation type is similar to the g format type):

        main_pres_type = 'fF'[(pres_type or 'g').isupper()]

        # The precision of the main parts must be adjusted so as
        # to take into account the special role of the decimal
        # point:
        if signif_limit is not None:  # If signif_limit is pertinent
            # The decimal point location is always included in the
            # printed digits (e.g., printing 3456 with only 2
            # significant digits requires to print at least four
            # digits, like in 3456 or 3500).
            #
            # The max() is important for example for
            # 1234567.89123+/-12345.678 with the f format: in this
            # case, signif_limit is +3 (2 significant digits necessary
            # for the error, as per the PDG rules), but the (Python
            # float formatting) precision to be used for the main
            # parts is 0 (all digits must be shown).
            #
            # The 1 for the None pres_type represents "at least one
            # digit past the decimal point" of Python
            # (https://docs.python.org/3.4/library/string.html#format-specification-mini-language). This
            # is only applied for null uncertainties.
            prec = max(-signif_limit,
                       1 if pres_type is None and not std_dev
                       else 0)
        ## print "PREC", prec

        ########################################

        # print (
        #     "FORMAT_NUM parameters: nom_val_mantissa={},"
        #     " std_dev_mantissa={}, common_exp={},"
        #     " match.groupdict()={}, prec={}, main_pres_type={},"
        #     " options={}".format(
        #     nom_val_mantissa, std_dev_mantissa, common_exp,
        #     match.groupdict(),
        #     prec,
        #     main_pres_type,
        #     options))

        # Final formatting:
        return format_num(nom_val_mantissa, std_dev_mantissa, common_exp,
                          match.groupdict(),
                          prec=prec,
                          main_pres_type=main_pres_type,
                          options=options)

    # Alternate name for __format__, for use with Python < 2.6 (and
    # other Python versions if the user so chooses: this helps moving
    # code from Python 2.6 to more recent versions):
    @set_doc("""
        Return the same result as self.__format__(format_spec), or
        equivalently as the format(self, format_spec) of Python 2.6+.

        This method is meant to be used for formatting numbers with
        uncertainties in Python < 2.6, with '... %s ...' %
        num.format('.2e').
        """)
    def format(*args, **kwargs):
        return args[0].__format__(*args[1:], **kwargs)

    def std_score(self, value):
        """
        Return 'value' - nominal value, in units of the standard
        deviation.

        Raises a ValueError exception if the standard deviation is zero.
        """
        try:
            # The ._nominal_value is a float: there is no integer division,
            # here:
            return (value - self._nominal_value) / self.std_dev
        except ZeroDivisionError:
            raise ValueError("The standard deviation is zero:"
                             " undefined result")

    def __deepcopy__(self, memo):
        """
        Hook for the standard copy module.

        The returned AffineScalarFunc is a completely fresh copy,
        which is fully independent of any variable defined so far.
        New variables are specially created for the returned
        AffineScalarFunc object.
        """
        return AffineScalarFunc(self._nominal_value,
                                copy.deepcopy(self._linear_part))

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
            if isinstance(slot_names, basestring):
                all_slots.add(slot_names)  # Single name
            else:
                all_slots.update(slot_names)

        # The slot values are stored:
        for name in all_slots:
            try:
                # !! It might happen that '__dict__' is itself a slot
                # name. In this case, its value is saved
                # again. Alternatively, the loop could be done on
                # all_slots - {'__dict__'}:
                all_attrs[name] = getattr(self, name)
            except AttributeError:
                pass  # Undefined slot attribute

        return all_attrs

    def __setstate__(self, data_dict):
        """
        Hook for the pickle module.
        """
        for (name, value) in data_dict.items():
            # Contrary to the default __setstate__(), this does not
            # necessarily save to the instance dictionary (because the
            # instance might contain slots):
            setattr(self, name, value)

# Nicer name, for users: isinstance(ufloat(...), UFloat) is
# True. Also: isinstance(..., UFloat) is the test for "is this a
# number with uncertainties from the uncertainties package?":
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
    Return operators with a reflection, along with their partial derivatives.

    Operators are things like +, /, etc. Those considered here have two
    arguments and can be called through Python's reflected methods __r…__ (e.g.
    __radd__).

    See the code for details.
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
        'sub': ("1.", "-1."),
        'truediv': ("1/y", "-x/y**2")
        }

    # Conversion to Python functions:
    ops_with_reflection = {}
    for (op, derivatives) in derivatives_list.items():
        ops_with_reflection[op] = [
            eval("lambda x, y: %s" % expr) for expr in derivatives ]

        ops_with_reflection["r"+op] = [
            eval("lambda y, x: %s" % expr) for expr in reversed(derivatives)]


    # The derivatives of pow() are more complicated:

    # The case x**y is constant one the line x = 0 and in y = 0;
    # the corresponding derivatives must be zero in these
    # cases. If the function is actually not defined (e.g. 0**-3),
    # then an exception will be raised when the nominal value is
    # calculated.  These derivatives are transformed to NaN if an
    # error happens during their calculation:

    def pow_deriv_0(x, y):
        if y == 0:
            return 0.
        elif x != 0 or y % 1 == 0:
            return y*x**(y-1)
        else:
            return float('nan')

    def pow_deriv_1(x, y):
        if x == 0 and y > 0:
            return 0.
        else:
            return log(x)*x**y

    ops_with_reflection['pow'] = [pow_deriv_0, pow_deriv_1]
    ops_with_reflection['rpow'] = [lambda y, x: pow_deriv_1(x, y),
                                   lambda y, x: pow_deriv_0(x, y)]

    # Undefined derivatives are converted to NaN when the function
    # itself can be calculated:
    for op in ['pow']:
        ops_with_reflection[op] = [
            nan_if_exception(func) for func in ops_with_reflection[op]]
        ops_with_reflection['r'+op] = [
            nan_if_exception(func) for func in ops_with_reflection['r'+op]]

    return ops_with_reflection

# Operators that have a reflection, along with their derivatives:
ops_with_reflection = get_ops_with_reflection()

# Some effectively modified operators (for the automated tests):
modified_operators = []
modified_ops_with_reflection = []

# Custom versions of some operators (instead of extending some float
# __*__ operators to AffineScalarFunc, the operators in custom_ops
# are used):
if sys.version_info < (3,):

    custom_ops = {}

else:

    # !!! This code is not run by the tests. It would be nice to have
    # it be tested.
    def no_complex_result(func):
        '''
        Return a function that does like func, but that raises a
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
    custom_ops = {
        'pow': no_complex_result(float.__pow__),
        'rpow': no_complex_result(float.__rpow__)
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

    def _simple_add_deriv(x):
        if x >= 0:
            return 1.
        else:
            return -1.

    # Single-argument operators that should be adapted from floats to
    # AffineScalarFunc objects, associated to their derivative:
    simple_numerical_operators_derivatives = {
        'abs': _simple_add_deriv,
        'neg': lambda x: -1.,
        'pos': lambda x: 1.,
        'trunc': lambda x: 0.
        }

    for (op, derivative) in (
        iter(simple_numerical_operators_derivatives.items())):

        attribute_name = "__%s__" % op

        # float objects don't exactly have the same attributes between
        # different versions of Python (for instance, __trunc__ was
        # introduced with Python 2.6):
        try:
            setattr(AffineScalarFunc, attribute_name,
                    wrap(getattr(float, attribute_name), [derivative]))
        except AttributeError:
            # Version of Python where floats don't have attribute_name:
            pass
        else:
            modified_operators.append(op)

    ########################################
    # Final definition of the operators for AffineScalarFunc objects:

    # Reversed versions (useful for float*AffineScalarFunc, for instance):
    for (op, derivatives) in ops_with_reflection.items():
        attribute_name = '__%s__' % op

        # float objects don't exactly have the same attributes between
        # different versions of Python (for instance, __div__ and
        # __rdiv__ were removed, in Python 3):

        # float objects don't exactly have the same attributes between
        # different versions of Python (for instance, __trunc__ was
        # introduced with Python 2.6):
        try:
            if op not in custom_ops:
                func_to_wrap = getattr(float, attribute_name)
            else:
                func_to_wrap = custom_ops[op]
        except AttributeError:
            # Version of Python with floats that don't have attribute_name:
            pass
        else:
            setattr(AffineScalarFunc, attribute_name,
                    wrap(func_to_wrap, derivatives))
            modified_ops_with_reflection.append(op)

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

class NegativeStdDev(Exception):
    '''Raise for a negative standard deviation'''
    pass

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
        are set.

        The value is converted to float.

        The standard deviation std_dev can be NaN. It should normally
        be a float or an integer.

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
        super(Variable, self).__init__(value, LinearCombination({self: 1.}))

        self.std_dev = std_dev  # Assignment through a Python property

        self.tag = tag

    @property
    def std_dev(self):
        return self._std_dev

    # Standard deviations can be modified (this is a feature).
    # AffineScalarFunc objects that depend on the Variable have their
    # std_dev automatically modified (recalculated with the new
    # std_dev of their Variables):
    @std_dev.setter
    def std_dev(self, std_dev):

        # We force the error to be float-like.  Since it is considered
        # as a standard deviation, it must be either positive or NaN:
        # (Note: if NaN < 0 is False, there is no need to test
        # separately for NaN. But this is not guaranteed, even if it
        # should work on most platforms.)
        if std_dev < 0 and not isinfinite(std_dev):
            raise NegativeStdDev("The standard deviation cannot be negative")

        self._std_dev = CallableStdDev(std_dev)

    # Support for legacy method:
    def set_std_dev(self, value):  # Obsolete
        deprecation('instead of set_std_dev(), please use'
                    ' .std_dev = ...')
        self.std_dev = value

    # The following method is overridden so that we can represent the tag:
    def __repr__(self):

        num_repr  = super(Variable, self).__repr__()

        if self.tag is None:
            return num_repr
        else:
            return "< %s = %s >" % (self.tag, num_repr)

    def __hash__(self):
        # All Variable objects are by definition independent
        # variables, so they never compare equal; therefore, their
        # id() are allowed to differ
        # (http://docs.python.org/reference/datamodel.html#object.__hash__):
        return id(self)

    def __copy__(self):
        """
        Hook for the standard copy module.
        """

        # !!!!!! The comment below might not be valid anymore now that
        # Variables do not contain derivatives anymore.

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
    Return the nominal value of x if it is a quantity with
    uncertainty (i.e., an AffineScalarFunc object); otherwise, returns
    x unchanged.

    This utility function is useful for transforming a series of
    numbers, when only some of them generally carry an uncertainty.
    """

    if isinstance(x, AffineScalarFunc):
        return x.nominal_value
    else:
        return x

def std_dev(x):
    """
    Return the standard deviation of x if it is a quantity with
    uncertainty (i.e., an AffineScalarFunc object); otherwise, returns
    the float 0.

    This utility function is useful for transforming a series of
    numbers, when only some of them generally carry an uncertainty.
    """

    if isinstance(x, AffineScalarFunc):
        return x.std_dev
    else:
        return 0.
def ufloat_fromstr(representation, tag=None):
    """
    Return a new random variable (Variable object) from a string.

    Strings 'representation' of the form '12.345+/-0.015',
    '12.345(15)', '12.3' or u'1.2±0.1' (Unicode string) are recognized
    (see more complete list below).  In the last case, an uncertainty
    of +/-1 is assigned to the last digit.

    Invalid representations raise a ValueError.

    This function tries to parse back most of the formats that are made
    available by this module. Examples of valid string representations:

        12.3e10+/-5e3
        (-3.1415 +/- 0.0001)e+02  # Factored exponent

        # Pretty-print notation (only with a unicode string):
        12.3e10 ± 5e3  # ± symbol
        (12.3 ± 5.0) × 10⁻¹²  # Times symbol, superscript
        12.3 ± 5e3  # Mixed notation (± symbol, but e exponent)

        # Double-exponent values:
        (-3.1415 +/- 1e-4)e+200
        (1e-20 +/- 3)e100

        0.29
        31.
        -31.
        31
        -3.1e10

        -1.23(3.4)
        -1.34(5)
        1(6)
        3(4.2)
        -9(2)
        1234567(1.2)
        12.345(15)
        -12.3456(78)e-6
        12.3(0.4)e-5
        169.0(7)
        169.1(15)
        .123(4)
        .1(.4)

        # NaN uncertainties:
        12.3(nan)
        12.3(NAN)
        3±nan

    Surrounding spaces are ignored.

    About the "shorthand" notation: 1.23(3) = 1.23 ± 0.03 but
    1.23(3.) = 1.23 ± 3.00. Thus, the presence of a decimal point in
    the uncertainty signals an absolute uncertainty (instead of an
    uncertainty on the last digits of the nominal value).
    """

    (nominal_value, std_dev) = str_to_number_with_uncert(
        representation.strip())

    return ufloat(nominal_value, std_dev, tag)

def ufloat_obsolete(representation, tag=None):
    '''
    Legacy version of ufloat(). Will eventually be removed.

    representation -- either a (nominal_value, std_dev) tuple, or a
    string representation of a number with uncertainty, in a format
    recognized by ufloat_fromstr().
    '''

    if isinstance(representation, tuple):
        return ufloat(representation[0], representation[1], tag)
    else:
        return ufloat_fromstr(representation, tag)

# The arguments are named for the new version, instead of bearing
# names that are closer to their obsolete use (e.g., std_dev could be
# instead std_dev_or_tag, since it can be the tag, in the obsolete
# ufloat((3, 0.14), "pi") form). This has the advantage of allowing
# new code to use keyword arguments as in ufloat(nominal_value=3,
# std_dev=0.14), without breaking when the obsolete form is not
# supported anymore.
def ufloat(nominal_value, std_dev=None, tag=None):
    """
    Return a new random variable (Variable object).

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

    nominal_value -- nominal value of the random variable. It is more
    meaningful to use a value close to the central value or to the
    mean. This value is propagated by mathematical operations as if it
    was a float.

    std_dev -- standard deviation of the random variable. The standard
    deviation must be convertible to a positive float, or be NaN.

    tag -- optional string tag for the variable.  Variables don't have
    to have distinct tags.  Tags are useful for tracing what values
    (and errors) enter in a given result (through the
    error_components() method).
    """

    try:
        # Standard case:
        return Variable(nominal_value, std_dev, tag=tag)
    # Exception types raised by, respectively: tuple or string that
    # can be converted through float() (case of a number with no
    # uncertainty), and string that cannot be converted through
    # float():
    except (TypeError, ValueError):

        if tag is not None:
            tag_arg = tag  # tag keyword used:
        else:
            tag_arg = std_dev  # 2 positional arguments form

        try:
            final_ufloat = ufloat_obsolete(nominal_value, tag_arg)
        except:  # The input is incorrect, not obsolete
            raise
        else:
            # Obsolete, two-argument call:
            deprecation(
                'either use ufloat(nominal_value, std_dev),'
                ' ufloat(nominal_value, std_dev, tag), or the'
                ' ufloat_fromstr() function, for string representations.')
            return final_ufloat
