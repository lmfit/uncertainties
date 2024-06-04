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
if sys.version_info < (3,):
     from past.builtins import basestring
else:
     # Avoid importing from past in Python 3 since it utilizes the builtin
     # 'imp' module, which is deprecated as of Python 3.4, see
     # https://docs.python.org/3/library/imp.html. The 2to3 tool replaces
     # basestring with str, so that's what we effectively do here as well:
     basestring = str

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

from uncertainties.formatting import format_ufloat, nrmlze_superscript

# The following restricts the local function getargspec() to the common
# features of inspect.getargspec() and inspect.getfullargspec():
if sys.version_info < (3,):  # !! Could be removed when moving to Python 3 only
    from inspect import getargspec
else:
    from inspect import getfullargspec as getargspec

# Attributes that are always exported (some other attributes are
# exported only if the NumPy module is available...):
__all__ = [

    # All sub-modules and packages are not imported by default,
    # in particular because NumPy might be unavailable.

    'ufloat',  # Main function: returns a number with uncertainty
    'ufloat_fromstr',  # Important function: returns a number with uncertainty

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
    'Variable',

    # Wrapper for allowing non-pure-Python function to handle
    # quantitities with uncertainties:
    'wrap'

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
FLOAT_LIKE_TYPES = (numbers.Number,)
CONSTANT_TYPES = FLOAT_LIKE_TYPES+(complex,)

###############################################################################
## Definitions that depend on the availability of NumPy:


try:
    import numpy
except ImportError:
    pass
else:

    # NumPy numbers do not depend on Variable objects:
    FLOAT_LIKE_TYPES += (numpy.generic,)
    CONSTANT_TYPES += FLOAT_LIKE_TYPES[-1:]

    # Entering variables as a block of correlated values.  Only available
    # if NumPy is installed.

    #! It would be possible to dispense with NumPy, but a routine should be
    # written for obtaining the eigenvectors of a symmetric matrix.  See
    # for instance Numerical Recipes: (1) reduction to tri-diagonal
    # [Givens or Householder]; (2) QR / QL decomposition.

    def correlated_values(nom_values, covariance_mat, tags=None):
        """
        Return numbers with uncertainties (AffineScalarFunc objects)
        that correctly reproduce the given covariance matrix, and have
        the given (float) values as their nominal value.

        The correlated_values_norm() function returns the same result,
        but takes a correlation matrix instead of a covariance matrix.

        The list of values and the covariance matrix must have the
        same length, and the matrix must be a square (symmetric) one.

        The numbers with uncertainties returned depend on newly
        created, independent variables (Variable objects).

        nom_values -- sequence with the nominal (real) values of the
        numbers with uncertainties to be returned.

        covariance_mat -- full covariance matrix of the returned numbers with
        uncertainties. For example, the first element of this matrix is the
        variance of the first number with uncertainty. This matrix must be a
        NumPy array-like (list of lists, NumPy array, etc.).

        tags -- if 'tags' is not None, it must list the tag of each new
        independent variable.
        """

        # !!! It would in principle be possible to handle 0 variance
        # variables by first selecting the sub-matrix that does not contain
        # such variables (with the help of numpy.ix_()), and creating
        # them separately.

        std_devs = numpy.sqrt(numpy.diag(covariance_mat))

        # For numerical stability reasons, we go through the correlation
        # matrix, because it is insensitive to any change of scale in the
        # quantities returned. However, care must be taken with 0 variance
        # variables: calculating the correlation matrix cannot be simply done
        # by dividing by standard deviations. We thus use specific
        # normalization values, with no null value:
        norm_vector = std_devs.copy()
        norm_vector[norm_vector==0] = 1

        return correlated_values_norm(
            # !! The following zip() is a bit suboptimal: correlated_values()
            # separates back the nominal values and the standard deviations:
            list(zip(nom_values, std_devs)),
            covariance_mat/norm_vector/norm_vector[:,numpy.newaxis],
            tags)

    __all__.append('correlated_values')

    def correlated_values_norm(values_with_std_dev, correlation_mat,
                               tags=None):
        '''
        Return correlated values like correlated_values(), but takes
        instead as input:

        - nominal (float) values along with their standard deviation, and
        - a correlation matrix (i.e. a normalized covariance matrix).

        values_with_std_dev -- sequence of (nominal value, standard
        deviation) pairs. The returned, correlated values have these
        nominal values and standard deviations.

        correlation_mat -- correlation matrix between the given values, except
        that any value with a 0 standard deviation must have its correlations
        set to 0, with a diagonal element set to an arbitrary value (something
        close to 0-1 is recommended, for a better numerical precision).  When
        no value has a 0 variance, this is the covariance matrix normalized by
        standard deviations, and thus a symmetric matrix with ones on its
        diagonal.  This matrix must be an NumPy array-like (list of lists,
        NumPy array, etc.).

        tags -- like for correlated_values().
        '''

        # If no tags were given, we prepare tags for the newly created
        # variables:
        if tags is None:
            tags = (None,) * len(values_with_std_dev)

        (nominal_values, std_devs) = numpy.transpose(values_with_std_dev)

        # We diagonalize the correlation matrix instead of the
        # covariance matrix, because this is generally more stable
        # numerically. In fact, the covariance matrix can have
        # coefficients with arbitrary values, through changes of units
        # of its input variables. This creates numerical instabilities.
        #
        # The covariance matrix is diagonalized in order to define
        # the independent variables that model the given values:
        (variances, transform) = numpy.linalg.eigh(correlation_mat)

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

        # The coordinates of each new uncertainty as a function of the
        # new variables must include the variable scale (standard deviation):
        transform *= std_devs[:, numpy.newaxis]

        # Representation of the initial correlated values:
        values_funcs = tuple(
            AffineScalarFunc(
                value,
                LinearCombination(dict(zip(variables, coords))))
            for (coords, value) in zip(transform, nominal_values))

        return values_funcs

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

# Step constant for numerical derivatives in
# partial_derivative(). Value chosen to as to get better numerical
# results:
STEP_SIZE = sqrt(sys.float_info.epsilon)

# !! It would be possible to split the partial derivative calculation
# into two functions: one for positional arguments (case of integer
# arg_ref) and one for keyword arguments (case of string
# arg_ref). However, this would either duplicate the code for the
# numerical differentiation, or require a call, which is probably more
# expensive in time than the tests done here.
def partial_derivative(f, arg_ref):
    """
    Return a function that numerically calculates the partial
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
        if change_kwargs:
            args_with_var = kwargs
        else:
            args_with_var = list(args)

        # The step is relative to the parameter being varied, so that
        # shifting it does not suffer from finite precision limitations:
        step = STEP_SIZE*abs(args_with_var[arg_ref])
        if not step:
            # Arbitrary, but "small" with respect to 1:
            step = STEP_SIZE

        args_with_var[arg_ref] += step

        if change_kwargs:
            shifted_f_plus = f(*args, **args_with_var)
        else:
            shifted_f_plus = f(*args_with_var, **kwargs)

        args_with_var[arg_ref] -= 2*step  # Optimization: only 1 list copy

        if change_kwargs:
            shifted_f_minus = f(*args, **args_with_var)
        else:
            shifted_f_minus = f(*args_with_var, **kwargs)

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
        Return the n-th numerical derivative of the function.
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

                value = next(self.iterable)

                if value is None:
                    value = self.none_converter(pos)

                returned_elements.append(value)

            return returned_elements[index]

    def __str__(self):
        return '<%s: [%s...]>' % (
            self.__class__.__name__,
            ', '.join(map(str, self.returned_elements)))

def wrap(f, derivatives_args=None, derivatives_kwargs=None):
    """Wrap a function f into one that accepts Variables.

    The function f must return a float or integer value.  The returned
    wrapped function will return values with both uncertainties and
    correlations, but can be used as a drop-in replacement for the
    original function.

    Arguments:
    ----------
    derivatives_args: list or iterable
           list or tupleof derivative functionss or None with respect to
           the positional arguments of `f`.  See Note 1.
    derivatives_kwargs: dictionary
           dict of derivative functionss or None with respect to the
           keyword arguments of `f`.  See Note 1.

    Notes:
    -------
    1.  Each function must be the partial derivative of f with respect to the
        corresponding positional parameters, and must have the same signature
        as ``f``. `derivative_args` hold derivitative functions for positional
        arguments (include `*varargs`), while  `derivative_kwargs` holds
        derivitative functions for keyword arguments (include `**kwargs`). If an
        entry is `None` or not supplied, and if the argument value isa numeric
        Variable, a numerical derivative will be used. Non-numeric are ignored.
    2.  If derivatives are meaningless or the function is not function is not
        differentiable, the derivative funcion should return NaN for values
        for which the the function is not differentiable.

    Example:
    --------
    To wrap `sin`, one could do
       >>> from uncertainties import wrap, umath
       >>> import math
       >>> usin_a = wrap(math.sin)   # uses numerical derivative
       >>> usin_b = wrap(math.sin, [math.cos])  # use analytic derivative
       >>> usin_c = umath.sin        # builtin, same as usin_2

    These will all give the same result.
    """
    if derivatives_args is None:
        derivatives_args = []
    if derivatives_kwargs is None:
        derivatives_kwargs = {}
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

########################################
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
        return float(sqrt(sum(
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
        return format_ufloat(self, format_spec)
    
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
    Representation of a float-like scalar Variable with its uncertainty.

    Variables are independent from each other, but correlations between them
    are handled through the AffineScalarFunc class.
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

        self._std_dev = float(std_dev)

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

def covariance_matrix(nums_with_uncert):
    """
    Return a matrix that contains the covariances between the given
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
    for (i1, expr1) in enumerate(nums_with_uncert, 1):
        derivatives1 = expr1.derivatives  # Optimization
        vars1 = set(derivatives1)  # !! Python 2.7+: viewkeys() would work
        coefs_expr1 = []

        for expr2 in nums_with_uncert[:i1]:
            derivatives2 = expr2.derivatives  # Optimization
            coefs_expr1.append(sum(
                ((derivatives1[var]*derivatives2[var]*var._std_dev**2)
                # var is a variable common to both numbers with
                # uncertainties:
                for var in vars1.intersection(derivatives2)),
                # The result is always a float (sum() with no terms
                # returns an integer):
                0.))

        covariance_matrix.append(coefs_expr1)

    # We symmetrize the matrix:
    for (i, covariance_coefs) in enumerate(covariance_matrix):
        covariance_coefs.extend([covariance_matrix[j][i]
                                 for j in range(i+1, len(covariance_matrix))])

    return covariance_matrix

try:
    import numpy
except ImportError:
    pass
else:
    def correlation_matrix(nums_with_uncert):
        '''
        Return the correlation matrix of the given sequence of
        numbers with uncertainties, as a NumPy array of floats.
        '''

        cov_mat = numpy.array(covariance_matrix(nums_with_uncert))

        std_devs = numpy.sqrt(cov_mat.diagonal())

        return cov_mat/std_devs/std_devs[numpy.newaxis].T

    __all__.append('correlation_matrix')

###############################################################################
# Parsing of values with uncertainties:

# Parsing of (part of) numbers. The reason why the decimal part is
# parsed (if any), instead of using the parsing built in float(), is
# that the presence (or not) of a decimal point does matter, in the
# semantics of some representations (e.g. .1(2.) = .1+/-2, whereas
# .1(2) = .1+/-0.2), so just getting the numerical value of the part
# in parentheses would not be sufficient.
POSITIVE_DECIMAL_UNSIGNED_OR_NON_FINITE = r'((\d*)(\.\d*)?|nan|NAN|inf|INF)'

# Regexp for a number with uncertainty (e.g., "-1.234(2)e-6"), where
# the uncertainty is optional (in which case the uncertainty is
# implicit). The uncertainty can also be nan or NAN:
#
# !! WARNING: in Python 2, the code relies on "… % <unicode string>" returning
# a Unicode string (even if the template is not Unicode):
NUMBER_WITH_UNCERT_RE_STR = u'''
    ([+-])?  # Sign
    %s  # Main number
    (?:\\(%s\\))?  # Optional uncertainty
    (?:
        (?:[eE]|\\s*×\\s*10)
        (.*)
    )?  # Optional exponent
    ''' % (POSITIVE_DECIMAL_UNSIGNED_OR_NON_FINITE,
           POSITIVE_DECIMAL_UNSIGNED_OR_NON_FINITE)

NUMBER_WITH_UNCERT_RE_MATCH = re.compile(
    u"%s$" % NUMBER_WITH_UNCERT_RE_STR, re.VERBOSE).match

# Number with uncertainty with a factored exponent (e.g., of the form
# (... +/- ...)e10): this is a loose matching, so as to accommodate
# for multiple formats:
NUMBER_WITH_UNCERT_GLOBAL_EXP_RE_MATCH = re.compile(u'''
    \\(
    (?P<simple_num_with_uncert>.*)
    \\)
    (?:[eE]|\\s*×\\s*10) (?P<exp_value>.*)
    $''', re.VERBOSE).match

class NotParenUncert(ValueError):
    '''
    Raised when a string representing an exact number or a number with
    an uncertainty indicated between parentheses was expected but not
    found.
    '''

def parse_error_in_parentheses(representation):
    # !!!! The code seems to handle superscript exponents, but the
    # docstring doesn't reflect this!?
    """
    Return (value, error) from a string representing a number with
    uncertainty like 12.34(5), 12.34(142), 12.5(3.4), 12.3(4.2)e3, or
    13.4(nan)e10.  If no parenthesis is given, an uncertainty of one
    on the last digit is assumed.

    The digits between parentheses correspond to the same number of digits
    at the end of the nominal value (the decimal point in the uncertainty
    is optional). Example: 12.34(142) = 12.34±1.42.

    Raises ValueError if the string cannot be parsed.
    """

    match = NUMBER_WITH_UNCERT_RE_MATCH(representation)

    if match:
        # The 'main' part is the nominal value, with 'int'eger part, and
        # 'dec'imal part.  The 'uncert'ainty is similarly broken into its
        # integer and decimal parts.
        (sign, main, _, main_dec, uncert, uncert_int, uncert_dec,
         exponent) = match.groups()
    else:
        raise NotParenUncert("Unparsable number representation: '%s'."
                             " See the documentation of ufloat_fromstr()."
                             % representation)

    # Global exponent:
    if exponent:
        factor = 10.**nrmlze_superscript(exponent)
    else:
        factor = 1

    # Nominal value:
    value = float((sign or '')+main)*factor

    if uncert is None:
        # No uncertainty was found: an uncertainty of 1 on the last
        # digit is assumed:
        uncert_int = '1'  # The other parts of the uncertainty are None

    # Do we have a fully explicit uncertainty?
    if uncert_dec is not None or uncert in {'nan', 'NAN', 'inf', 'INF'}:
        uncert_value = float(uncert)
    else:
        # uncert_int represents an uncertainty on the last digits:

        # The number of digits after the period defines the power of
        # 10 that must be applied to the provided uncertainty:
        if main_dec is None:
            num_digits_after_period = 0
        else:
            num_digits_after_period = len(main_dec)-1

        uncert_value = int(uncert_int)/10.**num_digits_after_period

    # We apply the exponent to the uncertainty as well:
    uncert_value *= factor

    return (value, uncert_value)

# Regexp for catching the two variable parts of -1.2×10⁻¹²:
PRETTY_PRINT_MATCH = re.compile(u'(.*?)\\s*×\\s*10(.*)').match

def to_float(value_str):
    '''
    Converts a string representing a float to a float.

    The usual valid Python float() representations are correctly
    parsed.

    In addition, the pretty-print notation -1.2×10⁻¹² is also
    converted.

    ValueError is raised if no float can be obtained.
    '''

    try:
        return float(value_str)
    except ValueError:
        pass

    # The pretty-print notation is tried:
    match = PRETTY_PRINT_MATCH(value_str)
    if match:
        try:
            return float(match.group(1))*10.**nrmlze_superscript(match.group(2))
        except ValueError:
            raise ValueError('Mantissa or exponent incorrect in pretty-print'
                             ' form %s' % value_str)
    else:
        raise ValueError('No valid Python float or pretty-print form'
                         ' recognized in %s' % value_str)


cannot_parse_ufloat_msg_pat = (
    'Cannot parse %s: see the documentation for ufloat_fromstr() for a'
    ' list of accepted formats')

# The following function is not exposed because it can in effect be
# obtained by doing x = ufloat_fromstr(representation) and reading
# x.nominal_value and x.std_dev:
def str_to_number_with_uncert(representation):
    """
    Given a string that represents a number with uncertainty, returns the
    nominal value and the uncertainty.

    See the documentation for ufloat_fromstr() for a list of accepted
    formats.

    When no numerical error is given, an uncertainty of 1 on the last
    digit is implied.

    Raises ValueError if the string cannot be parsed.

    representation -- string with no leading or trailing spaces.
    """

    # The "p" format can add parentheses around the whole printed result: we
    # remove them:
    if representation.startswith('(') and representation.endswith(')'):
        representation = representation[1:-1]

    match = NUMBER_WITH_UNCERT_GLOBAL_EXP_RE_MATCH(representation)

    # The representation is simplified, but the global factor is
    # calculated:

    if match:

        # We have a form with a factored exponent: (1.23 +/- 0.01)e10,
        # etc.

        exp_value_str = match.group('exp_value')

        try:
            exponent = nrmlze_superscript(exp_value_str)
        except ValueError:
            raise ValueError(cannot_parse_ufloat_msg_pat % representation)

        factor = 10.**exponent

        representation = match.group('simple_num_with_uncert')
    else:
        factor = 1  # No global exponential factor

    match = re.match(u'(.*)(?:\\+/-|±)(.*)', representation)
    if match:

        (nom_value, uncert) = match.groups()

        try:
            # Simple form 1234.45+/-1.2 or 1234.45±1.2, or 1.23e-10+/-1e-23
            # or -1.2×10⁻¹²±1e23:
            parsed_value = (to_float(nom_value)*factor,
                            to_float(uncert)*factor)
        except ValueError:
            raise ValueError(cannot_parse_ufloat_msg_pat % representation)

    else:
        # Form with error parentheses or no uncertainty:
        try:
            parsed_value = parse_error_in_parentheses(representation)
        except NotParenUncert:
            raise ValueError(cannot_parse_ufloat_msg_pat % representation)

    return parsed_value

def ufloat_fromstr(representation, tag=None):
    """
    Create an uncertainties Variable from a string representation.
    Several representation formats are supported.

    Arguments:
    ----------
    representation: string
        string representation of a value with uncertainty
    tag:   string or `None`
        optional tag for tracing and organizing Variables ['None']

    Returns:
    --------
    uncertainties Variable.

    Notes:
    --------
    1. Invalid representations raise a ValueError.

    2. Using the form "nominal(std)" where "std" is an integer creates
       a Variable with "std" giving the least significant digit(s).
       That is, "1.25(3)" is the same as `ufloat(1.25, 0.03)`,
       while "1.25(3.)" is the same as `ufloat(1.25, 3.)`

    Examples:
    -----------

    >>> x = ufloat_fromsstr("12.58+/-0.23")  # = ufloat(12.58, 0.23)
    >>> x = ufloat_fromsstr("12.58 ± 0.23")  # = ufloat(12.58, 0.23)
    >>> x = ufloat_fromsstr("3.85e5 +/- 2.3e4")  # = ufloat(3.8e5, 2.3e4)
    >>> x = ufloat_fromsstr("(38.5 +/- 2.3)e4")  # = ufloat(3.8e5, 2.3e4)

    >>> x = ufloat_fromsstr("72.1(2.2)")  # = ufloat(72.1, 2.2)
    >>> x = ufloat_fromsstr("72.15(4)")  # = ufloat(72.15, 0.04)
    >>> x = ufloat_fromstr("680(41)e-3")  # = ufloat(0.68, 0.041)
    >>> x = ufloat_fromstr("23.2")  # = ufloat(23.2, 0.1)
    >>> x = ufloat_fromstr("23.29")  # = ufloat(23.29, 0.01)

    >>> x = ufloat_fromstr("680.3(nan)") # = ufloat(680.3, numpy.nan)
    """
    (nom, std) = str_to_number_with_uncert(representation.strip())
    return ufloat(nom, std, tag)


def ufloat(nominal_value, std_dev=None, tag=None):
    """
    Create an uncertainties Variable

    Arguments:
    ----------
    nominal_value: float
        nominal value of Variable
    std_dev:   float or `None`
        standard error of Variable, or `None` if not available [`None`]
    tag:   string or `None`
        optional tag for tracing and organizing Variables ['None']

    Returns:
    --------
    uncertainties Variable

    Examples
    ----------
    >>> a = ufloat(5, 0.2)
    >>> b = ufloat(1000, 30, tag='kilo')


    Notes:
    --------
    1. `nominal_value` is typically interpreted as `mean` or `central value`
    2. `std_dev` is typically interpreted as `standard deviation` or the
        1-sigma level uncertainty.
    3. The returned Variable will have attributes `nominal_value`, `std_dev`,
       and `tag` which match the input values.
    """

    return Variable(nominal_value, std_dev, tag=tag)
