# This file contains code for AffineScalarFunc's arithmetic and comparative ops.

from math import sqrt, log  # Optimization: no attribute look-up
import sys
import itertools
from inspect import getfullargspec
import numbers

from uncertainties.ucombo import UCombo

# Some types known to not depend on Variable objects are put in
# CONSTANT_TYPES.  The most common types can be put in front, as this
# may slightly improve the execution speed.
FLOAT_LIKE_TYPES = (numbers.Number,)
CONSTANT_TYPES = FLOAT_LIKE_TYPES + (complex,)

try:
    import numpy
except ImportError:
    pass
else:
    # NumPy numbers do not depend on Variable objects:
    FLOAT_LIKE_TYPES += (numpy.generic,)
    CONSTANT_TYPES += FLOAT_LIKE_TYPES[-1:]


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
    """
    Wrapper around f(x, y) that let f return NaN when f raises one of
    a few numerical exceptions.
    """

    def wrapped_f(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except (ValueError, ZeroDivisionError, OverflowError):
            return float("nan")

    return wrapped_f


def pow_deriv_0(x, y):
    """
    The formula below works if x is positive or if y is an integer and x is negative
    of y is an integer, x is zero and y is greater than or equal to 1.
    """
    if x > 0 or (y % 1 == 0 and (x < 0 or y >= 1)):
        return y * x ** (y - 1)
    elif x == 0 and y == 0:
        return 0
    else:
        return float("nan")


def pow_deriv_1(x, y):
    if x > 0:
        return log(x) * x**y
    elif x == 0 and y > 0:
        return 0
    else:
        return float("nan")


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
        "add": ("1.", "1."),
        # 'div' is the '/' operator when __future__.division is not in
        # effect.  Since '/' is applied to
        # AffineScalarFunc._nominal_value numbers, it is applied on
        # floats, and is therefore the "usual" mathematical division.
        "div": ("1/y", "-x/y**2"),
        # The derivative wrt the 2nd arguments is something like (..., x//y),
        # but it is calculated numerically, for convenience:
        "mul": ("y", "x"),
        "sub": ("1.", "-1."),
        "truediv": ("1/y", "-x/y**2"),
    }

    # Conversion to Python functions:
    ops_with_reflection = {}
    for op, derivatives in derivatives_list.items():
        ops_with_reflection[op] = [
            eval("lambda x, y: %s" % expr) for expr in derivatives
        ]

        ops_with_reflection["r" + op] = [
            eval("lambda y, x: %s" % expr) for expr in reversed(derivatives)
        ]

    ops_with_reflection["pow"] = [pow_deriv_0, pow_deriv_1]
    ops_with_reflection["rpow"] = [
        lambda y, x: pow_deriv_1(x, y),
        lambda y, x: pow_deriv_0(x, y),
    ]

    # Undefined derivatives are converted to NaN when the function
    # itself can be calculated:
    for op in ["pow"]:
        ops_with_reflection[op] = [
            nan_if_exception(func) for func in ops_with_reflection[op]
        ]
        ops_with_reflection["r" + op] = [
            nan_if_exception(func) for func in ops_with_reflection["r" + op]
        ]

    return ops_with_reflection


# Operators that have a reflection, along with their derivatives:
ops_with_reflection = get_ops_with_reflection()

# Some effectively modified operators (for the automated tests):
modified_operators = []
modified_ops_with_reflection = []


# !!! This code is not run by the tests. It would be nice to have
# it be tested.
def no_complex_result(func):
    """
    Return a function that does like func, but that raises a
    ValueError if the result is complex.
    """

    def no_complex_func(*args, **kwargs):
        """
        Like %s, but raises a ValueError exception if the result
        is complex.
        """ % func.__name__

        value = func(*args, **kwargs)
        if isinstance(value, complex):
            raise ValueError(
                "The uncertainties module does not handle" " complex results"
            )
        else:
            return value

    return no_complex_func


# This module does not handle uncertainties on complex numbers:
# complex results for the nominal value of some operations cannot
# be calculated with an uncertainty:
custom_ops = {
    "pow": no_complex_result(float.__pow__),
    "rpow": no_complex_result(float.__rpow__),
}


def add_arithmetic_ops(cls):
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
            return 1.0
        else:
            return -1.0

    # Single-argument operators that should be adapted from floats to
    # AffineScalarFunc objects, associated to their derivative:
    simple_numerical_operators_derivatives = {
        "neg": lambda x: -1.0,
        "pos": lambda x: 1.0,
    }

    for op, derivative in iter(simple_numerical_operators_derivatives.items()):
        attribute_name = "__%s__" % op

        # float objects don't exactly have the same attributes between
        # different versions of Python (for instance, __trunc__ was
        # introduced with Python 2.6):
        try:
            setattr(
                cls,
                attribute_name,
                _wrap(cls, getattr(float, attribute_name), [derivative]),
            )
        except AttributeError:
            # Version of Python where floats don't have attribute_name:
            pass
        else:
            modified_operators.append(op)

    ########################################
    # Final definition of the operators for AffineScalarFunc objects:

    # Reversed versions (useful for float*AffineScalarFunc, for instance):
    for op, derivatives in ops_with_reflection.items():
        attribute_name = "__%s__" % op

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
            setattr(cls, attribute_name, _wrap(cls, func_to_wrap, derivatives))
            modified_ops_with_reflection.append(op)

    ########################################
    # Conversions to pure numbers are meaningless.  Note that the
    # behavior of float(1j) is similar.
    for coercion_type in ("complex", "int", "long", "float"):

        def raise_error(self):
            raise TypeError(
                "can't convert an affine function (%s)"
                " to %s; use x.nominal_value" % (self.__class__, coercion_type)
                # In case AffineScalarFunc is sub-classed:
            )

        setattr(cls, "__%s__" % coercion_type, raise_error)


class IndexableIter(object):
    """
    Iterable whose values can also be accessed through indexing.

    The input iterable values are cached.

    Some attributes:

    iterable -- iterable used for returning the elements one by one.

    returned_elements -- list with the elements directly accessible.
    through indexing. Additional elements are obtained from self.iterable.

    none_converter -- function that takes an index and returns the
    value to be returned when None is obtained form the iterable
    (instead of None).
    """

    def __init__(self, iterable, none_converter=lambda index: None):
        """
        iterable -- iterable whose values will be returned.

        none_converter -- function applied to None returned
        values. The value that replaces None is none_converter(index),
        where index is the index of the element.
        """
        self.iterable = iterable
        self.returned_elements = []
        self.none_converter = none_converter

    def __getitem__(self, index):
        returned_elements = self.returned_elements

        try:
            return returned_elements[index]

        except IndexError:  # Element not yet cached
            for pos in range(len(returned_elements), index + 1):
                value = next(self.iterable)

                if value is None:
                    value = self.none_converter(pos)

                returned_elements.append(value)

            return returned_elements[index]

    def __str__(self):
        return "<%s: [%s...]>" % (
            self.__class__.__name__,
            ", ".join(map(str, self.returned_elements)),
        )


def ufloat_from_uncertainty(cls, nominal_value: float, uncertainty: UCombo):
    # TODO: It is a hack that needs to be removed that this function uses a generic
    #   cls input. This issue stems from the monkey patching that connects core.py and
    #   ops.py.
    result = object.__new__(cls)
    result._nominal_value = nominal_value
    result._uncertainty = uncertainty
    return result


def _wrap(cls, f, derivatives_args=None, derivatives_kwargs=None):
    if derivatives_args is None:
        derivatives_args = []
    if derivatives_kwargs is None:
        derivatives_kwargs = {}
    derivatives_args_index = IndexableIter(
        # Automatic addition of numerical derivatives in case the
        # supplied derivatives_args is shorter than the number of
        # arguments in *args:
        itertools.chain(derivatives_args, itertools.repeat(None))
    )

    # Derivatives for keyword arguments (includes var-keyword
    # parameters **kwargs, but also var-or-keyword parameters, and
    # keyword-only parameters (Python 3):

    derivatives_all_kwargs = {}

    for name, derivative in derivatives_kwargs.items():
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
        argspec = getfullargspec(f)
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

        for index, name in enumerate(argspec.args):
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

    none_converter = lambda index: partial_derivative(f, index)  # noqa

    for index, derivative in enumerate(derivatives_args_index.returned_elements):
        if derivative is None:
            derivatives_args_index.returned_elements[index] = none_converter(index)

    # Future None values are also automatically converted:
    derivatives_args_index.none_converter = none_converter

    ## Wrapped function:

    #! Setting the doc string after "def f_with...()" does not
    # seem to work.  We define it explicitly:
    @set_doc(
        """\
    Version of %s(...) that returns an affine approximation
    (AffineScalarFunc object), if its result depends on variables
    (Variable objects).  Otherwise, returns a simple constant (when
    applied to constant arguments).

    Warning: arguments of the function that are not AffineScalarFunc
    objects must not depend on uncertainties.Variable objects in any
    way.  Otherwise, the dependence of the result in
    uncertainties.Variable objects will be incorrect.

    Original documentation:
    %s"""
        % (f.__name__, f.__doc__)
    )
    def f_with_affine_output(*args, **kwargs):
        ########################################
        # The involved random variables must first be gathered, so
        # that they can be independently updated.

        # The arguments that contain an uncertainty (AffineScalarFunc
        # objects) are gathered, as positions or names; they will be
        # replaced by their nominal value in order to calculate
        # the necessary derivatives of f.

        pos_w_uncert = [
            index for (index, value) in enumerate(args) if isinstance(value, cls)
        ]
        names_w_uncert = [
            key for (key, value) in kwargs.items() if isinstance(value, cls)
        ]

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
        uncertainty = UCombo(())

        for pos in pos_w_uncert:
            arg_uncertainty = args[pos]._uncertainty
            # if arg_uncertainty:
            derivative_val = derivatives_args_index[pos](*args_values, **kwargs)
            uncertainty += derivative_val * arg_uncertainty

        for name in names_w_uncert:
            # Optimization: caching of the automatic numerical
            # derivatives for keyword arguments that are
            # discovered. This gives a speedup when the original
            # function is called repeatedly with the same keyword
            # arguments:
            # if not kwargs_uncert_values[name].uncertainty:
            #     continue
            derivative_func = derivatives_all_kwargs.setdefault(
                name,
                # Derivative never needed before:
                partial_derivative(f, name),
            )
            derivative_val = derivative_func(*args_values, **kwargs)
            uncertainty += derivative_val * kwargs_uncert_values[name]._uncertainty

        # The function now returns the necessary linear approximation
        # to the function:
        return ufloat_from_uncertainty(cls, f_nominal_value, uncertainty)

    f_with_affine_output = set_doc(
        """\
    Version of %s(...) that returns an affine approximation
    (AffineScalarFunc object), if its result depends on variables
    (Variable objects).  Otherwise, returns a simple constant (when
    applied to constant arguments).

    Warning: arguments of the function that are not AffineScalarFunc
    objects must not depend on uncertainties.Variable objects in any
    way.  Otherwise, the dependence of the result in
    uncertainties.Variable objects will be incorrect.

    Original documentation:
    %s"""
        % (f.__name__, f.__doc__)
    )(f_with_affine_output)

    # It is easier to work with f_with_affine_output, which represents
    # a wrapped version of 'f', when it bears the same name as 'f':
    # ! __name__ is read-only, in Python 2.3:
    f_with_affine_output.name = f.__name__

    return f_with_affine_output


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
    change_kwargs = isinstance(arg_ref, str)

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
        step = STEP_SIZE * abs(args_with_var[arg_ref])
        if not step:
            # Arbitrary, but "small" with respect to 1:
            step = STEP_SIZE

        args_with_var[arg_ref] += step

        if change_kwargs:
            shifted_f_plus = f(*args, **args_with_var)
        else:
            shifted_f_plus = f(*args_with_var, **kwargs)

        args_with_var[arg_ref] -= 2 * step  # Optimization: only 1 list copy

        if change_kwargs:
            shifted_f_minus = f(*args, **args_with_var)
        else:
            shifted_f_minus = f(*args_with_var, **kwargs)

        return (shifted_f_plus - shifted_f_minus) / 2 / step

    return partial_derivative_of_f
