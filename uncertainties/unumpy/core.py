"""
Core functions used by unumpy and some of its submodules.

(c) 2010-2016 by Eric O. LEBIGOT (EOL).
"""

# The functions found in this module cannot be defined in unumpy or
# its submodule: this creates import loops, when unumpy explicitly
# imports one of the submodules in order to make it available to the
# user.

# Standard modules:
from builtins import zip
from functools import wraps
import sys
from typing import Callable, Union
from numbers import Real

# 3rd-party modules:
import numpy

# Local modules:
import uncertainties.umath_core as umath_core
import uncertainties.core as uncert_core
from uncertainties import UFloat
from uncertainties.ucombo import UCombo


__all__ = [
    # Factory functions:
    "uarray",
    "umatrix",
    # Utilities:
    "nominal_values",
    "std_devs",
    # Classes:
    "matrix",
]


def inject_to_args_kwargs(param, injected_arg, *args, **kwargs):
    if isinstance(param, int):
        new_kwargs = kwargs
        new_args = []
        for idx, arg in enumerate(args):
            if idx == param:
                new_args.append(injected_arg)
            else:
                new_args.append(arg)
    elif isinstance(param, str):
        new_args = args
        new_kwargs = kwargs
        new_kwargs[param] = injected_arg
    else:
        raise TypeError(f"{param} must be an int or str, not {type(param)}.")
    return new_args, new_kwargs


def get_args_kwargs_list(*args, **kwargs):
    args_kwargs_list = []
    for idx, arg in enumerate(args):
        args_kwargs_list.append((idx, arg))
    for key, arg in kwargs.items():
        args_kwargs_list.append((key, arg))
    return args_kwargs_list


SQRT_EPS = numpy.sqrt(sys.float_info.epsilon)


def array_numerical_partial_derivative(
    f: Callable[..., Real],
    target_param: Union[str, int],
    array_multi_index: tuple = None,
    *args,
    **kwargs,
) -> float:
    """
    Numerically calculate the partial derivative of a function f with respect to the
    target_param (string name or position number of the float parameter to f to be
    varied) holding all other arguments, *args and **kwargs, constant.
    """
    if isinstance(target_param, int):
        x = args[target_param]
    else:
        x = kwargs[target_param]

    if array_multi_index is None:
        dx = abs(x) * SQRT_EPS  # Numerical Recipes 3rd Edition, eq. 5.7.5
        x_lower = x - dx
        x_upper = x + dx
    else:
        dx = numpy.mean(numpy.abs(x)) * SQRT_EPS
        x_lower = numpy.copy(x)
        x_upper = numpy.copy(x)
        x_lower[array_multi_index] -= dx
        x_upper[array_multi_index] += dx

    lower_args, lower_kwargs = inject_to_args_kwargs(
        target_param,
        x_lower,
        *args,
        **kwargs,
    )
    upper_args, upper_kwargs = inject_to_args_kwargs(
        target_param,
        x_upper,
        *args,
        **kwargs,
    )

    lower_y = f(*lower_args, **lower_kwargs)
    upper_y = f(*upper_args, **upper_kwargs)

    derivative = (upper_y - lower_y) / (2 * dx)
    return derivative


def to_uarray_func(func, derivs=None):
    if derivs is None:
        derivs = {}

    @wraps(func)
    def wrapped(*args, **kwargs):
        return_uarray = False
        return_matrix = False

        float_args = []
        for arg in args:
            if isinstance(arg, UFloat):
                float_args.append(arg.nominal_value)
                return_uarray = True
            elif isinstance(arg, numpy.ndarray):
                if isinstance(arg, numpy.matrix):
                    return_matrix = True
                if isinstance(arg.flat[0], UFloat):
                    float_args.append(to_nominal_values(arg))
                    return_uarray = True
                else:
                    float_args.append(arg)
            else:
                float_args.append(arg)

        float_kwargs = {}
        for key, arg in kwargs.items():
            if isinstance(arg, UFloat):
                float_kwargs[key] = arg.nominal_value
                return_uarray = True
            elif isinstance(arg, numpy.ndarray):
                if isinstance(arg, numpy.matrix):
                    return_matrix = True
                if isinstance(arg.flat[0], UFloat):
                    float_kwargs[key] = to_nominal_values(arg)
                    return_uarray = True
                else:
                    float_kwargs[key] = arg
            else:
                float_kwargs[key] = arg

        new_nominal_array = func(*float_args, **float_kwargs)
        if not return_uarray:
            return new_nominal_array

        args_kwargs_list = get_args_kwargs_list(*args, **kwargs)

        ucombo_array = numpy.full(new_nominal_array.shape, UCombo(()))

        for label, arg in args_kwargs_list:
            if isinstance(arg, UFloat):
                if label in derivs and derivs[label] is not None:
                    deriv_func = derivs[label]
                    deriv_arr = deriv_func(
                        label,
                        None,
                        *float_args,
                        **float_kwargs,
                    )
                else:
                    deriv_arr = array_numerical_partial_derivative(
                        func, label, None, *float_args, **float_kwargs
                    )
                ucombo_array += deriv_arr * arg.uncertainty
            elif isinstance(arg, numpy.ndarray):
                if isinstance(arg.flat[0], UFloat):
                    it = numpy.nditer(arg, flags=["multi_index", "refs_ok"])
                    for sub_arg in it:
                        u_num = sub_arg.item()
                        if isinstance(u_num, UFloat):
                            multi_index = it.multi_index
                            if label in derivs and derivs[label] is not None:
                                deriv_func = derivs[label]
                                deriv_arr = deriv_func(
                                    label,
                                    multi_index,
                                    *float_args,
                                    **float_kwargs,
                                )
                            else:
                                deriv_arr = array_numerical_partial_derivative(
                                    func,
                                    label,
                                    multi_index,
                                    *float_args,
                                    **float_kwargs,
                                )
                            ucombo_array += numpy.array(deriv_arr) * u_num.uncertainty

        u_array = numpy.vectorize(UFloat)(new_nominal_array, ucombo_array)
        if return_matrix:
            u_array = numpy.matrix(u_array)
        return u_array

    return wrapped


###############################################################################
# Utilities:

# nominal_values() and std_devs() are defined as functions (instead of
# as additional methods of the unumpy.matrix class) because the user
# might well directly build arrays of numbers with uncertainties
# without going through the factory functions found in this module
# (uarray() and umatrix()).  Thus,
# numpy.array([uncert_core.ufloat((1, 0.1))]) would not
# have a nominal_values() method.  Adding such a method to, say,
# unumpy.matrix, would break the symmetry between NumPy arrays and
# matrices (no nominal_values() method), and objects defined in this
# module.

# ! Warning: the __doc__ is set, but help(nominal_values) does not
# display it, but instead displays the documentation for the type of
# nominal_values (i.e. the documentation of its class):

to_nominal_values = numpy.vectorize(
    uncert_core.nominal_value,
    otypes=[float],  # Because vectorize() has side effects (dtype setting)
    doc=(
        "Return the nominal value of the numbers with uncertainties contained"
        " in a NumPy (or unumpy) array (this includes matrices)."
    ),
)

to_std_devs = numpy.vectorize(
    uncert_core.std_dev,
    otypes=[float],  # Because vectorize() has side effects (dtype setting)
    doc=(
        "Return the standard deviation of the numbers with uncertainties"
        " contained in a NumPy array, or zero for other objects."
    ),
)

to_uncertainties = numpy.vectorize(
    uncert_core.uncertainty,
    otypes=[object],
)


def unumpy_to_numpy_matrix(arr):
    """
    If arr in a unumpy.matrix, it is converted to a numpy.matrix.
    Otherwise, it is returned unchanged.
    """

    if isinstance(arr, matrix):
        return arr.view(numpy.matrix)
    else:
        return arr


def nominal_values(arr):
    """
    Return the nominal values of the numbers in NumPy array arr.

    Elements that are not numbers with uncertainties (derived from a
    class from this module) are passed through untouched (because a
    numpy.array can contain numbers with uncertainties and pure floats
    simultaneously).

    If arr is of type unumpy.matrix, the returned array is a
    numpy.matrix, because the resulting matrix does not contain
    numbers with uncertainties.
    """

    return unumpy_to_numpy_matrix(to_nominal_values(arr))


def std_devs(arr):
    """
    Return the standard deviations of the numbers in NumPy array arr.

    Elements that are not numbers with uncertainties (derived from a
    class from this module) are passed through untouched (because a
    numpy.array can contain numbers with uncertainties and pure floats
    simultaneously).

    If arr is of type unumpy.matrix, the returned array is a
    numpy.matrix, because the resulting matrix does not contain
    numbers with uncertainties.
    """

    return unumpy_to_numpy_matrix(to_std_devs(arr))


###############################################################################
# Arrays


def uarray(nominal_values, std_devs=None):
    """
    Return a NumPy array of numbers with uncertainties
    initialized with the given nominal values and standard
    deviations.

    nominal_values, std_devs -- valid arguments for numpy.array, with
    identical shapes (list of numbers, list of lists, numpy.ndarray,
    etc.).

    std_devs=None is only used for supporting legacy code, where
    nominal_values can be the tuple of nominal values and standard
    deviations.
    """

    if std_devs is None:  # Obsolete, single tuple argument call
        raise TypeError("uarray() should be called with two arguments.")

    return numpy.vectorize(
        # ! Looking up uncert_core.Variable beforehand through
        # '_Variable = uncert_core.Variable' does not result in a
        # significant speed up:
        lambda v, s: uncert_core.UFloat(v, s),
        otypes=[object],
    )(nominal_values, std_devs)


###############################################################################


########## Matrix inverse


def inv_deriv(label, multi_index, *args, **kwargs):
    if isinstance(label, int):
        arr = args[label]
    elif isinstance(label, str):
        arr = kwargs[label]
    else:
        raise ValueError

    deriv_arr = numpy.zeros_like(arr)
    deriv_arr[multi_index] = 1

    inv_arr = numpy.linalg.inv(arr)
    return -inv_arr @ deriv_arr @ inv_arr


inv = to_uarray_func(numpy.linalg.inv, derivs={0: inv_deriv, "a": inv_deriv})
inv.__doc__ = (
    """\
    Version of numpy.linalg.inv that works with array-like objects
    that contain numbers with uncertainties.

    The result is a unumpy.matrix if numpy.linalg.pinv would return a
    matrix for the array of nominal values.

    Analytical formulas are used.

    Original documentation:
    %s
    """
    % numpy.linalg.inv.__doc__
)

########## Matrix pseudo-inverse


def pinv_deriv(label, multi_index, *args, **kwargs):
    if isinstance(label, int):
        A = args[label]
    elif isinstance(label, str):
        A = kwargs[label]
    else:
        raise ValueError

    """
    Analytical calculation of the derivative of the Pseudo-Inverse of matrix A with
    respect to an element of matrix A. This calculatinon is done according to
    formula (4.12) from

    G. H. Golub, V. Pereyra, "The Differentiation of Pseudo-Inverses and Nonlinear Least
    Squares Problems Whose Variables Separate", Journal on Numerical Analysis, Vol. 10,
    No. 2 (Apr., 1973), pp. 413-432

    See also
    http://mathoverflow.net/questions/25778/analytical-formula-for-numerical-derivative-of-the-matrix-pseudo-inverse
    """

    Aprime = numpy.zeros_like(A)
    Aprime[multi_index] = 1

    Aplus = numpy.linalg.pinv(*args, **kwargs)

    n, m = A.shape[-2:]
    eye_n = numpy.eye(n)
    eye_m = numpy.eye(m)

    ndim = len(A.shape)
    trans_axes = list(range(ndim))
    trans_axes[0], trans_axes[1] = trans_axes[1], trans_axes[0]

    AplusT = numpy.transpose(Aplus, axes=trans_axes)
    AprimeT = numpy.transpose(Aprime, axes=trans_axes)

    return (
        -Aplus @ Aprime @ Aplus
        + Aplus @ AplusT @ AprimeT @ (eye_n - A @ Aplus)
        + (eye_m - Aplus @ A) @ AprimeT @ AplusT @ Aplus
    )


pinv = to_uarray_func(numpy.linalg.pinv, derivs={0: pinv_deriv, "a": pinv_deriv})

pinv = uncert_core.set_doc(
    """
    Version of numpy.linalg.pinv that works with array-like objects
    that contain numbers with uncertainties.

    The result is a unumpy.matrix if numpy.linalg.pinv would return a
    matrix for the array of nominal values.

    Analytical formulas are used.

    Original documentation:
    %s
    """
    % numpy.linalg.pinv.__doc__
)(pinv)

########## Matrix class


class matrix(numpy.matrix):
    # The name of this class is the same as NumPy's, which is why it
    # does not follow PEP 8.
    """
    Class equivalent to numpy.matrix, but that behaves better when the
    matrix contains numbers with uncertainties.
    """

    def __rmul__(self, other):
        # ! NumPy's matrix __rmul__ uses an apparently restrictive
        # dot() function that cannot handle the multiplication of a
        # scalar and of a matrix containing objects (when the
        # arguments are given in this order).  We go around this
        # limitation:
        if numpy.isscalar(other):
            return numpy.dot(self, other)
        else:
            return numpy.dot(other, self)  # The order is important

    def getI(self):
        """Matrix inverse or pseudo-inverse."""
        m, n = self.shape
        return (inv if m == n else pinv)(self)

    I = numpy.matrix.I.getter(getI)  # noqa

    # !!! The following function is not in the official documentation
    # of the module. Maybe this is because arrays with uncertainties
    # do not have any equivalent in this module, and they should be
    # the first ones to have such methods?
    @property
    def nominal_values(self):
        """
        Nominal value of all the elements of the matrix.
        """
        return nominal_values(self)

    # !!! The following function is not in the official documentation
    # of the module. Maybe this is because arrays with uncertainties
    # do not have any equivalent in this module, and they should be
    # the first ones to have such methods?
    @property
    def std_devs(self):
        return numpy.matrix(std_devs(self))


def umatrix(nominal_values, std_devs=None):
    """
    Constructs a matrix that contains numbers with uncertainties.

    The arguments are the same as for uarray(...): nominal values, and
    standard deviations.

    The returned matrix can be inverted, thanks to the fact that it is
    a unumpy.matrix object instead of a numpy.matrix one.
    """

    if std_devs is None:  # Obsolete, single tuple argument call
        raise TypeError("umatrix() should be called with two arguments.")

    return uarray(nominal_values, std_devs).view(matrix)


###############################################################################


def define_vectorized_funcs():
    """
    Defines vectorized versions of functions from uncertainties.umath_core.

    Some functions have their name translated, so as to follow NumPy's
    convention (example: math.acos -> numpy.arccos).
    """

    this_module = sys.modules[__name__]
    # NumPy does not always use the same function names as the math
    # module:
    func_name_translations = dict(
        [
            (f_name, "arc" + f_name[1:])
            for f_name in ["acos", "acosh", "asin", "atan", "atan2", "atanh"]
        ]
    )

    new_func_names = [
        func_name_translations.get(function_name, function_name)
        # The functions from umath_core.non_std_wrapped_funcs
        # (available from umath) are normally not in
        # NumPy, so they are not included here:
        for function_name in umath_core.many_scalars_to_scalar_funcs
    ]

    for function_name, unumpy_name in zip(
        umath_core.many_scalars_to_scalar_funcs, new_func_names
    ):
        # ! The newly defined functions (uncertainties.unumpy.cos, etc.)
        # do not behave exactly like their NumPy equivalent (numpy.cos,
        # etc.): cos(0) gives an array() and not a
        # numpy.float... (equality tests succeed, though).
        func = getattr(umath_core, function_name)

        # Data type of the result of the unumpy function:
        otypes = (
            # It is much more convenient to preserve the type of
            # functions that return a number without
            # uncertainty. Thus, for example, unumpy.isnan() can
            # return an array with a boolean data type (instead of
            # object), which allows the result to be used with NumPy's
            # boolean indexing.
            {}
            if function_name in umath_core.locally_cst_funcs
            # If by any chance a function returns, in a particular
            # case, an integer instead of a number with uncertainty,
            # side-effects in vectorize() would fix the resulting
            # dtype to integer, which is not what is wanted (as
            # vectorize(), at least in NumPy around 2010 maybe,
            # decided about the output data type by looking at the
            # type of first element only).
            else {"otypes": [object]}
        )

        setattr(
            this_module,
            unumpy_name,
            #!!!! For umath_core.locally_cst_funcs, would it make sense
            # to optimize this by using instead the equivalent (? see
            # above) vectorized NumPy function on the nominal values?
            numpy.vectorize(
                func,
                doc="""\
Vectorized version of umath.%s.

Original documentation:
%s"""
                % (function_name, func.__doc__),
                **otypes,
            ),
        )

        __all__.append(unumpy_name)


define_vectorized_funcs()
