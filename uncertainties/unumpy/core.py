"""
Core functions used by unumpy and some of its submodules.

(c) 2010-2016 by Eric O. LEBIGOT (EOL).
"""

# The functions found in this module cannot be defined in unumpy or
# its submodule: this creates import loops, when unumpy explicitly
# imports one of the submodules in order to make it available to the
# user.

from __future__ import division

# Standard modules:
from builtins import next
from builtins import zip
from builtins import range
import sys
import inspect

# 3rd-party modules:
import numpy

# Local modules:
import uncertainties.umath_core as umath_core
import uncertainties.core as uncert_core

__all__ = [
    # Factory functions:
    'uarray', 'umatrix',

    # Utilities:
    'nominal_values', 'std_devs',

    # Classes:
    'matrix'
    ]

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
    doc=("Return the nominal value of the numbers with uncertainties contained"
         " in a NumPy (or unumpy) array (this includes matrices)."))

to_std_devs = numpy.vectorize(
    uncert_core.std_dev,
    otypes=[float],  # Because vectorize() has side effects (dtype setting)
    doc=("Return the standard deviation of the numbers with uncertainties"
         " contained in a NumPy array, or zero for other objects."))

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

def derivative(u, var):
    """
    Return the derivative of u along var, if u is an
    uncert_core.AffineScalarFunc instance, and if var is one of the
    variables on which it depends.  Otherwise, return 0.
    """
    if isinstance(u, uncert_core.AffineScalarFunc):
        try:
            return u.derivatives[var]
        except KeyError:
            return 0.
    else:
        return 0.

def wrap_array_func(func):
    # !!! This function is not used in the code, except in the tests.
    #
    # !!! The implementation seems superficially similar to
    # uncertainties.core.wrap(): is there code/logic duplication
    # (which should be removed)?
    """
    Return a version of the function func() that works even when
    func() is given a NumPy array that contains numbers with
    uncertainties, as first argument.

    This wrapper is similar to uncertainties.core.wrap(), except that
    it handles an array argument instead of float arguments, and that
    the result can be an array.

    However, the returned function is more restricted: the array
    argument cannot be given as a keyword argument with the name in
    the original function (it is not a drop-in replacement).

    func -- function whose first argument is a single NumPy array,
    and which returns a NumPy array.
    """

    @uncert_core.set_doc("""\
    Version of %s(...) that works even when its first argument is a NumPy
    array that contains numbers with uncertainties.

    Warning: elements of the first argument array that are not
    AffineScalarFunc objects must not depend on uncert_core.Variable
    objects in any way.  Otherwise, the dependence of the result in
    uncert_core.Variable objects will be incorrect.

    Original documentation:
    %s""" % (func.__name__, func.__doc__))
    def wrapped_func(arr, *args, **kwargs):
        # Nominal value:
        arr_nominal_value = nominal_values(arr)
        func_nominal_value = func(arr_nominal_value, *args, **kwargs)

        # The algorithm consists in numerically calculating the derivatives
        # of func:

        # Variables on which the array depends are collected:
        variables = set()
        for element in arr.flat:
            # floats, etc. might be present
            if isinstance(element, uncert_core.AffineScalarFunc):
                # !!!! The following forces an evaluation of the
                # derivatives!? Isn't this very slow, when
                # working with a large number of arrays?
                #
                # !! set() is only needed for Python 2 compatibility:
                variables |= set(element.derivatives.keys())

        # If the matrix has no variables, then the function value can be
        # directly returned:
        if not variables:
            return func_nominal_value

        # Calculation of the derivatives of each element with respect
        # to the variables.  Each element must be independent of the
        # others.  The derivatives have the same shape as the output
        # array (which might differ from the shape of the input array,
        # in the case of the pseudo-inverse).
        derivatives = numpy.vectorize(lambda _: {})(func_nominal_value)
        for var in variables:

            # A basic assumption of this package is that the user
            # guarantees that uncertainties cover a zone where
            # evaluated functions are linear enough.  Thus, numerical
            # estimates of the derivative should be good over the
            # standard deviation interval.  This is true for the
            # common case of a non-zero standard deviation of var.  If
            # the standard deviation of var is zero, then var has no
            # impact on the uncertainty of the function func being
            # calculated: an incorrect derivative has no impact.  One
            # scenario can give incorrect results, however, but it
            # should be extremely uncommon: the user defines a
            # variable x with 0 standard deviation, sets y = func(x)
            # through this routine, changes the standard deviation of
            # x, and prints y; in this case, the uncertainty on y
            # might be incorrect, because this program had no idea of
            # the scale on which func() is linear, when it calculated
            # the numerical derivative.

            # The standard deviation might be numerically too small
            # for the evaluation of the derivative, though: we set the
            # minimum variable shift.

            shift_var = max(var._std_dev/1e5, 1e-8*abs(var._nominal_value))
            # An exceptional case is that of var being exactly zero.
            # In this case, an arbitrary shift is used for the
            # numerical calculation of the derivative.  The resulting
            # derivative value might be quite incorrect, but this does
            # not matter as long as the uncertainty of var remains 0,
            # since it is, in this case, a constant.
            if not shift_var:
                shift_var = 1e-8

            # Shift of all the elements of arr when var changes by shift_var:
            shift_arr = array_derivative(arr, var)*shift_var

            # Origin value of array arr when var is shifted by shift_var:
            shifted_arr_values = arr_nominal_value + shift_arr
            func_shifted = func(shifted_arr_values, *args, **kwargs)
            numerical_deriv = (func_shifted-func_nominal_value)/shift_var

            # Update of the list of variables and associated
            # derivatives, for each element:
            for (derivative_dict, derivative_value) in (
                zip(derivatives.flat, numerical_deriv.flat)):

                if derivative_value:
                    derivative_dict[var] = derivative_value

        # numbers with uncertainties are built from the result:
        return numpy.vectorize(uncert_core.AffineScalarFunc)(
            func_nominal_value,
            numpy.vectorize(uncert_core.LinearCombination)(derivatives))

    wrapped_func = uncert_core.set_doc("""\
    Version of %s(...) that works even when its first argument is a NumPy
    array that contains numbers with uncertainties.

    Warning: elements of the first argument array that are not
    AffineScalarFunc objects must not depend on uncert_core.Variable
    objects in any way.  Otherwise, the dependence of the result in
    uncert_core.Variable objects will be incorrect.

    Original documentation:
    %s""" % (func.__name__, func.__doc__))(wrapped_func)

    # It is easier to work with wrapped_func, which represents a
    # wrapped version of 'func', when it bears the same name as
    # 'func' (the name is used by repr(wrapped_func)).
    wrapped_func.__name__ = func.__name__

    return wrapped_func

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
        raise TypeError('uarray() should be called with two arguments.')

    return (numpy.vectorize(
        # ! Looking up uncert_core.Variable beforehand through
        # '_Variable = uncert_core.Variable' does not result in a
        # significant speed up:
        lambda v, s: uncert_core.Variable(v, s), otypes=[object])
        (nominal_values, std_devs))

###############################################################################

def array_derivative(array_like, var):
    """
    Return the derivative of the given array with respect to the
    given variable.

    The returned derivative is a NumPy ndarray of the same shape as
    array_like, that contains floats.

    array_like -- array-like object (list, etc.)  that contains
    scalars or numbers with uncertainties.

    var -- Variable object.
    """
    return numpy.vectorize(lambda u: derivative(u, var),
                           # The type is set because an
                           # integer derivative should not
                           # set the output type of the
                           # array:
                           otypes=[float])(array_like)

def func_with_deriv_to_uncert_func(func_with_derivatives):
    # This function is used for instance for the calculation of the
    # inverse and pseudo-inverse of a matrix with uncertainties.
    """
    Return a function that can be applied to array-like objects that
    contain numbers with uncertainties (lists, lists of lists, NumPy
    arrays, etc.).

    func_with_derivatives -- defines a function that takes an
    array-like object containing scalars and returns an array.  Both
    the value and the derivatives of this function with respect to
    multiple scalar parameters are calculated by this
    func_with_derivatives() argument.

    func_with_derivatives(arr, input_type, derivatives, *args,
    **kwargs) must return an iterator.  The first element returned by
    this iterator is the value of the function at the n-dimensional
    array-like 'arr' (with the correct type).  The following elements
    are arrays that represent the derivative of the function for each
    derivative array from the iterator 'derivatives'.

    func_with_derivatives() takes the following arguments:

      arr -- NumPy ndarray of scalars where the function must be
      evaluated.

      input_type -- data type of the input array-like object.  This
      type is used for determining the type that the function should
      return.

      derivatives -- iterator that returns the derivatives of the
      argument of the function with respect to multiple scalar
      variables.  func_with_derivatives() returns the derivatives of
      the defined function with respect to these variables.

      args -- additional arguments that define the result (example:
      for the pseudo-inverse numpy.linalg.pinv: numerical cutoff).

    Examples of func_with_derivatives: inv_with_derivatives().
    """

    def wrapped_func(array_like, *args, **kwargs):
        """
        array_like -- n-dimensional array-like object that contains
        numbers with uncertainties (list, NumPy ndarray or matrix,
        etc.).

        args -- additional arguments that are passed directly to
        func_with_derivatives.
        """

        # The calculation below is not lazy, contrary to the linear
        # error propagation done in AffineScalarFunc. Making it lazy
        # in the same way would be quite a specific task: basically
        # this would amount to generalizing scalar coefficients in
        # core.LinearCombination to more general matrix
        # multiplications, and to replace Variable differentials by
        # full matrices of coefficients. This does not look very
        # efficient, as matrices are quite big, and since caching the
        # result of a few matrix functions that are not typically
        # stringed one after the other (unlike a big sum of numbers)
        # should not be needed.

        # So that .flat works even if array_like is a list:
        array_version = numpy.asanyarray(array_like)

        # Variables on which the array depends are collected:
        variables = set()
        for element in array_version.flat:
            # floats, etc. might be present
            if isinstance(element, uncert_core.AffineScalarFunc):
                # !!! set() is only needed for Python 2 compatibility:
                variables |= set(element.derivatives.keys())

        array_nominal = nominal_values(array_version)
        # Function value, then derivatives at array_nominal (the
        # derivatives are with respect to the variables contained in
        # array_like):
        func_then_derivs = func_with_derivatives(
            array_nominal,
            type(array_like),
            (array_derivative(array_version, var) for var in variables),
            *args, **kwargs)

        func_nominal_value = next(func_then_derivs)

        if not variables:
            return func_nominal_value

        # The result is built progressively, with the contribution of
        # each variable added in turn:

        # Calculation of the derivatives of the result with respect to
        # the variables.
        derivatives = (
            numpy.array(
                [{} for _ in range(func_nominal_value.size)], dtype=object)
            .reshape(func_nominal_value.shape))

        # Memory-efficient approach.  A memory-hungry approach would
        # be to calculate the matrix derivatives will respect to all
        # variables and then combine them into a matrix of
        # AffineScalarFunc objects.  The approach followed here is to
        # progressively build the matrix of derivatives, by
        # progressively adding the derivatives with respect to
        # successive variables.
        for (var, deriv_wrt_var) in zip(variables,
                                                   func_then_derivs):

            # Update of the list of variables and associated
            # derivatives, for each element:
            for (derivative_dict, derivative_value) in zip(
                derivatives.flat, deriv_wrt_var.flat):
                if derivative_value:
                    derivative_dict[var] = derivative_value

        # An array of numbers with uncertainties is built from the
        # result:
        result = numpy.vectorize(uncert_core.AffineScalarFunc)(
            func_nominal_value,
            numpy.vectorize(uncert_core.LinearCombination)(derivatives))

        # NumPy matrices that contain numbers with uncertainties are
        # better as unumpy matrices:
        if isinstance(result, numpy.matrix):
            result = result.view(matrix)

        return result

    return wrapped_func

########## Matrix inverse

def inv_with_derivatives(arr, input_type, derivatives):
    """
    Defines the matrix inverse and its derivatives.

    See the definition of func_with_deriv_to_uncert_func() for its
    detailed semantics.
    """

    inverse = numpy.linalg.inv(arr)
    # The inverse of a numpy.matrix is a numpy.matrix.  It is assumed
    # that numpy.linalg.inv is such that other types yield
    # numpy.ndarrays:
    if issubclass(input_type, numpy.matrix):
        inverse = inverse.view(numpy.matrix)
    yield inverse

    # It is mathematically convenient to work with matrices:
    inverse_mat = numpy.asmatrix(inverse)

    # Successive derivatives of the inverse:
    for derivative in derivatives:
        derivative_mat = numpy.asmatrix(derivative)
        yield -inverse_mat * derivative_mat * inverse_mat

inv = func_with_deriv_to_uncert_func(inv_with_derivatives)
inv.__doc__ = """\
    Version of numpy.linalg.inv that works with array-like objects
    that contain numbers with uncertainties.

    The result is a unumpy.matrix if numpy.linalg.pinv would return a
    matrix for the array of nominal values.

    Analytical formulas are used.

    Original documentation:
    %s
    """ % numpy.linalg.inv.__doc__

########## Matrix pseudo-inverse

def pinv_with_derivatives(arr, input_type, derivatives, rcond):
    """
    Defines the matrix pseudo-inverse and its derivatives.

    Works with real or complex matrices.

    See the definition of func_with_deriv_to_uncert_func() for its
    detailed semantics.
    """

    inverse = numpy.linalg.pinv(arr, rcond)
    # The pseudo-inverse of a numpy.matrix is a numpy.matrix.  It is
    # assumed that numpy.linalg.pinv is such that other types yield
    # numpy.ndarrays:
    if issubclass(input_type, numpy.matrix):
        inverse = inverse.view(numpy.matrix)
    yield inverse

    # It is mathematically convenient to work with matrices:
    inverse_mat = numpy.asmatrix(inverse)

    # Formula (4.12) from The Differentiation of Pseudo-Inverses and
    # Nonlinear Least Squares Problems Whose Variables
    # Separate. Author(s): G. H. Golub and V. Pereyra. Source: SIAM
    # Journal on Numerical Analysis, Vol. 10, No. 2 (Apr., 1973),
    # pp. 413-432

    # See also
    # http://mathoverflow.net/questions/25778/analytical-formula-for-numerical-derivative-of-the-matrix-pseudo-inverse

    # Shortcuts.  All the following factors should be numpy.matrix objects:
    PA = arr*inverse_mat
    AP = inverse_mat*arr
    factor21 = inverse_mat*inverse_mat.H
    factor22 = numpy.eye(arr.shape[0])-PA
    factor31 = numpy.eye(arr.shape[1])-AP
    factor32 = inverse_mat.H*inverse_mat

    # Successive derivatives of the inverse:
    for derivative in derivatives:
        derivative_mat = numpy.asmatrix(derivative)
        term1 = -inverse_mat*derivative_mat*inverse_mat
        derivative_mat_H = derivative_mat.H
        term2 = factor21*derivative_mat_H*factor22
        term3 = factor31*derivative_mat_H*factor32
        yield term1+term2+term3

# Default rcond argument for the generalization of numpy.linalg.pinv:
#
# Most common modern case first:
try: 
    pinv_default = (
        inspect.signature(numpy.linalg.pinv).parameters["rcond"].default)
except AttributeError:  # No inspect.signature() before Python 3.3
    try:
        # In numpy 1.17+, pinv is wrapped using a decorator which unfortunately
        # results in the metadata (argument defaults) being lost. However, we
        # can still get at the original function using the __wrapped__
        # attribute (which is what inspect.signature() does).
        pinv_default = numpy.linalg.pinv.__wrapped__.__defaults__[0]
    except AttributeError:  # Function not wrapped in NumPy < 1.17
        pinv_default = numpy.linalg.pinv.__defaults__[0]  # Python 1, 2.6+:

pinv_with_uncert = func_with_deriv_to_uncert_func(pinv_with_derivatives)

def pinv(array_like, rcond=pinv_default):
    return pinv_with_uncert(array_like, rcond)

pinv = uncert_core.set_doc("""
    Version of numpy.linalg.pinv that works with array-like objects
    that contain numbers with uncertainties.

    The result is a unumpy.matrix if numpy.linalg.pinv would return a
    matrix for the array of nominal values.

    Analytical formulas are used.

    Original documentation:
    %s
    """ % numpy.linalg.pinv.__doc__)(pinv)

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

    I = numpy.matrix.I.getter(getI)


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
        raise TypeError('umatrix() should be called with two arguments.')

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
    func_name_translations = dict([
        (f_name, 'arc'+f_name[1:])
        for f_name in ['acos', 'acosh', 'asin', 'atan', 'atan2', 'atanh']])

    new_func_names = [
        func_name_translations.get(function_name, function_name)
        # The functions from umath_core.non_std_wrapped_funcs
        # (available from umath) are normally not in
        # NumPy, so they are not included here:
        for function_name in umath_core.many_scalars_to_scalar_funcs]

    for (function_name, unumpy_name) in zip(
        umath_core.many_scalars_to_scalar_funcs, new_func_names):

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
            {} if function_name in umath_core.locally_cst_funcs
            # If by any chance a function returns, in a particular
            # case, an integer instead of a number with uncertainty,
            # side-effects in vectorize() would fix the resulting
            # dtype to integer, which is not what is wanted (as
            # vectorize(), at least in NumPy around 2010 maybe,
            # decided about the output data type by looking at the
            # type of first element only).
            else {'otypes': [object]})

        setattr(
            this_module, unumpy_name,
            #!!!! For umath_core.locally_cst_funcs, would it make sense
            # to optimize this by using instead the equivalent (? see
            # above) vectorized NumPy function on the nominal values?
            numpy.vectorize(func,
                            doc="""\
Vectorized version of umath.%s.

Original documentation:
%s""" % (function_name, func.__doc__),
                            **otypes))


        __all__.append(unumpy_name)

define_vectorized_funcs()
