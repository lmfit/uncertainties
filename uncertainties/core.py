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

from builtins import str, zip, range, object
from math import sqrt  # Optimization: no attribute look-up
from numbers import Real
from typing import Optional, Union
import copy

from uncertainties.formatting import format_ufloat
from uncertainties.parsing import str_to_number_with_uncert
from . import ops
from uncertainties.ops import (
    _wrap,
    partial_derivative,
    set_doc,
    nan_if_exception,
    modified_operators,
    modified_ops_with_reflection,
)
from uncertainties.ucombo import UCombo, UAtom

# Attributes that are always exported (some other attributes are
# exported only if the NumPy module is available...):
__all__ = [
    # All sub-modules and packages are not imported by default,
    # in particular because NumPy might be unavailable.
    "ufloat",  # Main function: returns a number with uncertainty
    "ufloat_fromstr",  # Important function: returns a number with uncertainty
    # Uniform access to nominal values and standard deviations:
    "nominal_value",
    "std_dev",
    # Utility functions (more are exported if NumPy is present):
    "covariance_matrix",
    # Class for testing whether an object is a number with
    # uncertainty.  Not usually created by users (except through the
    # Variable subclass), but possibly manipulated by external code
    # ['derivatives()' method, etc.].
    "UFloat",
    # Wrapper for allowing non-pure-Python function to handle
    # quantitities with uncertainties:
    "wrap",
    # used internally and will be removed by linter if not here
    "nan_if_exception",
    "modified_operators",
    "modified_ops_with_reflection",
]

###############################################################################
## Definitions that depend on the availability of NumPy:


try:
    import numpy
except ImportError:
    pass
else:
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
        norm_vector[norm_vector == 0] = 1

        return correlated_values_norm(
            # !! The following zip() is a bit suboptimal: correlated_values()
            # separates back the nominal values and the standard deviations:
            list(zip(nom_values, std_devs)),
            covariance_mat / norm_vector / norm_vector[:, numpy.newaxis],
            tags,
        )

    __all__.append("correlated_values")

    def correlated_values_norm(values_with_std_dev, correlation_mat, tags=None):
        """
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
        """

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
        variances[variances < 0] = 0.0

        # Creation of new, independent variables:

        # We use the fact that the eigenvectors in 'transform' are
        # special: 'transform' is unitary: its inverse is its transpose:

        # variables = tuple(
        #     # The variables represent "pure" uncertainties:
        #     Variable(0, sqrt(variance), tag)
        #     for (variance, tag) in zip(variances, tags)
        # )
        ind_vars = tuple(UCombo(((UAtom(), sqrt(variance)),)) for variance in variances)

        # The coordinates of each new uncertainty as a function of the
        # new variables must include the variable scale (standard deviation):
        transform *= std_devs[:, numpy.newaxis]

        corr_vars = []
        for sub_coords in transform:
            corr_var = sum(
                (ind_var * coord for ind_var, coord in zip(ind_vars, sub_coords)),
                UCombo(()),
            )
            corr_vars.append(corr_var)

        # Representation of the initial correlated values:
        values_funcs = tuple(
            AffineScalarFunc(value, corr_var)
            for (corr_var, value) in zip(corr_vars, nominal_values)
        )

        return values_funcs

    __all__.append("correlated_values_norm")

###############################################################################


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


########################################
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
    __slots__ = ("_nominal_value", "_uncertainty", "tag")

    # !! Fix for mean() in NumPy 1.8.0:
    class dtype(object):
        type = staticmethod(lambda value: value)

    #! The code could be modified in order to accommodate for non-float
    # nominal values.  This could for instance be done through
    # the operator module: instead of delegating operations to
    # float.__*__ operations, they could be delegated to
    # operator.__*__ functions (while taking care of properly handling
    # reverse operations: __radd__, etc.).

    def __init__(
        self,
        nominal_value: float,
        uncertainty: Union[UCombo, float],
        tag: Optional[str] = None,
    ):
        """
        nominal_value -- average value
        uncertainty -- linear combination of uncertain terms
        tag -- label
        """
        self._nominal_value = float(nominal_value)
        if isinstance(uncertainty, Real):
            uncertainty = UCombo(((UAtom(tag=tag), float(uncertainty)),))
        self._uncertainty = uncertainty
        self.tag = tag

    @property
    def nominal_value(self):
        """Nominal value of the random number."""
        return self._nominal_value

    # Abbreviation (for formulas, etc.):
    n = nominal_value

    @property
    def uncertainty(self):
        return self._uncertainty

    ############################################################

    @property
    def _linear_part(self):
        # Legacy
        return self.uncertainty

    def covariance(self, other):
        return self.uncertainty.covariance(other.uncertainty)

    ########################################

    # Uncertainties handling:

    def error_components(self):
        """
        The uncertainty is stored as a float-weighted linear combination of
        UAtom objects. Each UAtom is unique and independent of all other UAtoms.
        Each UAtom has a standard deviation of unity.

        This method returns a dictionary mapping each UAtom with which the
        AffineScalarFunc is correlated to its corresponding float weight.
        """
        return self.uncertainty.expanded

    @property
    def std_dev(self):
        """
        Standard deviation of the affine function.
        """
        return self.uncertainty.std_dev

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
            std_dev_str = "0"

        num_repr = "%r+/-%s" % (self.nominal_value, std_dev_str)
        if self.tag is not None:
            return "< %s = %s >" % (self.tag, num_repr)

        return num_repr

    def __str__(self):
        # An empty format string and str() usually return the same
        # string
        # (http://docs.python.org/2/library/string.html#format-specification-mini-language):
        return self.format("")

    @set_doc(format_ufloat.__doc__)
    def __format__(self, format_spec):
        return format_ufloat(self, format_spec)

    @set_doc("""
        Return the same result as self.__format__(format_spec), or
        equivalently as the format(self, format_spec) of Python 2.6+.

        This method is meant to be used for formatting numbers with
        uncertainties in Python < 2.6, with '... %s ...' %
        num.format('.2e').
        """)
    def format(self, format_spec):
        return format_ufloat(self, format_spec)

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
            raise ValueError("The standard deviation is zero:" " undefined result")

    def __deepcopy__(self, memo):
        """
        Hook for the standard copy module.

        The returned AffineScalarFunc is a completely fresh copy,
        which is fully independent of any variable defined so far.
        New variables are specially created for the returned
        AffineScalarFunc object.
        """
        return AffineScalarFunc(self._nominal_value, copy.deepcopy(self._linear_part))

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
            all_attrs["__dict__"] = self.__dict__
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
            slot_names = getattr(cls, "__slots__", ())

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
                # all_slots - {'__dict__'}:
                all_attrs[name] = getattr(self, name)
            except AttributeError:
                pass  # Undefined slot attribute

        return all_attrs

    def __setstate__(self, data_dict):
        """
        Hook for the pickle module.
        """
        for name, value in data_dict.items():
            # Contrary to the default __setstate__(), this does not
            # necessarily save to the instance dictionary (because the
            # instance might contain slots):
            setattr(self, name, value)


ops.add_arithmetic_ops(AffineScalarFunc)
ops.add_comparative_ops(AffineScalarFunc)
to_affine_scalar = AffineScalarFunc._to_affine_scalar

# Nicer name, for users: isinstance(ufloat(...), UFloat) is
# True. Also: isinstance(..., UFloat) is the test for "is this a
# number with uncertainties from the uncertainties package?":
UFloat = AffineScalarFunc


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
    return _wrap(
        AffineScalarFunc,
        f,
        derivatives_args=derivatives_args,
        derivatives_kwargs=derivatives_kwargs,
    )


###############################################################################


class NegativeStdDev(Exception):
    """Raise for a negative standard deviation"""

    pass


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
        return 0.0


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
    for i1, expr1 in enumerate(nums_with_uncert, 1):
        expanded_dict_1 = expr1.uncertainty.expanded_dict
        uatoms_1 = set(expanded_dict_1)
        coefs_expr1 = []

        for expr2 in nums_with_uncert[:i1]:
            expanded_dict_2 = expr2.uncertainty.expanded_dict
            uatom_2 = set(expanded_dict_2)
            coefs_expr1.append(
                sum(
                    (
                        (expanded_dict_1[uatom] * expanded_dict_2[uatom])
                        # uatom is common to both numbers with uncertainties:
                        for uatom in uatoms_1.intersection(uatom_2)
                    ),
                    # The result is always a float (sum() with no terms
                    # returns an integer):
                    0.0,
                )
            )

        covariance_matrix.append(coefs_expr1)

    # We symmetrize the matrix:
    for i, covariance_coefs in enumerate(covariance_matrix):
        covariance_coefs.extend(
            [covariance_matrix[j][i] for j in range(i + 1, len(covariance_matrix))]
        )

    return covariance_matrix


try:
    import numpy
except ImportError:
    pass
else:

    def correlation_matrix(nums_with_uncert):
        """
        Return the correlation matrix of the given sequence of
        numbers with uncertainties, as a NumPy array of floats.
        """

        cov_mat = numpy.array(covariance_matrix(nums_with_uncert))

        std_devs = numpy.sqrt(cov_mat.diagonal())

        return cov_mat / std_devs / std_devs[numpy.newaxis].T

    __all__.append("correlation_matrix")


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

    return AffineScalarFunc(nominal_value, std_dev, tag=tag)
