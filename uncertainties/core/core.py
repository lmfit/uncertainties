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


from . compat import (deprecation, FLOAT_LIKE_TYPES, CONSTANT_TYPES,
                      CallableStdDev,
                      )
from . util import (set_doc, covariance_matrix, partial_derivative,
                    NumericalDerivatives, IndexableIter, basestring,
                    )
from . parsing import (str_to_number_with_uncert,)
from . affinescalarfunc import (AffineScalarFunc, LinearCombination, 
                                AffineScalarFuncBase
                                )
from .affinescalarfunc.ops import nan_if_exception
###############################################################################

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
        norm_vector[norm_vector==0] = 1

        return correlated_values_norm(
            # !! The following zip() is a bit suboptimal: correlated_values()
            # separates back the nominal values and the standard deviations:
            list(zip(nom_values, std_devs)),
            covariance_mat/norm_vector/norm_vector[:,numpy.newaxis],
            tags)

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
        
    def correlation_matrix(nums_with_uncert):
        '''
        Return the correlation matrix of the given sequence of
        numbers with uncertainties, as a NumPy array of floats.
        '''

        cov_mat = numpy.array(covariance_matrix(nums_with_uncert))

        std_devs = numpy.sqrt(cov_mat.diagonal())

        return cov_mat/std_devs/std_devs[numpy.newaxis].T


# Nicer name, for users: isinstance(ufloat(...), UFloat) is
# True. Also: isinstance(..., UFloat) is the test for "is this a
# number with uncertainties from the uncertainties package?":
UFloat = AffineScalarFunc
to_affine_scalar = AffineScalarFunc._to_affine_scalar
wrap = UFloat.wrap

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

    if isinstance(x, AffineScalarFuncBase):
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

    if isinstance(x, AffineScalarFuncBase):
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
