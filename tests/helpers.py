from math import isnan, isinf

import uncertainties.core as uncert_core
from uncertainties.core import ufloat, UFloat


def get_single_uatom(num_with_uncertainty: UFloat):
    error_components = num_with_uncertainty.error_components
    if len(error_components) > 1:
        raise ValueError("UFloat has more than one error component.")
    return next(iter(error_components.keys()))


def power_all_cases(op):
    """
    Checks all cases for the value and derivatives of power-like
    operator op (op is typically the built-in pow(), or math.pow()).

    Checks only the details of special results like 0, 1 or NaN).

    Different cases for the value of x**p and its derivatives are
    tested by dividing the (x, p) plane with:

    - x < 0, x = 0, x > 0
    - p integer or not, p < 0, p = 0, p > 0

    (not all combinations are distinct: for instance x > 0 gives
    identical formulas for all p).
    """

    zero = ufloat(0, 0.1)
    zero2 = ufloat(0, 0.1)
    one = ufloat(1, 0.1)
    positive = ufloat(0.3, 0.01)
    positive2 = ufloat(0.3, 0.01)
    negative = ufloat(-0.3, 0.01)
    integer = ufloat(-3, 0)
    non_int_larger_than_one = ufloat(3.1, 0.01)
    positive_smaller_than_one = ufloat(0.3, 0.01)

    assert integer.uncertainty.ucombo_tuple == ()

    negative_uatom = get_single_uatom(negative)
    positive_uatom = get_single_uatom(positive)
    positive2_uatom = get_single_uatom(positive2)
    one_uatom = get_single_uatom(one)
    zero_uatom = get_single_uatom(zero)
    zero2_uatom = get_single_uatom(zero2)
    non_int_larger_than_one_uatom = get_single_uatom(non_int_larger_than_one)
    positive_smaller_than_one_uatom = get_single_uatom(positive_smaller_than_one)
    ## negative**integer

    result = op(negative, integer)
    assert not isnan(result.error_components[negative_uatom])

    # Limit cases:
    result = op(negative, one)
    assert result.error_components[negative_uatom] == negative.std_dev
    assert isnan(result.error_components[one_uatom])

    result = op(negative, zero)
    assert result.error_components[negative_uatom] == 0
    assert isnan(result.error_components[zero_uatom])

    ## negative**non-integer

    ## zero**...

    result = op(zero, non_int_larger_than_one)
    assert isnan(result.error_components[zero_uatom])
    assert result.error_components[non_int_larger_than_one_uatom] == 0

    # Special cases:
    result = op(zero, one)
    assert result.error_components[zero_uatom] == zero.std_dev
    assert result.error_components[one_uatom] == 0

    result = op(zero, 2 * one)
    assert result.error_components[zero_uatom] == 0
    assert result.error_components[one_uatom] == 0

    result = op(zero, positive_smaller_than_one)
    assert isnan(result.error_components[zero_uatom])
    assert result.error_components[positive_smaller_than_one_uatom] == 0

    result = op(zero, zero2)
    assert result.error_components[zero_uatom] == 0
    assert isnan(result.error_components[zero2_uatom])

    ## positive**...: this is a quite regular case where the value and
    ## the error_components are all defined.

    result = op(positive, positive2)
    assert not isnan(result.error_components[positive_uatom])
    assert not isnan(result.error_components[positive2_uatom])

    result = op(positive, zero)
    assert result.error_components[positive_uatom] == 0
    assert not isnan(result.error_components[zero_uatom])

    result = op(positive, negative)
    assert not isnan(result.error_components[positive_uatom])
    assert not isnan(result.error_components[negative_uatom])


def power_special_cases(op):
    """
    Checks special cases of the uncertainty power operator op (where
    op is typically the built-in pow or uncertainties.umath.pow).

    The values x = 0, x = 1 and x = NaN are special, as are null,
    integral and NaN values of p.
    """

    zero = ufloat(0, 0)
    one = ufloat(1, 0)
    p = ufloat(0.3, 0.01)

    assert op(0, p) == 0
    assert op(zero, p) == 0

    # The outcome of 1**nan and nan**0 was undefined before Python
    # 2.6 (http://docs.python.org/library/math.html#math.pow):
    assert op(float("nan"), zero) == 1.0
    assert op(one, float("nan")) == 1.0

    # …**0 == 1.0:
    assert op(p, 0) == 1.0
    assert op(zero, 0) == 1.0
    assert op((-p), 0) == 1.0
    # …**zero:
    assert op((-10.3), zero) == 1.0
    assert op(0, zero) == 1.0
    assert op(0.3, zero) == 1.0
    assert op((-p), zero) == 1.0
    assert op(zero, zero) == 1.0
    assert op(p, zero) == 1.0

    # one**… == 1.0
    assert op(one, -3) == 1.0
    assert op(one, -3.1) == 1.0
    assert op(one, 0) == 1.0
    assert op(one, 3) == 1.0
    assert op(one, 3.1) == 1.0

    # … with two numbers with uncertainties:
    assert op(one, (-p)) == 1.0
    assert op(one, zero) == 1.0
    assert op(one, p) == 1.0
    # 1**… == 1.0:
    assert op(1.0, (-p)) == 1.0
    assert op(1.0, zero) == 1.0
    assert op(1.0, p) == 1.0


def power_wrt_ref(op, ref_op):
    """
    Checks special cases of the uncertainty power operator op (where
    op is typically the built-in pow or uncertainties.umath.pow), by
    comparing its results to the reference power operator ref_op
    (which is typically the built-in pow or math.pow).
    """

    # Negative numbers with uncertainty can be exponentiated to an
    # integral power:
    assert op(ufloat(-1.1, 0.1), -9).nominal_value == ref_op(-1.1, -9)

    # Case of numbers with no uncertainty: should give the same result
    # as numbers with uncertainties:
    assert op(ufloat(-1, 0), 9) == ref_op(-1, 9)
    assert op(ufloat(-1.1, 0), 9) == ref_op(-1.1, 9)


###############################################################
# TODO: move to uncertainties/testing.py
###############################################################################

# Utilities for unit testing


def numbers_close(x, y, tolerance=1e-6):
    """
    Returns True if the given floats are close enough.

    The given tolerance is the relative difference allowed, or the absolute
    difference, if one of the numbers is 0.

    NaN is allowed: it is considered close to itself.
    """

    # !!! Python 3.5+ has math.isclose(): maybe it could be used here.

    # Instead of using a try and ZeroDivisionError, we do a test,
    # NaN could appear silently:

    if x != 0 and y != 0:
        if isinf(x):
            return isinf(y)
        elif isnan(x):
            return isnan(y)
        else:
            # Symmetric form of the test:
            return 2 * abs(x - y) / (abs(x) + abs(y)) < tolerance

    else:  # Either x or y is zero
        return abs(x or y) < tolerance


def ufloats_close(x, y, tolerance=1e-6):
    """
    Tests if two numbers with uncertainties are close, as random
    variables: this is stronger than testing whether their nominal
    value and standard deviation are close.

    The tolerance is applied to both the nominal value and the
    standard deviation of the difference between the numbers.
    """

    diff = x - y
    return numbers_close(diff.nominal_value, 0, tolerance) and numbers_close(
        diff.std_dev, 0, tolerance
    )


class DerivativesDiffer(Exception):
    pass


###############################################################################


try:
    import numpy  # noqa
except ImportError:
    pass
else:

    def uarrays_close(m1, m2, precision=1e-4):
        """
        Returns True iff m1 and m2 are almost equal, where elements
        can be either floats or AffineScalarFunc objects.

        Two independent AffineScalarFunc objects are deemed equal if
        both their nominal value and uncertainty are equal (up to the
        given precision).

        m1, m2 -- NumPy arrays.

        precision -- precision passed through to
        uncertainties.test_uncertainties.numbers_close().
        """

        # ! numpy.allclose() is similar to this function, but does not
        # work on arrays that contain numbers with uncertainties, because
        # of the isinf() function.

        for elmt1, elmt2 in zip(m1.flat, m2.flat):
            # For a simpler comparison, both elements are
            # converted to AffineScalarFunc objects:
            elmt1 = uncert_core.to_affine_scalar(elmt1)
            elmt2 = uncert_core.to_affine_scalar(elmt2)

            if not numbers_close(elmt1.nominal_value, elmt2.nominal_value, precision):
                return False

            if not numbers_close(elmt1.std_dev, elmt2.std_dev, precision):
                return False

        return True
