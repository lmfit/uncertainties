from math import isclose, isnan, isinf

import uncertainties.core as uncert_core


def nan_close(first, second):
    if isnan(first):
        return isnan(second)
    else:
        return isclose(first, second)


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
