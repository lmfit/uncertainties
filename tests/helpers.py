from math import isclose, isnan


def nan_close(first, second, *, rel_tol=1e-9, abs_tol=0.0):
    if isnan(first):
        return isnan(second)
    else:
        return isclose(first, second, rel_tol=rel_tol, abs_tol=abs_tol)


###############################################################
# TODO: move to uncertainties/testing.py
###############################################################################

# Utilities for unit testing


def ufloat_nan_close(x, y, tolerance=1e-6):
    """
    Tests if two numbers with uncertainties are close, as random
    variables: this is stronger than testing whether their nominal
    value and standard deviation are close.

    The tolerance is applied to both the nominal value and the
    standard deviation of the difference between the numbers.
    """

    diff = x - y
    nominal_values_close = nan_close(
        diff.n,
        0,
        rel_tol=tolerance,
        abs_tol=tolerance,
    )
    std_devs_close = nan_close(
        diff.s,
        0,
        rel_tol=tolerance,
        abs_tol=tolerance,
    )
    return nominal_values_close and std_devs_close


###############################################################################


try:
    import numpy as np  # noqa
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
        diff_arr = m1 - m2
        nominal_values_arr = np.vectorize(lambda x: x.n)(diff_arr)
        std_devs_arr = np.vectorize(lambda x: x.s)(diff_arr)
        nominal_values_close_arr = np.isclose(
            nominal_values_arr,
            0,
            rtol=precision,
            atol=precision,
        )
        std_devs_close_arr = np.isclose(
            std_devs_arr,
            0,
            rtol=precision,
            atol=precision,
        )
        return np.logical_and(nominal_values_close_arr, std_devs_close_arr)
