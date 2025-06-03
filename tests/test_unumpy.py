try:
    import numpy
except ImportError:
    import sys

    sys.exit()  # There is no reason to test the interface to NumPy

import uncertainties.umath
import uncertainties.core as uncert_core
from uncertainties import ufloat, unumpy
from uncertainties.unumpy import core
from helpers import nan_close, uarrays_close


def test_numpy():
    """
    Interaction with NumPy, including matrix inversion,
    correlated_values, and calculation of the mean.
    """

    arr = numpy.arange(3)
    num = ufloat(3.14, 0.01)

    # NumPy arrays can be multiplied by Variable objects,
    # whatever the order of the operands:
    prod1 = arr * num
    prod2 = num * arr
    # Additional check:
    assert (prod1 == prod2).all()

    # Operations with arrays work (they are first handled by NumPy,
    # then by this module):
    prod1 * prod2  # This should be calculable

    # The following does not completely work, because NumPy does not
    # implement numpy.exp on an array of general objects, apparently:
    assert numpy.exp(arr).all()  # All elements > 0
    # Equivalent with an array of AffineScalarFunc objects:
    try:
        numpy.exp(arr + ufloat(0, 0))
    except (AttributeError, TypeError):
        # In numpy<1.17, an AttributeError is raised in this situation. This was
        # considered a bug however, and in numpy 1.17 it was changed to a
        # TypeError (see PR #12700 in numpy repository)
        pass
    else:
        raise Exception("numpy.exp unexpectedly worked")

    # Calculation of the mean, global and with a specific axis:

    arr_floats = numpy.random.random((10, 3, 5))
    arr = unumpy.uarray(arr_floats, arr_floats / 100)
    assert arr.mean(axis=0).shape == (3, 5)
    assert arr.mean(axis=1).shape == (10, 5)
    arr.mean()  # Global mean


def derivatives_close(x, y):
    """
    Returns True iff the AffineScalarFunc objects x and y have
    derivatives that are close to each other (they must depend
    on the same variables).
    """

    # x and y must depend on the same variables:
    if set(x.derivatives) != set(y.derivatives):
        return False  # Not the same variables

    return all(
        nan_close(x.derivatives[var], y.derivatives[var]) for var in x.derivatives
    )


def test_inverse():
    "Tests of the matrix inverse"

    m = numpy.array([[ufloat(10, 1), -3.1], [0, ufloat(3, 0)]])
    m_nominal_values = unumpy.nominal_values(m)

    # "Regular" inverse matrix, when uncertainties are not taken
    # into account:
    m_no_uncert_inv = numpy.linalg.inv(m_nominal_values)

    # The matrix inversion should not yield numbers with uncertainties:
    assert m_no_uncert_inv.dtype == numpy.dtype(float)

    # Inverse with uncertainties:
    m_inv_uncert = core.inv(m)  # AffineScalarFunc elements
    # The inverse contains uncertainties: it must support custom
    # operations on matrices with uncertainties:
    assert isinstance(m_inv_uncert, numpy.ndarray)
    assert type(m_inv_uncert[0, 0]) == uncert_core.AffineScalarFunc

    # Checks of the numerical values: the diagonal elements of the
    # inverse should be the inverses of the diagonal elements of
    # m (because we started with a triangular matrix):
    assert nan_close(
        1 / m_nominal_values[0, 0], core.nominal_values(m_inv_uncert[0, 0])
    ), "Wrong value"

    assert nan_close(
        1 / m_nominal_values[1, 1], core.nominal_values(m_inv_uncert[1, 1])
    ), "Wrong value"

    ####################

    # Checks of the covariances between elements:
    x = ufloat(10, 1)
    m = numpy.array([[x, x], [0, 3 + 2 * x]])

    m_inverse = core.inv(m)

    # Check of the properties of the inverse:
    m_double_inverse = core.inv(m_inverse)
    # The initial matrix should be recovered, including its
    # derivatives, which define covariances:
    assert nan_close(m_double_inverse[0, 0].nominal_value, core.nominal_values(m[0, 0]))
    assert nan_close(m_double_inverse[0, 0].std_dev, core.std_devs(m[0, 0]))

    assert uarrays_close(m_double_inverse, m).all()

    # Partial test:
    assert derivatives_close(m_double_inverse[0, 0], m[0, 0])
    assert derivatives_close(m_double_inverse[1, 1], m[1, 1])

    ####################

    # Tests of covariances during the inversion:

    # There are correlations if both the next two derivatives are
    # not zero:
    assert m_inverse[0, 0].derivatives[x]
    assert m_inverse[0, 1].derivatives[x]

    # Correlations between m and m_inverse should create a perfect
    # inversion:
    assert uarrays_close(m @ m_inverse, numpy.eye(m.shape[0])).all()


def test_wrap_array_func():
    """
    Test of numpy.wrap_array_func(), with optional arguments and
    keyword arguments.
    """

    # Function that works with numbers with uncertainties in mat (if
    # mat is an uncertainties.unumpy.matrix):
    def f_unc(mat, *args, **kwargs):
        return core.pinv(mat) + args[0] * kwargs["factor"]

    # Test with optional arguments and keyword arguments:
    def f(mat, *args, **kwargs):
        # This function is wrapped: it should only be called with pure
        # numbers:
        assert not any(isinstance(v, uncert_core.UFloat) for v in mat.flat)
        return f_unc(mat, *args, **kwargs)

    # Wrapped function:
    f_wrapped = core.wrap_array_func(f)

    ##########
    # Full rank rectangular matrix:
    m = numpy.array([[ufloat(10, 1), -3.1], [0, ufloat(3, 0)], [1, -3.1]])

    # Numerical and package (analytical) pseudo-inverses: they must be
    # the same:
    m_f_wrapped = f_wrapped(m, 2, factor=10)
    m_f_unc = f_unc(m, 2, factor=10)

    assert uarrays_close(m_f_wrapped, m_f_unc).all()


def test_pseudo_inverse():
    "Tests of the pseudo-inverse"

    # Numerical version of the pseudo-inverse:
    pinv_num = core.wrap_array_func(numpy.linalg.pinv)

    ##########
    # Full rank rectangular matrix:
    m = numpy.array([[ufloat(10, 1), -3.1], [0, ufloat(3, 0)], [1, -3.1]])

    # Numerical and package (analytical) pseudo-inverses: they must be
    # the same:
    rcond = 1e-8  # Test of the second argument to pinv()
    m_pinv_num = pinv_num(m, rcond)
    m_pinv_package = core.pinv(m, rcond)
    assert uarrays_close(m_pinv_num, m_pinv_package).all()

    ##########
    # Example with a non-full rank rectangular matrix:
    vector = [ufloat(10, 1), -3.1, 11]
    m = numpy.array([vector, vector])
    m_pinv_num = pinv_num(m, rcond)
    m_pinv_package = core.pinv(m, rcond)
    assert uarrays_close(m_pinv_num, m_pinv_package).all()

    ##########
    # Example with a non-full-rank square matrix:
    m = numpy.array([[ufloat(10, 1), 0], [3, 0]])
    m_pinv_num = pinv_num(m, rcond)
    m_pinv_package = core.pinv(m, rcond)
    assert uarrays_close(m_pinv_num, m_pinv_package).all()


def test_broadcast_funcs():
    """
    Test of mathematical functions that work with NumPy arrays of
    numbers with uncertainties.
    """

    x = ufloat(0.2, 0.1)
    arr = numpy.array([x, 2 * x])
    assert unumpy.cos(arr)[1] == uncertainties.umath.cos(arr[1])

    # Some functions do not bear the same name in the math module and
    # in NumPy (acos instead of arccos, etc.):
    assert unumpy.arccos(arr)[1] == uncertainties.umath.acos(arr[1])

    # The acos() function should not exist in unumpy because the function
    # should have been renamed to arccos(). Starting with numpy 2 numpy.acos()
    # is an alias to numpy.arccos(). If similar aliases are added to unumpy,
    # the following tests can be removed.
    assert not hasattr(unumpy, "acos")

    # Test of the __all__ variable:
    assert "acos" not in unumpy.__all__


def test_array_creation():
    "Test of custom array creation"

    arr = unumpy.uarray([1, 2], [0.1, 0.2])

    assert core.nominal_values(arr)[1] == 2
    assert core.std_devs(arr)[1] == 0.2


def test_component_extraction():
    "Extracting the nominal values and standard deviations from an array"

    arr = unumpy.uarray([1, 2], [0.1, 0.2])

    assert numpy.all(unumpy.nominal_values(arr) == [1, 2])
    assert numpy.all(unumpy.std_devs(arr) == [0.1, 0.2])


def test_array_comparisons():
    "Test of array and matrix comparisons"

    arr = unumpy.uarray([1, 2], [1, 4])
    assert numpy.all((arr == [arr[0], 4]) == [True, False])
