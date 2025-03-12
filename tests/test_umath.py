import math
from math import isnan

import pytest

from uncertainties import ufloat
import uncertainties.core as uncert_core
import uncertainties.umath_core as umath_core

from helpers import (
    power_derivative_cases,
    power_float_result_cases,
    power_reference_cases,
    compare_derivatives,
    nan_close,
    numbers_close,
)
###############################################################################
# Unit tests


def test_fixed_derivatives_math_funcs():
    """
    Comparison between function derivatives and numerical derivatives.

    This comparison is useful for derivatives that are analytical.
    """

    for name in umath_core.many_scalars_to_scalar_funcs:
        func = getattr(umath_core, name)
        # Numerical derivatives of func: the nominal value of func() results
        # is used as the underlying function:
        numerical_derivatives = uncert_core.NumericalDerivatives(
            lambda *args: func(*args)
        )
        compare_derivatives(func, numerical_derivatives)

    # Functions that are not in umath_core.many_scalars_to_scalar_funcs:

    ##
    # modf(): returns a tuple:
    def frac_part_modf(x):
        return umath_core.modf(x)[0]

    def int_part_modf(x):
        return umath_core.modf(x)[1]

    compare_derivatives(
        frac_part_modf, uncert_core.NumericalDerivatives(lambda x: frac_part_modf(x))
    )
    compare_derivatives(
        int_part_modf, uncert_core.NumericalDerivatives(lambda x: int_part_modf(x))
    )

    ##
    # frexp(): returns a tuple:
    def mantissa_frexp(x):
        return umath_core.frexp(x)[0]

    def exponent_frexp(x):
        return umath_core.frexp(x)[1]

    compare_derivatives(
        mantissa_frexp, uncert_core.NumericalDerivatives(lambda x: mantissa_frexp(x))
    )
    compare_derivatives(
        exponent_frexp, uncert_core.NumericalDerivatives(lambda x: exponent_frexp(x))
    )


def test_compound_expression():
    """
    Test equality between different formulas.
    """

    x = ufloat(3, 0.1)

    # Prone to numerical errors (but not much more than floats):
    assert umath_core.tan(x) == umath_core.sin(x) / umath_core.cos(x)


def test_numerical_example():
    "Test specific numerical examples"

    x = ufloat(3.14, 0.01)
    result = umath_core.sin(x)
    # In order to prevent big errors such as a wrong, constant value
    # for all analytical and numerical derivatives, which would make
    # test_fixed_derivatives_math_funcs() succeed despite incorrect
    # calculations:
    assert (
        "%.6f +/- %.6f" % (result.nominal_value, result.std_dev)
        == "0.001593 +/- 0.010000"
    )

    # Regular calculations should still work:
    assert "%.11f" % umath_core.sin(3) == "0.14112000806"


def test_monte_carlo_comparison():
    """
    Full comparison to a Monte-Carlo calculation.

    Both the nominal values and the covariances are compared between
    the direct calculation performed in this module and a Monte-Carlo
    simulation.
    """

    try:
        import numpy
        import numpy.random
    except ImportError:
        import warnings

        warnings.warn("Test not performed because NumPy is not available")
        return

    # Works on numpy.arrays of Variable objects (whereas umath_core.sin()
    # does not):
    sin_uarray_uncert = numpy.vectorize(umath_core.sin, otypes=[object])

    # Example expression (with correlations, and multiple variables combined
    # in a non-linear way):
    def function(x, y):
        """
        Function that takes two NumPy arrays of the same size.
        """
        # The uncertainty due to x is about equal to the uncertainty
        # due to y:
        return 10 * x**2 - x * sin_uarray_uncert(y**3)

    x = ufloat(0.2, 0.01)
    y = ufloat(10, 0.001)
    function_result_this_module = function(x, y)
    nominal_value_this_module = function_result_this_module.nominal_value

    # Covariances "f*f", "f*x", "f*y":
    covariances_this_module = numpy.array(
        uncert_core.covariance_matrix((x, y, function_result_this_module))
    )

    def monte_carlo_calc(n_samples):
        """
        Calculate function(x, y) on n_samples samples and returns the
        median, and the covariances between (x, y, function(x, y)).
        """
        # Result of a Monte-Carlo simulation:
        x_samples = numpy.random.normal(x.nominal_value, x.std_dev, n_samples)
        y_samples = numpy.random.normal(y.nominal_value, y.std_dev, n_samples)

        # !! astype() is a fix for median() in NumPy 1.8.0:
        function_samples = function(x_samples, y_samples).astype(float)

        cov_mat = numpy.cov([x_samples, y_samples], function_samples)

        return (numpy.median(function_samples), cov_mat)

    (nominal_value_samples, covariances_samples) = monte_carlo_calc(1000000)

    ## Comparison between both results:

    # The covariance matrices must be close:

    # We rely on the fact that covariances_samples very rarely has
    # null elements:

    # !!! The test could be done directly with NumPy's comparison
    # tools, no? See assert_allclose, assert_array_almost_equal_nulp
    # or assert_array_max_ulp. This is relevant for all vectorized
    # occurrences of numbers_close.

    assert numpy.vectorize(numbers_close)(
        covariances_this_module, covariances_samples, 0.06
    ).all(), (
        "The covariance matrices do not coincide between"
        " the Monte-Carlo simulation and the direct calculation:\n"
        "* Monte-Carlo:\n%s\n* Direct calculation:\n%s"
        % (covariances_samples, covariances_this_module)
    )

    # The nominal values must be close:
    assert numbers_close(
        nominal_value_this_module,
        nominal_value_samples,
        # The scale of the comparison depends on the standard
        # deviation: the nominal values can differ by a fraction of
        # the standard deviation:
        math.sqrt(covariances_samples[2, 2]) / abs(nominal_value_samples) * 0.5,
    ), (
        "The nominal value (%f) does not coincide with that of"
        " the Monte-Carlo simulation (%f), for a standard deviation of %f."
        % (
            nominal_value_this_module,
            nominal_value_samples,
            math.sqrt(covariances_samples[2, 2]),
        )
    )


def test_math_module():
    "Operations with the math module"

    x = ufloat(-1.5, 0.1)

    # The exponent must not be differentiated, when calculating the
    # following (the partial derivative with respect to the exponent
    # is not defined):
    assert (x**2).nominal_value == 2.25

    # Regular operations are chosen to be unchanged:
    assert isinstance(umath_core.sin(3), float)

    # factorial() must not be "damaged" by the umath_core module, so as
    # to help make it a drop-in replacement for math (even though
    # factorial() does not work on numbers with uncertainties
    # because it is restricted to integers, as for
    # math.factorial()):
    assert umath_core.factorial(4) == 24

    # fsum is special because it does not take a fixed number of
    # variables:
    assert umath_core.fsum([x, x]).nominal_value == -3

    # Functions that give locally constant results are tested: they
    # should give the same result as their float equivalent:
    for name in umath_core.locally_cst_funcs:
        try:
            func = getattr(umath_core, name)
        except AttributeError:
            continue  # Not in the math module, so not in umath_core either

        assert func(x) == func(x.nominal_value)
        # The type should be left untouched. For example, isnan()
        # should always give a boolean:
        assert isinstance(func(x), type(func(x.nominal_value)))

    # The same exceptions should be generated when numbers with uncertainties
    # are used:

    # The type of the expected exception is first determined, because
    # it varies between versions of Python (OverflowError in Python
    # 2.6+, ValueError in Python 2.5,...):
    try:
        math.log(0)
    except Exception as err_math:
        # Python 3 does not make exceptions local variables: they are
        # restricted to their except block:
        err_math_args = err_math.args
        exception_class = err_math.__class__

    try:
        umath_core.log(0)
    except exception_class as err_ufloat:
        assert err_math_args == err_ufloat.args
    else:
        raise Exception("%s exception expected" % exception_class.__name__)
    try:
        umath_core.log(ufloat(0, 0))
    except exception_class as err_ufloat:
        assert err_math_args == err_ufloat.args
    else:
        raise Exception("%s exception expected" % exception_class.__name__)
    try:
        umath_core.log(ufloat(0, 1))
    except exception_class as err_ufloat:
        assert err_math_args == err_ufloat.args
    else:
        raise Exception("%s exception expected" % exception_class.__name__)


def test_hypot():
    """
    Special cases where derivatives cannot be calculated:
    """
    x = ufloat(0, 1)
    y = ufloat(0, 2)
    # Derivatives that cannot be calculated simply return NaN, with no
    # exception being raised, normally:
    result = umath_core.hypot(x, y)
    assert isnan(result.derivatives[x])
    assert isnan(result.derivatives[y])


@pytest.mark.parametrize(
    "first_ufloat, second_ufloat, first_der, second_der",
    power_derivative_cases,
)
def test_power_derivatives(first_ufloat, second_ufloat, first_der, second_der):
    result = umath_core.pow(first_ufloat, second_ufloat)
    first_der_result = result.derivatives[first_ufloat]
    second_der_result = result.derivatives[second_ufloat]
    assert nan_close(first_der_result, first_der)
    assert nan_close(second_der_result, second_der)


@pytest.mark.parametrize(
    "first_ufloat, second_ufloat, result_float",
    power_float_result_cases,
)
def test_power_float_result_cases(first_ufloat, second_ufloat, result_float):
    assert umath_core.pow(first_ufloat, second_ufloat) == result_float


zero = ufloat(0, 0)
positive = ufloat(0.3, 0.01)
negative = ufloat(-0.3, 0.01)
"""
math.pow raises ValueError in these cases, in contrast to pow which raises
ZeroDivisionError so these test cases are slightly different than those that appear for
test_power_exceptions in test_uncertainties.py.
"""
power_exception_cases = [
    (ufloat(0, 0), negative, ValueError),
    (ufloat(0, 0.1), negative, ValueError),
    (negative, positive, ValueError),
]


@pytest.mark.parametrize("first_ufloat, second_ufloat, exc_type", power_exception_cases)
def test_power_exceptions(first_ufloat, second_ufloat, exc_type):
    with pytest.raises(exc_type):
        umath_core.pow(first_ufloat, second_ufloat)


@pytest.mark.parametrize("first_ufloat, second_float", power_reference_cases)
def test_power_wrt_ref(first_ufloat, second_float):
    expected_result = math.pow(first_ufloat.n, second_float)
    actual_result = umath_core.pow(first_ufloat, second_float).n
    assert actual_result == expected_result
