import inspect
import itertools
import json
import math
from math import isnan
from pathlib import Path

import pytest

from helpers import get_single_uatom
from uncertainties import ufloat
import uncertainties.core as uncert_core
import uncertainties.umath_core as umath_core
from uncertainties import umath

from helpers import (
    numbers_close,
    ufloats_close,
)
###############################################################################
# Unit tests


single_input_data_path = Path(Path(__file__).parent, "cases", "single_inputs.json")
with open(single_input_data_path, "r") as f:
    single_input_dict = json.load(f)

double_input_data_path = Path(Path(__file__).parent, "cases", "double_inputs.json")
with open(double_input_data_path, "r") as f:
    double_input_dict = json.load(f)

real_single_input_funcs = (
    "asinh",
    "atan",
    "cos",
    "cosh",
    "degrees",
    "erf",
    "erfc",
    "exp",
    "expm1",
    "fabs",
    "gamma",
    "lgamma",
    "radians",
    "sin",
    "sinh",
    "tan",
    "tanh",
)
positive_single_input_funcs = (
    "log",
    "log1p",
    "log10",
    "sqrt",
)
minus_one_to_plus_one_single_input_funcs = (
    "acos",
    "asin",
    "atanh",
)
greater_than_one_single_input_funcs = ("acosh",)

real_single_input_cases = list(
    itertools.product(real_single_input_funcs, single_input_dict["real"])
)
positive_single_input_cases = list(
    itertools.product(positive_single_input_funcs, single_input_dict["positive"])
)
minus_one_to_plus_one_single_input_cases = list(
    itertools.product(
        minus_one_to_plus_one_single_input_funcs,
        single_input_dict["minus_one_to_plus_one"],
    )
)
greater_than_one_single_input_cases = list(
    itertools.product(
        greater_than_one_single_input_funcs,
        single_input_dict["greater_than_one"],
    )
)
single_input_cases = (
    real_single_input_cases
    + positive_single_input_cases
    + minus_one_to_plus_one_single_input_cases
    + greater_than_one_single_input_cases
)


@pytest.mark.parametrize("func, value", single_input_cases)
def test_single_input_func_derivatives(func, value):
    uval = ufloat(value, 1.0)
    uatom = get_single_uatom(uval)

    umath_func = getattr(umath, func)
    result = umath_func(uval)
    analytic_deriv = result.error_components[uatom]

    math_func = getattr(math, func)
    numerical_deriv = uncert_core.partial_derivative(math_func, 0)(value)

    assert math.isclose(analytic_deriv, numerical_deriv, rel_tol=1e-6, abs_tol=1e-6)


real_double_input_funcs = (
    "atan2",
    "copysign",
    "fmod",
    "pow",
)

positive_double_input_funcs = (
    "hypot",
    "log",
)

real_double_input_cases = list(
    itertools.product(real_double_input_funcs, *zip(*double_input_dict["real"]))
)
print(real_double_input_cases)
positive_double_input_cases = list(
    itertools.product(positive_double_input_funcs, *zip(*double_input_dict["positive"]))
)

double_input_cases = real_double_input_cases + positive_double_input_cases


@pytest.mark.parametrize("func, value_0, value_1", double_input_cases)
def test_double_input_func_derivatives(func, value_0, value_1):
    if func == "pow":
        value_0 = abs(value_0)

    uval_0 = ufloat(value_0, 1.0)
    uval_1 = ufloat(value_1, 1.0)

    uatom_0 = get_single_uatom(uval_0)
    uatom_1 = get_single_uatom(uval_1)

    math_func = getattr(math, func)
    umath_func = getattr(umath, func)

    result = umath_func(uval_0, uval_1)

    analytic_deriv_0 = result.error_components[uatom_0]
    analytic_deriv_1 = result.error_components[uatom_1]

    numerical_deriv_0 = uncert_core.partial_derivative(
        math_func,
        0,
    )(value_0, value_1)
    numerical_deriv_1 = uncert_core.partial_derivative(
        math_func,
        1,
    )(value_0, value_1)

    assert math.isclose(analytic_deriv_0, numerical_deriv_0, rel_tol=1e-6, abs_tol=1e-6)
    assert math.isclose(analytic_deriv_1, numerical_deriv_1, rel_tol=1e-6, abs_tol=1e-6)


def test_compound_expression():
    """
    Test equality between different formulas.
    """

    x = ufloat(3, 0.1)
    y1 = umath_core.tan(x)
    y2 = umath_core.sin(x) / umath_core.cos(x)
    assert ufloats_close(y1, y2)


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
    x_uatom = get_single_uatom(x)
    y_uatom = get_single_uatom(y)
    # Derivatives that cannot be calculated simply return NaN, with no
    # exception being raised, normally:
    result = umath_core.hypot(x, y)
    assert isnan(result.error_components[x_uatom])
    assert isnan(result.error_components[y_uatom])


@pytest.mark.parametrize("function_name", umath_core.deprecated_functions)
def test_deprecated_function(function_name):
    num_args = len(inspect.signature(getattr(math, function_name)).parameters)
    args = [ufloat(1, 0.1)]
    if num_args == 1:
        if function_name == "factorial":
            args[0] = 6
    else:
        if function_name == "ldexp":
            args.append(3)
        else:
            args.append(ufloat(-12, 2.4))
    with pytest.warns(FutureWarning, match="will be removed"):
        getattr(umath_core, function_name)(*args)
