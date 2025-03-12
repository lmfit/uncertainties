from math import pow as math_pow

import pytest

from uncertainties import ufloat
from uncertainties.umath_core import pow as umath_pow

from helpers import nan_close


zero = ufloat(0, 0.1)
zero2 = ufloat(0, 0.1)
one = ufloat(1, 0.1)
two = ufloat(2, 0.2)
positive = ufloat(0.3, 0.01)
positive2 = ufloat(0.3, 0.01)
negative = ufloat(-0.3, 0.01)
integer = ufloat(-3, 0)
non_int_larger_than_one = ufloat(3.1, 0.01)
positive_smaller_than_one = ufloat(0.3, 0.01)


power_derivative_cases = (
    (negative, integer, -370.37037037037044, float("nan")),
    (negative, one, 1.0, float("nan")),
    (negative, zero, 0.0, float("nan")),
    (zero, non_int_larger_than_one, float("nan"), 0.0),
    (zero, one, 1.0, 0.0),
    (zero, two, 0.0, 0.0),
    (zero, positive_smaller_than_one, float("nan"), 0.0),
    (zero, zero2, 0.0, float("nan")),
    (positive, positive2, 0.696845301935949, -0.8389827923531782),
    (positive, zero, 0.0, -1.2039728043259361),
    (positive, negative, -1.4350387341664474, -1.7277476090907193),
)


@pytest.mark.parametrize(
    "first_ufloat, second_ufloat, first_der, second_der",
    power_derivative_cases,
)
def test_power_derivatives(first_ufloat, second_ufloat, first_der, second_der):
    result = pow(first_ufloat, second_ufloat)
    first_der_result = result.derivatives[first_ufloat]
    second_der_result = result.derivatives[second_ufloat]
    assert nan_close(first_der_result, first_der)
    assert nan_close(second_der_result, second_der)

    result = umath_pow(first_ufloat, second_ufloat)
    first_der_result = result.derivatives[first_ufloat]
    second_der_result = result.derivatives[second_ufloat]
    assert nan_close(first_der_result, first_der)
    assert nan_close(second_der_result, second_der)


zero = ufloat(0, 0)
one = ufloat(1, 0)
p = ufloat(0.3, 0.01)

power_float_result_cases = [
    (0, p, 0),
    (zero, p, 0),
    (float("nan"), zero, 1),
    (one, float("nan"), 1),
    (p, 0, 1),
    (zero, 0, 1),
    (-p, 0, 1),
    (-10.3, zero, 1),
    (0, zero, 1),
    (0.3, zero, 1),
    (-p, zero, 1),
    (zero, zero, 1),
    (p, zero, 1),
    (one, -3, 1),
    (one, -3.1, 1),
    (one, 0, 1),
    (one, 3, 1),
    (one, 3.1, 1),
    (one, -p, 1),
    (one, zero, 1),
    (one, p, 1),
    (1, -p, 1),
    (1, zero, 1),
    (1, p, 1),
]


@pytest.mark.parametrize(
    "first_ufloat, second_ufloat, result_float",
    power_float_result_cases,
)
def test_power_float_result_cases(first_ufloat, second_ufloat, result_float):
    assert pow(first_ufloat, second_ufloat) == result_float
    assert umath_pow(first_ufloat, second_ufloat) == result_float


zero = ufloat(0, 0)
positive = ufloat(0.3, 0.01)
negative = ufloat(-0.3, 0.01)


power_reference_cases = [
    (ufloat(-1.1, 0.1), -9),
    (ufloat(-1, 0), 9),
    (ufloat(-1.1, 0), 9),
]


@pytest.mark.parametrize("first_ufloat, second_float", power_reference_cases)
def test_power_wrt_ref(first_ufloat, second_float):
    assert pow(first_ufloat, second_float).n == pow(first_ufloat.n, second_float)
    assert umath_pow(first_ufloat, second_float).n == math_pow(
        first_ufloat.n, second_float
    )


zero = ufloat(0, 0)
positive = ufloat(0.3, 0.01)
negative = ufloat(-0.3, 0.01)
power_exception_cases = [
    (ufloat(0, 0), negative, ZeroDivisionError),
    (ufloat(0, 0.1), negative, ZeroDivisionError),
    (negative, positive, ValueError),
]


@pytest.mark.parametrize("first_ufloat, second_ufloat, exc_type", power_exception_cases)
def test_power_exceptions(first_ufloat, second_ufloat, exc_type):
    with pytest.raises(exc_type):
        pow(first_ufloat, second_ufloat)


"""
math.pow raises ValueError in these cases, in contrast to pow which raises
ZeroDivisionError so these test cases are slightly different than those that appear for
test_power_exceptions in test_uncertainties.py.
"""
umath_power_exception_cases = [
    (ufloat(0, 0), negative, ValueError),
    (ufloat(0, 0.1), negative, ValueError),
    (negative, positive, ValueError),
]


@pytest.mark.parametrize(
    "first_ufloat, second_ufloat, exc_type",
    umath_power_exception_cases,
)
def test_umath_power_exceptions(first_ufloat, second_ufloat, exc_type):
    with pytest.raises(exc_type):
        umath_pow(first_ufloat, second_ufloat)
