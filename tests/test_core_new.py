from math import sqrt, sin, cos

import pytest

from uncertainties.core_new import UFloat, ToUFuncPositional

repr_cases = cases = [
        (UFloat(10, 1), 'UFloat(10.0, 1.0)'),
        (UFloat(20, 2), 'UFloat(20.0, 2.0)'),
        (UFloat(30, 3), 'UFloat(30.0, 3.0)'),
        (UFloat(-30, 3), 'UFloat(-30.0, 3.0)'),
        (UFloat(-30, float('nan')), 'UFloat(-30.0, nan)'),
    ]


@pytest.mark.parametrize("unum, expected_repr_str", repr_cases)
def test_repr(unum: UFloat, expected_repr_str: str):
    assert repr(unum) == expected_repr_str


x = UFloat(10, 1)
unary_cases = [
    (-x, -10, 1),
    (+x, 10, 1),
    (abs(x), 10, 1),
    (abs(-x), 10, 1),
]


@pytest.mark.parametrize(
    "unum, expected_val, expected_std_dev",
    unary_cases,
)
def test_unary(
        unum: UFloat,
        expected_val: float,
        expected_std_dev: float,
):
    assert unum.val == expected_val
    assert unum.std_dev == expected_std_dev


x = UFloat(10, 1)
y = UFloat(20, 2)
binary_cases = [
    (x + 20, 30, 1),
    (x - 20, -10, 1),
    (x * 20, 200, 20),
    (x / 20, 0.5, 0.05),
    (20 + x, 30, 1),
    (-20 + x, -10, 1),
    (20 * x, 200, 20),
    (x + y, 30, sqrt(2**2 + 1**2)),
    (x * y, 200, sqrt(20**2 + 20**2)),
    (x / y, 0.5, sqrt((1/20)**2 + (2*10/(20**2))**2)),
]


@pytest.mark.parametrize(
    "unum, expected_val, expected_std_dev",
    binary_cases,
)
def test_binary(
        unum: UFloat,
        expected_val: float,
        expected_std_dev: float,
):
    assert unum.val == expected_val
    assert unum.std_dev == expected_std_dev


u_zero = UFloat(0, 0)
x = UFloat(10, 1)
y = UFloat(10, 1)
equals_cases = [
    (x, x),
    (x-x, u_zero),
    (2*x - x, x),
    (x*0, u_zero),
    (x*0, y*0),
]


@pytest.mark.parametrize(
    "first, second",
    equals_cases,
)
def test_equals(first, second):
    assert first == second


u_zero = UFloat(0, 0)
x = UFloat(10, 1)
y = UFloat(10, 1)
not_equals_cases = [
    (x, y),
    (x-y, u_zero),
    (x, 10),
]


@pytest.mark.parametrize(
    "first, second",
    not_equals_cases,
)
def test_not_equals(first, second):
    assert first != second


usin = ToUFuncPositional((lambda x: cos(x),))(sin)
x = UFloat(10, 2)
sin_cases = [
    (usin(x), sin(10), 2 * cos(10))
]


@pytest.mark.parametrize(
    "unum, expected_val, expected_std_dev",
    binary_cases,
)
def test_sin(
        unum: UFloat,
        expected_val: float,
        expected_std_dev: float,
):
    assert unum.val == expected_val
    assert unum.std_dev == expected_std_dev


u_zero = UFloat(0, 0)
x = UFloat(10, 2)
y = UFloat(10, 2)
bool_val_cases = [
    (u_zero, False),
    (x, True),
    (y, True),
    (x-y, True),
    (x-x, False),
    (y-y, False),
    (0*x, False),
    (0*y, False),
]


@pytest.mark.parametrize(
    "unum, bool_val",
    bool_val_cases,
)
def test_bool(unum: UFloat, bool_val: bool):
    assert bool(unum) is bool_val


def test_negative_std():
    with pytest.raises(ValueError, match=r'Uncertainty must be non-negative'):
        unum = UFloat(-1.0, -1.0)
