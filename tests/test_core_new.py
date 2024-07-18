import math

import pytest

from uncertainties import umath_new
from uncertainties.core_new import UFloat, ToUFunc, ToUFuncPositional

from helpers import ufloats_close


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
    (x + y, 30, math.sqrt(2**2 + 1**2)),
    (x * y, 200, math.sqrt(20**2 + 20**2)),
    (x / y, 0.5, math.sqrt((1/20)**2 + (2*10/(20**2))**2)),
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


usin = ToUFuncPositional((lambda t: math.cos(t),))(math.sin)
sin_cases = [
    (
        usin(UFloat(10, 2)),
        math.sin(10),
        2 * math.cos(10),
    ),
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
        _ = UFloat(-1.0, -1.0)


func_derivs = ((k, v) for k, v in umath_new.deriv_dict.items())


@pytest.mark.parametrize("ufunc_name, ufunc_derivs", func_derivs)
def test_ufunc_analytic_numerical_partial(ufunc_name, ufunc_derivs):
    if ufunc_name == "acosh":
        # cosh returns values > 1
        args = (UFloat(1.1, 0.1),)
    elif ufunc_name == "atan2":
        # atan2 requires two arguments
        args = (UFloat(1.1, 0.1), UFloat(3.1, 0.2))
    else:
        args = (UFloat(0.1, 0.01),)
    ufunc = getattr(umath_new, ufunc_name)
    nfunc = ToUFunc(range(len(ufunc_derivs)))(getattr(math, ufunc_name))
    assert ufloats_close(ufunc(*args), nfunc(*args), tolerance=1e-6)
