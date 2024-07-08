from uncertainties import unumpy, umath, ufloat, UFloat
from helpers import (power_special_cases, power_all_cases, power_wrt_ref, numbers_close,
    ufloats_close, compare_derivatives, uarrays_close, nominal_and_std_dev_close)
import numpy as np
import pytest 

a = ufloat(1, 0.1)
b = ufloat(2, 0.2)


class TestArithmetic():
    @pytest.mark.parametrize(
        "first, second, expected",
        [
            (a, b, ufloat(3.0, 0.223606797749979)),
            (a, a, ufloat(2.0, 0.2)),
        ],
    )
    def test_add(self, first, second, expected):
        result = first + second
        assert nominal_and_std_dev_close(result, expected)

        result = np.add(first, second)
        assert nominal_and_std_dev_close(result, expected)

    @pytest.mark.parametrize(
        "first, second, expected",
        [
            (a, b, ufloat(-1.00, 0.223606797749979)),
            (a, a, ufloat(0.0, 0.0)),
        ],
    )
    def test_subtact(self, first, second, expected):
        result = first - second
        assert nominal_and_std_dev_close(result, expected)

        result = np.subtract(first, second)
        assert nominal_and_std_dev_close(result, expected)

    @pytest.mark.parametrize(
        "first, second, expected",
        [
            (a, b, ufloat(2.0, 0.28284271247461906)),
            (a, a, ufloat(1.0, 0.2)),
        ],
    )
    def test_multiply(self, first, second, expected):
        result = first * second
        assert nominal_and_std_dev_close(result, expected)

        result = np.multiply(first, second)
        assert nominal_and_std_dev_close(result, expected)

    @pytest.mark.parametrize(
        "first, second, expected",
        [
            (a, b, ufloat(0.5, 0.07071067811865477)),
            (a, a, ufloat(1.0, 0.0)),
        ],
    )
    def test_divide(self, first, second, expected):
        result = first / second
        assert nominal_and_std_dev_close(result, expected)

        result = np.divide(first, second)
        assert nominal_and_std_dev_close(result, expected)

        result = np.true_divide(first, second)
        assert nominal_and_std_dev_close(result, expected)

        
    @pytest.mark.parametrize(
        "first, second, expected",
        [
            (a, b, ufloat(0.0, 0.0)),
            (a, a, ufloat(1.0, 0.0)),
        ],
    )
    def test_floor_divide(self, first, second, expected):
        result = first // second
        assert nominal_and_std_dev_close(result, expected)

        result = np.floor_divide(first, second)
        assert nominal_and_std_dev_close(result, expected)
    

class TestComparative():
    @pytest.mark.parametrize(
        "first, second, expected",
        [
            (a, b, False),
            (a, a, True),
        ],
    )
    def test_equal(self, first, second, expected):
        result = first == second
        assert result == expected

        result = np.equal(first, second)
        assert result == expected

    @pytest.mark.parametrize(
        "first, second, expected",
        [
            (a, b, True),
            (a, a, False),
        ],
    )
    def test_not_equal(self, first, second, expected):
        result = first != second
        assert result == expected

        result = np.not_equal(first, second)
        assert result == expected

    @pytest.mark.parametrize(
        "first, second, expected",
        [
            (a, b, True),
            (a, a, False),
        ],
    )
    def test_less(self, first, second, expected):
        result = first < second
        assert result == expected

        result = np.less(first, second)
        assert result == expected

    @pytest.mark.parametrize(
        "first, second, expected",
        [
            (a, b, True),
            (a, a, True),
        ],
    )
    def test_less_equal(self, first, second, expected):
        result = first <= second
        assert result == expected

        result = np.less_equal(first, second)
        assert result == expected

    @pytest.mark.parametrize(
        "first, second, expected",
        [
            (a, b, False),
            (a, a, False),
        ],
    )
    def test_greater(self, first, second, expected):
        result = first > second
        assert result == expected

        result = np.greater(first, second)
        assert result == expected

    @pytest.mark.parametrize(
        "first, second, expected",
        [
            (a, b, False),
            (a, a, True),
        ],
    )
    def test_greater_equal(self, first, second, expected):
        result = first >= second
        assert result == expected

        result = np.greater_equal(first, second)
        assert result == expected


class TestUfuncs():
    zero = ufloat(0.0, 0.1)
    one = ufloat(1.0, 0.1)
    pi_4 = ufloat(0.7853981633974483, 0.1)  # pi/4
    pi_2 = ufloat(1.5707963267948966, 0.1)  # pi/2
    @pytest.mark.parametrize(
        "numpy_func, umath_func, arg, expected",
        [
            ('cos', 'cos',  zero, ufloat(1.0, 0.0)),
            ('cos', 'cos',  pi_4, ufloat(0.7071067811865476, 0.07071067811865477)),
            ('cos', 'cos',  pi_2, ufloat(6.123233995736766e-17, 0.1)),
            ('cosh', 'cosh',  zero, ufloat(1.0, 0.0)),
            ('cosh', 'cosh',  pi_4, ufloat(1.324609089252006, 0.08686709614860096)),
            ('cosh', 'cosh',  pi_2, ufloat(2.5091784786580567, 0.2301298902307295)),
            ('sin', 'sin',  zero, ufloat(0.0, 0.1)),
            ('sin', 'sin',  pi_4, ufloat(0.7071067811865476, 0.07071067811865477)),
            ('sin', 'sin',  pi_2, ufloat(1.0, 6.123233995736766e-18)),
            ('sinh', 'sinh',  zero, ufloat(0.0, 0.1)),
            ('sinh', 'sinh',  pi_4, ufloat(0.8686709614860095, 0.1324609089252006)),
            ('sinh', 'sinh',  pi_2, ufloat(2.3012989023072947, 0.2509178478658057)),
            ('tan', 'tan',  zero, ufloat(0.0, 0.1)),
            ('tan', 'tan',  pi_4, ufloat(0.9999999999999999, 0.19999999999999998)),
            ('tan', 'tan',  pi_2, ufloat(1.633123935319537e+16, 2.6670937881135717e+31)),
            ('tanh', 'tanh',  zero, ufloat(0.0, 0.1)),
            ('tanh', 'tanh',  pi_4, ufloat(0.6557942026326724, 0.05699339637933774)),
            ('tanh', 'tanh',  pi_2, ufloat(0.9171523356672744, 0.015883159318006324)),
            ('arccos', 'acos',  zero, ufloat(1.5707963267948966, 0.1)),
            ('arccos', 'acos',  one, ufloat(0.0, float("nan"))),
            ('arccosh', 'acosh',  one, ufloat(0.0, float("nan"))),
            ('arcsin', 'asin',  zero, ufloat(0.0, 0.1)),
            ('arcsin', 'asin',  one, ufloat(1.5707963267948966, float("nan"))),
            ('arcsinh', 'asinh',  zero, ufloat(0.0, 0.1)),
            ('arcsinh', 'asinh',  one, ufloat(0.8813735870195429, 0.07071067811865475)),
            ('arctan', 'atan',  zero, ufloat(0.0, 0.1)),
            ('arctan', 'atan',  one, ufloat(0.7853981633974483, 0.05)),
            ('arctanh', 'atanh',  zero, ufloat(0.0, 0.1)),
            ('exp', 'exp',  zero, ufloat(1.0, 0.1)),
            ('exp', 'exp',  one, ufloat(2.718281828459045, 0.27182818284590454)),
            ('exp2', None,  zero, ufloat(1.0, 0.06931471805599453)),
            ('exp2', None,  one, ufloat(2.0, 0.13862943611198905)),
            ('expm1', 'expm1',  zero, ufloat(0.0, 0.1)),
            ('expm1', 'expm1',  one, ufloat(1.718281828459045, 0.27182818284590454)),
            ('log10', 'log10',  one, ufloat(0.0, 0.04342944819032518)),
            ('log1p', 'log1p',  zero, ufloat(0.0, 0.1)),
            ('log1p', 'log1p',  one, ufloat(0.6931471805599453, 0.05)),
            ('degrees', 'degrees',  zero, ufloat(0.0, 5.729577951308233)),
            ('degrees', 'degrees',  one, ufloat(57.29577951308232, 5.729577951308233)),
            ('radians', 'radians',  zero, ufloat(0.0, 0.0017453292519943296)),
            ('radians', 'radians',  one, ufloat(0.017453292519943295, 0.0017453292519943296)),
            ('rad2deg', 'degrees',  zero, ufloat(0.0, 5.729577951308233)),
            ('rad2deg', 'degrees',  one, ufloat(57.29577951308232, 5.729577951308233)),
            ('deg2rad', 'radians',  zero, ufloat(0.0, 0.0017453292519943296)),
            ('deg2rad', 'radians',  one, ufloat(0.017453292519943295, 0.0017453292519943296)),
            ('sqrt', 'sqrt',  zero, ufloat(0.0, float("nan"))),
            ('sqrt', 'sqrt',  one, ufloat(1.0, 0.05)),
        ],
    )
    def test_single_arg(self, numpy_func, umath_func, arg, expected):
        func = getattr(np, numpy_func)
        result = func(arg)
        assert nominal_and_std_dev_close(result, expected)

        if umath_func:
            func = getattr(umath, umath_func)
            result = func(arg)
            assert nominal_and_std_dev_close(result, expected)
        