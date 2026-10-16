import math
import warnings
from decimal import Decimal, localcontext

import pytest

from uncertainties import ufloat
from uncertainties.umath import tanh


@pytest.mark.parametrize(
    "value", [-100, -50, -30, -20, -15, -10, -1, 0, 1, 10, 15, 20, 30, 50, 100]
)
def test_tanh_retains_uncertainty_at_saturated_inputs(value):
    with localcontext() as context:
        context.prec = 100
        x = Decimal(value)
        cosh = (x.exp() + (-x).exp()) / 2
        expected = float(Decimal.from_float(0.01) / cosh**2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        result = tanh(ufloat(value, 0.01))
        assert result.n == math.tanh(value)
        assert result.s == pytest.approx(expected, rel=2e-15, abs=0)


@pytest.mark.parametrize("value", [-math.inf, math.inf])
def test_tanh_infinite_input(value):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        result = tanh(ufloat(value, 0.01))
        assert result.n == math.tanh(value)
        assert result.s == 0
