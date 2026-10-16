import math
import warnings
from decimal import Decimal, localcontext

import pytest

from uncertainties import ufloat


@pytest.mark.parametrize("scale", [1e-300, 1e-200, 1e-160, 1, 1e160, 1e200, 1e300])
def test_std_dev_preserves_finite_component_norm(scale):
    x = ufloat(scale, scale * 0.01)
    y = ufloat(2 * scale, scale * 0.02)
    with localcontext() as context:
        context.prec = 80
        expected = float(
            (
                Decimal.from_float(scale * 0.01) ** 2
                + Decimal.from_float(scale * 0.02) ** 2
            ).sqrt()
        )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        result = x + y
        assert result.s == pytest.approx(expected, rel=1e-15, abs=0)
        assert (2 * x - x).s == x.std_dev
        assert (x - x).s == 0


@pytest.mark.parametrize("errors", [(math.inf,), (math.nan,), (math.inf, math.nan)])
def test_std_dev_preserves_nonfinite_propagation(errors):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        value = sum(ufloat(1, error) for error in errors)
        if any(math.isnan(error) for error in errors):
            assert math.isnan(value.s)
        else:
            assert value.s == math.inf
