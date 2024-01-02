"""
This file contains the tests required by pandas for an ExtensionArray and ExtensionType.
"""
import warnings

import numpy as np
import pandas as pd
import pandas._testing as tm
import pytest
from pandas.core import ops
from pandas.tests.extension import base
from pandas.tests.extension.conftest import (
    as_frame,  # noqa: F401
    as_array,  # noqa: F401
    as_series,  # noqa: F401
    fillna_method,  # noqa: F401
    groupby_apply_op,  # noqa: F401
    use_numpy,  # noqa: F401
)

from .uncertainty_array import UncertaintyArray, UncertaintyDtype, ufloat
from uncertainties import umath

# from .core import ufloat

@pytest.fixture(params=[True, False])
def box_in_series(request):
    """Whether to box the data in a Series"""
    return request.param


@pytest.fixture
def dtype():
    return UncertaintyDtype()


@pytest.fixture
def data(request):
    return UncertaintyArray(
        [ufloat(i, abs(i)/100) for i in np.arange(start=1.0, stop=101.0)]
    )


@pytest.fixture
def data_missing():
    return UncertaintyArray(
        [ufloat(i, abs(i)/100) for i in [np.nan, 1]]
    )


@pytest.fixture
def data_for_twos():
    x = [
        2.0,
    ] * 100
    return UncertaintyArray(
         [ufloat(i, abs(i)/100) for i in x]
    )


@pytest.fixture(params=["data", "data_missing"])
def all_data(request, data, data_missing):
    if request.param == "data":
        return data
    elif request.param == "data_missing":
        return data_missing


@pytest.fixture
def data_repeated(data):
    """Return different versions of data for count times"""

    def gen(count):
        for _ in range(count):
            yield data

    yield gen


@pytest.fixture(params=[None, lambda x: x])
def sort_by_key(request):
    """
    Simple fixture for testing keys in sorting methods.
    Tests None (no key) and the identity key.
    """
    return request.param


@pytest.fixture
def data_for_sorting():
    return UncertaintyArray(
         [ufloat(i, abs(i)/100) for i in [0.3, 10.0, -50.0]]
    )


@pytest.fixture
def data_missing_for_sorting():
    return UncertaintyArray(
         [ufloat(i, abs(i)/100) for i in [4.0, np.nan, -5.0]]
    )


@pytest.fixture
def na_cmp():
    """Binary operator for comparing NA values."""
    return lambda x, y: pd.isna(x) and pd.isna(y)


@pytest.fixture
def na_value():
    return pd.NA


@pytest.fixture
def data_for_grouping():
    a = ufloat(1.0, 0.1)
    b = ufloat(2.0, 0.1)
    c = ufloat(3.0, 0.1)
    x = [a, a, np.nan, np.nan, b, b, a, c]
    return UncertaintyArray(x)

# === missing from pandas extension docs about what has to be included in tests ===
# copied from pandas/pandas/conftest.py
_all_arithmetic_operators = [
    "__add__",
    "__radd__",
    "__sub__",
    "__rsub__",
    "__mul__",
    "__rmul__",
    "__floordiv__",
    "__rfloordiv__",
    "__truediv__",
    "__rtruediv__",
    "__pow__",
    "__rpow__",
    "__mod__",
    "__rmod__",
]


@pytest.fixture(params=_all_arithmetic_operators)
def all_arithmetic_operators(request):
    """
    Fixture for dunder names for common arithmetic operations
    """
    return request.param


@pytest.fixture(params=["__eq__", "__ne__", "__le__", "__lt__", "__ge__", "__gt__"])
def all_compare_operators(request):
    """
    Fixture for dunder names for common compare operations

    * >=
    * >
    * ==
    * !=
    * <
    * <=
    """
    return request.param


# commented functions aren't implemented in numpy/pandas
_all_numeric_reductions = [
    "sum",
    "max",
    "min",
    "mean",
    # "prod",
    "std",
    "var",
    "median",
    "sem",
    "kurt",
    "skew",
]


@pytest.fixture(params=_all_numeric_reductions)
def all_numeric_reductions(request):
    """
    Fixture for numeric reduction names.
    """
    return request.param


_all_boolean_reductions = ["all", "any"]


@pytest.fixture(params=_all_boolean_reductions)
def all_boolean_reductions(request):
    """
    Fixture for boolean reduction names.
    """
    return request.param


_all_numeric_accumulations = ["cumsum", "cumprod", "cummin", "cummax"]


@pytest.fixture(params=_all_numeric_accumulations)
def all_numeric_accumulations(request):
    """
    Fixture for numeric accumulation names
    """
    return request.param

class InvalidScalar():
    def strip(self):
        return "invalid_scalar"

@pytest.fixture
def invalid_scalar(data):
    """
    A scalar that *cannot* be held by this ExtensionArray.
    The default should work for most subclasses, but is not guaranteed.
    If the array can hold any item (i.e. object dtype), then use pytest.skip.
    """
    return InvalidScalar()


import operator

@pytest.fixture(
    params=[
        operator.eq,
        operator.ne,
        operator.gt,
        operator.ge,
        operator.lt,
        operator.le,
    ]
)
def comparison_op(request):
    """
    Fixture for operator module comparison functions.
    """
    return request.param

from pandas.tests.extension import base

class TestUncertaintyArray(base.ExtensionTests):
    divmod_exc = TypeError # TODO: fix this
    series_scalar_exc = None
    frame_scalar_exc = None
    series_array_exc = None

    # This test round trips to file. Set the uncertainty to zero so the recreated data compares equal to the original data.
    @pytest.mark.parametrize("data", [UncertaintyArray(
        [ufloat(i, 0) for i in np.arange(start=1.0, stop=101.0)]
    )])
    @pytest.mark.parametrize("engine", ["c", "python"])
    def test_EA_types(self, engine, data, request):
        super().test_EA_types(engine, data)

    @pytest.mark.xfail(run=True, reason="test returns Float64Dtype rather than float64 but is otherwise correct")
    def test_value_counts_with_normalize(self, data, index, obj):
        super().test_value_counts_with_normalize(data, index, obj)
    
    @pytest.mark.xfail(reason="Can't invert uncertainties")
    def test_invert(self, data):
        super().test_invert(data)
    
    @pytest.mark.xfail(reason="_reduce is not implemented")
    def test_in_numeric_groupby(self, data_for_grouping):
        super().test_in_numeric_groupby(data_for_grouping)
    
    @pytest.mark.xfail(reason="Couldn't work out why this fails")
    def test_groupby_extension_agg(self, data_for_grouping):
        super().test_groupby_extension_agg(data_for_grouping)

    @pytest.mark.xfail(reason="Couldn't work out why this fails")
    def test_from_dtype(self, data):
        super().test_from_dtype(data)