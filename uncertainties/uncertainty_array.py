from __future__ import annotations

import numpy
import pandas as pd
import pandas.api.extensions
from pandas.api.types import is_dtype_equal, is_list_like, is_scalar, pandas_dtype

from .core import UFloat, ufloat

import re
import numpy as np
from uncertainties import umath
from pandas.core.arrays.base import ExtensionArray

from pandas.core import (
    algorithms as algos,
    arraylike,
    missing,
    nanops,
    ops,
)

from pandas.core.construction import (
    array as pd_array,
    ensure_wrapped_if_datetimelike,
    extract_array,
)
@pandas.api.extensions.register_extension_dtype
class UncertaintyDtype(pandas.api.extensions.ExtensionDtype):
    """
    Extension dtype for uncertainty objects.
    """

    na_value = ufloat(np.nan, np.nan)
    kind = "O"
    names = None
    name = "UncertaintyDtype"
    type = UFloat
    
    @property
    def _is_numeric(self):
        return True

    @classmethod
    def construct_from_string(cls, name: str):
        if not isinstance(name, str):
            raise TypeError(
                f"'construct_from_string' expects a string, got {type(name)}"
            )

        if name != cls.name:
            raise TypeError(f"Cannot construct a '{cls.__name__}' from 'another_type'")

        return cls()

    def construct_array_type(self):
        return UncertaintyArray

    def __repr__(self):
        """
        Return a string representation for this object.

        Invoked by unicode(df) in py2 only. Yields a Unicode String in both
        py2/py3.
        """

        return self.name





class UncertaintyArray(
    pandas.core.arraylike.OpsMixin, pandas.core.arrays._mixins.NDArrayBackedExtensionArray
):
    dtype = UncertaintyDtype()
    name = "UncertaintyArray"
    # scalar used to denote NA value inside our self._ndarray, e.g. -1 for
    # Categorical, iNaT for Period. Outside of object dtype, self.isna() should
    # be exactly locations in self._ndarray with _internal_fill_value. See:
    # https://github.com/pandas-dev/pandas/blob/main/pandas/core/arrays/_mixins.py
    _internal_fill_value = ufloat(np.nan, np.nan)

    def __init__(self, values, dtype=None, copy: bool = False):
        if not (
            isinstance(values, numpy.ndarray) and values.dtype == numpy.dtype("object")
            and all(isinstance(v, UFloat) for v in values)
        ):
            values = self.__ndarray(values)
        elif copy:
            values = values.copy()

        super().__init__(values=values, dtype=values.dtype)

    @classmethod
    def __ndarray(cls, scalars):
        return numpy.array(
            [ufloat(scalar) if not isinstance(scalar, UFloat) else scalar for scalar in scalars],
            "object",
        )

    @classmethod
    def _from_sequence(cls, scalars, *, dtype=None, copy=False):
        if dtype is not None:
            assert dtype.__class__ is cls.dtype.__class__
        return cls(cls.__ndarray(scalars))

    _from_sequence_of_strings = _from_sequence

    def astype(self, dtype, copy=True):
        dtype = pandas_dtype(dtype)
        if is_dtype_equal(dtype, self.dtype):
            if not copy:
                return self
            else:
                return self.copy()

        return super().astype(dtype, copy=copy)

    def _cmp_method(self, other, op):
        """Compare array values, for use in OpsMixin."""

        if is_scalar(other) and (
            pandas.isna(other) or isinstance(other, self.dtype.type)
        ):
            other = type(self)([other])

        if type(other) is not type(self):
            return NotImplemented

        oshape = getattr(other, "shape", None)
        if oshape != self.shape and oshape != (1,) and self.shape != (1,):
            raise TypeError(
                "Can't compare arrays with different shapes", self.shape, oshape
            )
        return op(self._ndarray, other._ndarray)

    def _from_factorized(self, unique, original):
        return self.__class__(unique)

    def isna(self):
        return np.array([umath.isnan(x) for x in self._ndarray], dtype=bool)

    def _validate_scalar(self, value):
        """
        Validate and convert a scalar value to datetime64[ns] for storage in
        backing NumPy array.
        """
        return self._ufloat(value)

    def _validate_searchsorted_value(self, value):
        """
        Convert a value for use in searching for a value in the backing numpy array.

        TODO: With pandas 2.0, this may be unnecessary. https://github.com/pandas-dev/pandas/pull/45544#issuecomment-1052809232
        """
        return self._validate_setitem_value(value)

    def __setitem__(self, index, value):
        if is_list_like(value) and is_scalar(index):
                raise ValueError("Length of index must match length of values")
        super().__setitem__(index, value)


    def _validate_setitem_value(self, value):
        """
        Convert a value for use in setting a value in the backing numpy array.
        """
        if is_list_like(value):
            _ufloat = self._ufloat
            return [_ufloat(v) for v in value]

        return self._ufloat(value)
    
    def _ufloat(self, value):
        if pd.isna(value):
            return self.dtype.na_value
        if not isinstance(value, UFloat):
            return ufloat(value)
        return value

        
    def _arith_method(self, other, op):
        res_name = ops.get_op_result_name(self, other)

        lvalues = self._ndarray
        rvalues = extract_array(other, extract_numpy=True, extract_range=True)
        rvalues = ops.maybe_prepare_scalar_for_op(rvalues, lvalues.shape)
        rvalues = ensure_wrapped_if_datetimelike(rvalues)
        if isinstance(rvalues, range):
            rvalues = np.arange(rvalues.start, rvalues.stop, rvalues.step)

        with np.errstate(all="ignore"):
            result = ops.arithmetic_op(lvalues, rvalues, op)

        return self._from_sequence(result, dtype=self.dtype, copy=False)

    _logical_method = _arith_method

    def __invert__(self) -> NumpyExtensionArray:
        return type(self)(~self._ndarray)

    def __neg__(self) -> NumpyExtensionArray:
        return type(self)(-self._ndarray)

    def __pos__(self) -> NumpyExtensionArray:
        return type(self)(+self._ndarray)

    def __abs__(self) -> NumpyExtensionArray:
        return type(self)(abs(self._ndarray))

