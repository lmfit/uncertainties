import math
from numbers import Real
import sys
from typing import Union

from uncertainties.core_new import UFloat, ToUFuncPositional


UReal = Union[Real, UFloat]


def acos(value: UReal) -> UReal: ...


def acosh(value: UReal) -> UReal: ...


def asinh(value: UReal) -> UReal: ...


def atan(value: UReal) -> UReal: ...


def atan2(x: UReal, y: UReal) -> UReal: ...


def atanh(value: UReal) -> UReal: ...


def cos(value: UReal) -> UReal: ...


def cosh(value: UReal) -> UReal: ...


def degrees(value: UReal) -> UReal: ...


def erf(value: UReal) -> UReal: ...


def erfc(value: UReal) -> UReal: ...


def exp(value: UReal) -> UReal: ...


def log(value: UReal) -> UReal: ...


def log10(value: UReal) -> UReal: ...


def radians(value: UReal) -> UReal: ...


def sin(value: UReal) -> UReal: ...


def sinh(value: UReal) -> UReal: ...


def sqrt(value: UReal) -> UReal: ...


def tan(value: UReal) -> UReal: ...


def tanh(value: UReal) -> UReal: ...


def log_der0(*args):
    """
    Derivative of math.log() with respect to its first argument.

    Works whether 1 or 2 arguments are given.
    """
    if len(args) == 1:
        return 1 / args[0]
    else:
        return 1 / args[0] / math.log(args[1])  # 2-argument form


deriv_dict = {
    # In alphabetical order, here:
    "acos": ("-1/math.sqrt(1-x**2)",),
    "acosh": ("1/math.sqrt(x**2-1)",),
    "asinh": ("1/math.sqrt(1+x**2)",),
    "atan": ("1/(1+x**2)",),
    "atan2": ('y/(x**2+y**2)', "-x/(x**2+y**2)"),
    "atanh": ("1/(1-x**2)",),
    "cos": ("-math.sin(x)",),
    "cosh": ("math.sinh(x)",),
    "degrees": ("math.degrees(1)",),
    "erf": ("(2/math.sqrt(math.pi))*math.exp(-(x**2))",),
    "erfc": ("-(2/math.sqrt(math.pi))*math.exp(-(x**2))",),
    "exp": ("math.exp(x)",),
    "log": (log_der0, "-math.log(x, y) / y / math.log(y)"),
    "log10": ("1/x/math.log(10)",),
    "radians": ("math.radians(1)",),
    "sin": ("math.cos(x)",),
    "sinh": ("math.cosh(x)",),
    "sqrt": ("0.5/math.sqrt(x)",),
    "tan": ("1 + math.tan(x)**2",),
    "tanh": ("1 - math.tanh(x)**2",),
}

this_module = sys.modules[__name__]

for func_name, deriv_funcs in deriv_dict.items():
    func = getattr(math, func_name)
    ufunc = ToUFuncPositional(deriv_funcs, eval_locals={"math": math})(func)
    setattr(this_module, func_name, ufunc)
