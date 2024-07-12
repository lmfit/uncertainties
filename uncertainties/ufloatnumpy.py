import numpy as np
import math
# ufuncs are listed at https://numpy.org/doc/stable/reference/ufuncs.html
from . import ops
# from .umath_core import log_der0,_deriv_copysign, _deriv_fabs, _deriv_pow_0, _deriv_pow_1

from .ops import nan_if_exception


def log_der0(*args):
    """
    Derivative of math.log() with respect to its first argument.

    Works whether 1 or 2 arguments are given.
    """
    if len(args) == 1:
        return 1/args[0]
    else:
        return 1/args[0]/math.log(args[1])  # 2-argument form

    # The following version goes about as fast:

    ## A 'try' is used for the most common case because it is fast when no
    ## exception is raised:
    #try:
    #    return log_1arg_der(*args)  # Argument number check
    #except TypeError:
    #    return 1/args[0]/math.log(args[1])  # 2-argument form

def _deriv_copysign(x,y):
    if x >= 0:
        return math.copysign(1, y)
    else:
        return -math.copysign(1, y)

def _deriv_fabs(x):
    if x >= 0:
        return 1
    else:
        return -1

def _deriv_pow_0(x, y):
    if y == 0:
        return  0.
    elif x != 0 or y % 1 == 0:
        return y*math.pow(x, y-1)
    else:
        return float('nan')

def _deriv_pow_1(x, y):
    if x == 0 and y > 0:
        return 0.
    else:
        return math.log(x) * math.pow(x, y)

def is_upcast_type(t):
    # This can be used to allow downstream modules to overide operations; see pint
    # TODO add upcast_type list or dict to a public interface
    return False

def implements(numpy_func_string, func_type):
    """Register an __array_function__/__array_ufunc__ implementation for UArray
    objects.
    """
    print(numpy_func_string, func_type)

    def decorator(func):
        if func_type == "function":
            HANDLED_FUNCTIONS[numpy_func_string] = func
        elif func_type == "ufunc":
            HANDLED_UFUNCS[numpy_func_string] = func
        else:
            raise ValueError(f"Invalid func_type {func_type}")
        return func

    return decorator

HANDLED_FUNCTIONS = {} # noqa
HANDLED_UFUNCS = {}


def apply_func_elementwise(func, inputs, kwargs, result_dtype="object"):
    if len(inputs) == 1:
        result = func(*inputs, **kwargs)
    elif isinstance(inputs[0], np.ndarray):
        result = np.empty_like(inputs[0], dtype=result_dtype)
        for index, x in np.ndenumerate(inputs[0]):
            inputs_ = [x if i == 0 else inputs[i] for i in range(len(inputs))]
            result[index] = func(*inputs_, **kwargs)
        if inputs[0].ndim == 0:
            result = result.item()
    elif isinstance(inputs[1], np.ndarray):
        result = np.empty_like(inputs[1], dtype=result_dtype)
        for index, x in np.ndenumerate(inputs[1]):
            inputs_ = [x if i == 1 else inputs[i] for i in range(len(inputs))]
            result[index] = func(*inputs_, **kwargs)
        if inputs[1].ndim == 0:
            result = result.item()
    else:
        result = func(*inputs, **kwargs)
    return result

def numpy_wrap(func_type, func, args, kwargs, types):
    """Return the result from a NumPy function/ufunc as wrapped by uncertainties."""

    if func_type == "function":
        handled = HANDLED_FUNCTIONS
        # Need to handle functions in submodules
        name = ".".join(func.__module__.split(".")[1:] + [func.__name__])
    elif func_type == "ufunc":
        handled = HANDLED_UFUNCS
        # ufuncs do not have func.__module__
        name = func.__name__
    else:
        raise ValueError(f"Invalid func_type {func_type}")

    if name not in handled or any(is_upcast_type(t) for t in types):
        print("NotImplemented L54")
        raise TypeError
        return NotImplemented
    return handled[name](*args, **kwargs)

class UFloatNumpy(object):
    # NumPy function/ufunc support
    __array_priority__ = 17

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if method != "__call__":
            # Only handle ufuncs as callables
            return NotImplemented

        # Replicate types from __array_function__
        types = {
            type(arg)
            for arg in list(inputs) + list(kwargs.values())
            if hasattr(arg, "__array_ufunc__")
        }

        return numpy_wrap("ufunc", ufunc, inputs, kwargs, types)

    def __array_function__(self, func, types, args, kwargs):
        return numpy_wrap("function", func, args, kwargs, types)

    # original code for _add_numpy_ufuncs. may be helpful for writing a generic wraps
    # this can be deleted:
    # @classmethod
    # def _add_numpy_ufuncs(cls):
    #     def implement_ufunc(func_str, derivatives):
    #         func = getattr(np, func_str)
    #         @implements(func_str, "ufunc")
    #         def implementation(*inputs, **kwargs):
    #             if isinstance(inputs[0], np.ndarray):
    #                 result = np.empty_like(inputs[0], dtype=object)
    #                 for index, x in np.ndenumerate(inputs[0]):
    #                     inputs_ = (x if i == 0 else inputs[i] for i in range(len(inputs)))
    #                     result[index] = cls.wrap(func, derivatives)(*inputs_, **kwargs)
    #             elif isinstance(inputs[1], np.ndarray):
    #                 result = np.empty_like(inputs[1], dtype=object)
    #                 for index, x in np.ndenumerate(inputs[1]):
    #                     inputs_ = (x if i == 1 else inputs[i] for i in range(len(inputs)))
    #                     result[index] = cls.wrap(func, derivatives)(*inputs_, **kwargs)
    #             else:
    #                 result = cls.wrap(func, derivatives)(*inputs, **kwargs)
    #             return result

    #         return implementation

    #     for func_str, derivatives in ufunc_derivatives.items():
    #         implement_ufunc(func_str, derivatives)


    @classmethod
    def _add_numpy_arithmetic_ufuncs(cls):
        def implement_ufunc(func_str, derivatives):
            func = getattr(np, func_str)
            @implements(func_str, "ufunc")
            def implementation(*inputs, **kwargs):
                return apply_func_elementwise(
                    ops._wrap(cls, func, derivatives), inputs, kwargs)
            return implementation

        ufunc_derivatives = {
            'add': [lambda x, y: 1., lambda x, y: 1.],
            'subtract': [lambda x, y: 1., lambda x, y: -1.],
            'multiply': [lambda x, y: y, lambda x, y: x],
            'divide': [lambda x, y: 1./y, lambda x, y: -x/y**2],
            'true_divide': [lambda x, y: 1./y, lambda x, y: -x/y**2],
            'floor_divide': [lambda x, y: 0., lambda x, y: 0.],
            
            'arccos': [nan_if_exception(lambda x: -1/math.sqrt(1-x**2))],
            'arccosh': [nan_if_exception(lambda x: 1/math.sqrt(x**2-1))],
            'arcsin': [nan_if_exception(lambda x: 1/math.sqrt(1-x**2))],
            'arcsinh': [nan_if_exception(lambda x: 1/math.sqrt(1+x**2))],
            'arctan': [nan_if_exception(lambda x: 1/(1+x**2))],
            'arctan2': [nan_if_exception(lambda y, x: x/(x**2+y**2)),  # Correct for x == 0
                        nan_if_exception(lambda y, x: -y/(x**2+y**2))],  # Correct for x == 0
            'arctanh': [nan_if_exception(lambda x: 1/(1-x**2))],
            'cos': [lambda x: -math.sin(x)],
            'cosh': [math.sinh],
            'sin': [math.cos],
            'sinh': [math.cosh],
            'tan': [nan_if_exception(lambda x: 1+math.tan(x)**2)],
            'tanh': [nan_if_exception(lambda x: 1-math.tanh(x)**2)],
            'exp': [math.exp],
            "exp2": [lambda y: _deriv_pow_1(2, y)],
            'expm1': [math.exp],
            'log10': [nan_if_exception(lambda x: 1/x/math.log(10))],
            'log1p': [nan_if_exception(lambda x: 1/(1+x))],
            'degrees': [lambda x: math.degrees(1)],
            'rad2deg': [lambda x: math.degrees(1)],
            'radians': [lambda x: math.radians(1)],
            'deg2rad': [lambda x: math.radians(1)],
            'power': [_deriv_pow_0, _deriv_pow_1],
            'sqrt': [nan_if_exception(lambda x: 0.5/math.sqrt(x))],
            'hypot': [lambda x, y: x/math.hypot(x, y),
                      lambda x, y: y/math.hypot(x, y)],
            }
        # TODO: test arctan2, power, hypot
        for func_str, derivatives in ufunc_derivatives.items():
            implement_ufunc(func_str, derivatives)

    @classmethod
    def _add_numpy_comparative_ufuncs(cls):
        def implement_ufunc(func_str, func):
            @implements(func_str, "ufunc")
            def implementation(*inputs, **kwargs):
                inputs = [cls._to_affine_scalar(i) for i in inputs]
                return apply_func_elementwise(func, inputs, kwargs, result_dtype=bool)

            return implementation

        ufunc_comparatives = {
            'equal': ops.eq_on_aff_funcs,
            'not_equal': ops.ne_on_aff_funcs,
            'less': ops.lt_on_aff_funcs,
            'less_equal': ops.le_on_aff_funcs,
            'greater': ops.gt_on_aff_funcs,
            'greater_equal': ops.ge_on_aff_funcs,
        }
        for func_str, func in ufunc_comparatives.items():
            implement_ufunc(func_str, func)

    # 'copysign': [_deriv_copysign,
    #              lambda x, y: 0],
    # 'degrees': [lambda x: math.degrees(1)],
    # 'erf': [lambda x: math.exp(-x**2)*erf_coef],
    # 'erfc': [lambda x: -math.exp(-x**2)*erf_coef],
    # 'fabs': [_deriv_fabs],
    # 'hypot': [lambda x, y: x/math.hypot(x, y),
    #           lambda x, y: y/math.hypot(x, y)],
    # 'log': [log_der0,
    #         lambda x, y: -math.log(x, y)/y/math.log(y)],
    # 'pow': [_deriv_pow_0, _deriv_pow_1],
    # 'radians': [lambda x: math.radians(1)],
    # 'sqrt': [lambda x: 0.5/math.sqrt(x)],
    # def _add_trig_ufuncs(cls):

    #         "cumprod": ("", ""),
    #         "arccos": ("", "radian"),
    #         "arcsin": ("", "radian"),
    #         "arctan": ("", "radian"),
    #         "arccosh": ("", "radian"),
    #         "arcsinh": ("", "radian"),
    #         "arctanh": ("", "radian"),
    #         "exp": ("", ""),
    #         "expm1": ("", ""),
    #         "exp2": ("", ""),
    #         "log": ("", ""),
    #         "log10": ("", ""),
    #         "log1p": ("", ""),
    #         "log2": ("", ""),
    #         "sin": ("radian", ""),
    #         "cos": ("radian", ""),
    #         "tan": ("radian", ""),
    #         "sinh": ("radian", ""),
    #         "cosh": ("radian", ""),
    #         "tanh": ("radian", ""),
    #         "radians": ("degree", "radian"),
    #         "degrees": ("radian", "degree"),
    #         "deg2rad": ("degree", "radian"),
    #         "rad2deg": ("radian", "degree"),
    #         "logaddexp": ("", ""),
    #         "logaddexp2": ("", ""),

            