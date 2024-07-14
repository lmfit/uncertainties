from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from functools import lru_cache, wraps
import inspect
from math import sqrt
import sys
from typing import Callable, Collection, Dict, Optional, Tuple, Union, TYPE_CHECKING
import uuid

from uncertainties.parsing import str_to_number_with_uncert

if TYPE_CHECKING:
    from inspect import Signature


@dataclass(frozen=True)
class UncertaintyAtom:
    """
    Custom class to keep track of "atoms" of uncertainty. Note e.g. that
    UncertaintyAtom(3) is UncertaintyAtom(3)
    returns False.
    """
    unc: float
    uuid: uuid.UUID = field(default_factory=uuid.uuid4, init=False)


UncertaintyCombo = Tuple[
    Tuple[
        Union[UncertaintyAtom, "UncertaintyCombo"],
        float
    ],
    ...
]
UncertaintyComboExpanded = Tuple[
    Tuple[
        UncertaintyAtom,
        float
    ],
    ...
]


@lru_cache
def get_expanded_combo(combo: UncertaintyCombo) -> UncertaintyComboExpanded:
    """
    Recursively expand a linear combination of uncertainties out into the base Atoms.
    It is a performance optimization to sometimes store unexpanded linear combinations.
    For example, there may be a long calculation involving many layers of UFloat
    manipulations. We need not expand the linear combination until the end when a
    calculation of a standard deviation on a UFloat is requested.
    """
    expanded_dict = defaultdict(float)
    for unc_combo_1, weight_1 in combo:
        if isinstance(unc_combo_1, UncertaintyAtom):
            expanded_dict[unc_combo_1] += weight_1
        else:
            expanded_sub_combo = get_expanded_combo(unc_combo_1)
            for unc_atom, weight_2 in expanded_sub_combo:
                expanded_dict[unc_atom] += weight_2 * weight_1

    return tuple((unc_atom, weight) for unc_atom, weight in expanded_dict.items())


Value = Union["UValue", float]


class UFloat:
    """
    Core class. Stores a mean value (val, nominal_value, n) and an uncertainty stored
    as a (possibly unexpanded) linear combination of uncertainty atoms. Two UFloat's
    which share non-zero weight for a certain uncertainty atom are correlated.

    UFloats can be combined using arithmetic and more sophisticated mathematical
    operations. The uncertainty is propagation using rules of linear uncertainty
    propagation.
    """
    def __init__(
            self,
            /,
            val: float,
            unc: Union[UncertaintyCombo, float] = (),
            tag: Optional[str] = None
    ):
        self._val = float(val)
        if isinstance(unc, (float, int)):
            unc_atom = UncertaintyAtom(float(unc))
            unc_combo = ((unc_atom, 1.0),)
            self.unc_linear_combo = unc_combo
        else:
            self.unc_linear_combo = unc
        self.tag = tag

    @property
    def val(self: "UFloat") -> float:
        return self._val

    @property
    def unc(self: "UFloat") -> float:
        expanded_combo = get_expanded_combo(self.unc_linear_combo)
        return float(sqrt(sum([(weight * unc_atom.unc)**2 for unc_atom, weight in expanded_combo])))

    @property
    def nominal_value(self: "UFloat") -> float:
        return self.val

    @property
    def n(self: "UFloat") -> float:
        return self.val

    @property
    def std_dev(self: "UFloat") -> float:
        return self.unc

    @property
    def s(self: "UFloat") -> float:
        return self.unc

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self.val}, {self.unc})'


SQRT_EPS = sqrt(sys.float_info.epsilon)


def get_param_name(sig: Signature, param: Union[int, str]):
    if isinstance(param, int):
        param_name = list(sig.parameters.keys())[param]
    else:
        param_name = param
    return param_name


def partial_derivative(
        f: Callable[..., float],
        target_param: Union[str, int],
        *args,
        **kwargs
) -> float:
    """
    Numerically calculate the partial derivative of a function f with respect to the
    target_param (string name or position number of the float parameter to f to be
    varied) holding all other arguments, *args and **kwargs, constant.
    """
    sig = inspect.signature(f)
    lower_bound_sig = sig.bind(*args, **kwargs)
    upper_bound_sig = sig.bind(*args, **kwargs)

    for arg, val in lower_bound_sig.arguments.items():
        if isinstance(val, UFloat):
            lower_bound_sig.arguments[arg] = val.val
            upper_bound_sig.arguments[arg] = val.val

    target_param_name = get_param_name(sig, target_param)

    x = lower_bound_sig.arguments[target_param_name]
    dx = abs(x) * SQRT_EPS  # Numerical Recipes 3rd Edition, eq. 5.7.5

    # Inject x - dx into target_param and evaluate f
    lower_bound_sig.arguments[target_param_name] = x - dx
    lower_y = f(*lower_bound_sig.args, **lower_bound_sig.kwargs)

    # Inject x + dx into target_param and evaluate f
    upper_bound_sig.arguments[target_param_name] = x + dx
    upper_y = f(*upper_bound_sig.args, **upper_bound_sig.kwargs)

    derivative = (upper_y - lower_y) / (2 * dx)
    return derivative


ParamSpecifier = Union[str, int]
DerivFuncDict = Optional[Dict[ParamSpecifier, Optional[Callable[..., float]]]]


class ToUFunc:
    """
    Decorator to convert a function which typically accepts float inputs into a function
    which accepts UFloat inputs.

    >>> @ToUFunc(('x', 'y'))
    >>> def multiply(x, y, print_str='print this string!', do_print=False):
    ...     if do_print:
    ...         print(print_str)
    ...     return x * y

    Pass in a list of parameter names which correspond to float inputs that should now
    accept UFloat inputs.

    To calculate the output nominal value the decorator replaces all float inputs with
    their respective nominal values and evaluates the function directly.

    To calculate the output uncertainty linear combination the decorator calculates the
    partial derivative of the function with respect to each UFloat entry and appends the
    uncertainty linear combination corresponding to that UFloat, weighted by the
    corresponding partial derivative.

    The partial derivative is evaluated numerically by default using the
    partial_derivative() function. However, the user can optionaly pass in
    deriv_func_dict which maps each u_float parameter to a function that will calculate
    the partial derivative given *args and **kwargs supplied to the converted function.
    This latter approach may provide performance optimizations when it is faster to
    use an analytic formula to evaluate the partial derivative than the numerical
    calculation.
    """
    def __init__(
            self,
            ufloat_params: Collection[ParamSpecifier],
            deriv_func_dict: DerivFuncDict = None,
    ):
        self.ufloat_params = ufloat_params
        if deriv_func_dict is None:
            deriv_func_dict = {
                ufloat_param: None for ufloat_param in self.ufloat_params
            }
        self.deriv_func_dict: DerivFuncDict = deriv_func_dict

    def __call__(self, f: Callable[..., float]):
        sig = inspect.signature(f)

        @wraps(f)
        def wrapped(*args, **kwargs):
            """
            Calculate the
            """
            unc_linear_combo = []
            bound = sig.bind(*args, **kwargs)
            float_bound = sig.bind(*args, **kwargs)

            return_u_val = False
            for param, param_val in float_bound.arguments.items():
                if isinstance(param_val, UFloat):
                    float_bound.arguments[param] = param_val.val
                    return_u_val = True
                elif isinstance(param_val, int):
                    float_bound.arguments[param] = float(param_val)

            new_val = f(*float_bound.args, **float_bound.kwargs)
            if not return_u_val:
                return new_val

            for u_float_param in self.ufloat_params:
                u_float_param_name = get_param_name(sig, u_float_param)
                arg = bound.arguments[u_float_param_name]
                if isinstance(arg, UFloat):
                    sub_unc_linear_combo = arg.unc_linear_combo
                    deriv_func = self.deriv_func_dict[u_float_param]
                    if deriv_func is None:
                        derivative = partial_derivative(
                            f,
                            u_float_param_name,
                            *args,
                            **kwargs,
                        )
                    else:
                        derivative = deriv_func(*float_bound.args, **float_bound.kwargs)

                    unc_linear_combo.append((sub_unc_linear_combo, derivative))

            unc_linear_combo = tuple(unc_linear_combo)
            return UFloat(new_val, unc_linear_combo)

        return wrapped


def func_str_to_positional_func(func_str, nargs):
    if nargs == 1:
        def pos_func(x):
            return eval(func_str)
    elif nargs == 2:
        def pos_func(x, y):
            return eval(func_str)
    else:
        raise ValueError(f'Only nargs=1 or nargs=2 is supported, not {nargs=}.')
    return pos_func


def deriv_func_dict_positional_helper(deriv_funcs):
    if not isinstance(deriv_funcs, tuple):
        raise ValueError(f'deriv_funcs must be a tuple, not \"{deriv_funcs}\".')

    nargs = len(deriv_funcs)
    deriv_func_dict = {}

    for arg_num, deriv_func in enumerate(deriv_funcs):
        if isinstance(deriv_func, str):
            deriv_func = func_str_to_positional_func(deriv_func, nargs)
        elif deriv_func is None:
            pass
        else:
            if not callable(deriv_func):
                raise ValueError(
                    f'Derivative functions must be callable or strings. Not '
                    f'{deriv_func}.'
                )
        deriv_func_dict[arg_num] = deriv_func
    return deriv_func_dict


class ToUFuncPositional(ToUFunc):
    """
    Helper decorator for decorating a function to be UFloat compatible when only
    positional arguments are being converted. Instead of passing a list of parameter
    specifiers (names or number of parameters) and a dict of
    parameter specifiers : derivative functions
    we just pass a list of derivative functions. Each derivative function can either be
    a callable of a function string like '-x/y**2'.
    """
    def __init__(self, deriv_funcs: tuple[Callable[..., float]]):
        ufloat_params = tuple(range(len(deriv_funcs)))
        deriv_func_dict = deriv_func_dict_positional_helper(deriv_funcs)
        super().__init__(ufloat_params, deriv_func_dict)


def add_float_funcs_to_uvalue():
    """
    Monkey-patch common float instance methods over to UFloat

    Here I use a notation involving x and y which is parsed by
    resolve_deriv_func_dict_from_func_str_list. This is a compact way to specify the
    formulas to calculate the partial derivatives of binary and unary functions.

    # TODO: There's a bit of complexity added by allowing analytic derivative function
    #   in addition to the default numerical derivative function. It would be
    #   interesting to see performance differences between the two methods. Is the
    #   added complexity *actually* buying performance?
    """
    float_funcs_dict = {
        '__abs__': ('abs(x)/x',),
        '__pos__': ('1',),
        '__neg__': ('-1',),
        '__trunc__': ('0',),
        '__add__': ('1', '1'),
        '__radd__': ('1', '1'),
        '__sub__': ('1', '-1'),
        '__rsub__': ('-1', '1'),  # Note reversed order
        '__mul__': ('y', 'x'),
        '__rmul__': ('x', 'y'),  # Note reversed order
        '__truediv__': ('1/y', '-x/y**2'),
        '__rtruediv__': ('-x/y**2', '1/y'),  # Note reversed order
        '__floordiv__': ('0', '0'),  # ?
        '__rfloordiv__': ('0', '0'),  # ?
        '__pow__': (None, None),  # TODO: add these, see `uncertainties` source
        '__rpow__': (None, None),
        '__mod__': (None, None),
        '__rmod__': (None, None),
    }

    for func_name, deriv_funcs in float_funcs_dict.items():
        float_func = getattr(float, func_name)
        ufloat_ufunc = ToUFuncPositional(deriv_funcs)(float_func)
        setattr(UFloat, func_name, ufloat_ufunc)


add_float_funcs_to_uvalue()


def ufloat(val, unc, tag=None):
    return UFloat(val, unc, tag)

def ufloat_fromstr(string, tag=None):
    (nom, std) = str_to_number_with_uncert(string.strip())
    return ufloat(nom, std, tag)

"""
^^^
End libary code
____
Begin sample test code
vvvv
"""
from math import sin

usin = ToUFunc((0,))(sin)

x = UFloat(10, 1)

y = UFloat(10, 1)

z = UFloat(20, 2)

print(f'{x=}')
print(f'{-x=}')
print(f'{3*x=}')
print(f'{x-x=}  # A UFloat is correlated with itself')

print(f'{y=}')
print(f'{x-y=}  # Two distinct UFloats are not correlated unless they have the same Uncertainty Atoms')

print(f'{z=}')

print(f'{x*z=}')
print(f'{x/z=}')
print(f'{x**z=}')

print(f'{usin(x)=}  # We can UFloat-ify complex functions')

# x=UFloat(10.0, 1.0)
# -x=UFloat(-10.0, 1.0)
# 3*x=UFloat(30.0, 3.0)
# x-x=UFloat(0.0, 0.0)  # A UFloat is correlated with itself
# y=UFloat(10.0, 1.0)
# x-y=UFloat(0.0, 1.4142135623730951)  # Two distinct UFloats are not correlated unless they have the same Uncertainty Atoms
# z=UFloat(20.0, 2.0)
# x*z=UFloat(200.0, 28.284271247461902)
# x/z=UFloat(0.5, 0.07071067811865477)
# x**z=UFloat(1e+20, 5.0207163276303525e+20)
# usin(x)=UFloat(-0.5440211108893698, 0.8390715289860964)  # We can UFloat-ify complex functions



