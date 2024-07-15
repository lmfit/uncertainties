from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from functools import lru_cache, wraps
import inspect
from math import sqrt
from numbers import Real
import sys
from typing import Callable, Collection, Dict, Optional, Tuple, Union, TYPE_CHECKING
import uuid

from uncertainties.parsing import str_to_number_with_uncert

if TYPE_CHECKING:
    from inspect import Signature


@dataclass(frozen=True)
class UncertaintyAtom:
    """
    Custom class to keep track of "atoms" of uncertainty. Two UncertaintyAtoms are
    always uncorrelated.
    """
    std_dev: float
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
    Recursively expand a linear combination of uncertainties out into the base atoms.
    It is a performance optimization to sometimes store unexpanded linear combinations.
    For example, there may be a long calculation involving many layers of UFloat
    manipulations. We need not expand the linear combination until the end when a
    calculation of a standard deviation on a UFloat is requested.
    """
    expanded_dict = defaultdict(float)
    for combo, combo_weight in combo:
        if isinstance(combo, UncertaintyAtom):
            expanded_dict[combo] += combo_weight
        else:
            expanded_combo = get_expanded_combo(combo)
            for atom, atom_weight in expanded_combo:
                expanded_dict[atom] += atom_weight * combo_weight

    return tuple((atom, weight) for atom, weight in expanded_dict.items())


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
            value: Real,
            uncertainty: Union[UncertaintyCombo, Real] = (),
            tag: Optional[str] = None
    ):
        self._val = float(value)
        if isinstance(uncertainty, Real):
            atom = UncertaintyAtom(float(uncertainty))
            uncertainty_combo = ((atom, 1.0),)
            self.uncertainty_lin_combo = uncertainty_combo
        else:
            self.uncertainty_lin_combo = uncertainty
        self.tag = tag

    @property
    def val(self: "UFloat") -> float:
        return self._val

    @property
    def std_dev(self: "UFloat") -> float:
        # TODO: It would be interesting to memoize/cache this result. However, if we
        #   stored this result as an instance attribute that would qualify as a mutation
        #   of the object and have implications for hashability. For example, two UFloat
        #   objects might have different uncertainty_lin_combo, but when expanded
        #   they're the same so that the std_dev and even correlations with other UFloat
        #   are the same. Should these two have the same hash? My opinion is no.
        #   I think a good path forward could be to cache this as an instance attribute
        #   nonetheless, but to not include the std_dev in the hash. Also equality would
        #   be based on equality of uncertainty_lin_combo, not equality of std_dev.
        expanded_lin_combo = get_expanded_combo(self.uncertainty_lin_combo)
        list_of_squares = [
            (weight * atom.std_dev)**2 for atom, weight in expanded_lin_combo
        ]
        std_dev = sqrt(sum(list_of_squares))
        return std_dev

    def __eq__(self: "UFloat", other: "UFloat") -> bool:
        if not isinstance(other, UFloat):
            return False
        val_eq = self.val == other.val
        uncertainty_eq = self.uncertainty_lin_combo == other.uncertainty_lin_combo
        return val_eq and uncertainty_eq

    # def __gt__(self, other):

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self.val}, {self.std_dev})'

    # Aliases
    @property
    def nominal_value(self: "UFloat") -> float:
        return self.val

    @property
    def n(self: "UFloat") -> float:
        return self.val

    @property
    def s(self: "UFloat") -> float:
        return self.std_dev


SQRT_EPS = sqrt(sys.float_info.epsilon)


def get_param_name(sig: Signature, param: Union[int, str]):
    if isinstance(param, int):
        param_name = list(sig.parameters.keys())[param]
    else:
        param_name = param
    return param_name


def numerical_partial_derivative(
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

    for param, arg in lower_bound_sig.arguments.items():
        if isinstance(arg, UFloat):
            lower_bound_sig.arguments[param] = arg.val
            upper_bound_sig.arguments[param] = arg.val

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
            deriv_func_dict = {}
        for ufloat_param in ufloat_params:
            if ufloat_param not in deriv_func_dict:
                deriv_func_dict[ufloat_param] = None
        self.deriv_func_dict: DerivFuncDict = deriv_func_dict

    def __call__(self, f: Callable[..., float]):
        sig = inspect.signature(f)

        @wraps(f)
        def wrapped(*args, **kwargs):
            float_bound = sig.bind(*args, **kwargs)

            return_u_val = False
            for param, param_val in float_bound.arguments.items():
                if isinstance(param_val, UFloat):
                    float_bound.arguments[param] = param_val.val
                    return_u_val = True
                elif isinstance(param_val, Real):
                    float_bound.arguments[param] = float(param_val)

            new_val = f(*float_bound.args, **float_bound.kwargs)
            if not return_u_val:
                return new_val

            ufloat_bound = sig.bind(*args, **kwargs)
            new_uncertainty_lin_combo = []
            for u_float_param in self.ufloat_params:
                u_float_param_name = get_param_name(sig, u_float_param)
                arg = ufloat_bound.arguments[u_float_param_name]
                if isinstance(arg, UFloat):
                    deriv_func = self.deriv_func_dict[u_float_param]
                    if deriv_func is None:
                        derivative = numerical_partial_derivative(
                            f,
                            u_float_param_name,
                            *args,
                            **kwargs,
                        )
                    else:
                        derivative = deriv_func(*float_bound.args, **float_bound.kwargs)

                    new_uncertainty_lin_combo.append(
                        (arg.uncertainty_lin_combo, derivative)
                    )

            unc_linear_combo = tuple(new_uncertainty_lin_combo)
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


PositionalDerivFunc = Union[Callable[..., float], str]


def deriv_func_dict_positional_helper(
        deriv_funcs: Tuple[Optional[PositionalDerivFunc]],
):
    nargs = len(deriv_funcs)
    deriv_func_dict = {}

    for arg_num, deriv_func in enumerate(deriv_funcs):
        if deriv_func is None:
            pass
        elif callable(deriv_func):
            pass
        elif isinstance(deriv_func, str):
            deriv_func = func_str_to_positional_func(deriv_func, nargs)
        else:
            raise ValueError(
                f'Invalid deriv_func: {deriv_func}. Must be None, callable, or a '
                f'string.'
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
    def __init__(self, deriv_funcs: Tuple[Optional[PositionalDerivFunc]]):
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



import numpy as np
arr = np.array([x, y, z])
print(np.mean(arr))