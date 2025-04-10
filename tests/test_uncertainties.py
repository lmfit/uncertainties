import copy
import inspect
import math
import gc
import pickle
import random

import pytest

import uncertainties.core as uncert_core
from uncertainties.core import (
    ufloat,
    ufloat_fromstr,
    deprecated_methods,
)
from uncertainties import (
    umath,
    UFloat,
    correlated_values,
    correlated_values_norm,
    correlation_matrix,
)
from uncertainties.ops import partial_derivative
from helpers import (
    get_single_uatom,
    numbers_close,
    ufloats_close,
)


try:
    import numpy as np
except ImportError:
    np = None


def test_UFloat_class_construction():
    """Test creating UFloat directly."""
    x = UFloat(3, 0.14)
    assert x.nominal_value == 3
    assert x.std_dev == 0.14
    assert get_single_uatom(x).tag is None

    x = UFloat(3, 0.14, "pi")
    assert x.nominal_value == 3
    assert x.std_dev == 0.14
    assert get_single_uatom(x).tag == "pi"

    x = UFloat(3, 0.14, tag="pi")
    assert x.nominal_value == 3
    assert x.std_dev == 0.14
    assert get_single_uatom(x).tag == "pi"

    with pytest.raises(uncert_core.NegativeStdDev):
        _ = UFloat(3, -0.1)

    with pytest.raises(TypeError):
        UFloat(1)


def test_ufloat_function_construction():
    """Test creating UFloat via ufloat() function."""
    x = ufloat(3, 0.14)
    assert x.nominal_value == 3
    assert x.std_dev == 0.14
    assert get_single_uatom(x).tag is None

    x = ufloat(3, 0.14, "pi")
    assert x.nominal_value == 3
    assert x.std_dev == 0.14
    assert get_single_uatom(x).tag == "pi"

    # ... with tag keyword:
    x = ufloat(3, 0.14, tag="pi")
    assert x.nominal_value == 3
    assert x.std_dev == 0.14
    assert get_single_uatom(x).tag == "pi"

    with pytest.raises(uncert_core.NegativeStdDev):
        _ = ufloat(3, -0.1)

    with pytest.raises(TypeError):
        ufloat(1)


ufloat_from_str_cases = [
    ("-1.23(3.4)", -1.23, 3.4),
    ("  -1.23(3.4)  ", -1.23, 3.4),  # Test that leading and trailing spaces are ignored
    ("-1.34(5)", -1.34, 0.05),
    ("1(6)", 1, 6),
    ("3(4.2)", 3, 4.2),
    ("-9(2)", -9, 2),
    ("1234567(1.2)", 1234567, 1.2),
    ("12.345(15)", 12.345, 0.015),
    ("-12.3456(78)e-6", -12.3456e-6, 0.0078e-6),
    ("0.29", 0.29, 0.01),
    ("31.", 31, 1),
    ("-31.", -31, 1),
    # The following tests that the ufloat() routine does
    # not consider '31' like the tuple ('3', '1'), which would
    # make it expect two numbers (instead of 2 1-character
    # strings):
    ("31", 31, 1),
    ("-3.1e10", -3.1e10, 0.1e10),
    ("169.0(7)", 169, 0.7),
    ("-0.1+/-1", -0.1, 1),
    ("-13e-2+/-1e2", -13e-2, 1e2),
    ("-14.(15)", -14, 15),
    ("-100.0(15)", -100, 1.5),
    ("14.(15)", 14, 15),
    # Global exponent:
    ("(3.141+/-0.001)E+02", 314.1, 0.1),
    ## Pretty-print notation:
    # ± sign, global exponent (not pretty-printed):
    ("(3.141±0.001)E+02", 314.1, 0.1),
    # ± sign, individual exponent:
    ("3.141E+02±0.001e2", 314.1, 0.1),
    # ± sign, times symbol, superscript (= full pretty-print):
    ("(3.141 ± 0.001) × 10²", 314.1, 0.1),
    ## Others
    # Forced parentheses:
    ("(2 +/- 0.1)", 2, 0.1),
    # NaN uncertainty:
    ("(3.141±nan)E+02", 314.1, float("nan")),
    ("3.141e+02+/-nan", 314.1, float("nan")),
    ("3.4(nan)e10", 3.4e10, float("nan")),
    # NaN value:
    ("nan+/-3.14e2", float("nan"), 314),
    # Special float representation:
    ("-3(0.)", -3, 0),
]


@pytest.mark.parametrize("input_str,nominal_value,std_dev", ufloat_from_str_cases)
def test_ufloat_fromstr(input_str, nominal_value, std_dev):
    num = ufloat_fromstr(input_str)
    assert numbers_close(num.nominal_value, nominal_value)
    assert numbers_close(num.std_dev, std_dev)
    if std_dev != 0:
        assert get_single_uatom(num).tag is None
    else:
        assert num.uncertainty.ucombo_tuple == ()

    # With a tag as positional argument:
    num = ufloat_fromstr(input_str, "test variable")
    assert numbers_close(num.nominal_value, nominal_value)
    assert numbers_close(num.std_dev, std_dev)
    if std_dev != 0:
        assert get_single_uatom(num).tag == "test variable"
    else:
        assert num.uncertainty.ucombo_tuple == ()

    # With a tag as keyword argument:
    num = ufloat_fromstr(input_str, tag="test variable")
    assert numbers_close(num.nominal_value, nominal_value)
    assert numbers_close(num.std_dev, std_dev)
    if std_dev != 0:
        assert get_single_uatom(num).tag == "test variable"
    else:
        assert num.uncertainty.ucombo_tuple == ()


# Randomly generated but static test values.
deriv_propagation_cases = [
    ("__abs__", (1.1964838601545966,), 0.047308407404731856),
    ("__pos__", (1.5635699242286414,), 0.38219529954774223),
    ("__neg__", (-0.4520304708235554,), 0.8442835926901457),
    ("__trunc__", (0.4622631416873926,), 0.6540076679531033),
    ("__add__", (-0.7581877519537352, 1.6579645792821753), 0.5083165826806606),
    ("__radd__", (-0.976869259500134, 1.1542019729184076), 0.732839320238539),
    ("__sub__", (1.0233545960703134, 0.029354693323845993), 0.7475621525040559),
    ("__rsub__", (0.49861518245313663, -0.9927317702800833), 0.5421488555485847),
    ("__mul__", (0.0654070362874073, 1.9216078105121919), 0.6331001122119122),
    ("__rmul__", (-0.4006772142682373, 0.19628658198222926), 0.3300416314362784),
    ("__truediv__", (-0.5573378968194893, 0.28646277014641486), 0.42933306560556384),
    ("__rtruediv__", (1.7663869752268884, -0.1619387546963642), 0.6951025849642374),
    ("__floordiv__", (0.11750026664733992, -1.0120567560937617), 0.9557126076209381),
    ("__rfloordiv__", (-1.2872736512072698, -1.4416464249395973), 0.28262518984780205),
    ("__pow__", (0.34371967038364515, -0.8313605840956209), 0.6267147080961244),
    ("__rpow__", (1.593375683248082, 1.9890969272006154), 0.7171353266792271),
    ("__mod__", (0.7478106873313131, 1.2522332955942628), 0.5682413634363304),
    ("__rmod__", (1.5227432102303133, -0.5177923078991333), 0.25752786270795935),
]


@pytest.mark.parametrize("func_name, args, std_dev", deriv_propagation_cases)
def test_deriv_propagation(func_name, args, std_dev):
    func = getattr(UFloat, func_name)
    ufloat_args = [UFloat(arg, std_dev) for arg in args]
    output = func(*ufloat_args)

    for idx, _ in enumerate(ufloat_args):
        deriv = partial_derivative(func, idx)(*args)
        for atom, input_weight in ufloat_args[idx].error_components.items():
            output_weight = output.error_components[atom]
            assert numbers_close(output_weight, deriv * input_weight)


def test_copy():
    """Standard copy module integration."""
    x = ufloat(3, 0.1)
    assert x == x

    x_uatom = get_single_uatom(x)

    y = copy.copy(x)
    assert y == x
    assert get_single_uatom(y) == x_uatom

    z = copy.deepcopy(x)
    assert z == x
    assert get_single_uatom(z) == x_uatom

    t = x + 2 * z
    assert x_uatom in t.error_components

    t_copy = copy.copy(t)
    assert x_uatom in t_copy.error_components
    assert uncert_core.covariance_matrix([t, z]) == uncert_core.covariance_matrix(
        [t_copy, z]
    )

    t_deepcopy = copy.deepcopy(t)
    assert x_uatom in t_deepcopy.error_components
    assert uncert_core.covariance_matrix([t, z]) == uncert_core.covariance_matrix(
        [t_deepcopy, z]
    )

    del x
    gc.collect()
    assert x_uatom in y.error_components


"""
Classes to test pickling with different types of __slots__ inheritance
"""


class UFloatDict(UFloat):
    pass


class UFloatSlotsTuple(UFloat):
    __slots__ = ("new_attr",)


class UFloatSlotsStr(UFloat):
    __slots__ = "new_attr"


def test_pickling():
    """Standard pickle module integration."""
    x = UFloat(2, 0.1)
    x_unpickled = pickle.loads(pickle.dumps(x))

    assert x_unpickled == x

    f = 2 * x
    assert isinstance(f, UFloat)

    f_unpickled, x_unpickled2 = pickle.loads(pickle.dumps((f, x)))
    assert f_unpickled - 2 * x_unpickled2 == 0

    for subclass in (UFloatDict, UFloatSlotsTuple, UFloatSlotsStr):
        x = subclass(3, 0.14)

        # Pickling test with possibly uninitialized slots:
        assert pickle.loads(pickle.dumps(x)) == x

        # Unpickling test:
        x.new_attr = "New attr value"
        x_unpickled = pickle.loads(pickle.dumps(x))
        assert x_unpickled == x
        assert x_unpickled.new_attr == "New attr value"

    x = uncert_core.UCombo(())
    assert pickle.loads(pickle.dumps(x)).ucombo_tuple == ()


def test_int_div():
    "Integer division"
    # We perform all operations on floats, because derivatives can
    # otherwise be meaningless:
    x = ufloat(3.9, 2) // 2
    assert x.nominal_value == 1.0
    # All errors are supposed to be small, so the ufloat()
    # in x violates the assumption.  Therefore, the following is
    # correct:
    assert x.std_dev == 0.0


def test_comparison_ops():
    "Test of comparison operators"

    # Operations on quantities equivalent to Python numbers must still
    # be correct:
    a = ufloat(-3, 0)
    b = ufloat(10, 0)
    c = ufloat(10, 0)
    assert a < b
    assert a < 3
    assert 3 < b  # This is first given to int.__lt__()
    assert b == c

    x = ufloat(3, 0.1)

    # One constraint is that usual Python code for inequality testing
    # still work in a reasonable way (for instance, it is generally
    # desirable that functions defined by different formulas on
    # different intervals can still do "if 0 < x < 1:...".  This
    # supposes again that errors are "small" (as for the estimate of
    # the standard error).
    assert x > 1

    # The limit case is not obvious:
    assert not (x >= 3)
    assert not (x < 3)

    assert x == x
    # Comparaison between Variable and AffineScalarFunc:
    assert x == x + 0
    # Comparaison between 2 _different_ AffineScalarFunc objects
    # representing the same value:
    assert x / 2 == x / 2
    # With uncorrelated result that have the same behavior (value and
    # standard error):
    assert 2 * ufloat(1, 0.1) != ufloat(2, 0.2)
    # Comparaison between 2 _different_ Variable objects
    # that are uncorrelated:
    assert x != ufloat(3, 0.1)

    assert x != ufloat(3, 0.2)

    # Comparison to other types should work:
    assert x is not None  # Not comparable
    assert x - x == 0  # Comparable, even though the types are different
    assert x != [1, 2]

    ####################

    # Checks of the semantics of logical operations: they return True
    # iff they are always True when the parameters vary in an
    # infinitesimal interval inside sigma (sigma == 0 is a special
    # case):

    def test_all_comparison_ops(x, y):
        """
        Takes two Variable objects.

        Fails if any comparison operation fails to follow the proper
        semantics: a comparison only returns True if the correspond float
        comparison results are True for all the float values taken by
        the variables (of x and y) when they vary in an infinitesimal
        neighborhood within their uncertainty.

        This test is stochastic: it may, exceptionally, fail for
        correctly implemented comparison operators.
        """

        def random_float(var):
            """
            Returns a random value for Variable var, in an
            infinitesimal interval withing its uncertainty.  The case
            of a zero uncertainty is special.
            """
            return (random.random() - 0.5) * min(var.std_dev, 1e-5) + var.nominal_value

        # All operations are tested:
        for op in ["__%s__" % name for name in ("ne", "eq", "lt", "le", "gt", "ge")]:
            try:
                float_func = getattr(float, op)
            except AttributeError:  # Python 2.3's floats don't have __ne__
                continue

            # Determination of the correct truth value of func(x, y):

            sampled_results = []

            # The "main" value is an important particular case, and
            # the starting value for the final result
            # (correct_result):

            sampled_results.append(float_func(x.nominal_value, y.nominal_value))

            for check_num in range(50):  # Many points checked
                sampled_results.append(float_func(random_float(x), random_float(y)))

            min_result = min(sampled_results)
            max_result = max(sampled_results)

            if min_result == max_result:
                correct_result = min_result
            else:
                # Almost all results must be True, for the final value
                # to be True:
                num_min_result = sampled_results.count(min_result)

                # 1 exception is considered OK:
                correct_result = num_min_result == 1

            try:
                assert correct_result == getattr(x, op)(y)
            except AssertionError:
                print("Sampling results:", sampled_results)
                raise Exception(
                    "Semantic value of %s %s (%s) %s not"
                    " correctly reproduced." % (x, op, y, correct_result)
                )

    # With different numbers:
    test_all_comparison_ops(ufloat(3, 0.1), ufloat(-2, 0.1))
    test_all_comparison_ops(
        ufloat(0, 0),  # Special number
        ufloat(1, 1),
    )
    test_all_comparison_ops(
        ufloat(0, 0),  # Special number
        ufloat(0, 0.1),
    )
    # With identical numbers:
    test_all_comparison_ops(ufloat(0, 0), ufloat(0, 0))
    test_all_comparison_ops(ufloat(1, 1), ufloat(1, 1))


def test_logic():
    "Boolean logic: __nonzero__, bool."

    x = ufloat(3, 0)
    y = ufloat(0, 0)
    z = ufloat(0, 0.1)
    t = ufloat(-1, 2)

    assert bool(x)
    assert not bool(y)
    assert bool(z)
    assert bool(t)  # Only infinitseimal neighborhood are used


def test_basic_access_to_data():
    "Access to data from Variable and AffineScalarFunc objects."

    x = ufloat(3.14, 0.01, "x var")
    assert get_single_uatom(x).tag == "x var"
    assert x.nominal_value == 3.14
    assert x.std_dev == 0.01

    # Case of AffineScalarFunc objects:
    y = x + 0
    assert type(y) == UFloat
    assert y.nominal_value == 3.14
    assert y.std_dev == 0.01

    # Details on the sources of error:
    a = ufloat(-1, 0.001)
    y = 2 * x + 3 * x + 2 + a
    error_sources = y.error_components
    assert len(error_sources) == 2
    # 'a' and 'x'
    assert y.covariance(x) == 0.05 * 0.01
    assert y.covariance(a) == 0.001 * 0.001

    with pytest.raises(AttributeError):
        # std_dev cannot be modified
        x.std_dev = 1

    # Calculated values with uncertainties should not have a settable
    # standard deviation:
    y = 2 * x
    with pytest.raises(AttributeError):
        y.std_dev = 1

    # Calculation of deviations in units of the standard deviations:
    assert 10 / x.std_dev == x.std_score(10 + x.nominal_value)

    # "In units of the standard deviation" is not always meaningful:
    x = ufloat(1, 0)
    with pytest.raises(ValueError):
        x.std_score(1)


def test_correlations():
    "Correlations between variables"

    a = ufloat(1, 0)
    x = ufloat(4, 0.1)
    y = x * 2 + a
    # Correlations cancel "naive" additions of uncertainties:
    assert y.std_dev != 0
    normally_zero = y - (x * 2 + 1)
    assert normally_zero.nominal_value == 0
    assert normally_zero.std_dev == 0


def test_no_coercion():
    """
    Coercion of Variable object to a simple float.

    The coercion should be impossible, like for complex numbers.
    """

    x = ufloat(4, 1)
    try:
        assert float(x) == 4
    except TypeError:
        pass
    else:
        raise Exception("Conversion to float() should fail with TypeError")


def test_wrapped_func_no_args_no_kwargs():
    """
    Wrap a function that takes only positional-or-keyword parameters.
    """

    def f_auto_unc(x, y):
        return 2 * x + umath.sin(y)

    # Like f_auto_unc, but does not accept numbers with uncertainties:
    def f(x, y):
        assert not isinstance(x, uncert_core.UFloat)
        assert not isinstance(y, uncert_core.UFloat)
        return f_auto_unc(x, y)

    x = uncert_core.ufloat(1, 0.1)
    y = uncert_core.ufloat(10, 2)

    ### Automatic numerical derivatives:

    ## Fully automatic numerical derivatives:
    f_wrapped = uncert_core.wrap(f)
    assert ufloats_close(f_auto_unc(x, y), f_wrapped(x, y))

    # Call with keyword arguments:
    assert ufloats_close(f_auto_unc(y=y, x=x), f_wrapped(y=y, x=x))

    ## Automatic additional derivatives for non-defined derivatives,
    ## and explicit None derivative:
    f_wrapped = uncert_core.wrap(f, [None])  # No derivative for y
    assert ufloats_close(f_auto_unc(x, y), f_wrapped(x, y))

    # Call with keyword arguments:
    assert ufloats_close(f_auto_unc(y=y, x=x), f_wrapped(y=y, x=x))

    ### Explicit derivatives:

    ## Fully defined derivatives:
    f_wrapped = uncert_core.wrap(f, [lambda x, y: 2, lambda x, y: math.cos(y)])

    assert ufloats_close(f_auto_unc(x, y), f_wrapped(x, y))

    # Call with keyword arguments:
    assert ufloats_close(f_auto_unc(y=y, x=x), f_wrapped(y=y, x=x))

    ## Automatic additional derivatives for non-defined derivatives:
    f_wrapped = uncert_core.wrap(f, [lambda x, y: 2])  # No derivative for y
    assert ufloats_close(f_auto_unc(x, y), f_wrapped(x, y))

    # Call with keyword arguments:
    assert ufloats_close(f_auto_unc(y=y, x=x), f_wrapped(y=y, x=x))


def test_wrapped_func_args_no_kwargs():
    """
    Wrap a function that takes only positional-or-keyword and
    var-positional parameters.
    """

    def f_auto_unc(x, y, *args):
        return 2 * x + umath.sin(y) + 3 * args[1]

    # Like f_auto_unc, but does not accept numbers with uncertainties:
    def f(x, y, *args):
        assert not any(
            isinstance(value, uncert_core.UFloat) for value in [x, y] + list(args)
        )
        return f_auto_unc(x, y, *args)

    x = uncert_core.ufloat(1, 0.1)
    y = uncert_core.ufloat(10, 2)
    s = "string arg"
    z = uncert_core.ufloat(100, 3)

    args = [s, z, s]  # var-positional parameters

    ### Automatic numerical derivatives:

    ## Fully automatic numerical derivatives:
    f_wrapped = uncert_core.wrap(f)
    assert ufloats_close(f_auto_unc(x, y, *args), f_wrapped(x, y, *args))

    ## Automatic additional derivatives for non-defined derivatives,
    ## and explicit None derivative:
    f_wrapped = uncert_core.wrap(f, [None])  # No derivative for y
    assert ufloats_close(f_auto_unc(x, y, *args), f_wrapped(x, y, *args))

    ### Explicit derivatives:

    ## Fully defined derivatives:
    f_wrapped = uncert_core.wrap(
        f,
        [
            lambda x, y, *args: 2,
            lambda x, y, *args: math.cos(y),
            None,
            lambda x, y, *args: 3,
        ],
    )

    assert ufloats_close(f_auto_unc(x, y, *args), f_wrapped(x, y, *args))

    ## Automatic additional derivatives for non-defined derivatives:

    # No derivative for y:
    f_wrapped = uncert_core.wrap(f, [lambda x, y, *args: 2])
    assert ufloats_close(f_auto_unc(x, y, *args), f_wrapped(x, y, *args))


def test_wrapped_func_no_args_kwargs():
    """
    Wrap a function that takes only positional-or-keyword and
    var-keyword parameters.
    """

    def f_auto_unc(x, y, **kwargs):
        return 2 * x + umath.sin(y) + 3 * kwargs["z"]

    # Like f_auto_unc, but does not accept numbers with uncertainties:
    def f(x, y, **kwargs):
        assert not any(
            isinstance(value, uncert_core.UFloat)
            for value in [x, y] + list(kwargs.values())
        )
        return f_auto_unc(x, y, **kwargs)

    x = uncert_core.ufloat(1, 0.1)
    y = uncert_core.ufloat(10, 2)
    s = "string arg"
    z = uncert_core.ufloat(100, 3)

    kwargs = {"s": s, "z": z}  # Arguments not in signature

    ### Automatic numerical derivatives:

    ## Fully automatic numerical derivatives:
    f_wrapped = uncert_core.wrap(f)
    assert ufloats_close(f_auto_unc(x, y, **kwargs), f_wrapped(x, y, **kwargs))

    # Call with keyword arguments:
    assert ufloats_close(f_auto_unc(y=y, x=x, **kwargs), f_wrapped(y=y, x=x, **kwargs))

    ## Automatic additional derivatives for non-defined derivatives,
    ## and explicit None derivative:

    # No derivative for positional-or-keyword parameter y, no
    # derivative for optional-keyword parameter z:
    f_wrapped = uncert_core.wrap(f, [None])
    assert ufloats_close(f_auto_unc(x, y, **kwargs), f_wrapped(x, y, **kwargs))

    # Call with keyword arguments:
    assert ufloats_close(f_auto_unc(y=y, x=x, **kwargs), f_wrapped(y=y, x=x, **kwargs))

    # No derivative for positional-or-keyword parameter y, no
    # derivative for optional-keyword parameter z:
    f_wrapped = uncert_core.wrap(f, [None], {"z": None})
    assert ufloats_close(f_auto_unc(x, y, **kwargs), f_wrapped(x, y, **kwargs))

    # Call with keyword arguments:
    assert ufloats_close(f_auto_unc(y=y, x=x, **kwargs), f_wrapped(y=y, x=x, **kwargs))

    # No derivative for positional-or-keyword parameter y, derivative
    # for optional-keyword parameter z:
    f_wrapped = uncert_core.wrap(f, [None], {"z": lambda x, y, **kwargs: 3})
    assert ufloats_close(f_auto_unc(x, y, **kwargs), f_wrapped(x, y, **kwargs))

    # Call with keyword arguments:
    assert ufloats_close(f_auto_unc(y=y, x=x, **kwargs), f_wrapped(y=y, x=x, **kwargs))

    ### Explicit derivatives:

    ## Fully defined derivatives:
    f_wrapped = uncert_core.wrap(
        f,
        [lambda x, y, **kwargs: 2, lambda x, y, **kwargs: math.cos(y)],
        {"z:": lambda x, y, **kwargs: 3},
    )

    assert ufloats_close(f_auto_unc(x, y, **kwargs), f_wrapped(x, y, **kwargs))
    # Call with keyword arguments:
    assert ufloats_close(f_auto_unc(y=y, x=x, **kwargs), f_wrapped(y=y, x=x, **kwargs))

    ## Automatic additional derivatives for non-defined derivatives:

    # No derivative for y or z:
    f_wrapped = uncert_core.wrap(f, [lambda x, y, **kwargs: 2])
    assert ufloats_close(f_auto_unc(x, y, **kwargs), f_wrapped(x, y, **kwargs))

    # Call with keyword arguments:
    assert ufloats_close(f_auto_unc(y=y, x=x, **kwargs), f_wrapped(y=y, x=x, **kwargs))


def test_wrapped_func_args_kwargs():
    """
    Wrap a function that takes positional-or-keyword, var-positional
    and var-keyword parameters.
    """

    def f_auto_unc(x, y, *args, **kwargs):
        return 2 * x + umath.sin(y) + 4 * args[1] + 3 * kwargs["z"]

    # Like f_auto_unc, but does not accept numbers with uncertainties:
    def f(x, y, *args, **kwargs):
        assert not any(
            isinstance(value, uncert_core.UFloat)
            for value in [x, y] + list(args) + list(kwargs.values())
        )
        return f_auto_unc(x, y, *args, **kwargs)

    x = uncert_core.ufloat(1, 0.1)
    y = uncert_core.ufloat(10, 2)
    t = uncert_core.ufloat(1000, 4)
    s = "string arg"
    z = uncert_core.ufloat(100, 3)

    args = [s, t, s]
    kwargs = {"u": s, "z": z}  # Arguments not in signature

    ### Automatic numerical derivatives:

    ## Fully automatic numerical derivatives:
    f_wrapped = uncert_core.wrap(f)

    assert ufloats_close(
        f_auto_unc(x, y, *args, **kwargs),
        f_wrapped(x, y, *args, **kwargs),
        tolerance=1e-5,
    )

    ## Automatic additional derivatives for non-defined derivatives,
    ## and explicit None derivative:

    # No derivative for positional-or-keyword parameter y, no
    # derivative for optional-keyword parameter z:
    f_wrapped = uncert_core.wrap(f, [None, None, None, lambda x, y, *args, **kwargs: 4])
    assert ufloats_close(
        f_auto_unc(x, y, *args, **kwargs),
        f_wrapped(x, y, *args, **kwargs),
        tolerance=1e-5,
    )

    # No derivative for positional-or-keyword parameter y, no
    # derivative for optional-keyword parameter z:
    f_wrapped = uncert_core.wrap(f, [None], {"z": None})
    assert ufloats_close(
        f_auto_unc(x, y, *args, **kwargs),
        f_wrapped(x, y, *args, **kwargs),
        tolerance=1e-5,
    )

    # No derivative for positional-or-keyword parameter y, derivative
    # for optional-keyword parameter z:
    f_wrapped = uncert_core.wrap(f, [None], {"z": lambda x, y, *args, **kwargs: 3})
    assert ufloats_close(
        f_auto_unc(x, y, *args, **kwargs),
        f_wrapped(x, y, *args, **kwargs),
        tolerance=1e-5,
    )

    ### Explicit derivatives:

    ## Fully defined derivatives:
    f_wrapped = uncert_core.wrap(
        f,
        [lambda x, y, *args, **kwargs: 2, lambda x, y, *args, **kwargs: math.cos(y)],
        {"z:": lambda x, y, *args, **kwargs: 3},
    )

    assert ufloats_close(
        f_auto_unc(x, y, *args, **kwargs),
        f_wrapped(x, y, *args, **kwargs),
        tolerance=1e-5,
    )

    ## Automatic additional derivatives for non-defined derivatives:

    # No derivative for y or z:
    f_wrapped = uncert_core.wrap(f, [lambda x, y, *args, **kwargs: 2])
    assert ufloats_close(
        f_auto_unc(x, y, *args, **kwargs),
        f_wrapped(x, y, *args, **kwargs),
        tolerance=1e-5,
    )


def test_wrapped_func():
    """
    Test uncertainty-aware functions obtained through wrapping.
    """

    ########################################

    # Function which can automatically handle numbers with
    # uncertainties:
    def f_auto_unc(angle, *list_var):
        return umath.cos(angle) + sum(list_var)

    def f(angle, *list_var):
        # We make sure that this function is only ever called with
        # numbers with no uncertainty (since it is wrapped):
        assert not isinstance(angle, uncert_core.UFloat)
        assert not any(isinstance(arg, uncert_core.UFloat) for arg in list_var)
        return f_auto_unc(angle, *list_var)

    f_wrapped = uncert_core.wrap(f)

    my_list = [1, 2, 3]

    ########################################
    # Test of a wrapped function that only calls the original
    # function: it should obtain the exact same result:
    assert f_wrapped(0, *my_list) == f(0, *my_list)
    # 1 == 1 +/- 0, so the type must be checked too:
    assert isinstance(f_wrapped(0, *my_list), type(f(0, *my_list)))

    ########################################
    # Call with uncertainties:

    angle = uncert_core.ufloat(1, 0.1)
    list_value = uncert_core.ufloat(3, 0.2)

    # The random variables must be the same (full correlation):

    assert ufloats_close(f_wrapped(angle, *[1, angle]), f_auto_unc(angle, *[1, angle]))

    assert ufloats_close(
        f_wrapped(angle, *[list_value, angle]), f_auto_unc(angle, *[list_value, angle])
    )

    ########################################
    # Non-numerical arguments, and  explicit and implicit derivatives:
    def f(x, y, z, t, u):
        return x + 2 * z + 3 * t + 4 * u

    f_wrapped = uncert_core.wrap(
        f, [lambda *args: 1, None, lambda *args: 2, None]
    )  # No deriv. for u

    assert f_wrapped(10, "string argument", 1, 0, 0) == 12

    x = uncert_core.ufloat(10, 1)

    assert numbers_close(
        f_wrapped(x, "string argument", x, x, x).std_dev, (1 + 2 + 3 + 4) * x.std_dev
    )


def test_wrap_with_kwargs():
    """
    Tests wrap() on functions with keyword arguments.

    Includes both wrapping a function that takes optional keyword
    arguments and calling a wrapped function with keyword arguments
    (optional or not).
    """

    # Version of f() that automatically works with numbers with
    # uncertainties:
    def f_auto_unc(x, y, *args, **kwargs):
        return x + umath.sin(y) + 2 * args[0] + 3 * kwargs["t"]

    # We also add keyword arguments in the function which is wrapped:
    def f(x, y, *args, **kwargs):
        # We make sure that f is not called directly with a number with
        # uncertainty:

        for value in [x, y] + list(args) + list(kwargs.values()):
            assert not isinstance(value, uncert_core.UFloat)

        return f_auto_unc(x, y, *args, **kwargs)

    f_wrapped = uncert_core.wrap(f)

    x = ufloat(1, 0.1)
    y = ufloat(10, 0.11)
    z = ufloat(100, 0.111)
    t = ufloat(0.1, 0.1111)

    z_uatom = get_single_uatom(z)
    t_uatom = get_single_uatom(t)

    assert ufloats_close(
        f_wrapped(x, y, z, t=t), f_auto_unc(x, y, z, t=t), tolerance=1e-5
    )

    ########################################

    # We make sure that analytical derivatives are indeed used. We
    # also test the automatic handling of additional *args arguments
    # beyond the number of supplied derivatives.

    f_wrapped2 = uncert_core.wrap(f, [None, lambda x, y, *args, **kwargs: math.cos(y)])

    # The derivatives must be perfectly identical:

    # The *args parameter of f() is given as a keyword argument, so as
    # to try to confuse the code:

    assert (
        f_wrapped2(x, y, z, t=t).error_components[z_uatom]
        == f_auto_unc(x, y, z, t=t).error_components[z_uatom]
    )

    # Derivatives supplied through the keyword-parameter dictionary of
    # derivatives, and also derivatives supplied for the
    # var-positional arguments (*args[0]):

    f_wrapped3 = uncert_core.wrap(
        f,
        [None, None, lambda x, y, *args, **kwargs: 2],
        {"t": lambda x, y, *args, **kwargs: 3},
    )

    # The derivatives should be exactly the same, because they are
    # obtained with the exact same analytic formula:
    assert (
        f_wrapped3(x, y, z, t=t).error_components[z_uatom]
        == f_auto_unc(x, y, z, t=t).error_components[z_uatom]
    )
    assert (
        f_wrapped3(x, y, z, t=t).error_components[t_uatom]
        == f_auto_unc(x, y, z, t=t).error_components[t_uatom]
    )

    ########################################
    # Making sure that user-supplied derivatives are indeed called:

    class FunctionCalled(Exception):
        """
        Raised to signal that a function is indeed called.
        """

        pass

    def failing_func(x, y, *args, **kwargs):
        raise FunctionCalled

    f_wrapped4 = uncert_core.wrap(f, [None, failing_func], {"t": failing_func})

    try:
        f_wrapped4(x, 3.14, z, t=t)
    except FunctionCalled:
        pass
    else:
        raise Exception("User-supplied derivative should be called")

    try:
        f_wrapped4(x, y, z, t=3.14)
    except FunctionCalled:
        pass
    else:
        raise Exception("User-supplied derivative should be called")

    try:
        f_wrapped4(x, 3.14, z, t=3.14)
    except FunctionCalled:
        raise Exception("User-supplied derivative should *not* be called")


###############################################################################


def test_access_to_std_dev():
    "Uniform access to the standard deviation"

    x = ufloat(1, 0.1)
    y = 2 * x

    # std_dev for Variable and AffineScalarFunc objects:
    assert uncert_core.std_dev(x) == x.std_dev
    assert uncert_core.std_dev(y) == y.std_dev

    # std_dev for other objects:
    assert uncert_core.std_dev([]) == 0
    assert uncert_core.std_dev(None) == 0


###############################################################################


def test_covariances():
    "Covariance matrix"

    x = ufloat(1, 0.1)
    y = -2 * x + 10
    z = -3 * x
    covs = uncert_core.covariance_matrix([x, y, z])
    # Diagonal elements are simple:
    assert numbers_close(covs[0][0], 0.01)
    assert numbers_close(covs[1][1], 0.04)
    assert numbers_close(covs[2][2], 0.09)
    # Non-diagonal elements:
    assert numbers_close(covs[0][1], -0.02)


###############################################################################

# The tests below require NumPy, which is an optional package:
try:
    import numpy
    from helpers import uarrays_close
except ImportError:
    pass
else:

    def test_numpy_comparison():
        "Comparison with a NumPy array."

        x = ufloat(1, 0.1)

        # Comparison with a different type:
        assert x != [x, x]

        # NumPy arrays can be compared, through element-wise
        # comparisons.  Numbers with uncertainties should yield the
        # same kind of results as pure floats (i.e., a NumPy array,
        # etc.).

        # We test the comparison operators both for the uncertainties
        # package *and* the NumPy package:

        # Equalities, etc.:
        assert len(x == numpy.arange(10)) == 10
        assert len(numpy.arange(10) == x) == 10
        assert len(x != numpy.arange(10)) == 10
        assert len(numpy.arange(10) != x) == 10
        assert len(x == numpy.array([x, x, x])) == 3
        assert len(numpy.array([x, x, x]) == x) == 3
        assert numpy.all(x == numpy.array([x, x, x]))

        # Inequalities:
        assert len(x < numpy.arange(10)) == 10
        assert len(numpy.arange(10) > x) == 10
        assert len(x <= numpy.arange(10)) == 10
        assert len(numpy.arange(10) >= x) == 10
        assert len(x > numpy.arange(10)) == 10
        assert len(numpy.arange(10) < x) == 10
        assert len(x >= numpy.arange(10)) == 10
        assert len(numpy.arange(10) <= x) == 10

        # More detailed test, that shows that the comparisons are
        # meaningful (x >= 0, but not x <= 1):
        assert numpy.all((x >= numpy.arange(3)) == [True, False, False])

    def test_correlated_values():
        """
        Correlated variables.
        Test through the input of the (full) covariance matrix.
        """

        u = uncert_core.ufloat(1, 0.1)
        cov = uncert_core.covariance_matrix([u])
        # "1" is used instead of u.nominal_value because
        # u.nominal_value might return a float.  The idea is to force
        # the new variable u2 to be defined through an integer nominal
        # value:
        (u2,) = uncert_core.correlated_values([1], cov)
        expr = 2 * u2  # Calculations with u2 should be possible, like with u # noqa

        ####################

        # Covariances between output and input variables:

        x = ufloat(1, 0.1)
        y = ufloat(2, 0.3)
        z = -3 * x + y

        covs = uncert_core.covariance_matrix([x, y, z])

        # Test of the diagonal covariance elements:
        assert uarrays_close(
            numpy.array([v.std_dev**2 for v in (x, y, z)]), numpy.array(covs).diagonal()
        )

        # "Inversion" of the covariance matrix: creation of new
        # variables:
        (x_new, y_new, z_new) = uncert_core.correlated_values(
            [x.nominal_value, y.nominal_value, z.nominal_value],
            covs,
            tags=["x", "y", "z"],
        )

        # Even the uncertainties should be correctly reconstructed:
        assert uarrays_close(numpy.array((x, y, z)), numpy.array((x_new, y_new, z_new)))

        # ... and the covariances too:
        assert uarrays_close(
            numpy.array(covs),
            numpy.array(uncert_core.covariance_matrix([x_new, y_new, z_new])),
        )

        assert uarrays_close(numpy.array([z_new]), numpy.array([-3 * x_new + y_new]))

        ####################

        # ... as well as functional relations:

        u = ufloat(1, 0.05)
        v = ufloat(10, 0.1)
        sum_value = u + 2 * v

        # Covariance matrices:
        cov_matrix = uncert_core.covariance_matrix([u, v, sum_value])

        # Correlated variables can be constructed from a covariance
        # matrix, if NumPy is available:
        (u2, v2, sum2) = uncert_core.correlated_values(
            [x.nominal_value for x in [u, v, sum_value]], cov_matrix
        )

        # uarrays_close() is used instead of numbers_close() because
        # it compares uncertainties too:
        assert uarrays_close(numpy.array([u]), numpy.array([u2]))
        assert uarrays_close(numpy.array([v]), numpy.array([v2]))
        assert uarrays_close(numpy.array([sum_value]), numpy.array([sum2]))
        assert uarrays_close(numpy.array([0]), numpy.array([sum2 - (u2 + 2 * v2)]))

        # Spot checks of the correlation matrix:
        corr_matrix = uncert_core.correlation_matrix([u, v, sum_value])
        assert numbers_close(corr_matrix[0, 0], 1)
        assert numbers_close(corr_matrix[1, 2], 2 * v.std_dev / sum_value.std_dev)

        ####################

        # Test of numerical robustness despite wildly different
        # orders of magnitude (see
        # https://github.com/lebigot/uncertainties/issues/95):
        cov = numpy.diag([1e-70, 1e-70, 1e10])
        cov[0, 1] = cov[1, 0] = 0.9e-70
        cov[[0, 1], 2] = -3e-34
        cov[2, [0, 1]] = -3e-34
        variables = uncert_core.correlated_values([0] * 3, cov)

        # Since the numbers are very small, we need to compare them
        # in a stricter way, that handles the case of a 0 variance
        # in `variables`:
        assert numbers_close(
            1e66 * cov[0, 0], 1e66 * variables[0].s ** 2, tolerance=1e-5
        )
        assert numbers_close(
            1e66 * cov[1, 1], 1e66 * variables[1].s ** 2, tolerance=1e-5
        )

        ####################

        # 0 variances are a bit special, since the correlation matrix
        # cannot be calculated naively, so we test that there is no
        # specific problem in this case:

        cov = numpy.diag([0, 0, 10])
        nom_values = [1, 2, 3]
        variables = uncert_core.correlated_values(nom_values, cov)

        for variable, nom_value, variance in zip(variables, nom_values, cov.diagonal()):
            assert numbers_close(variable.n, nom_value)
            assert numbers_close(variable.s**2, variance)

        assert uarrays_close(cov, numpy.array(uncert_core.covariance_matrix(variables)))

    def test_correlated_values_correlation_mat():
        """
        Tests the input of correlated value.

        Test through their correlation matrix (instead of the
        covariance matrix).
        """

        x = ufloat(1, 0.1)
        y = ufloat(2, 0.3)
        z = -3 * x + y

        cov_mat = uncert_core.covariance_matrix([x, y, z])

        std_devs = numpy.sqrt(numpy.array(cov_mat).diagonal())

        corr_mat = cov_mat / std_devs / std_devs[numpy.newaxis].T

        # We make sure that the correlation matrix is indeed diagonal:
        assert (corr_mat - corr_mat.T).max() <= 1e-15
        # We make sure that there are indeed ones on the diagonal:
        assert (corr_mat.diagonal() - 1).max() <= 1e-15

        # We try to recover the correlated variables through the
        # correlation matrix (not through the covariance matrix):

        nominal_values = [v.nominal_value for v in (x, y, z)]
        std_devs = [v.std_dev for v in (x, y, z)]
        x2, y2, z2 = uncert_core.correlated_values_norm(
            list(zip(nominal_values, std_devs)), corr_mat
        )

        # uarrays_close() is used instead of numbers_close() because
        # it compares uncertainties too:

        # Test of individual variables:
        assert uarrays_close(numpy.array([x]), numpy.array([x2]))
        assert uarrays_close(numpy.array([y]), numpy.array([y2]))
        assert uarrays_close(numpy.array([z]), numpy.array([z2]))

        # Partial correlation test:
        assert uarrays_close(numpy.array([0]), numpy.array([z2 - (-3 * x2 + y2)]))

        # Test of the full covariance matrix:
        assert uarrays_close(
            numpy.array(cov_mat),
            numpy.array(uncert_core.covariance_matrix([x2, y2, z2])),
        )


@pytest.mark.skipif(
    np is not None,
    reason="This test is only run when numpy is not installed.",
)
def test_no_numpy():
    nom_values = [1, 2, 3]
    std_devs = [0.1, 0.2, 0.3]
    cov = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ]

    with pytest.raises(
        NotImplementedError,
        match="not able to import numpy",
    ):
        _ = correlated_values(nom_values, cov)

    with pytest.raises(
        NotImplementedError,
        match="not able to import numpy",
    ):
        _ = correlated_values_norm(
            list(zip(nom_values, std_devs)),
            cov,
        )

    x = ufloat(1, 0.1)
    y = ufloat(2, 0.2)
    z = ufloat(3, 0.3)

    with pytest.raises(
        NotImplementedError,
        match="not able to import numpy",
    ):
        _ = correlation_matrix([x, y, z])


@pytest.mark.parametrize("method_name", deprecated_methods)
def test_deprecated_method(method_name):
    x = ufloat(1, 0.1)
    y = ufloat(-12, 2.4)
    num_args = len(inspect.signature(getattr(float, method_name)).parameters)
    with pytest.warns(FutureWarning, match="will be removed"):
        if num_args == 1:
            getattr(x, method_name)()
        else:
            getattr(x, method_name)(y)
