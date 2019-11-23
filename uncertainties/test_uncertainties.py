# coding=utf-8

"""
Tests of the code in uncertainties/__init__.py.

These tests can be run through the Nose testing framework.

(c) 2010-2016 by Eric O. LEBIGOT (EOL).
"""




# Standard modules
import copy
import weakref
import math
from math import isnan, isinf
import random
import sys

# 3rd-party modules
# import nose.tools

# Local modules

import uncertainties.core as uncert_core
from uncertainties.core import ufloat, AffineScalarFunc, ufloat_fromstr
from uncertainties import umath

# The following information is useful for making sure that the right
# version of Python is running the tests (for instance with the Travis
# Continuous Integration system):
print("Testing with Python", sys.version)

###############################################################################

# Utilities for unit testing

def numbers_close(x, y, tolerance=1e-6):
    """
    Returns True if the given floats are close enough.

    The given tolerance is the relative difference allowed, or the absolute
    difference, if one of the numbers is 0.

    NaN is allowed: it is considered close to itself.
    """

    # !!! Python 3.5+ has math.isclose(): maybe it could be used here.

    # Instead of using a try and ZeroDivisionError, we do a test,
    # NaN could appear silently:

    if x != 0 and y != 0:
        if isinf(x):
            return isinf(y)
        elif isnan(x):
            return isnan(y)
        else:
            # Symmetric form of the test:
            return 2*abs(x-y)/(abs(x)+abs(y)) < tolerance

    else:  # Either x or y is zero
        return abs(x or y) < tolerance

def ufloats_close(x, y, tolerance=1e-6):
    '''
    Tests if two numbers with uncertainties are close, as random
    variables: this is stronger than testing whether their nominal
    value and standard deviation are close.

    The tolerance is applied to both the nominal value and the
    standard deviation of the difference between the numbers.
    '''

    diff = x-y
    return (numbers_close(diff.nominal_value, 0, tolerance)
            and numbers_close(diff.std_dev, 0, tolerance))

class DerivativesDiffer(Exception):
    pass


def compare_derivatives(func, numerical_derivatives,
                         num_args_list=None):
    """
    Checks the derivatives of a function 'func' (as returned by the
    wrap() wrapper), by comparing them to the
    'numerical_derivatives' functions.

    Raises a DerivativesDiffer exception in case of problem.

    These functions all take the number of arguments listed in
    num_args_list.  If num_args is None, it is automatically obtained.

    Tests are done on random arguments.
    """

    try:
        funcname = func.name
    except AttributeError:
        funcname = func.__name__

    # print "Testing", func.__name__

    if not num_args_list:

        # Detecting automatically the correct number of arguments is not
        # always easy (because not all values are allowed, etc.):

        num_args_table = {
            'atanh': [1],
            'log': [1, 2]  # Both numbers of arguments are tested
            }
        if funcname in num_args_table:
            num_args_list = num_args_table[funcname]
        else:

            num_args_list = []

            # We loop until we find reasonable function arguments:
            # We get the number of arguments by trial and error:
            for num_args in range(10):
                try:
                    #! Giving integer arguments is good for preventing
                    # certain functions from failing even though num_args
                    # is their correct number of arguments
                    # (e.g. math.ldexp(x, i), where i must be an integer)
                    func(*(1,)*num_args)
                except TypeError:
                    pass  # Not the right number of arguments
                else:  # No error
                    # num_args is a good number of arguments for func:
                    num_args_list.append(num_args)

            if not num_args_list:
                raise Exception("Can't find a reasonable number of arguments"
                                " for function '%s'." % funcname)

    for num_args in num_args_list:

        # Argument numbers that will have a random integer value:
        integer_arg_nums = set()

        if funcname == 'ldexp':
            # The second argument must be an integer:
            integer_arg_nums.add(1)

        while True:
            try:

                # We include negative numbers, for more thorough tests:
                args = []
                for arg_num in range(num_args):
                    if arg_num in integer_arg_nums:
                        args.append(random.choice(list(range(-10, 10))))
                    else:
                        args.append(
                            uncert_core.Variable(random.random()*4-2, 0))

                # 'args', but as scalar values:
                args_scalar = [uncert_core.nominal_value(v)
                               for v in args]

                func_approx = func(*args)

                # Some functions yield simple Python constants, after
                # wrapping in wrap(): no test has to be performed.
                # Some functions also yield tuples...
                if isinstance(func_approx, AffineScalarFunc):

                    # We compare all derivatives:
                    for (arg_num, (arg, numerical_deriv)) in (
                        enumerate(zip(args, numerical_derivatives))):

                        # Some arguments might not be differentiable:
                        if isinstance(arg, int):
                            continue

                        fixed_deriv_value = func_approx.derivatives[arg]

                        num_deriv_value = numerical_deriv(*args_scalar)

                        # This message is useful: the user can see that
                        # tests are really performed (instead of not being
                        # performed, silently):
                        print("Testing derivative #%d of %s at %s" % (
                            arg_num, funcname, args_scalar))

                        if not numbers_close(fixed_deriv_value,
                                             num_deriv_value, 1e-4):

                            # It is possible that the result is NaN:
                            if not isnan(func_approx):
                                raise DerivativesDiffer(
                                    "Derivative #%d of function '%s' may be"
                                    " wrong: at args = %s,"
                                    " value obtained = %.16f,"
                                    " while numerical approximation = %.16f."
                                    % (arg_num, funcname, args,
                                       fixed_deriv_value, num_deriv_value))

            except ValueError as err:  # Arguments out of range, or of wrong type
                # Factorial(real) lands here:
                if str(err).startswith('factorial'):
                    integer_arg_nums = set([0])
                continue  # We try with different arguments
            # Some arguments might have to be integers, for instance:
            except TypeError as err:
                if len(integer_arg_nums) == num_args:
                    raise Exception("Incorrect testing procedure: unable to "
                                    "find correct argument values for %s: %s"
                                    % (funcname, err))

                # Another argument might be forced to be an integer:
                integer_arg_nums.add(random.choice(list(range(num_args))))
            else:
                # We have found reasonable arguments, and the test passed:
                break

###############################################################################

def test_value_construction():
    '''
    Tests the various means of constructing a constant number with
    uncertainty *without a string* (see test_ufloat_fromstr(), for this).
    '''

    ## Simple construction:
    x = ufloat(3, 0.14)
    assert x.nominal_value == 3
    assert x.std_dev == 0.14
    assert x.tag is None

    # ... with tag as positional argument:
    x = ufloat(3, 0.14, 'pi')
    assert x.nominal_value == 3
    assert x.std_dev == 0.14
    assert x.tag == 'pi'

    # ... with tag keyword:
    x = ufloat(3, 0.14, tag='pi')
    assert x.nominal_value == 3
    assert x.std_dev == 0.14
    assert x.tag == 'pi'

    ## Comparison with the obsolete tuple form:

    # The following tuple is stored in a variable instead of being
    # repeated in the calls below, so that the automatic code update
    # does not replace ufloat((3, 0.14)) by ufloat(3, 14): the goal
    # here is to make sure that the obsolete form gives the same
    # result as the new form.

    representation = (3, 0.14)  # Obsolete representation

    x = ufloat(3, 0.14)
    x2 = ufloat(representation)  # Obsolete
    assert x.nominal_value == x2.nominal_value
    assert x.std_dev == x2.std_dev
    assert x.tag is None
    assert x2.tag is None

    # With tag as positional argument:
    x = ufloat(3, 0.14, "pi")
    x2 = ufloat(representation, "pi")  # Obsolete
    assert x.nominal_value == x2.nominal_value
    assert x.std_dev == x2.std_dev
    assert x.tag == 'pi'
    assert x2.tag == 'pi'

    # With tag keyword:
    x = ufloat(3, 0.14, tag="pi")
    x2 = ufloat(representation, tag="pi")  # Obsolete
    assert x.nominal_value == x2.nominal_value
    assert x.std_dev == x2.std_dev
    assert x.tag == 'pi'
    assert x2.tag == 'pi'

    # Negative standard deviations should be caught in a nice way
    # (with the right exception):
    try:
        x = ufloat(3, -0.1)
    except uncert_core.NegativeStdDev:
        pass

    try:
        # Obsolete form:
        x = ufloat((3, -0.1))
    except uncert_core.NegativeStdDev:
        pass

    ## Incorrect forms should not raise any deprecation warning, but
    ## raise an exception:

    try:
        ufloat(1)  # Form that has never been allowed
    except:
        pass
    else:
        raise Exception("An exception should be raised")

def test_ufloat_fromstr():
    "Input of numbers with uncertainties as a string"

    # String representation, and numerical values:
    tests = {
        "-1.23(3.4)": (-1.23, 3.4),  # (Nominal value, error)
        "  -1.23(3.4)  ": (-1.23, 3.4),  # Spaces ignored
        "-1.34(5)": (-1.34, 0.05),
        "1(6)": (1, 6),
        "3(4.2)": (3, 4.2),
        "-9(2)": (-9, 2),
        "1234567(1.2)": (1234567, 1.2),
        "12.345(15)": (12.345, 0.015),
        "-12.3456(78)e-6": (-12.3456e-6, 0.0078e-6),
        "0.29": (0.29, 0.01),
        "31.": (31, 1),
        "-31.": (-31, 1),
        # The following tests that the ufloat() routine does
        # not consider '31' like the tuple ('3', '1'), which would
        # make it expect two numbers (instead of 2 1-character
        # strings):
        "31": (31, 1),
        "-3.1e10": (-3.1e10, 0.1e10),
        "169.0(7)": (169, 0.7),
        "-0.1+/-1": (-0.1, 1),
        "-13e-2+/-1e2": (-13e-2, 1e2),
        '-14.(15)': (-14, 15),
        '-100.0(15)': (-100, 1.5),
        '14.(15)': (14, 15),
        # Global exponent:
        '(3.141+/-0.001)E+02': (314.1, 0.1),


        ## Pretty-print notation:

        # ± sign, global exponent (not pretty-printed):
        '(3.141±0.001)E+02': (314.1, 0.1),
        # ± sign, individual exponent:
        '3.141E+02±0.001e2': (314.1, 0.1),

        # ± sign, times symbol, superscript (= full pretty-print):
        '(3.141 ± 0.001) × 10²': (314.1, 0.1),

        # NaN uncertainty:
        '(3.141±nan)E+02': (314.1, float('nan')),
        '3.141e+02+/-nan': (314.1, float('nan')),
        '3.4(nan)e10': (3.4e10, float('nan')),
        # NaN value:
        'nan+/-3.14e2': (float('nan'), 314),
        # "Double-floats"
        '(-3.1415 +/- 1e-4)e+200': (-3.1415e200, 1e196),
        '(-3.1415e-10 +/- 1e-4)e+200': (-3.1415e190, 1e196),
        # Special float representation:
        '-3(0.)': (-3, 0)
        }

    for (representation, values) in tests.items():

        # Without tag:
        num = ufloat_fromstr(representation)
        assert numbers_close(num.nominal_value, values[0])
        assert numbers_close(num.std_dev, values[1])
        assert num.tag is None

        # With a tag as positional argument:
        num = ufloat_fromstr(representation, 'test variable')
        assert numbers_close(num.nominal_value, values[0])
        assert numbers_close(num.std_dev, values[1])
        assert num.tag == 'test variable'

        # With a tag as keyword argument:
        num = ufloat_fromstr(representation, tag='test variable')
        assert numbers_close(num.nominal_value, values[0])
        assert numbers_close(num.std_dev, values[1])
        assert num.tag == 'test variable'

        ## Obsolete forms

        num = ufloat(representation)  # Obsolete
        assert numbers_close(num.nominal_value, values[0])
        assert numbers_close(num.std_dev, values[1])
        assert num.tag is None

        # Call with a tag list argument:
        num = ufloat(representation, 'test variable')  # Obsolete
        assert numbers_close(num.nominal_value, values[0])
        assert numbers_close(num.std_dev, values[1])
        assert num.tag == 'test variable'

        # Call with a tag keyword argument:
        num = ufloat(representation, tag='test variable')  # Obsolete
        assert numbers_close(num.nominal_value, values[0])
        assert numbers_close(num.std_dev, values[1])
        assert num.tag == 'test variable'

###############################################################################

# Test of correctness of the fixed (usually analytical) derivatives:
def test_fixed_derivatives_basic_funcs():
    """
    Pre-calculated derivatives for operations on AffineScalarFunc.
    """

    def check_op(op, num_args):
        """
        Makes sure that the derivatives for function '__op__' of class
        AffineScalarFunc, which takes num_args arguments, are correct.

        If num_args is None, a correct value is calculated.
        """

        op_string = "__%s__" % op
        func = getattr(AffineScalarFunc, op_string)
        numerical_derivatives = uncert_core.NumericalDerivatives(
            # The __neg__ etc. methods of AffineScalarFunc only apply,
            # by definition, to AffineScalarFunc objects: we first map
            # possible scalar arguments (used for calculating
            # derivatives) to AffineScalarFunc objects:
            lambda *args: func(*list(map(uncert_core.to_affine_scalar, args))))
        compare_derivatives(func, numerical_derivatives, [num_args])

    # Operators that take 1 value:
    for op in uncert_core.modified_operators:
        check_op(op, 1)

    # Operators that take 2 values:
    for op in uncert_core.modified_ops_with_reflection:
        check_op(op, 2)

# Additional, more complex checks, for use with the nose unit testing
# framework.

def test_copy():
    "Standard copy module integration"
    import gc

    x = ufloat(3, 0.1)
    assert x == x

    y = copy.copy(x)
    assert x != y
    assert not(x == y)
    assert y in list(y.derivatives.keys())  # y must not copy the dependence on x

    z = copy.deepcopy(x)
    assert x != z

    # Copy tests on expressions:
    t = x + 2*z
    # t depends on x:
    assert x in t.derivatives

    # The relationship between the copy of an expression and the
    # original variables should be preserved:
    t_copy = copy.copy(t)
    # Shallow copy: the variables on which t depends are not copied:
    assert x in t_copy.derivatives
    assert (uncert_core.covariance_matrix([t, z]) ==
            uncert_core.covariance_matrix([t_copy, z]))

    # However, the relationship between a deep copy and the original
    # variables should be broken, since the deep copy created new,
    # independent variables:
    t_deepcopy = copy.deepcopy(t)
    assert x not in t_deepcopy.derivatives
    assert (uncert_core.covariance_matrix([t, z]) !=
            uncert_core.covariance_matrix([t_deepcopy, z]))

    # Test of implementations with weak references:

    # Weak references: destroying a variable should never destroy the
    # integrity of its copies (which would happen if the copy keeps a
    # weak reference to the original, in its derivatives member: the
    # weak reference to the original would become invalid):
    del x

    gc.collect()

    assert y in list(y.derivatives.keys())

## Classes for the pickling tests (put at the module level, so that
## they can be unpickled):

# Subclass without slots:
class NewVariable_dict(uncert_core.Variable):
    pass

# Subclass with slots defined by a tuple:
class NewVariable_slots_tuple(uncert_core.Variable):
    __slots__ = ('new_attr',)

# Subclass with slots defined by a string:
class NewVariable_slots_str(uncert_core.Variable):
    __slots__ = 'new_attr'

def test_pickling():
    "Standard pickle module integration."

    import pickle

    x = ufloat(2, 0.1)

    x_unpickled = pickle.loads(pickle.dumps(x))

    assert x != x_unpickled  # Pickling creates copies

    ## Tests with correlations and AffineScalarFunc objects:
    f = 2*x
    assert isinstance(f, AffineScalarFunc)
    (f_unpickled, x_unpickled2) = pickle.loads(pickle.dumps((f, x)))
    # Correlations must be preserved:
    assert f_unpickled - x_unpickled2 - x_unpickled2 == 0

    ## Tests with subclasses:

    for subclass in (NewVariable_dict, NewVariable_slots_tuple,
                     NewVariable_slots_str):

        x = subclass(3, 0.14)

        # Pickling test with possibly uninitialized slots:
        pickle.loads(pickle.dumps(x))

        # Unpickling test:
        x.new_attr = 'New attr value'
        x_unpickled = pickle.loads(pickle.dumps(x))
        # Must exist (From the slots of the parent class):
        x_unpickled.nominal_value
        x_unpickled.new_attr  # Must exist

    ##

    # Corner case test: when an attribute is present both in __slots__
    # and in __dict__, it is first looked up from the slots
    # (references:
    # http://docs.python.org/2/reference/datamodel.html#invoking-descriptors,
    # http://stackoverflow.com/a/15139208/42973). As a consequence,
    # the pickling process must pickle the correct value (i.e., not
    # the value from __dict__):
    x = NewVariable_dict(3, 0.14)
    x._nominal_value = 'in slots'
    # Corner case: __dict__ key which is also a slot name (it is
    # shadowed by the corresponding slot, so this is very unusual,
    # though):
    x.__dict__['_nominal_value'] = 'in dict'
    # Additional __dict__ attribute:
    x.dict_attr = 'dict attribute'

    x_unpickled = pickle.loads(pickle.dumps(x))
    # We make sure that the data is still there and untouched:
    assert x_unpickled._nominal_value == 'in slots'
    assert x_unpickled.__dict__ == x.__dict__

    ##

    # Corner case that should have no impact on the code but which is
    # not prevented by the documentation: case of constant linear
    # terms (the potential gotcha is that if the linear_combo
    # attribute is empty, __getstate__()'s result could be false, and
    # so __setstate__() would not be called and the original empty
    # linear combination would not be set in linear_combo.
    x = uncert_core.LinearCombination({})
    assert pickle.loads(pickle.dumps(x)).linear_combo == {}

def test_int_div():
    "Integer division"
    # We perform all operations on floats, because derivatives can
    # otherwise be meaningless:
    x = ufloat(3.9, 2)//2
    assert x.nominal_value == 1.
    # All errors are supposed to be small, so the ufloat()
    # in x violates the assumption.  Therefore, the following is
    # correct:
    assert x.std_dev == 0.0

def test_comparison_ops():
    "Test of comparison operators"

    import random

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
    assert not(x >= 3)
    assert not(x < 3)

    assert x == x
    # Comparaison between Variable and AffineScalarFunc:
    assert x == x + 0
    # Comparaison between 2 _different_ AffineScalarFunc objects
    # representing the same value:
    assert x/2 == x/2
    # With uncorrelated result that have the same behavior (value and
    # standard error):
    assert 2*ufloat(1, 0.1) != ufloat(2, 0.2)
    # Comparaison between 2 _different_ Variable objects
    # that are uncorrelated:
    assert x != ufloat(3, 0.1)

    assert x != ufloat(3, 0.2)

    # Comparison to other types should work:
    assert x != None  # Not comparable
    assert x-x == 0  # Comparable, even though the types are different
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

        import random

        def random_float(var):
            """
            Returns a random value for Variable var, in an
            infinitesimal interval withing its uncertainty.  The case
            of a zero uncertainty is special.
            """
            return ((random.random()-0.5) * min(var.std_dev, 1e-5)
                    + var.nominal_value)

        # All operations are tested:
        for op in ["__%s__" % name
                   for name in('ne', 'eq', 'lt', 'le', 'gt', 'ge')]:

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
                sampled_results.append(float_func(random_float(x),
                                                  random_float(y)))

            min_result = min(sampled_results)
            max_result = max(sampled_results)

            if min_result == max_result:
                correct_result = min_result
            else:

                # Almost all results must be True, for the final value
                # to be True:
                num_min_result = sampled_results.count(min_result)

                # 1 exception is considered OK:
                correct_result = (num_min_result == 1)

            try:
                assert correct_result == getattr(x, op)(y)
            except AssertionError:
                print("Sampling results:", sampled_results)
                raise Exception("Semantic value of %s %s (%s) %s not"
                                " correctly reproduced."
                                % (x, op, y, correct_result))

    # With different numbers:
    test_all_comparison_ops(ufloat(3, 0.1),
                            ufloat(-2, 0.1))
    test_all_comparison_ops(ufloat(0, 0),  # Special number
                            ufloat(1, 1))
    test_all_comparison_ops(ufloat(0, 0),  # Special number
                            ufloat(0, 0.1))
    # With identical numbers:
    test_all_comparison_ops(ufloat(0, 0),
                            ufloat(0, 0))
    test_all_comparison_ops(ufloat(1, 1),
                            ufloat(1, 1))


def test_logic():
    "Boolean logic: __nonzero__, bool."

    x = ufloat(3, 0)
    y = ufloat(0, 0)
    z = ufloat(0, 0.1)
    t = ufloat(-1, 2)

    assert bool(x) == True
    assert bool(y) == False
    assert bool(z) == True
    assert bool(t) == True  # Only infinitseimal neighborhood are used

def test_obsolete():
    'Tests some obsolete creation of number with uncertainties'
    x = ufloat(3, 0.1)
    # Obsolete function, protected against automatic modification:
    x.set_std_dev.__call__(0.2)  # Obsolete

    x_std_dev = x.std_dev
    assert x_std_dev() == 0.2  # Obsolete call

def test_basic_access_to_data():
    "Access to data from Variable and AffineScalarFunc objects."

    x = ufloat(3.14, 0.01, "x var")
    assert x.tag == "x var"
    assert x.nominal_value == 3.14
    assert x.std_dev == 0.01

    # Case of AffineScalarFunc objects:
    y = x + 0
    assert type(y) == AffineScalarFunc
    assert y.nominal_value == 3.14
    assert y.std_dev == 0.01

    # Details on the sources of error:
    a = ufloat(-1, 0.001)
    y = 2*x + 3*x + 2 + a
    error_sources = y.error_components()
    assert len(error_sources) == 2  # 'a' and 'x'
    assert error_sources[x] == 0.05
    assert error_sources[a] == 0.001

    # Derivative values should be available:
    assert y.derivatives[x] == 5

    # Modification of the standard deviation of variables:
    x.std_dev = 1
    assert y.error_components()[x] == 5  # New error contribution!

    # Calculated values with uncertainties should not have a settable
    # standard deviation:
    y = 2*x
    try:
        y.std_dev = 1
    except AttributeError:
        pass
    else:
        raise Exception(
            "std_dev should not be settable for calculated results")

    # Calculation of deviations in units of the standard deviations:
    assert 10/x.std_dev == x.std_score(10 + x.nominal_value)

    # "In units of the standard deviation" is not always meaningful:
    x.std_dev = 0
    try:
        x.std_score(1)
    except ValueError:
        pass  # Normal behavior

def test_correlations():
    "Correlations between variables"

    a = ufloat(1, 0)
    x = ufloat(4, 0.1)
    y = x*2 + a
    # Correlations cancel "naive" additions of uncertainties:
    assert y.std_dev != 0
    normally_zero = y - (x*2 + 1)
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
    '''
    Wrap a function that takes only positional-or-keyword parameters.
    '''

    def f_auto_unc(x, y):
        return 2*x+umath.sin(y)

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
    f_wrapped = uncert_core.wrap(f, [lambda x, y: 2,
                                     lambda x, y: math.cos(y)])

    assert ufloats_close(f_auto_unc(x, y), f_wrapped(x, y))

    # Call with keyword arguments:
    assert ufloats_close(f_auto_unc(y=y, x=x), f_wrapped(y=y, x=x))

    ## Automatic additional derivatives for non-defined derivatives:
    f_wrapped = uncert_core.wrap(f, [lambda x, y: 2])  # No derivative for y
    assert ufloats_close(f_auto_unc(x, y), f_wrapped(x, y))

    # Call with keyword arguments:
    assert ufloats_close(f_auto_unc(y=y, x=x), f_wrapped(y=y, x=x))

def test_wrapped_func_args_no_kwargs():
    '''
    Wrap a function that takes only positional-or-keyword and
    var-positional parameters.
    '''

    def f_auto_unc(x, y, *args):
        return 2*x+umath.sin(y)+3*args[1]

    # Like f_auto_unc, but does not accept numbers with uncertainties:
    def f(x, y, *args):
        assert not any(isinstance(value, uncert_core.UFloat)
                       for value in [x, y] + list(args))
        return f_auto_unc(x, y, *args)

    x = uncert_core.ufloat(1, 0.1)
    y = uncert_core.ufloat(10, 2)
    s = 'string arg'
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
    f_wrapped = uncert_core.wrap(f, [lambda x, y, *args: 2,
                                     lambda x, y, *args: math.cos(y),
                                     None,
                                     lambda x, y, *args: 3])

    assert ufloats_close(f_auto_unc(x, y, *args), f_wrapped(x, y, *args))

    ## Automatic additional derivatives for non-defined derivatives:

    # No derivative for y:
    f_wrapped = uncert_core.wrap(f, [lambda x, y, *args: 2])
    assert ufloats_close(f_auto_unc(x, y, *args), f_wrapped(x, y, *args))

def test_wrapped_func_no_args_kwargs():
    '''
    Wrap a function that takes only positional-or-keyword and
    var-keyword parameters.
    '''

    def f_auto_unc(x, y, **kwargs):
        return 2*x+umath.sin(y)+3*kwargs['z']

    # Like f_auto_unc, but does not accept numbers with uncertainties:
    def f(x, y, **kwargs):
        assert not any(isinstance(value, uncert_core.UFloat)
                       for value in [x, y] + list(kwargs.values()))
        return f_auto_unc(x, y, **kwargs)

    x = uncert_core.ufloat(1, 0.1)
    y = uncert_core.ufloat(10, 2)
    s = 'string arg'
    z = uncert_core.ufloat(100, 3)

    kwargs = {'s': s, 'z': z}  # Arguments not in signature

    ### Automatic numerical derivatives:

    ## Fully automatic numerical derivatives:
    f_wrapped = uncert_core.wrap(f)
    assert ufloats_close(f_auto_unc(x, y, **kwargs),
                         f_wrapped(x, y, **kwargs))

    # Call with keyword arguments:
    assert ufloats_close(f_auto_unc(y=y, x=x, **kwargs),
                         f_wrapped(y=y, x=x, **kwargs))

    ## Automatic additional derivatives for non-defined derivatives,
    ## and explicit None derivative:

    # No derivative for positional-or-keyword parameter y, no
    # derivative for optional-keyword parameter z:
    f_wrapped = uncert_core.wrap(f, [None])
    assert ufloats_close(f_auto_unc(x, y, **kwargs),
                         f_wrapped(x, y, **kwargs))

    # Call with keyword arguments:
    assert ufloats_close(f_auto_unc(y=y, x=x, **kwargs),
                         f_wrapped(y=y, x=x, **kwargs))

    # No derivative for positional-or-keyword parameter y, no
    # derivative for optional-keyword parameter z:
    f_wrapped = uncert_core.wrap(f, [None], {'z': None})
    assert ufloats_close(f_auto_unc(x, y, **kwargs),
                         f_wrapped(x, y, **kwargs))

    # Call with keyword arguments:
    assert ufloats_close(f_auto_unc(y=y, x=x, **kwargs),
                         f_wrapped(y=y, x=x, **kwargs))

    # No derivative for positional-or-keyword parameter y, derivative
    # for optional-keyword parameter z:
    f_wrapped = uncert_core.wrap(f, [None],
                                 {'z': lambda x, y, **kwargs: 3})
    assert ufloats_close(f_auto_unc(x, y, **kwargs),
                         f_wrapped(x, y, **kwargs))

    # Call with keyword arguments:
    assert ufloats_close(f_auto_unc(y=y, x=x, **kwargs),
                         f_wrapped(y=y, x=x, **kwargs))

    ### Explicit derivatives:

    ## Fully defined derivatives:
    f_wrapped = uncert_core.wrap(
        f,
        [lambda x, y, **kwargs: 2, lambda x, y, **kwargs: math.cos(y)],
        {'z:': lambda x, y, **kwargs: 3})

    assert ufloats_close(f_auto_unc(x, y, **kwargs),
                          f_wrapped(x, y, **kwargs))
    # Call with keyword arguments:
    assert ufloats_close(f_auto_unc(y=y, x=x, **kwargs),
                          f_wrapped(y=y, x=x, **kwargs))

    ## Automatic additional derivatives for non-defined derivatives:

    # No derivative for y or z:
    f_wrapped = uncert_core.wrap(f, [lambda x, y, **kwargs: 2])
    assert ufloats_close(f_auto_unc(x, y, **kwargs),
                          f_wrapped(x, y, **kwargs))

    # Call with keyword arguments:
    assert ufloats_close(f_auto_unc(y=y, x=x, **kwargs),
                          f_wrapped(y=y, x=x, **kwargs))

def test_wrapped_func_args_kwargs():
    '''
    Wrap a function that takes positional-or-keyword, var-positional
    and var-keyword parameters.
    '''

    def f_auto_unc(x, y, *args, **kwargs):
        return 2*x+umath.sin(y)+4*args[1]+3*kwargs['z']

    # Like f_auto_unc, but does not accept numbers with uncertainties:
    def f(x, y, *args, **kwargs):
        assert not any(isinstance(value, uncert_core.UFloat)
                       for value in [x, y]+list(args)+list(kwargs.values()))
        return f_auto_unc(x, y, *args, **kwargs)

    x = uncert_core.ufloat(1, 0.1)
    y = uncert_core.ufloat(10, 2)
    t = uncert_core.ufloat(1000, 4)
    s = 'string arg'
    z = uncert_core.ufloat(100, 3)

    args = [s, t, s]
    kwargs = {'u': s, 'z': z}  # Arguments not in signature

    ### Automatic numerical derivatives:

    ## Fully automatic numerical derivatives:
    f_wrapped = uncert_core.wrap(f)

    assert ufloats_close(f_auto_unc(x, y, *args, **kwargs),
                          f_wrapped(x, y, *args, **kwargs), tolerance=1e-5)

    ## Automatic additional derivatives for non-defined derivatives,
    ## and explicit None derivative:

    # No derivative for positional-or-keyword parameter y, no
    # derivative for optional-keyword parameter z:
    f_wrapped = uncert_core.wrap(f, [None, None, None,
                                     lambda x, y, *args, **kwargs: 4])
    assert ufloats_close(f_auto_unc(x, y, *args, **kwargs),
                         f_wrapped(x, y, *args, **kwargs), tolerance=1e-5)

    # No derivative for positional-or-keyword parameter y, no
    # derivative for optional-keyword parameter z:
    f_wrapped = uncert_core.wrap(f, [None], {'z': None})
    assert ufloats_close(f_auto_unc(x, y, *args, **kwargs),
                         f_wrapped(x, y, *args, **kwargs), tolerance=1e-5)

    # No derivative for positional-or-keyword parameter y, derivative
    # for optional-keyword parameter z:
    f_wrapped = uncert_core.wrap(f, [None],
                                 {'z': lambda x, y, *args, **kwargs: 3})
    assert ufloats_close(f_auto_unc(x, y, *args, **kwargs),
                         f_wrapped(x, y, *args, **kwargs), tolerance=1e-5)

    ### Explicit derivatives:

    ## Fully defined derivatives:
    f_wrapped = uncert_core.wrap(
        f,
        [lambda x, y, *args, **kwargs: 2,
         lambda x, y, *args, **kwargs: math.cos(y)],
        {'z:': lambda x, y, *args, **kwargs: 3})

    assert ufloats_close(f_auto_unc(x, y, *args, **kwargs),
                          f_wrapped(x, y, *args, **kwargs), tolerance=1e-5)

    ## Automatic additional derivatives for non-defined derivatives:

    # No derivative for y or z:
    f_wrapped = uncert_core.wrap(f, [lambda x, y, *args, **kwargs: 2])
    assert ufloats_close(f_auto_unc(x, y, *args, **kwargs),
                         f_wrapped(x, y, *args, **kwargs), tolerance=1e-5)


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
        assert not any(isinstance(arg, uncert_core.UFloat)
                       for arg in list_var)
        return f_auto_unc(angle, *list_var)

    f_wrapped = uncert_core.wrap(f)


    my_list = [1, 2, 3]

    ########################################
    # Test of a wrapped function that only calls the original
    # function: it should obtain the exact same result:
    assert f_wrapped(0, *my_list) == f(0, *my_list)
    # 1 == 1 +/- 0, so the type must be checked too:
    assert type(f_wrapped(0, *my_list)) == type(f(0, *my_list))

    ########################################
    # Call with uncertainties:

    angle = uncert_core.ufloat(1, 0.1)
    list_value = uncert_core.ufloat(3, 0.2)

    # The random variables must be the same (full correlation):

    assert ufloats_close(f_wrapped(angle, *[1, angle]),
                          f_auto_unc(angle, *[1, angle]))

    assert ufloats_close(f_wrapped(angle, *[list_value, angle]),
                          f_auto_unc(angle, *[list_value, angle]))

    ########################################
    # Non-numerical arguments, and  explicit and implicit derivatives:
    def f(x, y, z, t, u):
        return x+2*z+3*t+4*u

    f_wrapped = uncert_core.wrap(
        f, [lambda *args: 1, None, lambda *args:2, None])  # No deriv. for u

    assert f_wrapped(10, 'string argument', 1, 0, 0) == 12

    x = uncert_core.ufloat(10, 1)

    assert numbers_close(f_wrapped(x, 'string argument', x, x, x).std_dev,
                          (1+2+3+4)*x.std_dev)

def test_wrap_with_kwargs():
    '''
    Tests wrap() on functions with keyword arguments.

    Includes both wrapping a function that takes optional keyword
    arguments and calling a wrapped function with keyword arguments
    (optional or not).
    '''

    # Version of f() that automatically works with numbers with
    # uncertainties:
    def f_auto_unc(x, y, *args, **kwargs):
        return x + umath.sin(y) + 2*args[0] + 3*kwargs['t']

    # We also add keyword arguments in the function which is wrapped:
    def f(x, y, *args, **kwargs):
        # We make sure that f is not called directly with a number with
        # uncertainty:

        for value in [x, y]+list(args)+list(kwargs.values()):
            assert not isinstance(value, uncert_core.UFloat)

        return f_auto_unc(x, y, *args, **kwargs)

    f_wrapped = uncert_core.wrap(f)


    x = ufloat(1, 0.1)
    y = ufloat(10, 0.11)
    z = ufloat(100, 0.111)
    t = ufloat(0.1, 0.1111)

    assert ufloats_close(f_wrapped(x, y, z, t=t),
                          f_auto_unc(x, y, z, t=t), tolerance=1e-5)

    ########################################

    # We make sure that analytical derivatives are indeed used. We
    # also test the automatic handling of additional *args arguments
    # beyond the number of supplied derivatives.

    f_wrapped2 = uncert_core.wrap(
        f, [None, lambda x, y, *args, **kwargs: math.cos(y)])

    # The derivatives must be perfectly identical:

    # The *args parameter of f() is given as a keyword argument, so as
    # to try to confuse the code:

    assert (f_wrapped2(x, y, z, t=t).derivatives[y]
            == f_auto_unc(x, y, z, t=t).derivatives[y])

    # Derivatives supplied through the keyword-parameter dictionary of
    # derivatives, and also derivatives supplied for the
    # var-positional arguments (*args[0]):

    f_wrapped3 = uncert_core.wrap(
        f,
        [None, None, lambda x, y, *args, **kwargs: 2],
        {'t': lambda x, y, *args, **kwargs: 3})

    # The derivatives should be exactly the same, because they are
    # obtained with the exact same analytic formula:
    assert (f_wrapped3(x, y, z, t=t).derivatives[z]
            == f_auto_unc(x, y, z, t=t).derivatives[z])
    assert (f_wrapped3(x, y, z, t=t).derivatives[t]
            == f_auto_unc(x, y, z, t=t).derivatives[t])

    ########################################
    # Making sure that user-supplied derivatives are indeed called:

    class FunctionCalled(Exception):
        '''
        Raised to signal that a function is indeed called.
        '''
        pass

    def failing_func(x, y, *args, **kwargs):
        raise FunctionCalled

    f_wrapped4 = uncert_core.wrap(
        f,
        [None, failing_func],
        {'t': failing_func})

    try:
        f_wrapped4(x, 3.14, z, t=t)
    except FunctionCalled:
        pass
    else:
        raise Exception('User-supplied derivative should be called')

    try:
        f_wrapped4(x, y, z, t=3.14)
    except FunctionCalled:
        pass
    else:
        raise Exception('User-supplied derivative should be called')

    try:
        f_wrapped4(x, 3.14, z, t=3.14)
    except FunctionCalled:
        raise Exception('User-supplied derivative should *not* be called')

###############################################################################

def test_access_to_std_dev():
    "Uniform access to the standard deviation"

    x = ufloat(1, 0.1)
    y = 2*x

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
    y = -2*x+10
    z = -3*x
    covs = uncert_core.covariance_matrix([x, y, z])
    # Diagonal elements are simple:
    assert numbers_close(covs[0][0], 0.01)
    assert numbers_close(covs[1][1], 0.04)
    assert numbers_close(covs[2][2], 0.09)
    # Non-diagonal elements:
    assert numbers_close(covs[0][1], -0.02)


###############################################################################
def test_power_all_cases():
    '''
    Checks all cases for the value and derivatives of x**p.
    '''

    power_all_cases(pow)

def power_all_cases(op):
    '''
    Checks all cases for the value and derivatives of power-like
    operator op (op is typically the built-in pow(), or math.pow()).

    Checks only the details of special results like 0, 1 or NaN).

    Different cases for the value of x**p and its derivatives are
    tested by dividing the (x, p) plane with:

    - x < 0, x = 0, x > 0
    - p integer or not, p < 0, p = 0, p > 0

    (not all combinations are distinct: for instance x > 0 gives
    identical formulas for all p).
    '''

    zero = ufloat(0, 0.1)
    zero2 = ufloat(0, 0.1)
    one = ufloat(1, 0.1)
    positive = ufloat(0.3, 0.01)
    positive2 = ufloat(0.3, 0.01)
    negative = ufloat(-0.3, 0.01)
    integer = ufloat(-3, 0)
    non_int_larger_than_one = ufloat(3.1, 0.01)
    positive_smaller_than_one = ufloat(0.3, 0.01)

    ## negative**integer

    result = op(negative, integer)
    assert not isnan(result.derivatives[negative])
    assert isnan(result.derivatives[integer])

    # Limit cases:
    result = op(negative, one)
    assert result.derivatives[negative] == 1
    assert isnan(result.derivatives[one])

    result = op(negative, zero)
    assert result.derivatives[negative] == 0
    assert isnan(result.derivatives[zero])

    ## negative**non-integer

    ## zero**...

    result = op(zero, non_int_larger_than_one)
    assert isnan(result.derivatives[zero])
    assert result.derivatives[non_int_larger_than_one] == 0

    # Special cases:
    result = op(zero, one)
    assert result.derivatives[zero] == 1
    assert result.derivatives[one] == 0

    result = op(zero, 2*one)
    assert result.derivatives[zero] == 0
    assert result.derivatives[one] == 0

    result = op(zero, positive_smaller_than_one)
    assert isnan(result.derivatives[zero])
    assert result.derivatives[positive_smaller_than_one] == 0

    result = op(zero, zero2)
    assert result.derivatives[zero] == 0
    assert isnan(result.derivatives[zero2])

    ## positive**...: this is a quite regular case where the value and
    ## the derivatives are all defined.

    result = op(positive, positive2)
    assert not isnan(result.derivatives[positive])
    assert not isnan(result.derivatives[positive2])

    result = op(positive, zero)
    assert result.derivatives[positive] == 0
    assert not isnan(result.derivatives[zero])

    result = op(positive, negative)
    assert not isnan(result.derivatives[positive])
    assert not isnan(result.derivatives[negative])


###############################################################################

def test_power_special_cases():
    '''
    Checks special cases of x**p.
    '''
    power_special_cases(pow)

    # We want the same behavior for numbers with uncertainties and for
    # math.pow() at their nominal values:

    positive = ufloat(0.3, 0.01)
    negative = ufloat(-0.3, 0.01)

    # http://stackoverflow.com/questions/10282674/difference-between-the-built-in-pow-and-math-pow-for-floats-in-python

    try:
        pow(ufloat(0, 0), negative)
    except ZeroDivisionError:
        pass
    else:
        raise Exception("A proper exception should have been raised")

    try:
        pow(ufloat(0, 0.1), negative)
    except ZeroDivisionError:
        pass
    else:
        raise Exception('A proper exception should have been raised')

    try:
        result = pow(negative, positive)
    except ValueError:
        # The reason why it should also fail in Python 3 is that the
        # result of Python 3 is a complex number, which uncertainties
        # does not handle (no uncertainties on complex numbers). In
        # Python 2, this should always fail, since Python 2 does not
        # know how to calculate it.
        pass
    else:
        raise Exception('A proper exception should have been raised')

def power_special_cases(op):
    '''
    Checks special cases of the uncertainty power operator op (where
    op is typically the built-in pow or uncertainties.umath.pow).

    The values x = 0, x = 1 and x = NaN are special, as are null,
    integral and NaN values of p.
    '''

    zero = ufloat(0, 0)
    one = ufloat(1, 0)
    p = ufloat(0.3, 0.01)

    assert op(0, p) == 0
    assert op(zero, p) == 0

    # The outcome of 1**nan and nan**0 was undefined before Python
    # 2.6 (http://docs.python.org/library/math.html#math.pow):
    assert op(float('nan'), zero) == 1.0
    assert op(one, float('nan')) == 1.0

    # …**0 == 1.0:
    assert op(p, 0) == 1.0
    assert op(zero, 0) == 1.0
    assert op((-p), 0) == 1.0
    # …**zero:
    assert op((-10.3), zero) == 1.0
    assert op(0, zero) == 1.0
    assert op(0.3, zero) == 1.0
    assert op((-p), zero) == 1.0
    assert op(zero, zero) == 1.0
    assert op(p, zero) == 1.0

    # one**… == 1.0
    assert op(one, -3) == 1.0
    assert op(one, -3.1) == 1.0
    assert op(one, 0) == 1.0
    assert op(one, 3) == 1.0
    assert op(one, 3.1) == 1.0

    # … with two numbers with uncertainties:
    assert op(one, (-p)) == 1.0
    assert op(one, zero) == 1.0
    assert op(one, p) == 1.0
    # 1**… == 1.0:
    assert op(1., (-p)) == 1.0
    assert op(1., zero) == 1.0
    assert op(1., p) == 1.0


def test_power_wrt_ref():
    '''
    Checks special cases of the built-in pow() power operator.
    '''
    power_wrt_ref(pow, pow)

def power_wrt_ref(op, ref_op):
    '''
    Checks special cases of the uncertainty power operator op (where
    op is typically the built-in pow or uncertainties.umath.pow), by
    comparing its results to the reference power operator ref_op
    (which is typically the built-in pow or math.pow).
    '''

    # Negative numbers with uncertainty can be exponentiated to an
    # integral power:
    assert op(ufloat(-1.1, 0.1), -9).nominal_value == ref_op(-1.1, -9)

    # Case of numbers with no uncertainty: should give the same result
    # as numbers with uncertainties:
    assert op(ufloat(-1, 0), 9) == ref_op(-1, 9)
    assert op(ufloat(-1.1, 0), 9) == ref_op(-1.1, 9)


###############################################################################

def test_PDG_precision():
    '''
    Test of the calculation of the number of significant digits for
    the uncertainty.
    '''

    # The 3 cases of the rounding rules are covered in each case:
    tests = {
        # Very big floats:
        1.7976931348623157e308: (2, 1.7976931348623157e308),
        0.5e308: (1, 0.5e308),
        0.9976931348623157e+308: (2, 1e308),
        # Very small floats:
        1.3e-323: (2, 1.3e-323),
        5e-324: (1, 5e-324),
        9.99e-324: (2, 1e-323)
        }

    for (std_dev, result) in tests.items():
        assert uncert_core.PDG_precision(std_dev) == result

def test_repr():
    '''Test the representation of numbers with uncertainty.'''

    # The uncertainty is a power of 2, so that it can be exactly
    # represented:
    x = ufloat(3.14159265358979, 0.25)
    assert repr(x) == '3.14159265358979+/-0.25'

    x = ufloat(3.14159265358979, 0)
    assert repr(x) == '3.14159265358979+/-0'

    # Tagging:
    x = ufloat(3, 1, "length")
    assert repr(x) == '< length = 3.0+/-1.0 >'

def test_format():
    '''Test the formatting of numbers with uncertainty.'''

    # The way NaN is formatted with F, E and G depends on the version
    # of Python (NAN for Python 2.5+ at least):
    NaN_EFG = '%F' % float('nan')

    # !! The way NaN is formatted with F, E and G might depend on the
    # version of Python, if it is like NaN (could be tested with
    # Python 2.3 or 2.4 vs Python 2.7):
    Inf_EFG = '%F' % float('inf')

    # Tests of each point of the docstring of
    # AffineScalarFunc.__format__() in turn, mostly in the same order.

    # The LaTeX tests do not use the customization of
    # uncert_core.GROUP_SYMBOLS and uncert_core.EXP_PRINT: this
    # way, problems in the customization themselves are caught.

    tests = {  # (Nominal value, uncertainty): {format: result,...}

        # Usual float formatting, and individual widths, etc.:
        (3.1415, 0.0001): {
            '*^+7.2f': '*+3.14*+/-*0.00**',
            '+07.2f': '+003.14+/-0000.00',  # 0 fill
            '>10f': '  3.141500+/-  0.000100',  # Width and align
            '11.3e': '  3.142e+00+/-  0.000e+00',  # Duplicated exponent
            '0.4e': '3.1415e+00+/-0.0000e+00'  # Forced double exponent
        },

        # Full generalization of float formatting:
        (3.1415, 0.0001): {
            '+09.2uf': '+03.14150+/-000.00010',
            # Alignment is not available with the % formatting
            # operator of Python < 2.6:
            '*^+9.2uf': '+3.14150*+/-*0.00010*',
            '>9f': '  3.14150+/-  0.00010'  # Width and align
        },

        # Number of digits of the uncertainty fixed:
        (123.456789, 0.00123): {
            '.1uf': '123.457+/-0.001',
            '.2uf': '123.4568+/-0.0012',
            '.3uf': '123.45679+/-0.00123',
            '.2ue': '(1.234568+/-0.000012)e+02'
        },
        # Sign handling:
        (-123.456789, 0.00123): {
            '.1uf': '-123.457+/-0.001',
            '.2uf': '-123.4568+/-0.0012',
            '.3uf': '-123.45679+/-0.00123',
            '.2ue': '(-1.234568+/-0.000012)e+02'
        },
        # Uncertainty larger than the nominal value:
        (12.3, 456.78): {
            '': '12+/-457',
            '.1uf': '12+/-457',
            '.4uf': '12.3+/-456.8'
        },
        # ... Same thing, but with an exponent:
        (12.3, 456.78): {
            '.1ue': '(0+/-5)e+02',
            '.4ue': '(0.123+/-4.568)e+02',
            '.4ueS': '0.123(4.568)e+02'
        },

        (23456.789123, 1234.56789123): {
            '.6gS': '23456.8(1234.6)'
        },

        # Test of the various float formats: the nominal value should
        # have a similar representation as if it were directly
        # represented as a float:
        (1234567.89, 0.1): {
            '.0e': '(1+/-0)e+06',
            'e': '(1.23456789+/-0.00000010)e+06',
            'E': '(1.23456789+/-0.00000010)E+06',
            'f': '1234567.89+/-0.10',
            'F': '1234567.89+/-0.10',
            'g': '1234567.89+/-0.10',
            'G': '1234567.89+/-0.10',
            '%': '(123456789+/-10)%'
        },
        (1234567.89, 4.3): {
            'g': '1234568+/-4'
        },
        (1234567.89, 43): {  # Case where g triggers the exponent notation
            'g': '(1.23457+/-0.00004)e+06',
            'G': '(1.23457+/-0.00004)E+06'
        },


        (3.1415, 0.0001): {
            '+09.2uf': '+03.14150+/-000.00010'
            },

        (1234.56789, 0.1): {
            '.0f': '(1234+/-0.)',  # Approximate error indicated with "."
            'e': '(1.23456+/-0.00010)e+03',
            'E': '(1.23456+/-0.00010)E+03',
            'f': '1234.57+/-0.10',
            'F': '1234.57+/-0.10',
            'f': '1234.57+/-0.10',
            'F': '1234.57+/-0.10',
            '%': '123457+/-10%'
        },

        # Percent notation:
        (0.42, 0.0055): {
            # Because '%' does 0.0055*100, the value
            # 0.5499999999999999 is obtained, which rounds to 0.5. The
            # original rounded value is 0.006. The same behavior is
            # found in Python 2.7: '{:.1%}'.format(0.0055) is '0.5%'.
            '.1u%': '(42.0+/-0.5)%',
            '.1u%S': '42.0(5)%',
            '%P': '(42.0±0.5)%'
        },

        # Particle Data Group automatic convention, including limit cases:
        (1.2345678, 0.354): {'': '1.23+/-0.35'},
        (1.2345678, 0.3549): {'': '1.23+/-0.35'},
        (1.2345678, 0.355): {'': '1.2+/-0.4'},
        (1.5678, 0.355): {'': '1.6+/-0.4'},
        (1.2345678, 0.09499): {'': '1.23+/-0.09'},
        (1.2345678, 0.095): {'': '1.23+/-0.10'},

        # Automatic extension of the uncertainty up to the decimal
        # point:
        (1000, 123): {
            '.1uf': '1000+/-123',
            # The nominal value has 1 <= mantissa < 10. The precision
            # is the number of significant digits of the uncertainty:
            '.1ue': '(1.0+/-0.1)e+03'
        },

        # Spectroscopic notation:
        (-1.23, 3.4): {
            'S': '-1.2(3.4)',
            '.2ufS': '-1.2(3.4)',
            '.3ufS': '-1.23(3.40)',
        },
        (-123.456, 0.123): {
            'S': '-123.46(12)',
            '.1ufS': '-123.5(1)',
            '.2ufS': '-123.46(12)',
            '.3ufS': '-123.456(123)',
        },
        (-123.456, 0.567): {
            'S': '-123.5(6)',
            '.1ufS': '-123.5(6)',
            '.2ufS': '-123.46(57)',
            '.3ufS': '-123.456(567)',
        },
        (-123.456, 0.004): {
            # The decimal point shows that the uncertainty is not
            # exact:
            '.2fS': '-123.46(0.00)'
        },

        # LaTeX notation:
        #
        (1234.56789, 0.1): {
            'eL': r'\left(1.23457 \pm 0.00010\right) \times 10^{3}',
            'EL': r'\left(1.23457 \pm 0.00010\right) \times 10^{3}',
            'fL': '1234.57 \pm 0.10',
            'FL': '1234.57 \pm 0.10',
            'fL': '1234.57 \pm 0.10',
            'FL': '1234.57 \pm 0.10',
            '%L': r'\left(123457 \pm 10\right) \%'
        },
        #
        # ... combined with the spectroscopic notation:
        (-1.23, 3.4): {
            'SL': '-1.2(3.4)',
            'LS': '-1.2(3.4)',
            '.2ufSL': '-1.2(3.4)',
            '.2ufLS': '-1.2(3.4)'
        },

        # Special cases for the uncertainty (0, nan) and format
        # strings (extension S, L, U,..., global width, etc.).
        #
        # Python 3.2 and 3.3 give 1.4e-12*1e+12 = 1.4000000000000001
        # instead of 1.4 for Python 3.1. The problem does not appear
        # with 1.2, so 1.2 is used.
        (-1.2e-12, 0): {
            '12.2gPL': r'  -1.2×10⁻¹²±           0',
            # Pure "width" formats are not accepted by the % operator,
            # and only %-compatible formats are accepted, for Python <
            # 2.6:
            '13S': '  -1.2(0)e-12',
            '10P': '-1.2×10⁻¹²±         0',
            'L': r'\left(-1.2 \pm 0\right) \times 10^{-12}',
            # No factored exponent, LaTeX
            '1L': r'-1.2 \times 10^{-12} \pm 0',
            'SL': r'-1.2(0) \times 10^{-12}',
            'SP': r'-1.2(0)×10⁻¹²'
        },

        # Python 3.2 and 3.3 give 1.4e-12*1e+12 = 1.4000000000000001
        # instead of 1.4 for Python 3.1. The problem does not appear
        # with 1.2, so 1.2 is used.
        (-1.2e-12, float('nan')): {
            '.2uG': '(-1.2+/-%s)E-12' % NaN_EFG,  # u ignored, format used
            '15GS': '  -1.2(%s)E-12' % NaN_EFG,
            'SL': r'-1.2(\mathrm{nan}) \times 10^{-12}',  # LaTeX NaN
            # Pretty-print priority, but not for NaN:
            'PSL': '-1.2(\mathrm{nan})×10⁻¹²',
            'L': r'\left(-1.2 \pm \mathrm{nan}\right) \times 10^{-12}',
            # Uppercase NaN and LaTeX:
            '.1EL': (r'\left(-1.2 \pm \mathrm{%s}\right) \times 10^{-12}'
                     % NaN_EFG),
            '10': '  -1.2e-12+/-       nan',
            '15S': '  -1.2(nan)e-12'
        },

        (3.14e-10, 0.01e-10): {
            # Character (Unicode) strings:
            'P': '(3.140±0.010)×10⁻¹⁰',  # PDG rules: 2 digits
            'PL': r'(3.140±0.010)×10⁻¹⁰',  # Pretty-print has higher priority
            # Truncated non-zero uncertainty:
            '.1e': '(3.1+/-0.0)e-10',
            '.1eS': '3.1(0.0)e-10'
        },

        # Some special cases:
        (1, float('nan')): {
            'g': '1+/-nan',
            'G': '1+/-%s' % NaN_EFG,
            '%': '(100.000000+/-nan)%',  # The % format type is like f
            # Should be the same as '+05', for floats, but is not, in
            # Python 2.7:
            '+05g': '+0001+/-00nan',
            # 5 is the *minimal* width, 6 is the default number of
            # digits after the decimal point:
            '+05%': '(+100.000000+/-00nan)%',
            # There is a difference between '{}'.format(1.) and
            # '{:g}'.format(1.), which is not fully obvious in the
            # documentation, which indicates that a None format type
            # is like g. The reason is that the empty format string is
            # actually interpreted as str(), and that str() does not
            # have to behave like g ('{}'.format(1.234567890123456789)
            # and '{:g}'.format(1.234567890123456789) are different).
            '': '1.0+/-nan',
            # This is ugly, but consistent with
            # '{:+05}'.format(float('nan')) and format(1.) (which
            # differs from format(1)!):
            '+05': '+01.0+/-00nan'
        },

        (9.9, 0.1): {
            '.1ue': '(9.9+/-0.1)e+00',
            '.0fS': '10(0.)'
        },
        (9.99, 0.1): {
             # The precision has an effect on the exponent, like for
             # floats:
            '.2ue': '(9.99+/-0.10)e+00',  # Same exponent as for 9.99 alone
            '.1ue': '(1.00+/-0.01)e+01'  # Same exponent as for 9.99 alone
        },
        # 0 uncertainty: nominal value displayed like a float:
        (1.2345, 0): {
            '.2ue': '(1.23+/-0)e+00',
            '1.2ue': '1.23e+00+/-0',  # No factored exponent
            '.2uf': '1.23+/-0',
            '.2ufS': '1.23(0)',
            '.2fS': '1.23(0)',
            'g': '1.2345+/-0',
            '': '1.2345+/-0'
        },

        # Alignment and filling characters (supported in Python 2.6+):
        (3.1415e10, 0): {
            '<15': '31415000000.0  +/-0              ',
            '<20S': '31415000000.0(0)    ',
            # Trying to trip the format parsing with a fill character
            # which is an alignment character:
            '=>15': '==31415000000.0+/-==============0'
        },

        (1234.56789, 0): {
            '1.2ue': '1.23e+03+/-0',  # u ignored
            '1.2e': '1.23e+03+/-0',
            # Default precision = 6
            'eL': r'\left(1.234568 \pm 0\right) \times 10^{3}',
            'EL': r'\left(1.234568 \pm 0\right) \times 10^{3}',
            'fL': '1234.567890 \pm 0',
            'FL': '1234.567890 \pm 0',
            '%L': r'\left(123456.789000 \pm 0\right) \%'
        },

        (1e5, 0): {
            'g': '100000+/-0'
        },
        (1e6, 0): {
            # A default precision of 6 is used because the uncertainty
            # cannot be used for defining a default precision (it does
            # not have a magnitude):
            'g': '(1+/-0)e+06'
        },
        (1e6+10, 0): {
            # A default precision of 6 is used because the uncertainty
            # cannot be used for defining a default precision (it does
            # not have a magnitude):
            'g': '(1.00001+/-0)e+06'
        },
        # Rounding of the uncertainty that "changes" the number of
        # significant digits:
        (1, 0.994): {
            '.3uf': '1.000+/-0.994',
            '.2uf': '1.00+/-0.99',
            '.1uf': '1+/-1'  # Discontinuity in the number of digits
        },
        (12.3, 2.3): {
            '.2ufS': '12.3(2.3)'  # Decimal point on the uncertainty
        },
        (12.3, 2.3): {
            '.1ufS': '12(2)'  # No decimal point on the uncertainty
        },
        (0, 0): {  # Make defining the first significant digit problematic
            '.1f': '0.0+/-0',  # Simple float formatting
            'g': '0+/-0'
        },
        (1.2e-34, 5e-67): {
            '.6g': '(1.20000+/-0.00000)e-34',
            '13.6g': '  1.20000e-34+/-  0.00000e-34',
            '13.6G': '  1.20000E-34+/-  0.00000E-34',
            '.6GL': r'\left(1.20000 \pm 0.00000\right) \times 10^{-34}'
        },

        (float('nan'), 100): {  # NaN *nominal value*
            '': 'nan+/-100.0',  # Like '{}'.format(100.)
            'g': 'nan+/-100',  # Like '{:g}'.format(100.)
            '.1e': '(nan+/-1.0)e+02',  # Similar to 1±nan
            '.1E': '(%s+/-1.0)E+02' % NaN_EFG,
            '.1ue': '(nan+/-1)e+02',
            '10.1e': '       nan+/-   1.0e+02'
        },
        (float('nan'), 1e8): {  # NaN *nominal value*
            '': 'nan+/-100000000.0',  # Like '{}'.format(1e8)
            'g': '(nan+/-1)e+08',  # Like '{:g}'.format(1e8)
            '.1e': '(nan+/-1.0)e+08',
            '.1E': '(%s+/-1.0)E+08' % NaN_EFG,
            '.1ue': '(nan+/-1)e+08',
            '10.1e': '       nan+/-   1.0e+08'  # 'nane+08' would be strange
        },
        (float('nan'), 123456789): {  # NaN *nominal value*
            '': 'nan+/-123456789.0',  # Similar to '{}'.format(123456789.)
            'g': '(nan+/-1.23457)e+08',  # Similar to '{:g}'.format(123456789.)
            '.1e': '(nan+/-1.2)e+08',
            '.1E': '(%s+/-1.2)E+08' % NaN_EFG,
            '.1ue': '(nan+/-1)e+08',
            '.1ueL': r'\left(\mathrm{nan} \pm 1\right) \times 10^{8}',
            '10.1e': '       nan+/-   1.2e+08',
            '10.1eL': r'\mathrm{nan} \pm 1.2 \times 10^{8}'
        },
        (float('nan'), float('nan')): {  # *Double* NaN
            '': 'nan+/-nan',
            '.1e': 'nan+/-nan',
            '.1E': '%s+/-%s' % (NaN_EFG, NaN_EFG),
            '.1ue': 'nan+/-nan',
            'EL': r'\mathrm{%s} \pm \mathrm{%s}' % (NaN_EFG, NaN_EFG)
        },

        (float('inf'), 100): {  # Inf *nominal value*
            '': 'inf+/-100.0',  # Like '{}'.format(100.)
            'g': 'inf+/-100',  # Like '{:g}'.format(100.)
            '.1e': '(inf+/-1.0)e+02',  # Similar to 1±inf
            '.1E': '(%s+/-1.0)E+02' % Inf_EFG,
            '.1ue': '(inf+/-1)e+02',
            '10.1e': '       inf+/-   1.0e+02'
        },
        (float('inf'), 1e8): {  # Inf *nominal value*
            '': 'inf+/-100000000.0',  # Like '{}'.format(1e8)
            'g': '(inf+/-1)e+08',  # Like '{:g}'.format(1e8)
            '.1e': '(inf+/-1.0)e+08',
            '.1E': '(%s+/-1.0)E+08' % Inf_EFG,
            '.1ue': '(inf+/-1)e+08',
            '10.1e': '       inf+/-   1.0e+08'  # 'infe+08' would be strange
        },
        (float('inf'), 123456789): {  # Inf *nominal value*
            '': 'inf+/-123456789.0',  # Similar to '{}'.format(123456789.)
            'g': '(inf+/-1.23457)e+08',  # Similar to '{:g}'.format(123456789.)
            '.1e': '(inf+/-1.2)e+08',
            '.1E': '(%s+/-1.2)E+08' % Inf_EFG,
            '.1ue': '(inf+/-1)e+08',
            '.1ueL': r'\left(\infty \pm 1\right) \times 10^{8}',
            '10.1e': '       inf+/-   1.2e+08',
            '10.1eL': r'    \infty \pm 1.2 \times 10^{8}'
        },
        (float('inf'), float('inf')): {  # *Double* Inf
            '': 'inf+/-inf',
            '.1e': 'inf+/-inf',
            '.1E': '%s+/-%s' % (Inf_EFG, Inf_EFG),
            '.1ue': 'inf+/-inf',
            'EL': r'\infty \pm \infty'
        },

        # Like the tests for +infinity, but for -infinity:
        (float('-inf'), 100): {  # Inf *nominal value*
            '': '-inf+/-100.0',  # Like '{}'.format(100.)
            'g': '-inf+/-100',  # Like '{:g}'.format(100.)
            '.1e': '(-inf+/-1.0)e+02',  # Similar to 1±inf
            '.1E': '(-%s+/-1.0)E+02' % Inf_EFG,
            '.1ue': '(-inf+/-1)e+02',
            '10.1e': '      -inf+/-   1.0e+02'
        },
        (float('-inf'), 1e8): {  # Inf *nominal value*
            '': '-inf+/-100000000.0',  # Like '{}'.format(1e8)
            'g': '(-inf+/-1)e+08',  # Like '{:g}'.format(1e8)
            '.1e': '(-inf+/-1.0)e+08',
            '.1E': '(-%s+/-1.0)E+08' % Inf_EFG,
            '.1ue': '(-inf+/-1)e+08',
            '10.1e': '      -inf+/-   1.0e+08'  # 'infe+08' would be strange
        },
        (float('-inf'), 123456789): {  # Inf *nominal value*
            '': '-inf+/-123456789.0',  # Similar to '{}'.format(123456789.)
            'g': '(-inf+/-1.23457)e+08',  # Similar to '{:g}'.format(123456789.)
            '.1e': '(-inf+/-1.2)e+08',
            '.1E': '(-%s+/-1.2)E+08' % Inf_EFG,
            '.1ue': '(-inf+/-1)e+08',
            '.1ueL': r'\left(-\infty \pm 1\right) \times 10^{8}',
            '10.1e': '      -inf+/-   1.2e+08',
            '10.1eL': r'   -\infty \pm 1.2 \times 10^{8}'
        },
        (float('-inf'), float('inf')): {  # *Double* Inf
            '': '-inf+/-inf',
            '.1e': '-inf+/-inf',
            '.1E': '-%s+/-%s' % (Inf_EFG, Inf_EFG),
            '.1ue': '-inf+/-inf',
            'EL': r'-\infty \pm \infty'
        },

        # The Particle Data Group convention trumps the "at least one
        # digit past the decimal point" for Python floats, but only
        # with a non-zero uncertainty:
        (724.2, 26.4): {
            '': '724+/-26'
        },
        (724, 0): {
            '': '724.0+/-0'
        },

        # More NaN and infinity, in particular with LaTeX and various
        # options:
        (float('-inf'), float('inf')): {
            'S': '-inf(inf)',
            'LS': '-\infty(\infty)',
            'L': '-\infty \pm \infty',
            'LP': '-\infty±\infty',
            # The following is consistent with Python's own
            # formatting, which depends on the version of Python:
            # formatting float("-inf") with format(..., "020") gives
            # '-0000000000000000inf' with Python 2.7, but
            # '-00000000000000.0inf' with Python 2.6. However, Python
            # 2.6 gives the better, Python 2.7 form when format()ting
            # with "020g" instead, so this formatting would be better,
            # in principle, and similarly for "%020g" % ... Thus,
            # Python's format() breaks the official rule according to
            # which no format type is equivalent to "g", for
            # floats. If the better behavior was needed, internal
            # formatting could in principle force the "g" formatting
            # type when none is given; however, Python does not
            # actually fully treat the none format type in the same
            # was as the "g" format, so this solution cannot be used,
            # as it would break other formatting behaviors in this
            # code. It is thus best to mimic the native behavior of
            # none type formatting (even if it does not look so good
            # in Python 2.6).
            '020S': format(float("-inf"), '015')+'(inf)'
        },
        (-float('nan'), float('inf')): {
            'S': 'nan(inf)',
            'LS': '\mathrm{nan}(\infty)',
            'L': '\mathrm{nan} \pm \infty',
            'LP': '\mathrm{nan}±\infty'
        },

        # Leading zeroes in the shorthand notation:
        (-2, 3): {
            "020S": "-000000000002.0(3.0)"
        }

    }

    # ',' format option: introduced in Python 2.7
    if sys.version_info >= (2, 7):

        tests.update({
            (1234.56789, 0.012): {
                ',.1uf': '1,234.57+/-0.01'
                },

            (123456.789123, 1234.5678): {
                ',f': '123,457+/-1,235',  # Particle Data Group convention
                ',.4f': '123,456.7891+/-1,234.5678'
                }
        })

    # True if we can detect that the Jython interpreter is running this code:
    try:
        jython_detected = sys.subversion[0] == 'Jython'
    except AttributeError:
        jython_detected = False

    for (values, representations) in tests.items():

        value = ufloat(*values)

        for (format_spec, result) in representations.items():

            # print "FORMATTING {} WITH '{}'".format(repr(value), format_spec)

            # Jython 2.5.2 does not always represent NaN as nan or NAN
            # in the CPython way: for example, '%.2g' % float('nan')
            # is '\ufffd'. The test is skipped, in this case:
            if jython_detected and (
                    isnan(value.std_dev) or isnan(value.nominal_value)):
                continue

            # Call that works with Python < 2.6 too:
            representation = value.format(format_spec)

            assert representation == result, (
                # The representation is used, for terminal that do not
                # support some characters like ±, and superscripts:
                'Incorrect representation %r for format %r of %r:'
                ' %r expected.'
                % (representation, format_spec, value, result))

            # An empty format string is like calling str()
            # (http://docs.python.org/2/library/string.html#formatspec):
            if not format_spec:
                assert representation == str(value), (
                    'Empty format should give the same thing as str():'
                    ' %s obtained instead of %s'
                    % (representation, str(value)))

            # Parsing back into a number with uncertainty (unless the
            # LaTeX or comma notation is used):
            if (not set(format_spec).intersection('L,*%')  # * = fill with *
                # "0nan"
                and '0nan' not in representation.lower()
                # "0inf"
                and '0inf' not in representation.lower()
                # Specific case:
                and '=====' not in representation):

                value_back = ufloat_fromstr(representation)

                # The original number and the new one should be consistent
                # with each other:
                try:

                    # The nominal value can be rounded to 0 when the
                    # uncertainty is larger (because p digits on the
                    # uncertainty can still show 0.00... for the
                    # nominal value). The relative error is infinite,
                    # so this should not cause an error:
                    if value_back.nominal_value:
                        assert numbers_close(value.nominal_value,
                                              value_back.nominal_value, 2.4e-1)

                    # If the uncertainty is zero, then the relative
                    # change can be large:
                    assert numbers_close(value.std_dev,
                                         value_back.std_dev, 3e-1)

                except AssertionError:
                    # !! The following string formatting requires
                    # str() to work (to not raise an exception) on the
                    # values (which have a non-standard class):
                    raise AssertionError(
                        'Original value %s and value %s parsed from %r'
                        ' (obtained through format specification %r)'
                        ' are not close enough'
                        % (value, value_back, representation, format_spec))

def test_unicode_format():
    '''Test of the unicode formatting of numbers with uncertainties'''

    x = ufloat(3.14159265358979, 0.25)

    assert isinstance('Résultat = %s' % x.format(''), str)
    assert isinstance('Résultat = %s' % x.format('P'), str)

###############################################################################

# The tests below require NumPy, which is an optional package:
try:
    import numpy
except ImportError:
    pass
else:

    def arrays_close(m1, m2, precision=1e-4):
        """
        Returns True iff m1 and m2 are almost equal, where elements
        can be either floats or AffineScalarFunc objects.

        Two independent AffineScalarFunc objects are deemed equal if
        both their nominal value and uncertainty are equal (up to the
        given precision).

        m1, m2 -- NumPy arrays.

        precision -- precision passed through to
        uncertainties.test_uncertainties.numbers_close().
        """

        # ! numpy.allclose() is similar to this function, but does not
        # work on arrays that contain numbers with uncertainties, because
        # of the isinf() function.

        for (elmt1, elmt2) in zip(m1.flat, m2.flat):

            # For a simpler comparison, both elements are
            # converted to AffineScalarFunc objects:
            elmt1 = uncert_core.to_affine_scalar(elmt1)
            elmt2 = uncert_core.to_affine_scalar(elmt2)

            if not numbers_close(elmt1.nominal_value,
                                 elmt2.nominal_value, precision):
                return False

            if not numbers_close(elmt1.std_dev,
                                 elmt2.std_dev, precision):
                return False
        
        return True


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
        u2, = uncert_core.correlated_values([1], cov)
        expr = 2*u2  # Calculations with u2 should be possible, like with u

        ####################

        # Covariances between output and input variables:

        x = ufloat(1, 0.1)
        y = ufloat(2, 0.3)
        z = -3*x+y

        covs = uncert_core.covariance_matrix([x, y, z])

        # Test of the diagonal covariance elements:
        assert arrays_close(
            numpy.array([v.std_dev**2 for v in (x, y, z)]),
            numpy.array(covs).diagonal())

        # "Inversion" of the covariance matrix: creation of new
        # variables:
        (x_new, y_new, z_new) = uncert_core.correlated_values(
            [x.nominal_value, y.nominal_value, z.nominal_value],
            covs,
            tags = ['x', 'y', 'z'])

        # Even the uncertainties should be correctly reconstructed:
        assert arrays_close(numpy.array((x, y, z)),
                              numpy.array((x_new, y_new, z_new)))

        # ... and the covariances too:
        assert arrays_close(
            numpy.array(covs),
            numpy.array(uncert_core.covariance_matrix([x_new, y_new, z_new])))

        assert arrays_close(
            numpy.array([z_new]), numpy.array([-3*x_new+y_new]))

        ####################

        # ... as well as functional relations:

        u = ufloat(1, 0.05)
        v = ufloat(10, 0.1)
        sum_value = u+2*v

        # Covariance matrices:
        cov_matrix = uncert_core.covariance_matrix([u, v, sum_value])

        # Correlated variables can be constructed from a covariance
        # matrix, if NumPy is available:
        (u2, v2, sum2) = uncert_core.correlated_values(
            [x.nominal_value for x in [u, v, sum_value]],
            cov_matrix)

        # arrays_close() is used instead of numbers_close() because
        # it compares uncertainties too:
        assert arrays_close(numpy.array([u]), numpy.array([u2]))
        assert arrays_close(numpy.array([v]), numpy.array([v2]))
        assert arrays_close(numpy.array([sum_value]), numpy.array([sum2]))
        assert arrays_close(numpy.array([0]),
                            numpy.array([sum2-(u2+2*v2)]))


        # Spot checks of the correlation matrix:
        corr_matrix = uncert_core.correlation_matrix([u, v, sum_value])
        assert numbers_close(corr_matrix[0,0], 1)
        assert numbers_close(corr_matrix[1,2], 2*v.std_dev/sum_value.std_dev)

        ####################

        # Test of numerical robustness despite wildly different
        # orders of magnitude (see
        # https://github.com/lebigot/uncertainties/issues/95):
        cov = numpy.diag([1e-70, 1e-70, 1e10])
        cov[0, 1] = cov[1, 0] = 0.9e-70
        cov[[0, 1], 2] = -3e-34
        cov[2, [0, 1]] = -3e-34
        variables = uncert_core.correlated_values([0]*3, cov)

        # Since the numbers are very small, we need to compare them
        # in a stricter way, that handles the case of a 0 variance
        # in `variables`:
        assert numbers_close(
                1e66*cov[0,0], 1e66*variables[0].s**2, tolerance=1e-5)
        assert numbers_close(
                1e66*cov[1,1], 1e66*variables[1].s**2, tolerance=1e-5)

        ####################

        # 0 variances are a bit special, since the correlation matrix
        # cannot be calculated naively, so we test that there is no
        # specific problem in this case:

        cov = numpy.diag([0, 0, 10])
        nom_values = [1, 2, 3]
        variables = uncert_core.correlated_values(nom_values, cov)

        for (variable, nom_value, variance) in zip(
            variables, nom_values, cov.diagonal()):
            
            assert numbers_close(variable.n, nom_value)
            assert numbers_close(variable.s**2, variance) 
        
        assert arrays_close(
            cov,
            numpy.array(uncert_core.covariance_matrix(variables)))

    def test_correlated_values_correlation_mat():
        '''
        Tests the input of correlated value.

        Test through their correlation matrix (instead of the
        covariance matrix).
        '''

        x = ufloat(1, 0.1)
        y = ufloat(2, 0.3)
        z = -3*x+y

        cov_mat = uncert_core.covariance_matrix([x, y, z])

        std_devs = numpy.sqrt(numpy.array(cov_mat).diagonal())

        corr_mat = cov_mat/std_devs/std_devs[numpy.newaxis].T

        # We make sure that the correlation matrix is indeed diagonal:
        assert (corr_mat-corr_mat.T).max() <= 1e-15
        # We make sure that there are indeed ones on the diagonal:
        assert (corr_mat.diagonal()-1).max() <= 1e-15

        # We try to recover the correlated variables through the
        # correlation matrix (not through the covariance matrix):

        nominal_values = [v.nominal_value for v in (x, y, z)]
        std_devs = [v.std_dev for v in (x, y, z)]
        x2, y2, z2 = uncert_core.correlated_values_norm(
            list(zip(nominal_values, std_devs)), corr_mat)

        # arrays_close() is used instead of numbers_close() because
        # it compares uncertainties too:

        # Test of individual variables:
        assert arrays_close(numpy.array([x]), numpy.array([x2]))
        assert arrays_close(numpy.array([y]), numpy.array([y2]))
        assert arrays_close(numpy.array([z]), numpy.array([z2]))

        # Partial correlation test:
        assert arrays_close(numpy.array([0]), numpy.array([z2-(-3*x2+y2)]))

        # Test of the full covariance matrix:
        assert arrays_close(
            numpy.array(cov_mat),
            numpy.array(uncert_core.covariance_matrix([x2, y2, z2])))
