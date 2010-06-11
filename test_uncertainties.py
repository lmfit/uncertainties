"""
Tests of the code in uncertainties/__init__.py.

These tests can be run through the Nose testing framework.

(c) 2010 by Eric O. LEBIGOT (EOL).
"""

from __future__ import division

# Standard modules
import copy
import weakref
import math

# 3rd-party modules
import nose.tools

# Local modules

import uncertainties
from uncertainties import ufloat, AffineScalarFunc

from uncertainties import __author__

###############################################################################

# Utilities for unit testing

def _numbers_close(x, y, tolerance=1e-6):
    """
    Returns True if the numbers are close enough.

    The given tolerance is the relative difference allowed, or the absolute
    difference, if one of the numbers is 0.
    """

    # Instead of using a try and ZeroDivisionError, we do a test,
    # NaN could appear silently:

    if x != 0 and y != 0:
        return abs(1-y/x) < tolerance
    else:
        if x == 0:
            return abs(y) < tolerance
        else:
            return abs(x) < tolerance

class DerivativesDiffer(Exception):
    pass

def _compare_derivatives(func, numerical_derivatives,
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

    import random

    # print "Testing", func.__name__

    if not num_args_list: 

        # Detecting automatically the correct number of arguments is not
        # always easy (because not all values are allowed, etc.):

        num_args_table = {
            'atanh': [1],
            'log': [1, 2]  # Both numbers of arguments are tested
            }
        if func.__name__ in num_args_table:
            num_args_list = num_args_table[func.__name__]
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
                                " for function '%s'." % func.__name__)

    for num_args in num_args_list:

        # Argument numbers that will have a random integer value:
        integer_arg_nums = set()

        if func.__name__ == 'ldexp':
            # The second argument must be an integer:
            integer_arg_nums.add(1)

        while True:
            try:

                # We include negative numbers, for more thorough tests:
                args = [
                    random.choice(range(-10, 10))
                    if arg_num in integer_arg_nums
                    else uncertainties.Variable(random.random()*10-5, 0)
                    for arg_num in range(num_args)]

                func_approx = func(*args)

                # Some functions are simple Python constants, after
                # wrapping in wrap(): no test has to be performed:
                if isinstance(func_approx, AffineScalarFunc):

                    # We compare all derivatives:
                    for (arg_num, (arg, numerical_deriv)) in \
                        enumerate(zip(args, numerical_derivatives)):

                        # Some arguments might not be differentiable:
                        if isinstance(arg, int):
                            continue

                        fixed_deriv_value = func_approx.derivatives[arg]
                        num_deriv_value = numerical_deriv(*args)

                        # This message is useful: the user can see that
                        # tests are really performed (instead of not being
                        # performed, silently):
                        print "Testing %s at %s, arg #%d" % (
                            func.__name__, args, arg_num)

                        if not _numbers_close(fixed_deriv_value,
                                              num_deriv_value, 1e-4):

                            # It is possible that the result is NaN:

                            # ! Python 2.6+: this would be
                            # not math.isnan(func_approx):
                            if func_approx == func_approx:
                                raise DerivativesDiffer(
                                    "Derivative #%d of function '%s' may be"
                                    " wrong: at args = %s,"
                                    " value obtained = %f,"
                                    " while numerical approximation = %f."
                                    % (arg_num, func.__name__, args,
                                       fixed_deriv_value, num_deriv_value))

            except ValueError, err:  # Arguments out of range, or of wrong type
                # Factorial(real) lands here:
                if str(err).startswith('factorial'):
                    integer_arg_nums = set([0])
                continue  # We try with different arguments
            # Some arguments might have to be integers, for instance:
            except TypeError:
                if len(integer_arg_nums) == num_args:
                    raise Exception("Incorrect testing procedure: unable to "
                                    "find correct argument values for %s."
                                    % func.__name__)

                # Another argument might be forced to be an integer:
                integer_arg_nums.add(random.choice(range(num_args)))
            else:
                # We have found reasonable arguments, and the test passed:
                break

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
        # print "Checking %s..." % op_string
        func = getattr(AffineScalarFunc, op_string)
        numerical_derivatives = uncertainties.NumericalDerivatives(
            lambda *args: func(*args).nominal_value)
        _compare_derivatives(func, numerical_derivatives, [num_args])

    # Operators that take 1 value:
    for op in uncertainties._modified_operators:
        check_op(op, 1)

    # Operators that take 2 values:
    for op in uncertainties._ops_with_reflection.keys():
        check_op(op, 2)

# Additional, more complex checks, for use with the nose unit testing
# framework.

def test_copy():
    "Standard copy module integration"
    import gc
    
    x = ufloat((3, 0.1))
    assert x == x
    
    y = copy.copy(x)
    assert x != y
    assert not(x == y)
    assert y in y.derivatives.keys()  # y must not copy the dependence on x
    
    z = copy.deepcopy(x)
    assert x != z

    
    # Weak references: destroying a variable should never destroy the
    # integrity of its copies (which would happen if the copy keeps a
    # weak reference to the original, in its derivatives member: the
    # weak reference to the original would become invalid):
    del x

    gc.collect()

    assert y in y.derivatives.keys()
    
def test_pickling():
    "Standard pickle module integration."

    import pickle

    x = ufloat((2, 0.1))

    x_unpickled = pickle.loads(pickle.dumps(x))

    assert x != x_unpickled  # Pickling creates copies

    ## Tests with correlations and AffineScalarFunc objects:
    f = 2*x
    assert isinstance(f, AffineScalarFunc)
    (f_unpickled, x_unpickled2) = pickle.loads(pickle.dumps((f, x)))
    # Correlations must be preserved:
    assert f_unpickled - x_unpickled2 - x_unpickled2 == 0
    
    
def test_int_div():
    "Integer division"
    # We perform all operations on floats, because derivatives can
    # otherwise be meaningless:
    x = ufloat((3.9, 2))//2
    assert x.nominal_value == 1.
    # All errors are supposed to be small, so the ufloat()
    # in x violates the assumption.  Therefore, the following is
    # correct:
    assert x.std_dev() == 0.0

def test_comparison_ops():
    "Test of comparison operators"

    import random
    
    # Operations on quantities equivalent to Python numbers must still
    # be correct:
    a = ufloat((-3, 0))
    b = ufloat((10, 0))
    c = ufloat((10, 0))
    assert a < b
    assert a < 3
    assert 3 < b  # This is first given to int.__lt__()
    assert b == c

    x = ufloat((3, 0.1))
    
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
    assert 2*ufloat((1, 0.1)) != ufloat((2, 0.2))    
    # Comparaison between 2 _different_ Variable objects
    # that are uncorrelated:
    assert x != ufloat((3, 0.1))
    
    assert x != ufloat((3, 0.2))

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
            return ((random.random()-0.5) * min(var.std_dev(), 1e-5)
                    + var.nominal_value)

        # All operations are tested:
        for op in ("__%s__" % name
                   for name in('ne', 'eq', 'lt', 'le', 'gt', 'ge')):

            float_func = getattr(float, op)
            
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
                print "Sampling results:", sampled_results
                raise Exception("Semantic value of %s %s (%s) %s not"
                                " correctly reproduced."
                                % (x, op, y, correct_result))

    # With different numbers:
    test_all_comparison_ops(ufloat((3, 0.1)),
                            ufloat((-2, 0.1)))
    test_all_comparison_ops(ufloat((0, 0)),  # Special number
                            ufloat((1, 1)))
    test_all_comparison_ops(ufloat((0, 0)),  # Special number
                            ufloat((0, 0.1)))
    # With identical numbers:
    test_all_comparison_ops(ufloat((0, 0)),
                            ufloat((0, 0)))
    test_all_comparison_ops(ufloat((1, 1)),
                            ufloat((1, 1)))

    
def test_logic():
    "Boolean logic: __nonzero__, bool."

    x = ufloat((3, 0))
    y = ufloat((0, 0))
    z = ufloat((0, 0.1))
    t = ufloat((-1, 2))

    assert bool(x) == True
    assert bool(y) == False
    assert bool(z) == True
    assert bool(t) == True  # Only infinitseimal neighborhood are used

        
    
def test_basic_access_to_data():
    "Access to data from Variable and AffineScalarFunc objects."

    x = ufloat((3.14, 0.01), "x var")
    assert x.tag == "x var"
    assert x.nominal_value == 3.14
    assert x.std_dev() == 0.01

    # Case of AffineScalarFunc objects:
    y = x + 0
    assert type(y) == AffineScalarFunc
    assert y.nominal_value == 3.14
    assert y.std_dev() == 0.01

    # Details on the sources of error:
    a = ufloat((-1, 0.001))
    y = 2*x + 3*x + 2 + a
    error_sources = y.error_components()
    assert len(error_sources) == 2  # 'a' and 'x'
    assert error_sources[x] == 0.05
    assert error_sources[a] == 0.001

    # Derivative values should be available:
    assert y.derivatives[x] == 5

    # Modification of the standard deviation of variables:
    x.set_std_dev(1)
    assert y.error_components()[x] == 5  # New error contribution!

    # Calculation of deviations in units of the standard deviations:
    assert 10/x.std_dev() == x.position_in_sigmas(10 + x.nominal_value)

    # "In units of the standard deviation" is not always meaningfull:
    x.set_std_dev(0)
    try:
        x.position_in_sigmas(1)
    except ValueError:
        pass  # Normal behavior

def test_correlations():
    "Correlations between variables"

    a = ufloat((1, 0))
    x = ufloat((4, 0.1))
    y = x*2 + a
    # Correlations cancel "naive" additions of uncertainties:
    assert y.std_dev() != 0
    normally_zero = y - (x*2 + 1)
    assert normally_zero.nominal_value == 0
    assert normally_zero.std_dev() == 0

def test_str_input():

    "Input of numbers with uncertainties as a string"

    # String representation, and numerical values:
    tests = {
        "-1.23(3.4)": (-1.23, 3.4),  # (Nominal value, error)
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
        "-13e-2+/-1e2": (-13e-2, 1e2)
        }
          
    for (representation, values) in tests.iteritems():
        
        num = ufloat(representation)

        assert _numbers_close(num.nominal_value, values[0])
        assert _numbers_close(num.std_dev(), values[1])

    
def test_no_coercion():
    """
    Coercion of Variable object to a simple float.

    The coercion should be impossible, like for complex numbers.
    """

    x = ufloat((4, 1))
    try:
        assert float(x) == 4
    except TypeError:
        pass
    else:
        raise Exception("Conversion to float() should fail with TypeError")

def test_wrapped_func():
    """
    Test uncertainty-aware functions obtained through wrapping.
    """

    # This function can be wrapped so that it works when 'angle' has
    # an uncertainty (math.cos does not handle numbers with
    # uncertainties):
    def f(angle, list_var):
        return math.cos(angle) + sum(list_var)

    f_wrapped = uncertainties.wrap(f)
    my_list = [1, 2, 3]

    # Test of a wrapped function that only calls the original function:
    assert f_wrapped(0, my_list) == 1 + sum(my_list)

    # As a precaution, the wrapped function does not venture into
    # calculating f with uncertainties when one of the argument is not
    # a simple number, because this argument might contain variables:
    angle = ufloat((0, 0.1))

    assert f_wrapped(angle, [angle, angle]) == NotImplemented
    assert f_wrapped(angle, my_list) == NotImplemented
    

###############################################################################

# The tests below require NumPy:
try:
    import numpy
except ImportError:
    pass
else:
    def test_numpy_comparison():
        "Comparison with a Numpy array"

        x = ufloat((1, 0.1))
        
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

###############################################################################
        
def test_access_to_std_dev():
    "Uniform access to the standard deviation"

    x = ufloat((1, 0.1))
    y = 2*x

    # std_dev for Variable and AffineScalarFunc objects:
    assert uncertainties.std_dev(x) == x.std_dev()
    assert uncertainties.std_dev(y) == y.std_dev()

    # std_dev for other objects:
    assert uncertainties.std_dev([]) == 0
    assert uncertainties.std_dev(None) == 0
    
def test_legacy():
    "Test of legacy code"
    x = (1, 0.1)

    old1 = uncertainties.NumberWithUncert(x)  # Warnings are normal
    old2 = uncertainties.num_with_uncert(x)  # Warnings are normal
    new = uncertainties.ufloat(x)
    
    assert old1.nominal_value == new.nominal_value
    assert old2.nominal_value == new.nominal_value
    assert old1.std_dev() == new.std_dev()
    assert old2.std_dev() == new.std_dev()
