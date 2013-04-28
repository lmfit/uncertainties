# coding=utf-8

"""
Tests of the code in uncertainties/__init__.py.

These tests can be run through the Nose testing framework.

(c) 2010-2013 by Eric O. LEBIGOT (EOL).
"""

from __future__ import division

# Standard modules
import copy
import weakref
import math
import random
import sys

# 3rd-party modules
# import nose.tools

# Local modules

import uncertainties
from uncertainties import ufloat, AffineScalarFunc, ufloat_fromstr, isnan
from uncertainties import umath

from uncertainties import __author__

# The following information is useful for making sure that the right
# version of Python is running the tests (for instance with the Travis
# Continuous Integration system):
print "Testing with Python", sys.version

###############################################################################

# Utilities for unit testing

def _numbers_close(x, y, tolerance=1e-6):
    """
    Returns True if the given (real) numbers are close enough.

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

def _ufloats_close(x, y, tolerance=1e-6):
    '''
    Tests if two numbers with uncertainties (or floats) are close, as
    random variables: this is stronger than testing whether their
    nominal value and standard deviation are close.

    The tolerance is applied to both the nominal value and the
    standard deviation of the difference between the numbers.
    '''

    diff = x-y
    return (_numbers_close(diff.nominal_value, 0, tolerance)
            and _numbers_close(diff.std_dev, 0, tolerance))
    
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
                    else uncertainties.Variable(random.random()*4-2, 0)
                    for arg_num in range(num_args)]

                # 'args', but as scalar values:
                args_scalar = [uncertainties.nominal_value(v)
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
                                    " value obtained = %.16f,"
                                    " while numerical approximation = %.16f."
                                    % (arg_num, func.__name__, args,
                                       fixed_deriv_value, num_deriv_value))

            except ValueError, err:  # Arguments out of range, or of wrong type
                # Factorial(real) lands here:
                if str(err).startswith('factorial'):
                    integer_arg_nums = set([0])
                continue  # We try with different arguments
            # Some arguments might have to be integers, for instance:
            except TypeError, err:
                if len(integer_arg_nums) == num_args:
                    raise Exception("Incorrect testing procedure: unable to "
                                    "find correct argument values for %s: %s"
                                    % (func.__name__, err))

                # Another argument might be forced to be an integer:
                integer_arg_nums.add(random.choice(range(num_args)))
            else:
                # We have found reasonable arguments, and the test passed:
                break

###############################################################################

def test_value_construction():
    '''
    Tests the various means of constructing a constant number with
    uncertainty *without a string* (see test_str_input(), for this).
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
        "-13e-2+/-1e2": (-13e-2, 1e2),
        '-14.(15)': (-14, 15),
        '-100.0(15)': (-100, 1.5),
        '14.(15)': (14, 15)
        }
          
    for (representation, values) in tests.iteritems():

        # Without tag:
        num = ufloat_fromstr(representation)
        assert _numbers_close(num.nominal_value, values[0])
        assert _numbers_close(num.std_dev, values[1])
        assert num.tag is None
        
        # With a tag as positional argument:
        num = ufloat_fromstr(representation, 'test variable')
        assert _numbers_close(num.nominal_value, values[0])
        assert _numbers_close(num.std_dev, values[1])
        assert num.tag == 'test variable'

        # With a tag as keyword argument:
        num = ufloat_fromstr(representation, tag='test variable')
        assert _numbers_close(num.nominal_value, values[0])
        assert _numbers_close(num.std_dev, values[1])
        assert num.tag == 'test variable'
        
        ## Obsolete forms

        num = ufloat(representation)  # Obsolete
        assert _numbers_close(num.nominal_value, values[0])
        assert _numbers_close(num.std_dev, values[1])
        assert num.tag is None
        
        # Call with a tag list argument:
        num = ufloat(representation, 'test variable')  # Obsolete
        assert _numbers_close(num.nominal_value, values[0])
        assert _numbers_close(num.std_dev, values[1])
        assert num.tag == 'test variable'
        
        # Call with a tag keyword argument:
        num = ufloat(representation, tag='test variable')  # Obsolete
        assert _numbers_close(num.nominal_value, values[0])
        assert _numbers_close(num.std_dev, values[1])
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
        numerical_derivatives = uncertainties.NumericalDerivatives(
            # The __neg__ etc. methods of AffineScalarFunc only apply,
            # by definition, to AffineScalarFunc objects: we first map
            # possible scalar arguments (used for calculating
            # derivatives) to AffineScalarFunc objects:
            lambda *args: func(*map(uncertainties.to_affine_scalar, args)))
        _compare_derivatives(func, numerical_derivatives, [num_args])

    # Operators that take 1 value:
    for op in uncertainties._modified_operators:
        check_op(op, 1)

    # Operators that take 2 values:
    for op in uncertainties._modified_ops_with_reflection:
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
    assert y in y.derivatives.keys()  # y must not copy the dependence on x
    
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
    assert (uncertainties.covariance_matrix([t, z]) ==
            uncertainties.covariance_matrix([t_copy, z]))

    # However, the relationship between a deep copy and the original
    # variables should be broken, since the deep copy created new,
    # independent variables:
    t_deepcopy = copy.deepcopy(t)
    assert x not in t_deepcopy.derivatives    
    assert (uncertainties.covariance_matrix([t, z]) !=
            uncertainties.covariance_matrix([t_deepcopy, z]))

    # Test of implementations with weak references:

    # Weak references: destroying a variable should never destroy the
    # integrity of its copies (which would happen if the copy keeps a
    # weak reference to the original, in its derivatives member: the
    # weak reference to the original would become invalid):
    del x

    gc.collect()

    assert y in y.derivatives.keys()

## Classes for the pickling tests (put at the module level, so that
## they can be unpickled):
    
# Subclass without slots:
class NewVariable_dict(uncertainties.Variable):
    pass

# Subclass with slots defined by a tuple:
class NewVariable_slots_tuple(uncertainties.Variable):
    __slots__ = ('new_attr',)

# Subclass with slots defined by a string:
class NewVariable_slots_str(uncertainties.Variable):
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
    Wraps a function that takes only positional-or-keyword parameters.
    '''
    
    def f_auto_unc(x, y):
        return 2*x+umath.sin(y)

    # Like f_auto_unc, but does not accept numbers with uncertainties:
    def f(x, y):
        assert not isinstance(x, uncertainties.UFloat)
        assert not isinstance(y, uncertainties.UFloat)
        return f_auto_unc(x, y)

    x = uncertainties.ufloat(1, 0.1)
    y = uncertainties.ufloat(10, 2)

    ### Automatic numerical derivatives:
    
    ## Fully automatic numerical derivatives:
    f_wrapped = uncertainties.wrap(f)
    assert _ufloats_close(f_auto_unc(x, y), f_wrapped(x, y))

    # Call with keyword arguments:
    assert _ufloats_close(f_auto_unc(y=y, x=x), f_wrapped(y=y, x=x))

    ## Automatic additional derivatives for non-defined derivatives:
    f_wrapped = uncertainties.wrap(f, [None])  # No derivative for y
    assert _ufloats_close(f_auto_unc(x, y), f_wrapped(x, y))

    # Call with keyword arguments:
    assert _ufloats_close(f_auto_unc(y=y, x=x), f_wrapped(y=y, x=x))

    ### Explicit derivatives:

    ## Fully defined derivatives:
    f_wrapped = uncertainties.wrap(f, [lambda x, y: 2,
                                       lambda x, y: math.cos(y)])
    
    assert _ufloats_close(f_auto_unc(x, y), f_wrapped(x, y))

    # Call with keyword arguments:
    assert _ufloats_close(f_auto_unc(y=y, x=x), f_wrapped(y=y, x=x))

    ## Automatic additional derivatives for non-defined derivatives:
    f_wrapped = uncertainties.wrap(f, [lambda x, y: 2])  # No derivative for y
    assert _ufloats_close(f_auto_unc(x, y), f_wrapped(x, y))

    # Call with keyword arguments:
    assert _ufloats_close(f_auto_unc(y=y, x=x), f_wrapped(y=y, x=x))

def test_wrapped_func_no_args_no_kwargs():
    '''
    Wraps a function that takes only positional-or-keyword parameters.
    '''
    
    def f_auto_unc(x, y):
        return 2*x+umath.sin(y)

    # Like f_auto_unc, but does not accept numbers with uncertainties:
    def f(x, y):
        assert not isinstance(x, uncertainties.UFloat)
        assert not isinstance(y, uncertainties.UFloat)
        return f_auto_unc(x, y)

    x = uncertainties.ufloat(1, 0.1)
    y = uncertainties.ufloat(10, 2)

    ### Automatic numerical derivatives:
    
    ## Fully automatic numerical derivatives:
    f_wrapped = uncertainties.wrap(f)
    assert _ufloats_close(f_auto_unc(x, y), f_wrapped(x, y))

    # Call with keyword arguments:
    assert _ufloats_close(f_auto_unc(y=y, x=x), f_wrapped(y=y, x=x))

    ## Automatic additional derivatives for non-defined derivatives,
    ## and explicit None derivative:
    f_wrapped = uncertainties.wrap(f, [None])  # No derivative for y
    assert _ufloats_close(f_auto_unc(x, y), f_wrapped(x, y))

    # Call with keyword arguments:
    assert _ufloats_close(f_auto_unc(y=y, x=x), f_wrapped(y=y, x=x))

    ### Explicit derivatives:

    ## Fully defined derivatives:
    f_wrapped = uncertainties.wrap(f, [lambda x, y: 2,
                                       lambda x, y: math.cos(y)])
    
    assert _ufloats_close(f_auto_unc(x, y), f_wrapped(x, y))

    # Call with keyword arguments:
    assert _ufloats_close(f_auto_unc(y=y, x=x), f_wrapped(y=y, x=x))

    ## Automatic additional derivatives for non-defined derivatives:
    f_wrapped = uncertainties.wrap(f, [lambda x, y: 2])  # No derivative for y
    assert _ufloats_close(f_auto_unc(x, y), f_wrapped(x, y))

    # Call with keyword arguments:
    assert _ufloats_close(f_auto_unc(y=y, x=x), f_wrapped(y=y, x=x))

def test_wrapped_func_args_no_kwargs():
    '''
    Wraps a function that takes only positional-or-keyword and
    var-positional parameters.
    '''
    
    def f_auto_unc(x, y, *args):
        return 2*x+umath.sin(y)+3*args[1]

    # Like f_auto_unc, but does not accept numbers with uncertainties:
    def f(x, y, *args):
        assert not any(isinstance(value, uncertainties.UFloat)
                       for value in [x, y] + list(args))
        return f_auto_unc(x, y, *args)

    x = uncertainties.ufloat(1, 0.1)
    y = uncertainties.ufloat(10, 2)
    s = 'string arg'
    z = uncertainties.ufloat(100, 3)

    args = [s, z, s]  # var-positional parameters
    
    ### Automatic numerical derivatives:
    
    ## Fully automatic numerical derivatives:
    f_wrapped = uncertainties.wrap(f)
    assert _ufloats_close(f_auto_unc(x, y, *args), f_wrapped(x, y, *args))

    ## Automatic additional derivatives for non-defined derivatives,
    ## and explicit None derivative:
    f_wrapped = uncertainties.wrap(f, [None])  # No derivative for y
    assert _ufloats_close(f_auto_unc(x, y, *args), f_wrapped(x, y, *args))

    ### Explicit derivatives:

    ## Fully defined derivatives:
    f_wrapped = uncertainties.wrap(f, [lambda x, y, *args: 2,
                                       lambda x, y, *args: math.cos(y),
                                       None,
                                       lambda x, y, *args: 3])
    
    assert _ufloats_close(f_auto_unc(x, y, *args), f_wrapped(x, y, *args))

    ## Automatic additional derivatives for non-defined derivatives:
    
    # No derivative for y:    
    f_wrapped = uncertainties.wrap(f, [lambda x, y, *args: 2])
    assert _ufloats_close(f_auto_unc(x, y, *args), f_wrapped(x, y, *args))

def test_wrapped_func_no_args_kwargs():
    '''
    Wraps a function that takes only positional-or-keyword and
    var-keyword parameters.
    '''
    
    def f_auto_unc(x, y, **kwargs):
        return 2*x+umath.sin(y)+3*kwargs['z']

    # Like f_auto_unc, but does not accept numbers with uncertainties:
    def f(x, y, **kwargs):
        assert not any(isinstance(value, uncertainties.UFloat)
                       for value in [x, y] + kwargs.values())
        return f_auto_unc(x, y, **kwargs)

    x = uncertainties.ufloat(1, 0.1)
    y = uncertainties.ufloat(10, 2)
    s = 'string arg'
    z = uncertainties.ufloat(100, 3)

    kwargs = {'s': s, 'z': z}  # Arguments not in signature

    ### Automatic numerical derivatives:
    
    ## Fully automatic numerical derivatives:
    f_wrapped = uncertainties.wrap(f)
    assert _ufloats_close(f_auto_unc(x, y, **kwargs),
                          f_wrapped(x, y, **kwargs))

    # Call with keyword arguments:
    assert _ufloats_close(f_auto_unc(y=y, x=x, **kwargs),
                          f_wrapped(y=y, x=x, **kwargs))
    
    ## Automatic additional derivatives for non-defined derivatives,
    ## and explicit None derivative:

    # No derivative for positional-or-keyword parameter y, no
    # derivative for optional-keyword parameter z:
    f_wrapped = uncertainties.wrap(f, [None])
    assert _ufloats_close(f_auto_unc(x, y, **kwargs),
                          f_wrapped(x, y, **kwargs))

    # Call with keyword arguments:
    assert _ufloats_close(f_auto_unc(y=y, x=x, **kwargs),
                          f_wrapped(y=y, x=x, **kwargs))

    # No derivative for positional-or-keyword parameter y, no
    # derivative for optional-keyword parameter z:
    f_wrapped = uncertainties.wrap(f, [None], {'z': None})
    assert _ufloats_close(f_auto_unc(x, y, **kwargs),
                          f_wrapped(x, y, **kwargs))

    # Call with keyword arguments:
    assert _ufloats_close(f_auto_unc(y=y, x=x, **kwargs),
                          f_wrapped(y=y, x=x, **kwargs))
    
    # No derivative for positional-or-keyword parameter y, derivative
    # for optional-keyword parameter z:
    f_wrapped = uncertainties.wrap(f, [None],
                                   {'z': lambda x, y, **kwargs: 3})
    assert _ufloats_close(f_auto_unc(x, y, **kwargs),
                          f_wrapped(x, y, **kwargs))

    # Call with keyword arguments:
    assert _ufloats_close(f_auto_unc(y=y, x=x, **kwargs),
                          f_wrapped(y=y, x=x, **kwargs))
    
    ### Explicit derivatives:

    ## Fully defined derivatives:
    f_wrapped = uncertainties.wrap(
        f,
        [lambda x, y, **kwargs: 2, lambda x, y, **kwargs: math.cos(y)],
        {'z:': lambda x, y, **kwargs: 3})
    
    assert _ufloats_close(f_auto_unc(x, y, **kwargs),
                          f_wrapped(x, y, **kwargs))
    # Call with keyword arguments:
    assert _ufloats_close(f_auto_unc(y=y, x=x, **kwargs),
                          f_wrapped(y=y, x=x, **kwargs))
    
    ## Automatic additional derivatives for non-defined derivatives:
    
    # No derivative for y or z:    
    f_wrapped = uncertainties.wrap(f, [lambda x, y, **kwargs: 2])
    assert _ufloats_close(f_auto_unc(x, y, **kwargs),
                          f_wrapped(x, y, **kwargs))

    # Call with keyword arguments:
    assert _ufloats_close(f_auto_unc(y=y, x=x, **kwargs),
                          f_wrapped(y=y, x=x, **kwargs))

def test_wrapped_func_args_kwargs():
    '''
    Wraps a function that takes positional-or-keyword, var-positional
    and var-keyword parameters.
    '''
    
    def f_auto_unc(x, y, *args, **kwargs):
        return 2*x+umath.sin(y)+4*args[1]+3*kwargs['z']

    # Like f_auto_unc, but does not accept numbers with uncertainties:
    def f(x, y, *args, **kwargs):
        assert not any(isinstance(value, uncertainties.UFloat)
                       for value in [x, y]+list(args)+kwargs.values())
        return f_auto_unc(x, y, *args, **kwargs)

    x = uncertainties.ufloat(1, 0.1)
    y = uncertainties.ufloat(10, 2)
    t = uncertainties.ufloat(1000, 4)
    s = 'string arg'
    z = uncertainties.ufloat(100, 3)

    args = [s, t, s]
    kwargs = {'u': s, 'z': z}  # Arguments not in signature

    ### Automatic numerical derivatives:
    
    ## Fully automatic numerical derivatives:
    f_wrapped = uncertainties.wrap(f)
    
    assert _ufloats_close(f_auto_unc(x, y, *args, **kwargs),
                          f_wrapped(x, y, *args, **kwargs), tolerance=1e-5)

    ## Automatic additional derivatives for non-defined derivatives,
    ## and explicit None derivative:

    # No derivative for positional-or-keyword parameter y, no
    # derivative for optional-keyword parameter z:
    f_wrapped = uncertainties.wrap(f, [None, None, None,
                                       lambda x, y, *args, **kwargs: 4])
    assert _ufloats_close(f_auto_unc(x, y, *args, **kwargs),
                          f_wrapped(x, y, *args, **kwargs), tolerance=1e-5)

    # No derivative for positional-or-keyword parameter y, no
    # derivative for optional-keyword parameter z:
    f_wrapped = uncertainties.wrap(f, [None], {'z': None})
    assert _ufloats_close(f_auto_unc(x, y, *args, **kwargs),
                          f_wrapped(x, y, *args, **kwargs), tolerance=1e-5)
    
    # No derivative for positional-or-keyword parameter y, derivative
    # for optional-keyword parameter z:
    f_wrapped = uncertainties.wrap(f, [None],
                                   {'z': lambda x, y, *args, **kwargs: 3})
    assert _ufloats_close(f_auto_unc(x, y, *args, **kwargs),
                          f_wrapped(x, y, *args, **kwargs), tolerance=1e-5)
    
    ### Explicit derivatives:

    ## Fully defined derivatives:
    f_wrapped = uncertainties.wrap(
        f,
        [lambda x, y, *args, **kwargs: 2,
         lambda x, y, *args, **kwargs: math.cos(y)],
        {'z:': lambda x, y, *args, **kwargs: 3})
    
    assert _ufloats_close(f_auto_unc(x, y, *args, **kwargs),
                          f_wrapped(x, y, *args, **kwargs), tolerance=1e-5)
    
    ## Automatic additional derivatives for non-defined derivatives:
    
    # No derivative for y or z:    
    f_wrapped = uncertainties.wrap(f, [lambda x, y, *args, **kwargs: 2])
    assert _ufloats_close(f_auto_unc(x, y, *args, **kwargs),
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
        assert not isinstance(angle, uncertainties.UFloat)
        assert not any(isinstance(arg, uncertainties.UFloat)
                       for arg in list_var)
        return f_auto_unc(angle, *list_var)
    
    f_wrapped = uncertainties.wrap(f)


    my_list = [1, 2, 3]

    ########################################
    # Test of a wrapped function that only calls the original
    # function: it should obtain the exact same result:
    assert f_wrapped(0, *my_list) == f(0, *my_list)
    # 1 == 1 +/- 0, so the type must be checked too:
    assert type(f_wrapped(0, *my_list)) == type(f(0, *my_list))

    ########################################
    # Call with uncertainties:
    
    angle = uncertainties.ufloat(1, 0.1)
    list_value = uncertainties.ufloat(3, 0.2)

    # The random variables must be the same (full correlation):

    assert _ufloats_close(f_wrapped(angle, *[1, angle]),
                          f_auto_unc(angle, *[1, angle]))
    
    assert _ufloats_close(f_wrapped(angle, *[list_value, angle]),
                          f_auto_unc(angle, *[list_value, angle]))
    
    ########################################
    # Non-numerical arguments, and  explicit and implicit derivatives:
    def f(x, y, z, t, u):
        return x+2*z+3*t+4*u
    
    f_wrapped = uncertainties.wrap(
        f, [lambda *args: 1, None, lambda *args:2, None])  # No deriv. for u

    assert f_wrapped(10, 'string argument', 1, 0, 0) == 12

    x = uncertainties.ufloat(10, 1)

    assert _numbers_close(f_wrapped(x, 'string argument', x, x, x).std_dev,
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

        for value in [x, y]+list(args)+kwargs.values():
            assert not isinstance(value, uncertainties.UFloat)
        
        return f_auto_unc(x, y, *args, **kwargs)
    
    f_wrapped = uncertainties.wrap(f)


    x = ufloat(1, 0.1)
    y = ufloat(10, 0.11)
    z = ufloat(100, 0.111)
    t = ufloat(0.1, 0.1111)
        
    assert _ufloats_close(f_wrapped(x, y, z, t=t),
                          f_auto_unc(x, y, z, t=t), tolerance=1e-5)

    ########################################

    # We make sure that analytical derivatives are indeed used. We
    # also test the automatic handling of additional *args arguments
    # beyond the number of supplied derivatives.

    f_wrapped2 = uncertainties.wrap(
        f, [None, lambda x, y, *args, **kwargs: math.cos(y)])

    # The derivatives must be perfectly identical:

    # The *args parameter of f() is given as a keyword argument, so as
    # to try to confuse the code:
    
    assert (f_wrapped2(x, y, z, t=t).derivatives[y]
            == f_auto_unc(x, y, z, t=t).derivatives[y])
    
    # Derivatives supplied through the keyword-parameter dictionary of
    # derivatives, and also derivatives supplied for the
    # var-positional arguments (*args[0]):

    f_wrapped3 = uncertainties.wrap(
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

    f_wrapped4 = uncertainties.wrap(
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
    assert uncertainties.std_dev(x) == x.std_dev
    assert uncertainties.std_dev(y) == y.std_dev

    # std_dev for other objects:
    assert uncertainties.std_dev([]) == 0
    assert uncertainties.std_dev(None) == 0
    
###############################################################################

def test_covariances():
    "Covariance matrix"

    x = ufloat(1, 0.1)
    y = -2*x+10
    z = -3*x
    covs = uncertainties.covariance_matrix([x, y, z])
    # Diagonal elements are simple:
    assert _numbers_close(covs[0][0], 0.01)
    assert _numbers_close(covs[1][1], 0.04)
    assert _numbers_close(covs[2][2], 0.09)
    # Non-diagonal elements:
    assert _numbers_close(covs[0][1], -0.02)

    
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
    if sys.version_info >= (2, 6):
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
        
        m1, m2 -- NumPy matrices.
        precision -- precision passed through to
        uncertainties.test_uncertainties._numbers_close().
        """

        # ! numpy.allclose() is similar to this function, but does not
        # work on arrays that contain numbers with uncertainties, because
        # of the isinf() function.

        for (elmt1, elmt2) in zip(m1.flat, m2.flat):

            # For a simpler comparison, both elements are
            # converted to AffineScalarFunc objects:
            elmt1 = uncertainties.to_affine_scalar(elmt1)
            elmt2 = uncertainties.to_affine_scalar(elmt2)

            if not _numbers_close(elmt1.nominal_value,
                                  elmt2.nominal_value, precision):
                return False

            if not _numbers_close(elmt1.std_dev,
                                  elmt2.std_dev, precision):
                return False
        return True

    
    def test_numpy_comparison():
        "Comparison with a Numpy array."

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

        u = uncertainties.ufloat(1, 0.1)
        cov = uncertainties.covariance_matrix([u])
        # "1" is used instead of u.nominal_value because
        # u.nominal_value might return a float.  The idea is to force
        # the new variable u2 to be defined through an integer nominal
        # value:
        u2, = uncertainties.correlated_values([1], cov)
        expr = 2*u2  # Calculations with u2 should be possible, like with u

        ####################    

        # Covariances between output and input variables:

        x = ufloat(1, 0.1)
        y = ufloat(2, 0.3)
        z = -3*x+y

        covs = uncertainties.covariance_matrix([x, y, z])

        # Test of the diagonal covariance elements:
        assert arrays_close(
            numpy.array([v.std_dev**2 for v in (x, y, z)]),
            numpy.array(covs).diagonal())
        
        # "Inversion" of the covariance matrix: creation of new
        # variables:
        (x_new, y_new, z_new) = uncertainties.correlated_values(
            [x.nominal_value, y.nominal_value, z.nominal_value],
            covs,
            tags = ['x', 'y', 'z'])

        # Even the uncertainties should be correctly reconstructed:
        assert arrays_close(numpy.array((x, y, z)),
                              numpy.array((x_new, y_new, z_new)))

        # ... and the covariances too:
        assert arrays_close(
            numpy.array(covs),
            numpy.array(uncertainties.covariance_matrix([x_new, y_new, z_new])))

        assert arrays_close(
            numpy.array([z_new]), numpy.array([-3*x_new+y_new]))

        ####################

        # ... as well as functional relations:

        u = ufloat(1, 0.05)
        v = ufloat(10, 0.1)
        sum_value = u+2*v

        # Covariance matrices:
        cov_matrix = uncertainties.covariance_matrix([u, v, sum_value])

        # Correlated variables can be constructed from a covariance
        # matrix, if NumPy is available:
        (u2, v2, sum2) = uncertainties.correlated_values(
            [x.nominal_value for x in [u, v, sum_value]],
            cov_matrix)

        # arrays_close() is used instead of _numbers_close() because
        # it compares uncertainties too:
        assert arrays_close(numpy.array([u]), numpy.array([u2]))
        assert arrays_close(numpy.array([v]), numpy.array([v2]))
        assert arrays_close(numpy.array([sum_value]), numpy.array([sum2]))
        assert arrays_close(numpy.array([0]),
                              numpy.array([sum2-(u2+2*v2)]))


    def test_correlated_values_correlation_mat():
        '''
        Tests the input of correlated value.

        Test through their correlation matrix (instead of the
        covariance matrix).
        '''
        
        x = ufloat(1, 0.1)
        y = ufloat(2, 0.3)
        z = -3*x+y

        cov_mat = uncertainties.covariance_matrix([x, y, z])

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
        x2, y2, z2 = uncertainties.correlated_values_norm(
            zip(nominal_values, std_devs), corr_mat)
        
        # arrays_close() is used instead of _numbers_close() because
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
            numpy.array(uncertainties.covariance_matrix([x2, y2, z2])))

