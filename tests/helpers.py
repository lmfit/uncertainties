import random
from math import isclose, isnan, isinf

import uncertainties.core as uncert_core
from uncertainties.core import ufloat, AffineScalarFunc


power_derivative_cases = (
    ((-0.3, 0.01), (-3.0, 0.0), -370.37037037037044, float("nan")),
    ((-0.3, 0.01), (1.0, 0.1), 1.0, float("nan")),
    ((-0.3, 0.01), (0.0, 0.1), 0.0, float("nan")),
    ((0.0, 0.1), (3.1, 0.01), float("nan"), 0.0),
    ((0.0, 0.1), (1.0, 0.1), 1.0, 0.0),
    ((0.0, 0.1), (2.0, 0.2), 0.0, 0.0),
    ((0.0, 0.1), (0.3, 0.01), float("nan"), 0.0),
    ((0.0, 0.1), (0.0, 0.1), 0.0, float("nan")),
    ((0.3, 0.01), (0.3, 0.01), 0.696845301935949, -0.8389827923531782),
    ((0.3, 0.01), (0.0, 0.1), 0.0, -1.2039728043259361),
    ((0.3, 0.01), (-0.3, 0.01), -1.4350387341664474, -1.7277476090907193),
)


def power_all_cases(op):
    for (
        (first_val, first_std),
        (second_val, second_std),
        first_der,
        second_der,
    ) in power_derivative_cases:
        first_ufloat = ufloat(first_val, first_std)
        second_ufloat = ufloat(second_val, second_std)
        result = op(first_ufloat, second_ufloat)
        first_der_result = result.derivatives[first_ufloat]
        second_der_result = result.derivatives[second_ufloat]
        if isnan(first_der):
            assert isnan(first_der_result)
        else:
            assert isclose(first_der_result, first_der)
        if isnan(second_der):
            assert isnan(second_der_result)
        else:
            assert isclose(second_der_result, second_der)


zero = ufloat(0, 0)
one = ufloat(1, 0)
p = ufloat(0.3, 0.01)

power_float_result_cases = [
    (0, p, 0),
    (zero, p, 0),
    (float("nan"), zero, 1),
    (one, float("nan"), 1),
    (p, 0, 1),
    (zero, 0, 1),
    (-p, 0, 1),
    (-10.3, zero, 1),
    (0, zero, 1),
    (0.3, zero, 1),
    (-p, zero, 1),
    (zero, zero, 1),
    (p, zero, 1),
    (one, -3, 1),
    (one, -3.1, 1),
    (one, 0, 1),
    (one, 3, 1),
    (one, 3.1, 1),
    (one, -p, 1),
    (one, zero, 1),
    (one, p, 1),
    (1, -p, 1),
    (1, zero, 1),
    (1, p, 1),
]


def power_special_cases(op):
    for first, second, result in power_float_result_cases:
        assert op(first, second) == result


def power_wrt_ref(op, ref_op):
    """
    Checks special cases of the uncertainty power operator op (where
    op is typically the built-in pow or uncertainties.umath.pow), by
    comparing its results to the reference power operator ref_op
    (which is typically the built-in pow or math.pow).
    """

    # Negative numbers with uncertainty can be exponentiated to an
    # integral power:
    assert op(ufloat(-1.1, 0.1), -9).nominal_value == ref_op(-1.1, -9)

    # Case of numbers with no uncertainty: should give the same result
    # as numbers with uncertainties:
    assert op(ufloat(-1, 0), 9) == ref_op(-1, 9)
    assert op(ufloat(-1.1, 0), 9) == ref_op(-1.1, 9)


###############################################################
# TODO: move to uncertainties/testing.py
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
            return 2 * abs(x - y) / (abs(x) + abs(y)) < tolerance

    else:  # Either x or y is zero
        return abs(x or y) < tolerance


def ufloats_close(x, y, tolerance=1e-6):
    """
    Tests if two numbers with uncertainties are close, as random
    variables: this is stronger than testing whether their nominal
    value and standard deviation are close.

    The tolerance is applied to both the nominal value and the
    standard deviation of the difference between the numbers.
    """

    diff = x - y
    return numbers_close(diff.nominal_value, 0, tolerance) and numbers_close(
        diff.std_dev, 0, tolerance
    )


class DerivativesDiffer(Exception):
    pass


def compare_derivatives(func, numerical_derivatives, num_args_list=None):
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
            "atanh": [1],
            "log": [1, 2],  # Both numbers of arguments are tested
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
                    func(*(1,) * num_args)
                except TypeError:
                    pass  # Not the right number of arguments
                else:  # No error
                    # num_args is a good number of arguments for func:
                    num_args_list.append(num_args)

            if not num_args_list:
                raise Exception(
                    "Can't find a reasonable number of arguments"
                    " for function '%s'." % funcname
                )

    for num_args in num_args_list:
        # Argument numbers that will have a random integer value:
        integer_arg_nums = set()

        if funcname == "ldexp":
            # The second argument must be an integer:
            integer_arg_nums.add(1)

        while True:
            try:
                # We include negative numbers, for more thorough tests:
                args = []
                for arg_num in range(num_args):
                    if arg_num in integer_arg_nums:
                        args.append(random.choice(range(-10, 10)))
                    else:
                        args.append(uncert_core.Variable(random.random() * 4 - 2, 0))

                # 'args', but as scalar values:
                args_scalar = [uncert_core.nominal_value(v) for v in args]

                func_approx = func(*args)

                # Some functions yield simple Python constants, after
                # wrapping in wrap(): no test has to be performed.
                # Some functions also yield tuples...
                if isinstance(func_approx, AffineScalarFunc):
                    # We compare all derivatives:
                    for arg_num, (arg, numerical_deriv) in enumerate(
                        zip(args, numerical_derivatives)
                    ):
                        # Some arguments might not be differentiable:
                        if isinstance(arg, int):
                            continue

                        fixed_deriv_value = func_approx.derivatives[arg]

                        num_deriv_value = numerical_deriv(*args_scalar)

                        # This message is useful: the user can see that
                        # tests are really performed (instead of not being
                        # performed, silently):
                        print(
                            "Testing derivative #%d of %s at %s"
                            % (arg_num, funcname, args_scalar)
                        )

                        if not numbers_close(fixed_deriv_value, num_deriv_value, 1e-4):
                            # It is possible that the result is NaN:
                            if not isnan(func_approx):
                                raise DerivativesDiffer(
                                    "Derivative #%d of function '%s' may be"
                                    " wrong: at args = %s,"
                                    " value obtained = %.16f,"
                                    " while numerical approximation = %.16f."
                                    % (
                                        arg_num,
                                        funcname,
                                        args,
                                        fixed_deriv_value,
                                        num_deriv_value,
                                    )
                                )

            except ValueError as err:  # Arguments out of range, or of wrong type
                # Factorial(real) lands here:
                if str(err).startswith("factorial"):
                    integer_arg_nums = set([0])
                continue  # We try with different arguments
            # Some arguments might have to be integers, for instance:
            except TypeError as err:
                if len(integer_arg_nums) == num_args:
                    raise Exception(
                        "Incorrect testing procedure: unable to "
                        "find correct argument values for %s: %s" % (funcname, err)
                    )

                # Another argument might be forced to be an integer:
                integer_arg_nums.add(random.choice(range(num_args)))
            else:
                # We have found reasonable arguments, and the test passed:
                break


###############################################################################


try:
    import numpy  # noqa
except ImportError:
    pass
else:

    def uarrays_close(m1, m2, precision=1e-4):
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

        for elmt1, elmt2 in zip(m1.flat, m2.flat):
            # For a simpler comparison, both elements are
            # converted to AffineScalarFunc objects:
            elmt1 = uncert_core.to_affine_scalar(elmt1)
            elmt2 = uncert_core.to_affine_scalar(elmt2)

            if not numbers_close(elmt1.nominal_value, elmt2.nominal_value, precision):
                return False

            if not numbers_close(elmt1.std_dev, elmt2.std_dev, precision):
                return False

        return True
