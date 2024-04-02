from uncertainties.core import ufloat
from math import isnan

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

