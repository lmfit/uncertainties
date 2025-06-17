# Some tests are already performed in test_unumpy (unumpy contains a
# matrix inversion, for instance).  They are not repeated here.

try:
    import numpy
except ImportError:
    import sys

    sys.exit()  # There is no reason to test the interface to NumPy

from uncertainties import unumpy, ufloat
from helpers import uarrays_close


def test_list_inverse():
    "Test of the inversion of a square matrix"

    mat_list = [[1, 1], [1, 0]]

    # numpy.linalg.inv(mat_list) does calculate the inverse even
    # though mat_list is a list of lists (and not a matrix).  Can
    # ulinalg do the same?  Here is a test:
    mat_list_inv = unumpy.ulinalg.inv(mat_list)

    # unumpy.ulinalg should behave in the same way as numpy.linalg,
    # with respect to types:
    mat_list_inv_numpy = numpy.linalg.inv(mat_list)
    assert type(mat_list_inv) == type(mat_list_inv_numpy)

    # Individual element check:
    assert isinstance(mat_list_inv[1, 1], float)
    assert mat_list_inv[1, 1] == -1


def test_list_pseudo_inverse():
    "Test of the pseudo-inverse"

    x = ufloat(1, 0.1)
    y = ufloat(2, 0.1)
    mat = numpy.array([[x, x], [y, 0]])

    # Internal consistency: the inverse and the pseudo-inverse yield
    # the same result on square matrices:
    assert uarrays_close(unumpy.ulinalg.inv(mat), unumpy.ulinalg.pinv(mat), 1e-4).all()
    assert uarrays_close(
        unumpy.ulinalg.inv(mat),
        # Support for the optional pinv argument is
        # tested:
        unumpy.ulinalg.pinv(mat, 1e-15),
        1e-4,
    ).all()
