from math import isnan
import sys

from uncertainties import ufloat, ufloat_fromstr, formatting

from helpers import numbers_close


def test_PDG_precision():
    """
    Test of the calculation of the number of significant digits for
    the uncertainty.
    """

    # The 3 cases of the rounding rules are covered in each case:
    tests = {
        # Very big floats:
        1.7976931348623157e308: (2, 1.7976931348623157e308),
        0.5e308: (1, 0.5e308),
        0.9976931348623157e308: (2, 1e308),
        # Very small floats:
        1.3e-323: (2, 1.3e-323),
        5e-324: (1, 5e-324),
        9.99e-324: (2, 1e-323),
    }

    for std_dev, result in tests.items():
        assert formatting.PDG_precision(std_dev) == result


def test_repr():
    """Test the representation of numbers with uncertainty."""

    # The uncertainty is a power of 2, so that it can be exactly
    # represented:
    x = ufloat(3.14159265358979, 0.25)
    assert repr(x) == "3.14159265358979+/-0.25"

    x = ufloat(3.14159265358979, 0)
    assert repr(x) == "3.14159265358979+/-0"

    # Tagging:
    x = ufloat(3, 1, "length")
    assert repr(x) == "< length = 3.0+/-1.0 >"


def test_format():
    """Test the formatting of numbers with uncertainty."""

    # The way NaN is formatted with F, E and G depends on the version
    # of Python (NAN for Python 2.5+ at least):
    NaN_EFG = "%F" % float("nan")

    # !! The way NaN is formatted with F, E and G might depend on the
    # version of Python, if it is like NaN (could be tested with
    # Python 2.3 or 2.4 vs Python 2.7):
    Inf_EFG = "%F" % float("inf")

    # Tests of each point of the docstring of
    # AffineScalarFunc.__format__() in turn, mostly in the same order.

    # The LaTeX tests do not use the customization of
    # uncert_core.GROUP_SYMBOLS and uncert_core.EXP_PRINT: this
    # way, problems in the customization themselves are caught.

    tests = {  # (Nominal value, uncertainty): {format: result,...}
        # Usual float formatting, and individual widths, etc.:
        (3.1415, 0.0001): {
            "*^+7.2f": "*+3.14*+/-*0.00**",
            "+07.2f": "+003.14+/-0000.00",  # 0 fill
            ">10f": "  3.141500+/-  0.000100",  # Width and align
            "11.3e": "  3.142e+00+/-  0.000e+00",  # Duplicated exponent
            "0.4e": "3.1415e+00+/-0.0000e+00",  # Forced double exponent
        },
        # Full generalization of float formatting:
        (3.1415, 0.0001): {  # noqa
            "+09.2uf": "+03.14150+/-000.00010",
            # Alignment is not available with the % formatting
            # operator of Python < 2.6:
            "*^+9.2uf": "+3.14150*+/-*0.00010*",
            ">9f": "  3.14150+/-  0.00010",  # Width and align
        },
        # Number of digits of the uncertainty fixed:
        (123.456789, 0.00123): {
            ".1uf": "123.457+/-0.001",
            ".2uf": "123.4568+/-0.0012",
            ".3uf": "123.45679+/-0.00123",
            ".2ue": "(1.234568+/-0.000012)e+02",
        },
        # Sign handling:
        (-123.456789, 0.00123): {
            ".1uf": "-123.457+/-0.001",
            ".2uf": "-123.4568+/-0.0012",
            ".3uf": "-123.45679+/-0.00123",
            ".2ue": "(-1.234568+/-0.000012)e+02",
        },
        # Uncertainty larger than the nominal value:
        (12.3, 456.78): {"": "12+/-457", ".1uf": "12+/-457", ".4uf": "12.3+/-456.8"},
        # ... Same thing, but with an exponent:
        (12.3, 456.78): {  # noqa
            ".1ue": "(0+/-5)e+02",
            ".4ue": "(0.123+/-4.568)e+02",
            ".4ueS": "0.123(4.568)e+02",
        },
        (23456.789123, 1234.56789123): {".6gS": "23456.8(1234.6)"},
        # Test of the various float formats: the nominal value should
        # have a similar representation as if it were directly
        # represented as a float:
        (1234567.89, 0.1): {
            ".0e": "(1+/-0)e+06",
            "e": "(1.23456789+/-0.00000010)e+06",
            "E": "(1.23456789+/-0.00000010)E+06",
            "f": "1234567.89+/-0.10",
            "F": "1234567.89+/-0.10",
            "g": "1234567.89+/-0.10",
            "G": "1234567.89+/-0.10",
            "%": "(123456789+/-10)%",
        },
        (1234567.89, 4.3): {"g": "1234568+/-4"},
        (1234567.89, 43): {  # Case where g triggers the exponent notation
            "g": "(1.23457+/-0.00004)e+06",
            "G": "(1.23457+/-0.00004)E+06",
        },
        (3.1415, 0.0001): {"+09.2uf": "+03.14150+/-000.00010"},  # noqa
        (1234.56789, 0.1): {
            ".0f": "(1234+/-0.)",  # Approximate error indicated with "."
            "e": "(1.23456+/-0.00010)e+03",
            "E": "(1.23456+/-0.00010)E+03",
            "f": "1234.57+/-0.10",
            "F": "1234.57+/-0.10",
            "f": "1234.57+/-0.10",  # noqa
            "F": "1234.57+/-0.10",  # noqa
            "%": "123457+/-10%",
        },
        # Percent notation:
        (0.42, 0.0055): {
            # Because '%' does 0.0055*100, the value
            # 0.5499999999999999 is obtained, which rounds to 0.5. The
            # original rounded value is 0.006. The same behavior is
            # found in Python 2.7: '{:.1%}'.format(0.0055) is '0.5%'.
            ".1u%": "(42.0+/-0.5)%",
            ".1u%S": "42.0(5)%",
            "%P": "(42.0±0.5)%",
        },
        # Particle Data Group automatic convention, including limit cases:
        (1.2345678, 0.354): {"": "1.23+/-0.35"},
        (1.2345678, 0.3549): {"": "1.23+/-0.35"},
        (1.2345678, 0.355): {"": "1.2+/-0.4"},
        (1.5678, 0.355): {"": "1.6+/-0.4"},
        (1.2345678, 0.09499): {"": "1.23+/-0.09"},
        (1.2345678, 0.095): {"": "1.23+/-0.10"},
        # Automatic extension of the uncertainty up to the decimal
        # point:
        (1000, 123): {
            ".1uf": "1000+/-123",
            # The nominal value has 1 <= mantissa < 10. The precision
            # is the number of significant digits of the uncertainty:
            ".1ue": "(1.0+/-0.1)e+03",
        },
        # Spectroscopic notation:
        (-1.23, 3.4): {
            "S": "-1.2(3.4)",
            ".2ufS": "-1.2(3.4)",
            ".3ufS": "-1.23(3.40)",
        },
        (-123.456, 0.123): {
            "S": "-123.46(12)",
            ".1ufS": "-123.5(1)",
            ".2ufS": "-123.46(12)",
            ".3ufS": "-123.456(123)",
        },
        (-123.456, 0.567): {
            "S": "-123.5(6)",
            ".1ufS": "-123.5(6)",
            ".2ufS": "-123.46(57)",
            ".3ufS": "-123.456(567)",
        },
        (-123.456, 0.004): {
            # The decimal point shows that the uncertainty is not
            # exact:
            ".2fS": "-123.46(0.00)"
        },
        # LaTeX notation:
        #
        (1234.56789, 0.1): {  # noqa
            "eL": r"\left(1.23457 \pm 0.00010\right) \times 10^{3}",
            "EL": r"\left(1.23457 \pm 0.00010\right) \times 10^{3}",
            "fL": r"1234.57 \pm 0.10",
            "FL": r"1234.57 \pm 0.10",
            "fL": r"1234.57 \pm 0.10",  # noqa
            "FL": r"1234.57 \pm 0.10",  # noqa
            "%L": r"\left(123457 \pm 10\right) \%",
        },
        #
        # ... combined with the spectroscopic notation:
        (-1.23, 3.4): {  # noqa
            "SL": "-1.2(3.4)",
            "LS": "-1.2(3.4)",
            ".2ufSL": "-1.2(3.4)",
            ".2ufLS": "-1.2(3.4)",
        },
        # Special cases for the uncertainty (0, nan) and format
        # strings (extension S, L, U,..., global width, etc.).
        #
        # Python 3.2 and 3.3 give 1.4e-12*1e+12 = 1.4000000000000001
        # instead of 1.4 for Python 3.1. The problem does not appear
        # with 1.2, so 1.2 is used.
        (-1.2e-12, 0): {
            "12.2gPL": "  -1.2×10⁻¹²±           0",
            # Pure "width" formats are not accepted by the % operator,
            # and only %-compatible formats are accepted, for Python <
            # 2.6:
            "13S": "  -1.2(0)e-12",
            "10P": "-1.2×10⁻¹²±         0",
            "L": r"\left(-1.2 \pm 0\right) \times 10^{-12}",
            # No factored exponent, LaTeX
            "1L": r"-1.2 \times 10^{-12} \pm 0",
            "SL": r"-1.2(0) \times 10^{-12}",
            "SP": "-1.2(0)×10⁻¹²",
        },
        # Python 3.2 and 3.3 give 1.4e-12*1e+12 = 1.4000000000000001
        # instead of 1.4 for Python 3.1. The problem does not appear
        # with 1.2, so 1.2 is used.
        (-1.2e-12, float("nan")): {
            ".2uG": "(-1.2+/-%s)E-12" % NaN_EFG,  # u ignored, format used
            "15GS": "  -1.2(%s)E-12" % NaN_EFG,
            "SL": r"-1.2(\mathrm{nan}) \times 10^{-12}",  # LaTeX NaN
            # Pretty-print priority, but not for NaN:
            "PSL": r"-1.2(\mathrm{nan})×10⁻¹²",
            "L": r"\left(-1.2 \pm \mathrm{nan}\right) \times 10^{-12}",
            # Uppercase NaN and LaTeX:
            ".1EL": (r"\left(-1.2 \pm \mathrm{%s}\right) \times 10^{-12}" % NaN_EFG),
            "10": "  -1.2e-12+/-       nan",
            "15S": "  -1.2(nan)e-12",
        },
        (3.14e-10, 0.01e-10): {
            # Character (Unicode) strings:
            "P": "(3.140±0.010)×10⁻¹⁰",  # PDG rules: 2 digits
            "PL": "(3.140±0.010)×10⁻¹⁰",  # Pretty-print has higher priority
            # Truncated non-zero uncertainty:
            ".1e": "(3.1+/-0.0)e-10",
            ".1eS": "3.1(0.0)e-10",
        },
        # Some special cases:
        (1, float("nan")): {
            "g": "1+/-nan",
            "G": "1+/-%s" % NaN_EFG,
            "%": "(100.000000+/-nan)%",  # The % format type is like f
            # Should be the same as '+05', for floats, but is not, in
            # Python 2.7:
            "+05g": "+0001+/-00nan",
            # 5 is the *minimal* width, 6 is the default number of
            # digits after the decimal point:
            "+05%": "(+100.000000+/-00nan)%",
            # There is a difference between '{}'.format(1.) and
            # '{:g}'.format(1.), which is not fully obvious in the
            # documentation, which indicates that a None format type
            # is like g. The reason is that the empty format string is
            # actually interpreted as str(), and that str() does not
            # have to behave like g ('{}'.format(1.234567890123456789)
            # and '{:g}'.format(1.234567890123456789) are different).
            "": "1.0+/-nan",
            # This is ugly, but consistent with
            # '{:+05}'.format(float('nan')) and format(1.) (which
            # differs from format(1)!):
            "+05": "+01.0+/-00nan",
        },
        (9.9, 0.1): {".1ue": "(9.9+/-0.1)e+00", ".0fS": "10(0.)"},
        (9.99, 0.1): {
            # The precision has an effect on the exponent, like for
            # floats:
            ".2ue": "(9.99+/-0.10)e+00",  # Same exponent as for 9.99 alone
            ".1ue": "(1.00+/-0.01)e+01",  # Same exponent as for 9.99 alone
        },
        # 0 uncertainty: nominal value displayed like a float:
        (1.2345, 0): {
            ".2ue": "(1.23+/-0)e+00",
            "1.2ue": "1.23e+00+/-0",  # No factored exponent
            ".2uf": "1.23+/-0",
            ".2ufS": "1.23(0)",
            ".2fS": "1.23(0)",
            "g": "1.2345+/-0",
            "": "1.2345+/-0",
        },
        # Alignment and filling characters (supported in Python 2.6+):
        (3.1415e10, 0): {
            "<15": "31415000000.0  +/-0              ",
            "<20S": "31415000000.0(0)    ",
            # Trying to trip the format parsing with a fill character
            # which is an alignment character:
            "=>15": "==31415000000.0+/-==============0",
        },
        (1234.56789, 0): {
            "1.2ue": "1.23e+03+/-0",  # u ignored
            "1.2e": "1.23e+03+/-0",
            # Default precision = 6
            "eL": r"\left(1.234568 \pm 0\right) \times 10^{3}",
            "EL": r"\left(1.234568 \pm 0\right) \times 10^{3}",
            "fL": r"1234.567890 \pm 0",
            "FL": r"1234.567890 \pm 0",
            "%L": r"\left(123456.789000 \pm 0\right) \%",
        },
        (1e5, 0): {"g": "100000+/-0"},
        (1e6, 0): {
            # A default precision of 6 is used because the uncertainty
            # cannot be used for defining a default precision (it does
            # not have a magnitude):
            "g": "(1+/-0)e+06"
        },
        (1e6 + 10, 0): {
            # A default precision of 6 is used because the uncertainty
            # cannot be used for defining a default precision (it does
            # not have a magnitude):
            "g": "(1.00001+/-0)e+06"
        },
        # Rounding of the uncertainty that "changes" the number of
        # significant digits:
        (1, 0.994): {
            ".3uf": "1.000+/-0.994",
            ".2uf": "1.00+/-0.99",
            ".1uf": "1+/-1",  # Discontinuity in the number of digits
        },
        (12.3, 2.3): {
            ".2ufS": "12.3(2.3)"  # Decimal point on the uncertainty
        },
        (12.3, 2.3): {  # noqa
            ".1ufS": "12(2)"  # No decimal point on the uncertainty
        },
        (0, 0): {  # Make defining the first significant digit problematic
            ".1f": "0.0+/-0",  # Simple float formatting
            "g": "0+/-0",
        },
        (1.2e-34, 5e-67): {
            ".6g": "(1.20000+/-0.00000)e-34",
            "13.6g": "  1.20000e-34+/-  0.00000e-34",
            "13.6G": "  1.20000E-34+/-  0.00000E-34",
            ".6GL": r"\left(1.20000 \pm 0.00000\right) \times 10^{-34}",
            ".6GLp": r"\left(1.20000 \pm 0.00000\right) \times 10^{-34}",
        },
        (float("nan"), 100): {  # NaN *nominal value*
            "": "nan+/-100.0",  # Like '{}'.format(100.)
            "g": "nan+/-100",  # Like '{:g}'.format(100.)
            ".1e": "(nan+/-1.0)e+02",  # Similar to 1±nan
            ".1E": "(%s+/-1.0)E+02" % NaN_EFG,
            ".1ue": "(nan+/-1)e+02",
            "10.1e": "       nan+/-   1.0e+02",
        },
        (float("nan"), 1e8): {  # NaN *nominal value*
            "": "nan+/-100000000.0",  # Like '{}'.format(1e8)
            "g": "(nan+/-1)e+08",  # Like '{:g}'.format(1e8)
            ".1e": "(nan+/-1.0)e+08",
            ".1E": "(%s+/-1.0)E+08" % NaN_EFG,
            ".1ue": "(nan+/-1)e+08",
            "10.1e": "       nan+/-   1.0e+08",  # 'nane+08' would be strange
        },
        (float("nan"), 123456789): {  # NaN *nominal value*
            "": "nan+/-123456789.0",  # Similar to '{}'.format(123456789.)
            "g": "(nan+/-1.23457)e+08",  # Similar to '{:g}'.format(123456789.)
            ".1e": "(nan+/-1.2)e+08",
            ".1E": "(%s+/-1.2)E+08" % NaN_EFG,
            ".1ue": "(nan+/-1)e+08",
            ".1ueL": r"\left(\mathrm{nan} \pm 1\right) \times 10^{8}",
            "10.1e": "       nan+/-   1.2e+08",
            "10.1eL": r"\mathrm{nan} \pm 1.2 \times 10^{8}",
        },
        (float("nan"), float("nan")): {  # *Double* NaN
            "": "nan+/-nan",
            ".1e": "nan+/-nan",
            ".1E": "%s+/-%s" % (NaN_EFG, NaN_EFG),
            ".1ue": "nan+/-nan",
            "EL": r"\mathrm{%s} \pm \mathrm{%s}" % (NaN_EFG, NaN_EFG),
        },
        (float("inf"), 100): {  # Inf *nominal value*
            "": "inf+/-100.0",  # Like '{}'.format(100.)
            "g": "inf+/-100",  # Like '{:g}'.format(100.)
            ".1e": "(inf+/-1.0)e+02",  # Similar to 1±inf
            ".1E": "(%s+/-1.0)E+02" % Inf_EFG,
            ".1ue": "(inf+/-1)e+02",
            "10.1e": "       inf+/-   1.0e+02",
        },
        (float("inf"), 1e8): {  # Inf *nominal value*
            "": "inf+/-100000000.0",  # Like '{}'.format(1e8)
            "g": "(inf+/-1)e+08",  # Like '{:g}'.format(1e8)
            ".1e": "(inf+/-1.0)e+08",
            ".1E": "(%s+/-1.0)E+08" % Inf_EFG,
            ".1ue": "(inf+/-1)e+08",
            "10.1e": "       inf+/-   1.0e+08",  # 'infe+08' would be strange
        },
        (float("inf"), 123456789): {  # Inf *nominal value*
            "": "inf+/-123456789.0",  # Similar to '{}'.format(123456789.)
            "g": "(inf+/-1.23457)e+08",  # Similar to '{:g}'.format(123456789.)
            ".1e": "(inf+/-1.2)e+08",
            ".1ep": "(inf+/-1.2)e+08",
            ".1E": "(%s+/-1.2)E+08" % Inf_EFG,
            ".1ue": "(inf+/-1)e+08",
            ".1ueL": r"\left(\infty \pm 1\right) \times 10^{8}",
            ".1ueLp": r"\left(\infty \pm 1\right) \times 10^{8}",
            "10.1e": "       inf+/-   1.2e+08",
            "10.1eL": r"    \infty \pm 1.2 \times 10^{8}",
        },
        (float("inf"), float("inf")): {  # *Double* Inf
            "": "inf+/-inf",
            ".1e": "inf+/-inf",
            ".1E": "%s+/-%s" % (Inf_EFG, Inf_EFG),
            ".1ue": "inf+/-inf",
            "EL": r"\infty \pm \infty",
            "ELp": r"\left(\infty \pm \infty\right)",
        },
        # Like the tests for +infinity, but for -infinity:
        (float("-inf"), 100): {  # Inf *nominal value*
            "": "-inf+/-100.0",  # Like '{}'.format(100.)
            "g": "-inf+/-100",  # Like '{:g}'.format(100.)
            ".1e": "(-inf+/-1.0)e+02",  # Similar to 1±inf
            ".1E": "(-%s+/-1.0)E+02" % Inf_EFG,
            ".1ue": "(-inf+/-1)e+02",
            "10.1e": "      -inf+/-   1.0e+02",
        },
        (float("-inf"), 1e8): {  # Inf *nominal value*
            "": "-inf+/-100000000.0",  # Like '{}'.format(1e8)
            "g": "(-inf+/-1)e+08",  # Like '{:g}'.format(1e8)
            ".1e": "(-inf+/-1.0)e+08",
            ".1E": "(-%s+/-1.0)E+08" % Inf_EFG,
            ".1ue": "(-inf+/-1)e+08",
            "10.1e": "      -inf+/-   1.0e+08",  # 'infe+08' would be strange
        },
        (float("-inf"), 123456789): {  # Inf *nominal value*
            "": "-inf+/-123456789.0",  # Similar to '{}'.format(123456789.)
            "g": "(-inf+/-1.23457)e+08",  # Similar to '{:g}'.format(123456789.)
            ".1e": "(-inf+/-1.2)e+08",
            ".1E": "(-%s+/-1.2)E+08" % Inf_EFG,
            ".1ue": "(-inf+/-1)e+08",
            ".1ueL": r"\left(-\infty \pm 1\right) \times 10^{8}",
            "10.1e": "      -inf+/-   1.2e+08",
            "10.1eL": r"   -\infty \pm 1.2 \times 10^{8}",
        },
        (float("-inf"), float("inf")): {  # *Double* Inf
            "": "-inf+/-inf",
            ".1e": "-inf+/-inf",
            ".1E": "-%s+/-%s" % (Inf_EFG, Inf_EFG),
            ".1ue": "-inf+/-inf",
            "EL": r"-\infty \pm \infty",
        },
        # The Particle Data Group convention trumps the "at least one
        # digit past the decimal point" for Python floats, but only
        # with a non-zero uncertainty:
        (724.2, 26.4): {"": "724+/-26", "p": "(724+/-26)"},
        (724, 0): {"": "724.0+/-0"},
        # More NaN and infinity, in particular with LaTeX and various
        # options:
        (float("-inf"), float("inf")): {  # noqa
            "S": "-inf(inf)",
            "LS": r"-\infty(\infty)",
            "L": r"-\infty \pm \infty",
            "LP": r"-\infty±\infty",
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
            "020S": format(float("-inf"), "015") + "(inf)",
        },
        (-float("nan"), float("inf")): {
            "S": "nan(inf)",
            "LS": r"\mathrm{nan}(\infty)",
            "L": r"\mathrm{nan} \pm \infty",
            "LP": r"\mathrm{nan}±\infty",
        },
        # Leading zeroes in the shorthand notation:
        (-2, 3): {"020S": "-000000000002.0(3.0)"},
    }

    # ',' format option: introduced in Python 2.7
    if sys.version_info >= (2, 7):
        tests.update(
            {
                (1234.56789, 0.012): {",.1uf": "1,234.57+/-0.01"},
                (123456.789123, 1234.5678): {
                    ",f": "123,457+/-1,235",  # Particle Data Group convention
                    ",.4f": "123,456.7891+/-1,234.5678",
                },
            }
        )

    # True if we can detect that the Jython interpreter is running this code:
    try:
        jython_detected = sys.subversion[0] == "Jython"
    except AttributeError:
        jython_detected = False

    for values, representations in tests.items():
        value = ufloat(*values)

        for format_spec, result in representations.items():
            # print "FORMATTING {} WITH '{}'".format(repr(value), format_spec)

            # Jython 2.5.2 does not always represent NaN as nan or NAN
            # in the CPython way: for example, '%.2g' % float('nan')
            # is '\ufffd'. The test is skipped, in this case:
            if jython_detected and (isnan(value.std_dev) or isnan(value.nominal_value)):
                continue

            # Call that works with Python < 2.6 too:
            representation = value.format(format_spec)

            assert representation == result, (
                # The representation is used, for terminal that do not
                # support some characters like ±, and superscripts:
                "Incorrect representation %r for format %r of %r:"
                " %r expected." % (representation, format_spec, value, result)
            )

            # An empty format string is like calling str()
            # (http://docs.python.org/2/library/string.html#formatspec):
            if not format_spec:
                assert representation == str(value), (
                    "Empty format should give the same thing as str():"
                    " %s obtained instead of %s" % (representation, str(value))
                )

            # Parsing back into a number with uncertainty (unless the
            # LaTeX or comma notation is used):
            if (
                not set(format_spec).intersection("L,*%")  # * = fill with *
                # "0nan"
                and "0nan" not in representation.lower()
                # "0inf"
                and "0inf" not in representation.lower()
                # Specific case:
                and "=====" not in representation
            ):
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
                        assert numbers_close(
                            value.nominal_value, value_back.nominal_value, 2.4e-1
                        )

                    # If the uncertainty is zero, then the relative
                    # change can be large:
                    assert numbers_close(value.std_dev, value_back.std_dev, 3e-1)

                except AssertionError:
                    # !! The following string formatting requires
                    # str() to work (to not raise an exception) on the
                    # values (which have a non-standard class):
                    raise AssertionError(
                        "Original value %s and value %s parsed from %r"
                        " (obtained through format specification %r)"
                        " are not close enough"
                        % (value, value_back, representation, format_spec)
                    )


def test_unicode_format():
    """Test of the unicode formatting of numbers with uncertainties"""

    x = ufloat(3.14159265358979, 0.25)

    assert isinstance("Résultat = %s" % x.format(""), str)
    assert isinstance("Résultat = %s" % x.format("P"), str)


def test_custom_pretty_print_and_latex():
    """Test of the pretty-print and LaTeX format customizations"""

    x = ufloat(2, 0.1) * 1e-11

    # We will later restore the defaults:
    PREV_CUSTOMIZATIONS = {
        var: getattr(formatting, var).copy()
        for var in ["PM_SYMBOLS", "MULT_SYMBOLS", "GROUP_SYMBOLS"]
    }

    # Customizations:
    for format in ["pretty-print", "latex"]:
        formatting.PM_SYMBOLS[format] = " ± "
        formatting.MULT_SYMBOLS[format] = "⋅"
        formatting.GROUP_SYMBOLS[format] = ("[", "]")

    assert "{:P}".format(x) == "[2.00 ± 0.10]⋅10⁻¹¹"
    assert "{:L}".format(x) == "[2.00 ± 0.10] ⋅ 10^{-11}"

    # We restore the defaults:
    for var, setting in PREV_CUSTOMIZATIONS.items():
        setattr(formatting, var, setting)
