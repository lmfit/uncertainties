import re
from uncertainties.formatting import nrmlze_superscript

###############################################################################
# Parsing of values with uncertainties:

# Parsing of (part of) numbers. The reason why the decimal part is
# parsed (if any), instead of using the parsing built in float(), is
# that the presence (or not) of a decimal point does matter, in the
# semantics of some representations (e.g. .1(2.) = .1+/-2, whereas
# .1(2) = .1+/-0.2), so just getting the numerical value of the part
# in parentheses would not be sufficient.
POSITIVE_DECIMAL_UNSIGNED_OR_NON_FINITE = r"((\d*)(\.\d*)?|nan|NAN|inf|INF)"

# Regexp for a number with uncertainty (e.g., "-1.234(2)e-6"), where
# the uncertainty is optional (in which case the uncertainty is
# implicit). The uncertainty can also be nan or NAN:
#
# !! WARNING: in Python 2, the code relies on "… % <unicode string>" returning
# a Unicode string (even if the template is not Unicode):
NUMBER_WITH_UNCERT_RE_STR = """
    ([+-])?  # Sign
    %s  # Main number
    (?:\\(%s\\))?  # Optional uncertainty
    (?:
        (?:[eE]|\\s*×\\s*10)
        (.*)
    )?  # Optional exponent
    """ % (
    POSITIVE_DECIMAL_UNSIGNED_OR_NON_FINITE,
    POSITIVE_DECIMAL_UNSIGNED_OR_NON_FINITE,
)

NUMBER_WITH_UNCERT_RE_MATCH = re.compile(
    "%s$" % NUMBER_WITH_UNCERT_RE_STR, re.VERBOSE
).match

# Number with uncertainty with a factored exponent (e.g., of the form
# (... +/- ...)e10): this is a loose matching, so as to accommodate
# for multiple formats:
NUMBER_WITH_UNCERT_GLOBAL_EXP_RE_MATCH = re.compile(
    """
    \\(
    (?P<simple_num_with_uncert>.*)
    \\)
    (?:[eE]|\\s*×\\s*10) (?P<exp_value>.*)
    $""",
    re.VERBOSE,
).match


class NotParenUncert(ValueError):
    """
    Raised when a string representing an exact number or a number with
    an uncertainty indicated between parentheses was expected but not
    found.
    """


def parse_error_in_parentheses(representation):
    # !!!! The code seems to handle superscript exponents, but the
    # docstring doesn't reflect this!?
    """
    Return (value, error) from a string representing a number with
    uncertainty like 12.34(5), 12.34(142), 12.5(3.4), 12.3(4.2)e3, or
    13.4(nan)e10.  If no parenthesis is given, an uncertainty of one
    on the last digit is assumed.

    The digits between parentheses correspond to the same number of digits
    at the end of the nominal value (the decimal point in the uncertainty
    is optional). Example: 12.34(142) = 12.34±1.42.

    Raises ValueError if the string cannot be parsed.
    """

    match = NUMBER_WITH_UNCERT_RE_MATCH(representation)

    if match:
        # The 'main' part is the nominal value, with 'int'eger part, and
        # 'dec'imal part.  The 'uncert'ainty is similarly broken into its
        # integer and decimal parts.
        (sign, main, _, main_dec, uncert, uncert_int, uncert_dec, exponent) = (
            match.groups()
        )
    else:
        raise NotParenUncert(
            "Unparsable number representation: '%s'."
            " See the documentation of ufloat_fromstr()." % representation
        )

    # Global exponent:
    if exponent:
        factor = 10.0 ** nrmlze_superscript(exponent)
    else:
        factor = 1

    # Nominal value:
    value = float((sign or "") + main) * factor

    if uncert is None:
        # No uncertainty was found: an uncertainty of 1 on the last
        # digit is assumed:
        uncert_int = "1"  # The other parts of the uncertainty are None

    # Do we have a fully explicit uncertainty?
    if uncert_dec is not None or uncert in {"nan", "NAN", "inf", "INF"}:
        uncert_value = float(uncert)
    else:
        # uncert_int represents an uncertainty on the last digits:

        # The number of digits after the period defines the power of
        # 10 that must be applied to the provided uncertainty:
        if main_dec is None:
            num_digits_after_period = 0
        else:
            num_digits_after_period = len(main_dec) - 1

        uncert_value = int(uncert_int) / 10.0**num_digits_after_period

    # We apply the exponent to the uncertainty as well:
    uncert_value *= factor

    return (value, uncert_value)


# Regexp for catching the two variable parts of -1.2×10⁻¹²:
PRETTY_PRINT_MATCH = re.compile("(.*?)\\s*×\\s*10(.*)").match


def to_float(value_str):
    """
    Converts a string representing a float to a float.

    The usual valid Python float() representations are correctly
    parsed.

    In addition, the pretty-print notation -1.2×10⁻¹² is also
    converted.

    ValueError is raised if no float can be obtained.
    """

    try:
        return float(value_str)
    except ValueError:
        pass

    # The pretty-print notation is tried:
    match = PRETTY_PRINT_MATCH(value_str)
    if match:
        try:
            return float(match.group(1)) * 10.0 ** nrmlze_superscript(match.group(2))
        except ValueError:
            raise ValueError(
                "Mantissa or exponent incorrect in pretty-print" " form %s" % value_str
            )
    else:
        raise ValueError(
            "No valid Python float or pretty-print form" " recognized in %s" % value_str
        )


cannot_parse_ufloat_msg_pat = (
    "Cannot parse %s: see the documentation for ufloat_fromstr() for a"
    " list of accepted formats"
)


# The following function is not exposed because it can in effect be
# obtained by doing x = ufloat_fromstr(representation) and reading
# x.nominal_value and x.std_dev:
def str_to_number_with_uncert(representation):
    """
    Given a string that represents a number with uncertainty, returns the
    nominal value and the uncertainty.

    See the documentation for ufloat_fromstr() for a list of accepted
    formats.

    When no numerical error is given, an uncertainty of 1 on the last
    digit is implied.

    Raises ValueError if the string cannot be parsed.

    representation -- string with no leading or trailing spaces.
    """

    # The "p" format can add parentheses around the whole printed result: we
    # remove them:
    if representation.startswith("(") and representation.endswith(")"):
        representation = representation[1:-1]

    match = NUMBER_WITH_UNCERT_GLOBAL_EXP_RE_MATCH(representation)

    # The representation is simplified, but the global factor is
    # calculated:

    if match:
        # We have a form with a factored exponent: (1.23 +/- 0.01)e10,
        # etc.

        exp_value_str = match.group("exp_value")

        try:
            exponent = nrmlze_superscript(exp_value_str)
        except ValueError:
            raise ValueError(cannot_parse_ufloat_msg_pat % representation)

        factor = 10.0**exponent

        representation = match.group("simple_num_with_uncert")
    else:
        factor = 1  # No global exponential factor

    match = re.match("(.*)(?:\\+/-|±)(.*)", representation)
    if match:
        (nom_value, uncert) = match.groups()

        try:
            # Simple form 1234.45+/-1.2 or 1234.45±1.2, or 1.23e-10+/-1e-23
            # or -1.2×10⁻¹²±1e23:
            parsed_value = (to_float(nom_value) * factor, to_float(uncert) * factor)
        except ValueError:
            raise ValueError(cannot_parse_ufloat_msg_pat % representation)

    else:
        # Form with error parentheses or no uncertainty:
        try:
            parsed_value = parse_error_in_parentheses(representation)
        except NotParenUncert:
            raise ValueError(cannot_parse_ufloat_msg_pat % representation)

    return parsed_value
