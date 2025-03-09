from math import isinf, isnan, isfinite
import math
import re


def first_digit(value):
    """
    Return the first digit position of the given value, as an integer.

    0 is the digit just before the decimal point. Digits to the right
    of the decimal point have a negative position.

    Return 0 for a null value.
    """
    try:
        return int(math.floor(math.log10(abs(value))))
    except ValueError:  # Case of value == 0
        return 0


def PDG_precision(std_dev):
    """
    Return the number of significant digits to be used for the given
    standard deviation, according to the rounding rules of the
    Particle Data Group (2010)
    (http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf).

    Also returns the effective standard deviation to be used for
    display.
    """

    exponent = first_digit(std_dev)

    # The first three digits are what matters: we get them as an
    # integer number in [100; 999).
    #
    # In order to prevent underflow or overflow when calculating
    # 10**exponent, the exponent is slightly modified first and a
    # factor to be applied after "removing" the new exponent is
    # defined.
    #
    # Furthermore, 10**(-exponent) is not used because the exponent
    # range for very small and very big floats is generally different.
    if exponent >= 0:
        # The -2 here means "take two additional digits":
        (exponent, factor) = (exponent - 2, 1)
    else:
        (exponent, factor) = (exponent + 1, 1000)
    digits = int(std_dev / 10.0**exponent * factor)  # int rounds towards zero

    # Rules:
    if digits <= 354:
        return (2, std_dev)
    elif digits <= 949:
        return (1, std_dev)
    else:
        # The parentheses matter, for very small or very large
        # std_dev:
        return (2, 10.0**exponent * (1000 / factor))


# Definition of a basic (format specification only, no full-feature
# format string) formatting function that works whatever the version
# of Python. This function exists so that the more capable format() is
# used instead of the % formatting operator, if available:
robust_format = format

# Exponent letter: the keys are the possible main_fmt_type values of
# format_num():
EXP_LETTERS = {"f": "e", "F": "E"}


def robust_align(orig_str, fill_char, align_option, width):
    """
    Aligns the given string with the given fill character.

    orig_str -- string to be aligned (str or unicode object).

    fill_char -- if empty, space is used.

    align_option -- as accepted by format().

    wdith -- string that contains the width.
    """

    # print "ALIGNING", repr(orig_str), "WITH", fill_char+align_option,
    # print "WIDTH", width

    return format(orig_str, fill_char + align_option + width)


# Maps some Unicode code points ("-", "+", and digits) to their
# superscript version:
TO_SUPERSCRIPT = {
    0x2B: "⁺",
    0x2D: "⁻",
    0x30: "⁰",
    0x31: "¹",
    0x32: "²",
    0x33: "³",
    0x34: "⁴",
    0x35: "⁵",
    0x36: "⁶",
    0x37: "⁷",
    0x38: "⁸",
    0x39: "⁹",
}

# Inverted TO_SUPERSCRIPT table, for use with unicode.translate():
#
#! Python 2.7+ can use a dictionary comprehension instead:
FROM_SUPERSCRIPT = {ord(sup): normal for (normal, sup) in TO_SUPERSCRIPT.items()}


def to_superscript(value):
    """
    Return a (Unicode) string with the given value as superscript characters.

    The value is formatted with the %d %-operator format.

    value -- integer.
    """

    return ("%d" % value).translate(TO_SUPERSCRIPT)


def nrmlze_superscript(number_str):
    """
    Return a string with superscript digits transformed into regular digits.

    Non-superscript digits are not changed before the conversion. Thus, the
    string can also contain regular digits.

    ValueError is raised if the conversion cannot be done.

    number_str -- string to be converted (of type str, but also possibly, for
    Python 2, unicode, which allows this string to contain superscript digits).
    """
    # !! Python 3 doesn't need this str(), which is only here for giving the
    # .translate() method to str objects in Python 2 (this str() comes
    # from the builtins module of the future package and is therefore
    # a subclass of unicode, in Python 2):
    return int(str(number_str).translate(FROM_SUPERSCRIPT))


PM_SYMBOLS = {"pretty-print": "±", "latex": r" \pm ", "default": "+/-"}

# Multiplication symbol for pretty printing (so that pretty printing can
# be customized):
MULT_SYMBOLS = {"pretty-print": "×", "latex": r"\times"}

# Function that transforms a numerical exponent produced by format_num() into
# the corresponding string notation (for non-default modes):
EXP_PRINT = {
    "pretty-print": lambda common_exp: "%s10%s"
    % (MULT_SYMBOLS["pretty-print"], to_superscript(common_exp)),
    "latex": lambda common_exp: r" %s 10^{%d}" % (MULT_SYMBOLS["latex"], common_exp),
}

# Symbols used for grouping (typically between parentheses) in format_num():
GROUP_SYMBOLS = {
    "pretty-print": ("(", ")"),
    # Because of possibly exponents inside the parentheses (case of a
    # specified field width), it is better to use auto-adjusting
    # parentheses. This has the side effect of making the part between
    # the parentheses non-breakable (the text inside parentheses in a
    # LaTeX math expression $...$ can be broken).
    "latex": (r"\left(", r"\right)"),
    "default": ("(", ")"),  # Basic text mode
}


def format_num(
    nom_val_main, error_main, common_exp, fmt_parts, prec, main_pres_type, options
):
    """
    Return a formatted number with uncertainty.

    Null errors (error_main) are displayed as the integer 0, with
    no decimal point.

    The formatting can be customized globally through the PM_SYMBOLS,
    MULT_SYMBOLS, GROUP_SYMBOLS and EXP_PRINT dictionaries, which contain
    respectively the symbol for ±, for multiplication, for parentheses, and a
    function that maps an exponent to something like "×10²" (using
    MULT_SYMBOLS).

    Each of these dictionary has (at least) a 'pretty-print' and a 'latex' key,
    that define the symbols to be used for these two output formats (the
    PM_SYMBOLS and GROUP_SYMBOLS also have a 'default' key for the default
    output format). For example, the defaults for the 'pretty-print' format
    are:

    - PM_SYMBOLS['pretty-print'] = '±'
    - MULT_SYMBOLS['pretty-print'] = '×'
    - GROUP_SYMBOLS['pretty-print'] = ( '(', ')' )
    - EXP_PRINT['pretty-print']: see the source code.

    Arguments:

    nom_val_main, error_main -- nominal value and error, before using
    common_exp (e.g., "1.23e2" would have a main value of 1.23;
    similarly, "12.3+/-0.01" would have a main value of 12.3).

    common_exp -- common exponent to use. If None, no common exponent
    is used.

    fmt_parts -- mapping that contains at least the following parts of
    the format specification: fill, align, sign, zero, width, comma,
    type; the value are strings. These format specification parts are
    handled. The width is applied to each value, or, if the shorthand
    notation is used, globally. If the error is special (zero, NaN, inf),
    the parts are applied as much as possible to the nominal value.

    prec -- precision to use with the main_pres_type format type
    (see below).

    main_pres_type -- format presentation type, either "f" or
    "F". This defines how the mantissas, exponents and NaN/inf values are
    represented (in the same way as for float). None, the empty
    string, or "%" are not accepted.

    options -- options (as an object that support membership testing, like for
    instance a string). "P" is for pretty-printing ("±" between the nominal
    value and the error, superscript exponents, etc.). "L" is for a LaTeX
    output. "S" is for the shorthand notation 1.23(1). "p" is for making sure
    that the …±… part is surrounded by parentheses.  "%" adds a final percent
    sign, and parentheses if the shorthand notation is not used. Options can
    be combined. The P option has priority over the L option (if both are
    given). For details, see the documentation for
    AffineScalarFunction.__format__().
    """

    # print (nom_val_main, error_main, common_exp,
    #        fmt_parts, prec, main_pres_type, options)

    # If a decimal point were always present in zero rounded errors
    # that are not zero, the formatting would be difficult, in general
    # (because the formatting options are very general): an example
    # is'{:04.0f}'.format(0.1), which gives "0000" and would have to
    # give "000.". Another example is '{:<4.0f}'.format(0.1), which
    # gives "0 " but should give "0.  ". This is cumbersome to
    # implement in the general case, because no format prints "0."
    # for 0. Furthermore, using the .0f format already brings the same
    # kind of difficulty: non-zero numbers can appear as the exact
    # integer zero, after rounding. The problem is not larger, for
    # numbers with an error.
    #
    # That said, it is good to indicate null errors explicitly when
    # possible: printing 3.1±0 with the default format prints 3.1+/-0,
    # which shows that the uncertainty is exactly zero.

    # The suffix of the result is calculated first because it is
    # useful for the width handling of the shorthand notation.

    # Printing type for parts of the result (exponent, parentheses),
    # taking into account the priority of the pretty-print mode over
    # the LaTeX mode. This setting does not apply to everything: for
    # example, NaN is formatted as \mathrm{nan} (or NAN) if the LaTeX
    # mode is required.
    if "P" in options:
        print_type = "pretty-print"
    elif "L" in options:
        print_type = "latex"
    else:
        print_type = "default"

    # Exponent part:
    if common_exp is None:
        exp_str = ""
    elif print_type == "default":
        # Case of e or E. The same convention as Python 2.7
        # to 3.3 is used for the display of the exponent:
        exp_str = EXP_LETTERS[main_pres_type] + "%+03d" % common_exp
    else:
        exp_str = EXP_PRINT[print_type](common_exp)

    # Possible % sign:
    percent_str = ""
    if "%" in options:
        if "L" in options:
            # % is a special character, in LaTeX: it must be escaped.
            #
            # Using '\\' in the code instead of r'\' so as not to
            # confuse emacs's syntax highlighting:
            percent_str += " \\"
        percent_str += "%"

    ####################

    # Only true if the error should not have an exponent (has priority
    # over common_exp):
    special_error = not error_main or not isfinite(error_main)

    # Nicer representation of the main nominal part, with no trailing
    # zeros, when the error does not have a defined number of
    # significant digits:
    if special_error and fmt_parts["type"] in ("", "g", "G"):
        # The main part is between 1 and 10 because any possible
        # exponent is taken care of by common_exp, so it is
        # formatted without an exponent (otherwise, the exponent
        # would have to be handled for the LaTeX option):
        fmt_suffix_n = (fmt_parts["prec"] or "") + fmt_parts["type"]
    else:
        fmt_suffix_n = ".%d%s" % (prec, main_pres_type)

    # print "FMT_SUFFIX_N", fmt_suffix_n

    ####################

    # Calculation of the mostly final numerical part value_str (no %
    # sign, no global width applied).

    # Error formatting:

    if "S" in options:  # Shorthand notation:
        # Calculation of the uncertainty part, uncert_str:

        if error_main == 0:
            # The error is exactly zero
            uncert_str = "0"
        elif isnan(error_main):
            uncert_str = robust_format(error_main, main_pres_type)
            if "L" in options:
                uncert_str = r"\mathrm{%s}" % uncert_str
        elif isinf(error_main):
            if "L" in options:
                uncert_str = r"\infty"
            else:
                uncert_str = robust_format(error_main, main_pres_type)
        else:  #  Error with a meaningful first digit (not 0, and real number)
            uncert = round(error_main, prec)

            # The representation uncert_str of the uncertainty (which will
            # be put inside parentheses) is calculated:

            # The uncertainty might straddle the decimal point: we
            # keep it as it is, in this case (e.g. 1.2(3.4), as this
            # makes the result easier to read); the shorthand
            # notation then essentially coincides with the +/-
            # notation:
            if first_digit(uncert) >= 0 and prec > 0:
                # This case includes a zero rounded error with digits
                # after the decimal point:
                uncert_str = "%.*f" % (prec, uncert)

            else:
                if uncert:
                    # The round is important because 566.99999999 can
                    # first be obtained when 567 is wanted (%d prints the
                    # integer part, not the rounded value):
                    uncert_str = "%d" % round(uncert * 10.0**prec)
                else:
                    # The decimal point indicates a truncated float
                    # (this is easy to do, in this case, since
                    # fmt_prefix_e is ignored):
                    uncert_str = "0."

        # End of the final number representation (width and alignment
        # not included). This string is important for the handling of
        # the width:
        value_end = "(%s)%s%s" % (uncert_str, exp_str, percent_str)
        any_exp_factored = True  # Single exponent in the output

        ##########
        # Nominal value formatting:

        # Calculation of fmt_prefix_n (prefix for the format of the
        # main part of the nominal value):

        if fmt_parts["zero"] and fmt_parts["width"]:
            # Padding with zeros must be done on the nominal value alone:

            # Remaining width (for the nominal value):
            nom_val_width = max(int(fmt_parts["width"]) - len(value_end), 0)
            fmt_prefix_n = "%s%s%d%s" % (
                fmt_parts["sign"],
                fmt_parts["zero"],
                nom_val_width,
                fmt_parts["comma"],
            )

        else:
            # Any 'zero' part should not do anything: it is not
            # included
            fmt_prefix_n = fmt_parts["sign"] + fmt_parts["comma"]

        # print "FMT_PREFIX_N", fmt_prefix_n
        # print "FMT_SUFFIX_N", fmt_suffix_n

        nom_val_str = robust_format(nom_val_main, fmt_prefix_n + fmt_suffix_n)

        ##########
        # Overriding of nom_val_str for LaTeX,; possibly based on the
        # existing value (for NaN vs nan):
        if "L" in options:
            if isnan(nom_val_main):
                nom_val_str = r"\mathrm{%s}" % nom_val_str
            elif isinf(nom_val_main):
                # !! It is wasteful, in this case, to replace
                # nom_val_str: could this be avoided while avoiding to
                # duplicate the formula for nom_val_str for the common
                # case (robust_format(...))?
                nom_val_str = r"%s\infty" % ("-" if nom_val_main < 0 else "")

        value_str = nom_val_str + value_end

        # Global width, if any:

        if fmt_parts["width"]:  # An individual alignment is needed:
            # Default alignment, for numbers: to the right (if no
            # alignment is specified, a string is aligned to the
            # left):
            value_str = robust_align(
                value_str,
                fmt_parts["fill"],
                fmt_parts["align"] or ">",
                fmt_parts["width"],
            )

    else:  # +/- notation:
        # The common exponent is factored or not, depending on the
        # width. This gives nice columns for the nominal values and
        # the errors (no shift due to a varying exponent), when a need
        # is given:
        any_exp_factored = not fmt_parts["width"]

        # True when the error part has any exponent directly attached
        # (case of an individual exponent for both the nominal value
        # and the error, when the error is a non-0, real number).
        # The goal is to avoid the strange notation nane-10, and to
        # avoid the 0e10 notation for an exactly zero uncertainty,
        # because .0e can give this for a non-zero error (the goal is
        # to have a zero uncertainty be very explicit):
        error_has_exp = not any_exp_factored and not special_error

        # Like error_has_exp, but only for real number handling
        # (there is no special meaning to a zero nominal value):
        nom_has_exp = not any_exp_factored and isfinite(nom_val_main)

        # Prefix for the parts:
        if fmt_parts["width"]:  # Individual widths
            # If zeros are needed, then the width is taken into
            # account now (before the exponent is added):
            if fmt_parts["zero"]:
                width = int(fmt_parts["width"])

                # Remaining (minimum) width after including the
                # exponent:
                remaining_width = max(width - len(exp_str), 0)

                fmt_prefix_n = "%s%s%d%s" % (
                    fmt_parts["sign"],
                    fmt_parts["zero"],
                    remaining_width if nom_has_exp else width,
                    fmt_parts["comma"],
                )

                fmt_prefix_e = "%s%d%s" % (
                    fmt_parts["zero"],
                    remaining_width if error_has_exp else width,
                    fmt_parts["comma"],
                )

            else:
                fmt_prefix_n = fmt_parts["sign"] + fmt_parts["comma"]
                fmt_prefix_e = fmt_parts["comma"]

        else:  # Global width
            fmt_prefix_n = fmt_parts["sign"] + fmt_parts["comma"]
            fmt_prefix_e = fmt_parts["comma"]

        ## print "ANY_EXP_FACTORED", any_exp_factored
        ## print "ERROR_HAS_EXP", error_has_exp
        ## print "NOM_HAS_EXP", nom_has_exp

        ####################
        # Nominal value formatting:

        # !! The following fails with Python < 2.6 when the format is
        # not accepted by the % operator. This can happen when
        # special_error is true, as the format used for the nominal
        # value is essentially the format provided by the user, which
        # may be empty:

        # print "FMT_PREFIX_N", fmt_prefix_n
        # print "FMT_SUFFIX_N", fmt_suffix_n

        nom_val_str = robust_format(nom_val_main, fmt_prefix_n + fmt_suffix_n)

        # print "NOM_VAL_STR", nom_val_str

        ####################
        # Error formatting:

        # !! Note: .0f applied to a float has no decimal point, but
        # this does not appear to be documented
        # (http://docs.python.org/2/library/string.html#format-specification-mini-language). This
        # feature is used anyway, because it allows a possible comma
        # format parameter to be handled more conveniently than if the
        # 'd' format was used.
        #
        # The following uses a special integer representation of a
        # zero uncertainty:
        if error_main:
            # The handling of NaN/inf in the nominal value identical to
            # the handling of NaN/inf in the standard deviation:
            if (
                not isfinite(nom_val_main)
                # Only some formats have a nicer representation:
                and fmt_parts["type"] in ("", "g", "G")
            ):
                # The error can be formatted independently:
                fmt_suffix_e = (fmt_parts["prec"] or "") + fmt_parts["type"]
            else:
                fmt_suffix_e = ".%d%s" % (prec, main_pres_type)
        else:
            fmt_suffix_e = ".0%s" % main_pres_type

        error_str = robust_format(error_main, fmt_prefix_e + fmt_suffix_e)

        ##########
        # Overriding of nom_val_str and error_str for LaTeX:
        if "L" in options:
            if isnan(nom_val_main):
                nom_val_str = r"\mathrm{%s}" % nom_val_str
            elif isinf(nom_val_main):
                nom_val_str = r"%s\infty" % ("-" if nom_val_main < 0 else "")

            if isnan(error_main):
                error_str = r"\mathrm{%s}" % error_str
            elif isinf(error_main):
                error_str = r"\infty"

        if nom_has_exp:
            nom_val_str += exp_str
        if error_has_exp:
            error_str += exp_str

        ####################
        # Final alignment of each field, if needed:

        if fmt_parts["width"]:  # An individual alignment is needed:
            # Default alignment, for numbers: to the right (if no
            # alignment is specified, a string is aligned to the
            # left):
            effective_align = fmt_parts["align"] or ">"

            # robust_format() is used because it may handle alignment
            # options, where the % operator does not:

            nom_val_str = robust_align(
                nom_val_str, fmt_parts["fill"], effective_align, fmt_parts["width"]
            )

            error_str = robust_align(
                error_str, fmt_parts["fill"], effective_align, fmt_parts["width"]
            )

        ####################
        pm_symbol = PM_SYMBOLS[print_type]  # Shortcut

        ####################

        # Construction of the final value, value_str, possibly with
        # grouping (typically inside parentheses):

        (LEFT_GROUPING, RIGHT_GROUPING) = GROUP_SYMBOLS[print_type]

        # The nominal value and the error might have to be explicitly
        # grouped together with parentheses, so as to prevent an
        # ambiguous notation. This is done in parallel with the
        # percent sign handling because this sign may too need
        # parentheses.
        if any_exp_factored and common_exp is not None:  # Exponent
            value_str = "".join(
                (
                    LEFT_GROUPING,
                    nom_val_str,
                    pm_symbol,
                    error_str,
                    RIGHT_GROUPING,
                    exp_str,
                    percent_str,
                )
            )
        else:  # No exponent
            value_str = "".join([nom_val_str, pm_symbol, error_str])
            if percent_str:
                value_str = "".join(
                    (LEFT_GROUPING, value_str, RIGHT_GROUPING, percent_str)
                )
            elif "p" in options:
                value_str = "".join((LEFT_GROUPING, value_str, RIGHT_GROUPING))

    return value_str


def signif_dgt_to_limit(value, num_signif_d):
    """
    Return the precision limit necessary to display value with
    num_signif_d significant digits.

    The precision limit is given as -1 for 1 digit after the decimal
    point, 0 for integer rounding, etc. It can be positive.
    """

    fst_digit = first_digit(value)

    limit_no_rounding = fst_digit - num_signif_d + 1

    # The number of significant digits of the uncertainty, when
    # rounded at this limit_no_rounding level, can be too large by 1
    # (e.g., with num_signif_d = 1, 0.99 gives limit_no_rounding = -1, but
    # the rounded value at that limit is 1.0, i.e. has 2
    # significant digits instead of num_signif_d = 1). We correct for
    # this effect by adjusting limit if necessary:
    rounded = round(value, -limit_no_rounding)
    fst_digit_rounded = first_digit(rounded)

    if fst_digit_rounded > fst_digit:
        # The rounded limit is fst_digit_rounded-num_signif_d+1;
        # but this can only be 1 above the non-rounded limit:
        limit_no_rounding += 1

    return limit_no_rounding


def format_ufloat(ufloat_to_format, format_spec):
    """
    Formats a number with uncertainty.

    The format specification are the same as for format() for
    floats, as defined for Python 2.6+ (restricted to what the %
    operator accepts, if using an earlier version of Python),
    except that the n presentation type is not supported. In
    particular, the usual precision, alignment, sign flag,
    etc. can be used. The behavior of the various presentation
    types (e, f, g, none, etc.) is similar. Moreover, the format
    is extended: the number of digits of the uncertainty can be
    controlled, as is the way the uncertainty is indicated (with
    +/- or with the short-hand notation 3.14(1), in LaTeX or with
    a simple text string,...).

    Beyond the use of options at the end of the format
    specification, the main difference with floats is that a "u"
    just before the presentation type (f, e, g, none, etc.)
    activates the "uncertainty control" mode (e.g.: ".6u").  This
    mode is also activated when not using any explicit precision
    (e.g.: "g", "10f", "+010,e" format specifications).  If the
    uncertainty does not have a meaningful number of significant
    digits (0 and NaN uncertainties), this mode is automatically
    deactivated.

    The nominal value and the uncertainty always use the same
    precision. This implies trailing zeros, in general, even with
    the g format type (contrary to the float case). However, when
    the number of significant digits of the uncertainty is not
    defined (zero or NaN uncertainty), it has no precision, so
    there is no matching. In this case, the original format
    specification is used for the nominal value (any "u" is
    ignored).

    Any precision (".p", where p is a number) is interpreted (if
    meaningful), in the uncertainty control mode, as indicating
    the number p of significant digits of the displayed
    uncertainty. Example: .1uf will return a string with one
    significant digit in the uncertainty (and no exponent).

    If no precision is given, the rounding rules from the
    Particle Data Group are used, if possible
    (http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf). For
    example, the "f" format specification generally does not use
    the default 6 digits after the decimal point, but applies the
    PDG rules.

    A common exponent is used if an exponent is needed for the
    larger of the nominal value (in absolute value) and the
    standard deviation, unless this would result in a zero
    uncertainty being represented as 0e... or a NaN uncertainty as
    NaNe.... Thanks to this common exponent, the quantity that
    best describes the associated probability distribution has a
    mantissa in the usual 1-10 range. The common exponent is
    factored (as in "(1.2+/-0.1)e-5"). unless the format
    specification contains an explicit width (" 1.2e-5+/- 0.1e-5")
    (this allows numbers to be in a single column, when printing
    numbers over many lines). Specifying a minimum width of 1 is a
    way of forcing any common exponent to not be factored out.

    The fill, align, zero and width parameters of the format
    specification are applied individually to each of the nominal
    value and standard deviation or, if the shorthand notation is
    used, globally.

    The sign parameter of the format specification is only applied
    to the nominal value (since the standard deviation is
    positive).

    In the case of a non-LaTeX output, the returned string can
    normally be parsed back with ufloat_fromstr(). This however
    excludes cases where numbers use the "," thousands separator,
    for example.

    Options can be added, at the end of the format
    specification. Multiple options can be specified:

    - When "P" is present, the pretty-printing mode is activated: "±"
      separates the nominal value from the standard deviation, exponents
      use superscript characters, etc.
    - When "S" is present (like in .1uS), the short-hand notation 1.234(5)
      is used, indicating an uncertainty on the last digits; if the digits
      of the uncertainty straddle the decimal point, it uses a fixed-point
      notation, like in 12.3(4.5).
    - When "L" is present, the output is formatted with LaTeX.
    - "p" ensures that there are parentheses around the …±… part (no
      parentheses are added if some are already present, for instance
      because of an exponent or of a trailing % sign, etc.). This produces
      outputs like (1.0±0.2) or (1.0±0.2)e7, which can be useful for
      removing any ambiguity if physical units are added after the printed
      number.

    An uncertainty which is exactly zero is represented as the
    integer 0 (i.e. with no decimal point).

    The "%" format type forces the percent sign to be at the end
    of the returned string (it is not attached to each of the
    nominal value and the standard deviation).

    Some details of the formatting can be customized as described
    in format_num().
    """

    # Convention on limits "between" digits: 0 = exactly at the
    # decimal point, -1 = after the first decimal, 1 = before the
    # units digit, etc.

    # Convention on digits: 0 is units (10**0), 1 is tens, -1 is
    # tenths, etc.

    # This method does the format specification parsing, and
    # calculates the various parts of the displayed value
    # (mantissas, exponent, position of the last digit). The
    # formatting itself is delegated to format_num().

    ########################################

    # Format specification parsing:

    match = re.match(
        r"""
        (?P<fill>[^{}]??)(?P<align>[<>=^]?)  # fill cannot be { or }
        (?P<sign>[-+ ]?)
        (?P<zero>0?)
        (?P<width>\d*)
        (?P<comma>,?)
        (?:\.(?P<prec>\d+))?
        (?P<uncert_prec>u?)  # Precision for the uncertainty?
        # The type can be omitted. Options must not go here:
        (?P<type>[eEfFgG%]??)  # n not supported
        (?P<options>[PSLp]*)  # uncertainties-specific flags
        $""",
        format_spec,
        re.VERBOSE,
    )

    # Does the format specification look correct?
    if not match:
        raise ValueError(
            "Format specification %r cannot be used with object of type"
            " %r. Note that uncertainties-specific flags must be put at"
            " the end of the format string."
            # Sub-classes handled:
            % (format_spec, ufloat_to_format.__class__.__name__)
        )

    # Effective format presentation type: f, e, g, etc., or None,
    # like in
    # https://docs.python.org/3.4/library/string.html#format-specification-mini-language. Contrary
    # to what is written in the documentation, it is not true that
    # None is "the same as 'g'": "{}".format() and "{:g}" do not
    # give the same result, on 31415000000.0. None is thus kept as
    # is instead of being replaced by "g".
    pres_type = match.group("type") or None

    # Shortcut:
    fmt_prec = match.group("prec")  # Can be None

    ########################################

    # Since the '%' (percentage) format specification can change
    # the value to be displayed, this value must first be
    # calculated. Calculating the standard deviation is also an
    # optimization: the standard deviation is generally
    # calculated: it is calculated only once, here:
    nom_val = ufloat_to_format.nominal_value
    std_dev = ufloat_to_format.std_dev

    # 'options' is the options that must be given to format_num():
    options = set(match.group("options"))

    ########################################

    # The '%' format is treated internally as a display option: it
    # should not be applied individually to each part:
    if pres_type == "%":
        # Because '%' does 0.0055*100, the value
        # 0.5499999999999999 is obtained, which rounds to 0.5. The
        # original rounded value is 0.006. The same behavior is
        # found in Python 2.7: '{:.1%}'.format(0.0055) is '0.5%'.
        # If a different behavior is needed, a solution to this
        # problem would be to do the rounding before the
        # multiplication.
        std_dev *= 100
        nom_val *= 100
        pres_type = "f"
        options.add("%")

    # At this point, pres_type is in eEfFgG or None (not %).

    ########################################

    # Non-real values (nominal value or standard deviation) must
    # be handled in a specific way:
    real_values = [value for value in [abs(nom_val), std_dev] if isfinite(value)]

    # Calculation of digits_limit, which defines the precision of
    # the nominal value and of the standard deviation (it can be
    # None when it does not matter, like for NaN±NaN):

    # Reference value for the calculation of a possible exponent,
    # if needed:
    if pres_type in (None, "e", "E", "g", "G"):
        # Reference value for the exponent: the largest value
        # defines what the exponent will be (another convention
        # could have been chosen, like using the exponent of the
        # nominal value, irrespective of the standard deviation):
        try:
            exp_ref_value = max(real_values)
        except ValueError:  # No non-NaN value: NaN±NaN…
            # No meaningful common exponent can be obtained:
            pass
        ## else:
        ##     print "EXP_REF_VAL", exp_ref_value

    # Should the precision be interpreted like for a float, or
    # should the number of significant digits on the uncertainty
    # be controlled?
    if (
        (
            # Default behavior: number of significant digits on the
            # uncertainty controlled (if useful, i.e. only in
            # situations where the nominal value and the standard
            # error digits are truncated at the same place):
            (not fmt_prec and len(real_values) == 2) or match.group("uncert_prec")
        )  # Explicit control
        # The number of significant digits of the uncertainty must
        # be meaningful, otherwise the position of the significant
        # digits of the uncertainty does not have a clear
        # meaning. This gives us the *effective* uncertainty
        # control mode:
        and std_dev
        and isfinite(std_dev)
    ):
        # The number of significant digits on the uncertainty is
        # controlled.

        # The limit digits_limit on the digits of nom_val and std_dev
        # to be displayed is calculated. If the exponent notation is
        # used, this limit is generally different from the finally
        # displayed limit (e.g. 314.15+/-0.01 has digits_limit=-2, but
        # will be displayed with an exponent as (3.1415+/-0.0001)e+02,
        # which corresponds to 4 decimals after the decimal point, not
        # 2).

        # Number of significant digits to use:
        if fmt_prec:
            num_signif_d = int(fmt_prec)  # Can only be non-negative
            if not num_signif_d:
                raise ValueError(
                    "The number of significant digits"
                    " on the uncertainty should be positive"
                )
        else:
            (num_signif_d, std_dev) = PDG_precision(std_dev)

        digits_limit = signif_dgt_to_limit(std_dev, num_signif_d)

    else:
        # No control of the number of significant digits on the
        # uncertainty.

        ## print "PRECISION NOT BASED ON UNCERTAINTY"

        # The precision has the same meaning as for floats (it is
        # not the uncertainty that defines the number of digits).

        # The usual default precision is used (this is useful for
        # 3.141592±NaN with an "f" format specification, for
        # example):
        #
        # prec is the precision for the main parts of the final
        # format (in the sense of float formatting):
        #
        # https://docs.python.org/3.4/library/string.html#format-specification-mini-language
        if fmt_prec:
            prec = int(fmt_prec)
        elif pres_type is None:
            prec = 12
        else:
            prec = 6

        if pres_type in ("f", "F"):
            digits_limit = -prec

        else:  # Format type in None, eEgG
            # We first calculate the number of significant digits
            # to be displayed (if possible):

            if pres_type in ("e", "E"):
                # The precision is the number of significant
                # digits required - 1 (because there is a single
                # digit before the decimal point, which is not
                # included in the definition of the precision with
                # the e/E format type):
                num_signif_digits = prec + 1

            else:  # Presentation type in None, g, G
                # Effective format specification precision: the rule
                # of
                # http://docs.python.org/2.7/library/string.html#format-specification-mini-language
                # is used:

                # The final number of significant digits to be
                # displayed is not necessarily obvious: trailing
                # zeros are removed (with the gG presentation
                # type), so num_signif_digits is the number of
                # significant digits if trailing zeros were not
                # removed. This quantity is relevant for the
                # rounding implied by the exponent test of the g/G
                # format:

                # 0 is interpreted like 1 (as with floats with a
                # gG presentation type):
                num_signif_digits = prec or 1

            # The number of significant digits is important for
            # example for determining the exponent:

            ## print "NUM_SIGNIF_DIGITS", num_signif_digits

            digits_limit = (
                signif_dgt_to_limit(exp_ref_value, num_signif_digits)
                if real_values
                else None
            )

            ## print "DIGITS_LIMIT", digits_limit

    #######################################

    # Common exponent notation: should it be used? use_exp is set
    # accordingly. If a common exponent should be used (use_exp is
    # True), 'common_exp' is set to the exponent that should be
    # used.

    if pres_type in ("f", "F"):
        use_exp = False

    else:  # e, E, g, G, None
        # The rules from
        # https://docs.python.org/3.4/library/string.html#format-specification-mini-language
        # are applied.

        # Python's native formatting (whose result could be parsed
        # in order to determine whether a common exponent should
        # be used) is not used because there is shared information
        # between the nominal value and the standard error (same
        # last digit, common exponent) and extracting this
        # information from Python would entail parsing its
        # formatted string, which is in principle inefficient
        # (internally, Python performs calculations that yield a
        # string, and the string would be parsed back into
        # separate parts and numbers, which is in principle
        # unnecessary).

        # Should the scientific notation be used? The same rule as
        # for floats is used ("-4 <= exponent of rounded value <
        # p"), on the nominal value.

        if use_exp := real_values:
            # The number of significant digits of the reference value
            # rounded at digits_limit is exponent-digits_limit+1:
            common_exp = first_digit(round(exp_ref_value, -digits_limit))
            common_factor = 10.0**common_exp
            # cases where this doesn't apply: too many digits when expressed like this,
            # or the common factor is way too small.
            if pres_type not in ("e", "E"):  # g, G, None
                use_exp = not (-4 <= common_exp < common_exp - digits_limit + 1) and (
                    common_factor != 0.0
                )

    ########################################

    # Calculation of signif_limit (position of the significant
    # digits limit in the final fixed point representations; this
    # is either a non-positive number, or None), of
    # nom_val_mantissa ("mantissa" for the nominal value,
    # i.e. value possibly corrected for a factorized exponent),
    # and std_dev_mantissa (similarly for the standard
    # deviation). common_exp is also set to None if no common
    # exponent should be used.

    if use_exp:
        nom_val_mantissa = nom_val / common_factor
        std_dev_mantissa = std_dev / common_factor
        # Limit for the last digit of the mantissas:
        signif_limit = digits_limit - common_exp

    else:  # No common exponent
        common_exp = None

        nom_val_mantissa = nom_val
        std_dev_mantissa = std_dev
        signif_limit = digits_limit

    ## print "SIGNIF_LIMIT", signif_limit

    ########################################

    # Format of the main (i.e. with no exponent) parts (the None
    # presentation type is similar to the g format type):

    main_pres_type = "fF"[(pres_type or "g").isupper()]

    # The precision of the main parts must be adjusted so as
    # to take into account the special role of the decimal
    # point:
    if signif_limit is not None:  # If signif_limit is pertinent
        # The decimal point location is always included in the
        # printed digits (e.g., printing 3456 with only 2
        # significant digits requires to print at least four
        # digits, like in 3456 or 3500).
        #
        # The max() is important for example for
        # 1234567.89123+/-12345.678 with the f format: in this
        # case, signif_limit is +3 (2 significant digits necessary
        # for the error, as per the PDG rules), but the (Python
        # float formatting) precision to be used for the main
        # parts is 0 (all digits must be shown).
        #
        # The 1 for the None pres_type represents "at least one
        # digit past the decimal point" of Python
        # (https://docs.python.org/3.4/library/string.html#format-specification-mini-language). This
        # is only applied for null uncertainties.
        prec = max(-signif_limit, 1 if pres_type is None and not std_dev else 0)
    ## print "PREC", prec

    ########################################

    # print (
    #     "FORMAT_NUM parameters: nom_val_mantissa={},"
    #     " std_dev_mantissa={}, common_exp={},"
    #     " match.groupdict()={}, prec={}, main_pres_type={},"
    #     " options={}".format(
    #     nom_val_mantissa, std_dev_mantissa, common_exp,
    #     match.groupdict(),
    #     prec,
    #     main_pres_type,
    #     options))

    # Final formatting:
    return format_num(
        nom_val_mantissa,
        std_dev_mantissa,
        common_exp,
        match.groupdict(),
        prec=prec,
        main_pres_type=main_pres_type,
        options=options,
    )
