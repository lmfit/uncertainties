
import math
from math import sqrt, log, isnan, isinf  # Optimization: no attribute look-up
try:
    from math import isinfinite  # !! Python 3.2+
except ImportError:
    def isinfinite(x):
        return isinf(x) or isnan(x)

def first_digit(value):
    '''
    Return the first digit position of the given value, as an integer.

    0 is the digit just before the decimal point. Digits to the right
    of the decimal point have a negative position.

    Return 0 for a null value.
    '''
    try:
        return int(math.floor(math.log10(abs(value))))
    except ValueError:  # Case of value == 0
        return 0

def PDG_precision(std_dev):
    '''
    Return the number of significant digits to be used for the given
    standard deviation, according to the rounding rules of the
    Particle Data Group (2010)
    (http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf).

    Also returns the effective standard deviation to be used for
    display.
    '''

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
        (exponent, factor) = (exponent-2, 1)
    else:
        (exponent, factor) = (exponent+1, 1000)
    digits = int(std_dev/10.**exponent*factor)  # int rounds towards zero

    # Rules:
    if digits <= 354:
        return (2, std_dev)
    elif digits <= 949:
        return (1, std_dev)
    else:
        # The parentheses matter, for very small or very large
        # std_dev:
        return (2, 10.**exponent*(1000/factor))

# Definition of a basic (format specification only, no full-feature
# format string) formatting function that works whatever the version
# of Python. This function exists so that the more capable format() is
# used instead of the % formatting operator, if available:
robust_format = format



# Exponent letter: the keys are the possible main_fmt_type values of
# format_num():
EXP_LETTERS = {'f': 'e', 'F': 'E'}

def robust_align(orig_str, fill_char, align_option, width):
    '''
    Aligns the given string with the given fill character.

    orig_str -- string to be aligned (str or unicode object).

    fill_char -- if empty, space is used.

    align_option -- as accepted by format().

    wdith -- string that contains the width.
    '''

    # print "ALIGNING", repr(orig_str), "WITH", fill_char+align_option,
    # print "WIDTH", width

    return format(orig_str, fill_char+align_option+width)

# Maps some Unicode code points ("-", "+", and digits) to their
# superscript version:
TO_SUPERSCRIPT = {
    0x2b: u'⁺',
    0x2d: u'⁻',
    0x30: u'⁰',
    0x31: u'¹',
    0x32: u'²',
    0x33: u'³',
    0x34: u'⁴',
    0x35: u'⁵',
    0x36: u'⁶',
    0x37: u'⁷',
    0x38: u'⁸',
    0x39: u'⁹'
    }

# Inverted TO_SUPERSCRIPT table, for use with unicode.translate():
#
#! Python 2.7+ can use a dictionary comprehension instead:
FROM_SUPERSCRIPT = {
    ord(sup): normal for (normal, sup) in TO_SUPERSCRIPT.items()}

def to_superscript(value):
    '''
    Return a (Unicode) string with the given value as superscript characters.

    The value is formatted with the %d %-operator format.

    value -- integer.
    '''

    return (u'%d' % value).translate(TO_SUPERSCRIPT)

def nrmlze_superscript(number_str):
    '''
    Return a string with superscript digits transformed into regular digits.

    Non-superscript digits are not changed before the conversion. Thus, the
    string can also contain regular digits.

    ValueError is raised if the conversion cannot be done.

    number_str -- string to be converted (of type str, but also possibly, for 
    Python 2, unicode, which allows this string to contain superscript digits).
    '''
    # !! Python 3 doesn't need this str(), which is only here for giving the
    # .translate() method to str objects in Python 2 (this str() comes
    # from the builtins module of the future package and is therefore
    # a subclass of unicode, in Python 2):
    return int(str(number_str).translate(FROM_SUPERSCRIPT))

PM_SYMBOLS = {'pretty-print': u'±', 'latex': r' \pm ', 'default': '+/-'}

# Multiplication symbol for pretty printing (so that pretty printing can
# be customized):
MULT_SYMBOLS = {'pretty-print': u'×', 'latex': r'\times'}

# Function that transforms a numerical exponent produced by format_num() into
# the corresponding string notation (for non-default modes):
EXP_PRINT = {
    'pretty-print': lambda common_exp: u'%s10%s' % (
        MULT_SYMBOLS['pretty-print'], to_superscript(common_exp)),
    'latex': lambda common_exp: r' %s 10^{%d}' % (
        MULT_SYMBOLS['latex'], common_exp)}

# Symbols used for grouping (typically between parentheses) in format_num():
GROUP_SYMBOLS = {
    'pretty-print': ('(', ')'),
    # Because of possibly exponents inside the parentheses (case of a
    # specified field width), it is better to use auto-adjusting
    # parentheses. This has the side effect of making the part between
    # the parentheses non-breakable (the text inside parentheses in a
    # LaTeX math expression $...$ can be broken).
    'latex': (r'\left(', r'\right)'),
    'default': ('(', ')')  # Basic text mode
    }

def format_num(nom_val_main, error_main, common_exp,
               fmt_parts, prec, main_pres_type, options):
    u'''
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
    '''

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
    if 'P' in options:
        print_type = 'pretty-print'
    elif 'L' in options:
        print_type = 'latex'
    else:
        print_type = 'default'

    # Exponent part:
    if common_exp is None:
        exp_str = ''
    elif print_type == 'default':
        # Case of e or E. The same convention as Python 2.7
        # to 3.3 is used for the display of the exponent:
        exp_str = EXP_LETTERS[main_pres_type]+'%+03d' % common_exp
    else:
        exp_str = EXP_PRINT[print_type](common_exp)

    # Possible % sign:
    percent_str = ''
    if '%' in options:
        if 'L' in options:
            # % is a special character, in LaTeX: it must be escaped.
            #
            # Using '\\' in the code instead of r'\' so as not to
            # confuse emacs's syntax highlighting:
            percent_str += ' \\'
        percent_str += '%'

    ####################

    # Only true if the error should not have an exponent (has priority
    # over common_exp):
    special_error = not error_main or isinfinite(error_main)

    # Nicer representation of the main nominal part, with no trailing
    # zeros, when the error does not have a defined number of
    # significant digits:
    if special_error and fmt_parts['type'] in ('', 'g', 'G'):
        # The main part is between 1 and 10 because any possible
        # exponent is taken care of by common_exp, so it is
        # formatted without an exponent (otherwise, the exponent
        # would have to be handled for the LaTeX option):
        fmt_suffix_n = (fmt_parts['prec'] or '')+fmt_parts['type']
    else:
        fmt_suffix_n = '.%d%s' % (prec, main_pres_type)


    # print "FMT_SUFFIX_N", fmt_suffix_n

    ####################

    # Calculation of the mostly final numerical part value_str (no %
    # sign, no global width applied).

    # Error formatting:


    if 'S' in options:  # Shorthand notation:

        # Calculation of the uncertainty part, uncert_str:

        if error_main == 0:
            # The error is exactly zero
            uncert_str = '0'
        elif isnan(error_main):
            uncert_str = robust_format(error_main, main_pres_type)
            if 'L' in options:
                uncert_str = r'\mathrm{%s}' % uncert_str
        elif isinf(error_main):
            if 'L' in options:
                uncert_str = r'\infty'
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
                uncert_str = '%.*f' % (prec, uncert)

            else:
                if uncert:
                    # The round is important because 566.99999999 can
                    # first be obtained when 567 is wanted (%d prints the
                    # integer part, not the rounded value):
                    uncert_str = '%d' % round(uncert*10.**prec)
                else:
                    # The decimal point indicates a truncated float
                    # (this is easy to do, in this case, since
                    # fmt_prefix_e is ignored):
                    uncert_str = '0.'

        # End of the final number representation (width and alignment
        # not included). This string is important for the handling of
        # the width:
        value_end = '(%s)%s%s' % (uncert_str, exp_str, percent_str)
        any_exp_factored = True  # Single exponent in the output

        ##########
        # Nominal value formatting:

        # Calculation of fmt_prefix_n (prefix for the format of the
        # main part of the nominal value):

        if fmt_parts['zero'] and fmt_parts['width']:

            # Padding with zeros must be done on the nominal value alone:

            # Remaining width (for the nominal value):
            nom_val_width = max(int(fmt_parts['width']) - len(value_end), 0)
            fmt_prefix_n = '%s%s%d%s' % (
                fmt_parts['sign'], fmt_parts['zero'], nom_val_width,
                fmt_parts['comma'])

        else:
            # Any 'zero' part should not do anything: it is not
            # included
            fmt_prefix_n = fmt_parts['sign']+fmt_parts['comma']

        # print "FMT_PREFIX_N", fmt_prefix_n
        # print "FMT_SUFFIX_N", fmt_suffix_n

        nom_val_str = robust_format(nom_val_main, fmt_prefix_n+fmt_suffix_n)

        ##########
        # Overriding of nom_val_str for LaTeX,; possibly based on the
        # existing value (for NaN vs nan):
        if 'L' in options:

            if isnan(nom_val_main):
                nom_val_str = r'\mathrm{%s}' % nom_val_str
            elif isinf(nom_val_main):
                # !! It is wasteful, in this case, to replace
                # nom_val_str: could this be avoided while avoiding to
                # duplicate the formula for nom_val_str for the common
                # case (robust_format(...))?
                nom_val_str = r'%s\infty' % ('-' if nom_val_main < 0 else '')

        value_str = nom_val_str+value_end

        # Global width, if any:

        if fmt_parts['width']:  # An individual alignment is needed:

            # Default alignment, for numbers: to the right (if no
            # alignment is specified, a string is aligned to the
            # left):
            value_str = robust_align(
                value_str, fmt_parts['fill'], fmt_parts['align'] or '>',
                fmt_parts['width'])

    else:  # +/- notation:

        # The common exponent is factored or not, depending on the
        # width. This gives nice columns for the nominal values and
        # the errors (no shift due to a varying exponent), when a need
        # is given:
        any_exp_factored = not fmt_parts['width']

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
        nom_has_exp = not any_exp_factored and not isinfinite(nom_val_main)

        # Prefix for the parts:
        if fmt_parts['width']:  # Individual widths

            # If zeros are needed, then the width is taken into
            # account now (before the exponent is added):
            if fmt_parts['zero']:

                width = int(fmt_parts['width'])

                # Remaining (minimum) width after including the
                # exponent:
                remaining_width = max(width-len(exp_str), 0)

                fmt_prefix_n = '%s%s%d%s' % (
                    fmt_parts['sign'], fmt_parts['zero'],
                    remaining_width if nom_has_exp else width,
                    fmt_parts['comma'])

                fmt_prefix_e = '%s%d%s' % (
                    fmt_parts['zero'],
                    remaining_width if error_has_exp else width,
                    fmt_parts['comma'])

            else:
                fmt_prefix_n = fmt_parts['sign']+fmt_parts['comma']
                fmt_prefix_e = fmt_parts['comma']

        else:  # Global width
            fmt_prefix_n = fmt_parts['sign']+fmt_parts['comma']
            fmt_prefix_e = fmt_parts['comma']

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

        nom_val_str = robust_format(nom_val_main, fmt_prefix_n+fmt_suffix_n)

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
            if (isinfinite(nom_val_main)
                # Only some formats have a nicer representation:
                and fmt_parts['type'] in ('', 'g', 'G')):
                # The error can be formatted independently:
                fmt_suffix_e = (fmt_parts['prec'] or '')+fmt_parts['type']
            else:
                fmt_suffix_e = '.%d%s' % (prec, main_pres_type)
        else:
            fmt_suffix_e = '.0%s' % main_pres_type

        error_str = robust_format(error_main, fmt_prefix_e+fmt_suffix_e)

        ##########
        # Overriding of nom_val_str and error_str for LaTeX:
        if 'L' in options:

            if isnan(nom_val_main):
                nom_val_str = r'\mathrm{%s}' % nom_val_str
            elif isinf(nom_val_main):
                nom_val_str = r'%s\infty' % ('-' if nom_val_main < 0 else '')

            if isnan(error_main):
                error_str = r'\mathrm{%s}' % error_str
            elif isinf(error_main):
                error_str = r'\infty'

        if nom_has_exp:
            nom_val_str += exp_str
        if error_has_exp:
            error_str += exp_str

        ####################
        # Final alignment of each field, if needed:

        if fmt_parts['width']:  # An individual alignment is needed:

            # Default alignment, for numbers: to the right (if no
            # alignment is specified, a string is aligned to the
            # left):
            effective_align = fmt_parts['align'] or '>'

            # robust_format() is used because it may handle alignment
            # options, where the % operator does not:

            nom_val_str = robust_align(
                nom_val_str, fmt_parts['fill'], effective_align,
                fmt_parts['width'])

            error_str = robust_align(
                error_str, fmt_parts['fill'], effective_align,
                fmt_parts['width'])

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
            value_str = ''.join((
                LEFT_GROUPING,
                nom_val_str, pm_symbol, error_str,
                RIGHT_GROUPING,
                exp_str, percent_str))
        else:  # No exponent
            value_str = ''.join([nom_val_str, pm_symbol, error_str])
            if percent_str:
                value_str = ''.join((
                    LEFT_GROUPING, value_str, RIGHT_GROUPING, percent_str))
            elif 'p' in options:
                value_str = ''.join((LEFT_GROUPING, value_str, RIGHT_GROUPING))

    return value_str

def signif_dgt_to_limit(value, num_signif_d):
    '''
    Return the precision limit necessary to display value with
    num_signif_d significant digits.

    The precision limit is given as -1 for 1 digit after the decimal
    point, 0 for integer rounding, etc. It can be positive.
    '''

    fst_digit = first_digit(value)

    limit_no_rounding = fst_digit-num_signif_d+1

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




    