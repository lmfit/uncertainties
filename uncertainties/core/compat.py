import numbers
import warnings

# Some types known to not depend on Variable objects are put in
# CONSTANT_TYPES.  The most common types can be put in front, as this
# may slightly improve the execution speed.
FLOAT_LIKE_TYPES = (numbers.Number,)
CONSTANT_TYPES = FLOAT_LIKE_TYPES+(complex,)

try:
    import numpy
except ImportError:
    pass
else:

    # NumPy numbers do not depend on Variable objects:
    FLOAT_LIKE_TYPES += (numpy.generic,)
    CONSTANT_TYPES += FLOAT_LIKE_TYPES[-1:]


###############################################################################

# Utility for issuing deprecation warnings

def deprecation(message):
    '''
    Warn the user with the given message, by issuing a
    DeprecationWarning.
    '''

    # stacklevel = 3 points to the original user call (not to the
    # function from this module that called deprecation()).
    # DeprecationWarning is ignored by default: not used.

    warnings.warn('Obsolete: %s Code can be automatically updated with'
                  ' python -m uncertainties.1to2 -w ProgramDirectory.'
                  % message, stacklevel=3)


########################################
class CallableStdDev(float):
    '''
    Class for standard deviation results, which used to be
    callable. Provided for compatibility with old code. Issues an
    obsolescence warning upon call.
    '''

    # This class is a float. It must be set to the standard deviation
    # upon construction.

    def __call__ (self):
        deprecation('the std_dev attribute should not be called'
                    ' anymore: use .std_dev instead of .std_dev().')
        return self