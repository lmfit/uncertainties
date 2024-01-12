from . core import ufloat, ufloat_fromstr, nominal_value, std_dev, covariance_matrix, UFloat, wrap


__all__ = [

    # All sub-modules and packages are not imported by default,
    # in particular because NumPy might be unavailable.

    'ufloat',  # Main function: returns a number with uncertainty
    'ufloat_fromstr',  # Important function: returns a number with uncertainty

    # Uniform access to nominal values and standard deviations:
    'nominal_value',
    'std_dev',

    # Utility functions (more are exported if NumPy is present):
    'covariance_matrix',

    # Class for testing whether an object is a number with
    # uncertainty.  Not usually created by users (except through the
    # Variable subclass), but possibly manipulated by external code
    # ['derivatives()' method, etc.].
    'UFloat',

    # Wrapper for allowing non-pure-Python function to handle
    # quantitities with uncertainties:
    'wrap'

    ]

    
try:
    import numpy
    from . core import correlated_values, correlated_values_norm, correlation_matrix

    __all__ += [
    'correlated_values',
    'correlated_values_norm',
    'correlation_matrix'
    ]
except ImportError:
    pass