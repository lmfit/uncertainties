"""
Partial back-port of functions for older versions of Python.
"""

import sys

# Whatever the version of Python, a few backport.* names are defined
# (this allows programs to rely on backport.* names: tests can be run
# even from more modern versions of Python).

# For Python < 2.5:
if sys.version_info < (2, 5):
    
    def any(iterable):
        for element in iterable:
            if element:
                return True
            return False
        
    def all(iterable):
        for element in iterable:
            if not element:
                return False
            return True
        
    if sys.version_info < (2, 4):
        
        def reversed(sequence):
            return sequence[::-1]

        from sets import Set as set

    else:

        reversed = reversed
        set = set
        
else:

    any = any
    all = all
