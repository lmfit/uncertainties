"""
Partial back-port of functions for older versions of Python.

This module is intended to be imported as 'from backport import *'.
in fact, this does not redefine any, all, etc. if they already exist.
Furthermore, no program must depend on backport.* functions.
"""

import sys

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
