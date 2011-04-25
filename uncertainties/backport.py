"""
Partial back-port of functions for older versions of Python.
"""

import sys

# For Python < 2.5:
if sys.version_info[:2] < (2, 5):
    
    def any(iterable):
        for element in iterable:
            if element:
                return True
            return False
        
    if sys.version_info[:2] < (2, 4):
        
        def reversed(sequence):
            return sequence[::-1]

        from sets import Set as set


