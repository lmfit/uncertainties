#!/usr/bin/env python

'''
Fixes code like the 2to3 Python utility, but with fixers from the local fixes
directory.

(c) 2013 by Eric O. LEBIGOT (EOL).
'''

# Code inspired by the 2to3 Python code.

# This module can be parsed through the import of fixers:
if __name__ == '__main__':
    
    import lib2to3.main
    import sys
    sys.exit(lib2to3.main.main('uncertainties.lib1to2.fixes'))
