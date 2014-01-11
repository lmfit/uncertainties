#!/usr/bin/env python

'''
Fixes code like the 2to3 Python utility, but with fixers from the local fixes
directory.

(c) 2013 by Eric O. LEBIGOT (EOL).
'''

# Code inspired by the 2to3 Python code.

import sys

if sys.version_info < (2, 6):
    sys.exit("Please run this program with Python 2.6+.")
    
import lib2to3.main

sys.exit(lib2to3.main.main('uncertainties.lib1to2.fixes'))
