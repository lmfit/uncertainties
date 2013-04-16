#!/usr/bin/env python

'''
Fixes code like the 2to3 Python utility, but with fixers from the local fixes
directory.
'''

import lib2to3.main
import sys

sys.exit(lib2to3.main.main('uncertainties.1to2.fixes'))
