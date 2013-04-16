#!/usr/bin/env python

'''
Fixes code like the 2to3 Python utility, but with fixers from the local fixes
directory.
'''

from lib2to3.main import main
import sys

sys.exit(main('uncertainties.1to2.fixes'))  # !! Would 'fixes' be enough?