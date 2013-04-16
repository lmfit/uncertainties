#!/usr/bin/env python

'''
Fixes code like the 2to3 Python utility, but with fixers from the lcoal fixes
directory.
'''

from lib2to3.main import main
import sys

sys.exit(main('fixes'))