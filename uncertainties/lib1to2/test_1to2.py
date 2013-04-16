#!/usr/bin/env python

'''
Unit tests for the uncertainties.1to2 package.

(c) 2013 by Eric O. LEBIGOT (EOL).
'''

#!!!!!!! test

import lib2to3.refactor
print lib2to3.refactor.get_fixers_from_package('lib2to3.fixes')

import sys
sys.path.insert(0, '.')
import uncertainties.lib1to2.fixes

###############################################################################
# Code inspired by:
#
# - lib2to3.tests.test_fixers.py

import lib2to3.tests.support as support

refactor = support.get_refactorer(fixer_pkg='uncertainties.lib1to2')

before = support.reformat("oldname = 123")
expected = support.reformat("newname = 123")
new = refactor.refactor_string(before, '<string>')

print before
print expected
print new
assert expected == unicode(new)
