#!/usr/bin/env python

'''
Unit tests for the uncertainties.1to2 package.

Meant to be run through nosetests.

(c) 2013 by Eric O. LEBIGOT (EOL).
'''

# Code inspired by:
#
# - lib2to3.tests.test_fixers.py

import sys
import lib2to3.tests.support as support
import os

# The following lib1to2 module must refer to the *local* package (not
# to any other installed module):
sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

refactor = support.get_refactorer(fixer_pkg='lib1to2')

def test_fix1():
    
    before = support.reformat("oldname = 123")
    expected = support.reformat("newname = 123")
    new = refactor.refactor_string(before, '<string>')

    print before
    print expected
    print new
    assert expected == unicode(new)

