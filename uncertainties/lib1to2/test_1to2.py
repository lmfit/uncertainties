#!/usr/bin/env python

'''
Unit tests for the uncertainties.lib1to2 code update package.

Meant to be run through nosetests.

(c) 2013 by Eric O. LEBIGOT (EOL).
'''

# Code inspired by:
#
# - lib2to3.tests.test_fixers.py

import sys
import os
import lib2to3.tests.support as support

# The following lib1to2 module must refer to the *local* package (not
# to any other installed module):
sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

def check_refactor(refactorer, source, expected):
    """
    Raises an AssertionError if the given
    lib2to3.refactor.RefactoringTool does not refactor 'source' into
    'expected'.

    source, expected -- strings (typically with Python code).
    """
    
    new = unicode(
        refactorer.refactor_string(support.reformat(source), '<string>'))

    assert support.reformat(expected) == new
    
def test_fix_name1():
    
    refactorer = support.get_refactorer(fixer_pkg='lib1to2', fixers=['name1'])
    check_refactor(refactorer, "oldname = 123", "newname = 123")
