#!/usr/bin/env python

'''
Unit tests for the uncertainties.1to2 package.
'''

from lib2to3.tests.support import get_refactorer

refactor = get_refactorer(fixer_pkg='uncertainties.1to2')

before = "oldname = 123"
expected = "othername = 42"
tree = refactor.refactor_string(before, self.filename)
print expected
print tree
assert expected == tree
