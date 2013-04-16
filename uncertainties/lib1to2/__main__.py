#!/usr/bin/env python

'''
Fixes code like the 2to3 Python utility, but with fixers from the local fixes
directory.

(c) 2013 by Eric O. LEBIGOT (EOL).
'''

# Code inspired by the 2to3 Python code.

import sys
import os
import lib2to3.main

# The lib1to2 module must refer to the *local* package (not to any
# other installed module) (this is done through the __import__() used via
# support.get_refactorer()):
# sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

sys.exit(lib2to3.main.main('uncertainties.lib1to2.fixes'))
