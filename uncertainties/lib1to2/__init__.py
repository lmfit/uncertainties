'''
Adds a path to this file to the module access path.

This is useful for using lib2to3, which needs this in order to access
the fixers defined in fixers.

(c) 2013 by Eric O. LEBIGOT (EOL).
'''

import os
import sys

# The lib1to2 module must refer to the *local* package (not to any
# other installed module) (this is done through the __import__() used via
# support.get_refactorer()):
sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))
