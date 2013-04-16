'''
Sets the module path so as to include the directory where this file is.

This allows lib2to3 module to find the fixers from the local fixes/
package.

(c) 2013 by Eric O. LEBIGOT (EOL).
'''

import sys
import os

# The lib1to2 module must refer to the *local* package (not to any
# other installed module) (this is done through the __import__() used via
# support.get_refactorer()):
sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))
