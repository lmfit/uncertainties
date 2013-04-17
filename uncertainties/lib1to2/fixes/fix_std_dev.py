'''
Fixer for lib2to3.

Transforms .std_dev() calls into .std_dev attribute access.

(c) 2013 by Eric O. LEBIGOT.
'''

from lib2to3.fixer_base import BaseFix

class FixStdDev(BaseFix):
    
    PATTERN = "power< any trailer< '.' 'std_dev' > trailer< '(' ')' > >"
    
    def transform(self, node, results):

        # '.std_dev' is followed by a call with no argument: the call
        # is removed:
        node.children[2].remove()

