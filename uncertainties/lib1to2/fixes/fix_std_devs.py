'''
Fixer for lib2to3.

Transform .std_devs() calls into .std_devs attribute access.

(c) 2016 by Eric O. LEBIGOT.
'''

from lib2to3.fixer_base import BaseFix
from lib2to3.fixer_util import Name, Assign

class FixStdDevs(BaseFix):

    PATTERN = """
    power< any* trailer< '.' 'std_devs' > trailer< '(' ')' > >
    """

    def transform(self, node, results):

        # '.std_dev' is followed by a call with no argument: the call
        # is removed:
        node.children[-1].remove()
