'''
Fixer for lib2to3.

Transforms .std_dev() calls into .std_dev attribute access.

(c) 2013 by Eric O. LEBIGOT.
'''

from lib2to3.fixer_base import BaseFix
from lib2to3.fixer_util import Name, Assign

class FixStdDev(BaseFix):
    
    PATTERN = """
    power< any* trailer< '.' 'std_dev' > trailer< '(' ')' > >
    |
    power< any* trailer< '.' 'set_std_dev' > trailer< '(' set_arg=any ')' > >
    """
    
    def transform(self, node, results):

        if 'set_arg' in results:  # Case of .set_std_dev()

            # set_std_dev => std_dev
            attribute = node.children[-2]  # .set_std_dev
            attribute.children[1].replace(Name('std_dev'))

            # Call "(arg)": removed
            node.children[-1].remove()

            # Replacement by an assignment:
            node.replace(Assign(node.clone(), results['set_arg'].clone()))
            
        else:
            # '.std_dev' is followed by a call with no argument: the call
            # is removed:
            node.children[-1].remove()
        
