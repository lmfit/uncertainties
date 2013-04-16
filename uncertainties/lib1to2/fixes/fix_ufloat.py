'''
Fixer for lib2to3.

Transforms ufloat(tuple,...) and ufloat(string,...) into
ufloat(nominal_value, std_dev,...) and ufloat_from_str

(c) 2013 by Eric O. LEBIGOT.
'''

from lib2to3.fixer_base import BaseFix
from lib2to3.fixer_util import Call, Name, Comma

class FixUfloat(BaseFix):
    
    PATTERN = """
        power< 'ufloat' trailer< '('
        atom< '(' testlist_gexp< arg0=any ',' arg1=any > ')' >
        ')' > >"""
    
    def transform(self, node, results):

        # print "ARG0", repr(results['arg0'])
        # print "ARG0", repr(results['arg1'])
        
        # Replacement by a direct call with the arguments:
        node.replace(Call(Name('ufloat'), prefix=node.prefix,
                     args=(results['arg0'].clone(),
                           Comma(),
                           results['arg1'].clone())))
