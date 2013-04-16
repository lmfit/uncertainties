'''
Fixer for lib2to3.

Transforms ufloat(tuple,...) and ufloat(string,...) into
ufloat(nominal_value, std_dev,...) and ufloat_fromstr

(c) 2013 by Eric O. LEBIGOT.
'''

from lib2to3.fixer_base import BaseFix
from lib2to3.fixer_util import Call, Name, Comma
import lib2to3.pgen2.token as token

###############################################################################
# lib2to3 grammar parts.

#! Warning: indentation is meaningful!

# (tuple):
tuple_call = """
        trailer< '('
            atom< '(' testlist_gexp< arg0=any ',' arg1=any > ')' >
        ')' >"""

# (tuple, any):
tuple_any_call = """
        trailer< '('
            arglist<
                atom< '(' testlist_gexp< arg0=any ',' arg1=any > ')' >
                ',' tag=any
            >
        ')' >"""


###############################################################################

class FixUfloat(BaseFix):

    # Non dotted access, then dotted access.
    # Tuple call, then string call.
    # No-tag call, then tag call.
    PATTERN = """
        power< 'ufloat' {tuple_call} >
        |
        power< 'ufloat' {tuple_any_call} >
        |
        power< 'ufloat' trailer< '(' string=STRING ')' > >
        |
        power< 'ufloat' trailer< '('
            arglist<
                string=STRING
                ',' tag=any
            >
        ')' > >
        |
        power< object=NAME trailer< '.' 'ufloat' > {tuple_call} >
        |
        power< object=NAME trailer< '.' 'ufloat' > {tuple_any_call} >
        |
        power< object=NAME trailer< '.' 'ufloat' >
        trailer< '(' string=STRING ')' > >
        |
        power< object=NAME trailer< '.' 'ufloat' > trailer< '('
            arglist<
                string=STRING
                ',' tag=any
            >
        ')' > >

        """.format(tuple_call=tuple_call,
                   tuple_any_call=tuple_any_call)

    
    def transform(self, node, results):

        if 'string' in results:
            # Single argument form:

            func_name = 'ufloat_fromstr'
            
            # A constant string can be handled:

            # New arguments:
            args=[results['string'].clone()]

        else:
            # Tuple as first argument:

            func_name = 'ufloat'
            
            # New arguments:
            args = [results['arg0'].clone(), Comma(), results['arg1'].clone()]


        if 'object' in results:
            func_name = '{}.{}'.format(results['object'], func_name)
        
        # Call with a tag:
        if 'tag' in results:
            args.extend([Comma(), results['tag'].clone()])            

        # Replacement by a direct call with the arguments:
        node.replace(Call(Name(func_name), args=args, prefix=node.prefix))
            
            
# print FixUfloat.PATTERN
