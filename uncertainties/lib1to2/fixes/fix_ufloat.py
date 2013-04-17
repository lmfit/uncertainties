'''
Fixer for lib2to3.

Transforms ufloat(tuple,...) and ufloat(string,...) into
ufloat(nominal_value, std_dev,...) and ufloat_fromstr

(c) 2013 by Eric O. LEBIGOT.
'''

from lib2to3.fixer_base import BaseFix
from lib2to3.fixer_util import ArgList, Call, Comma, Name, syms

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
        power< 'ufloat' {tuple_call} any* >
        |
        power< 'ufloat' {tuple_any_call} any* >
        |
        power< 'ufloat' trailer< '(' string=STRING ')' > any* >
        |
        power< 'ufloat' trailer< '('
            arglist<
                string=STRING
                ',' tag=any
            >
        ')' > any* >
        |
        power< object=NAME trailer< '.' 'ufloat' > {tuple_call} any* >
        |
        power< object=NAME trailer< '.' 'ufloat' > {tuple_any_call} any* >
        |
        power< object=NAME trailer< '.' 'ufloat' >
        trailer< '(' string=STRING ')' >
        any* >
        |
        power< object=NAME trailer< '.' 'ufloat' >
        trailer< '(' arglist< string=STRING ',' tag=any > ')' >
        any* >
        """.format(tuple_call=tuple_call,
                   tuple_any_call=tuple_any_call)

    
    def transform(self, node, results):
        
        # Handling of the first argument:
        
        if 'string' in results:  # String as first argument
            
            new_func_name = 'ufloat_fromstr'

            # New arguments:
            new_args=[results['string'].clone()]

        else:  # Tuple as first argument

            new_func_name = 'ufloat'
            
            # New arguments:
            new_args = [results['arg0'].clone(),
                        Comma(), results['arg1'].clone()]

        # Handling of the second argument (call with a tag):
        if 'tag' in results:
            new_args.extend([Comma(), results['tag'].clone()])            

        if 'object' in results:  # If dotted access: unc.ufloat()
            func_name = node.children[1].children[1]
            args = node.children[2]
        else:
            func_name = node.children[0]
            args = node.children[1]
            
        # Function name update:
        func_name.value = new_func_name
        #! func_name.changed()  # Necessary when only .value is changed
        
        # Argument list update:
        args.replace(ArgList(new_args))
