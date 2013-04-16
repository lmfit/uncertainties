'''
Fixer for lib2to3.

Transforms ufloat(tuple,...) and ufloat(string,...) into
ufloat(nominal_value, std_dev,...) and ufloat_fromstr

(c) 2013 by Eric O. LEBIGOT.
'''

from lib2to3.fixer_base import BaseFix
from lib2to3.fixer_util import Call, Name, Comma
import lib2to3.pgen2.token as token

class FixUfloat(BaseFix):
    
    PATTERN = """
        power< 'ufloat' trailer< '('
            atom< '(' testlist_gexp< arg0=any ',' arg1=any > ')' >
        ')' > >
        |
        power< 'ufloat' trailer< '('
            arglist<
                atom< '(' testlist_gexp< arg0=any ',' arg1=any > ')' >
                ',' tag=any
            >
        ')' > >
        |
        power< 'ufloat' trailer< '(' string=STRING ')' > >
        |
        power< 'ufloat' trailer< '('
            arglist<
                string=STRING
                ',' tag=any
            >
        ')' > >
        """
    
    def transform(self, node, results):


        if 'string' in results:
            # Single argument form:
            
            # A constant string can be handled:
            
            string_leaf = results['string']

            # New arguments:
            args=[string_leaf.clone()]

            # Call with a tag:
            if 'tag' in results:
                args.extend([Comma(), results['tag'].clone()])            
            
            node.replace(Call(
                Name('ufloat_fromstr'), args = args, prefix=node.prefix))

        else:
            # Tuple as first argument:
            
            # New arguments:
            args = [results['arg0'].clone(), Comma(), results['arg1'].clone()]

            # Call with a tag:
            if 'tag' in results:
                args.extend([Comma(), results['tag'].clone()])

            # Replacement by a direct call with the arguments:
            node.replace(Call(Name('ufloat'), args=args, prefix=node.prefix))
