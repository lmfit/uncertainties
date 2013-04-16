from lib2to3.fixer_base import BaseFix
from lib2to3.pgen2 import token

class FixStdDev(BaseFix):
    
    _accept_type = token.NAME

    def match(self, node):
        return node.value == 'std_dev'
    
    def transform(self, node, results):

        # If 'std_dev' is followed by a call with no argument, the
        # call is removed:
        
        print "NODE", node
        print "TYPE NODE", type(node)
        print "DIR NODE", dir(node)
        print "RESULTS", results

        # !!!!!!!! should I work by catching std_dev or ()?
        print "NEXT SIBLING", repr(node.next_sibling)
        node.next_sibling.remove()
        
        #node.value = 'newname'
        #node.changed()
