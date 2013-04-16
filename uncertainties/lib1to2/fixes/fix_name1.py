from lib2to3.fixer_base import BaseFix
from lib2to3.pgen2 import token
class FixName1(BaseFix):
    
    _accept_type = token.NAME

    def match(self, node):
        if node.value == 'oldname':
            return True
        return False
    
    def transform(self, node, results):
        node.value = 'newname'
        node.changed()
