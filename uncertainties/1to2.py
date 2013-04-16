#!/usr/bin/env python

'''
Program that partially converts programs that use a version of the
uncertainties package prior to 2.0 to the API introduced in 2.0.

(c) 2013 by Eric O. LEBIGOT (EOL).
'''

import re
import tempfile
import shutil

def ufloat_tuple_fix(match):
    '''
    Converts a parsed ufloat call to the current ufloat() recommended
    call signature. Returns the code for this call.
    
    Takes a match object from the re module, whose groups give the
    nominal value, the standard deviation and (if not None), the tag
    argument of an ufloat call (along with the attached comma).
    '''

    return 'ufloat({}, {}{})'.format(
        match.group('nom_value'), match.group('std_dev'),
        match.group('tag') if match.group('tag') is not None else '')

def ufloat_str_fix(match):
    '''
    Converts a parsed ufloat call to the current ufloat() recommended
    call signature. Returns the code for this call.
    
    Takes a match object from the re module, whose groups give string
    and (if not None), the tag argument of an ufloat call (along with
    the attached comma).
    '''

    return 'ufloat_from_str({}{})'.format(
        match.group('string'),
        match.group('tag') if match.group('tag') is not None else '')

def update_line(line):
    '''
    Returns an updated version of the line.

    Not all cases are covered, only simple ones.

    It is possible for the code update to be incorrect, but this
    should only happen in very rare circumstances (a case like
    ufloat((3, 0.14), tag=")") might trip the simple parser used
    because of the tag, which looks superficially like the end of the
    function call).
    '''
    
    # No-call standard deviation std_dev:
    line = re.sub('std_dev\(\)', 'std_dev', line)

    # Update of ufloat(tuple,...):
    line = re.sub(
        '''
        ufloat \(

        \(  # Tuple start

           # The only allowed "," and "()" symbols are assumed
           # to delimit the tuple and function call:

           (?P<nom_value> [^,()]+ )  # Nominal value
           ,
           (?P<std_dev> [^,()]+ )  # Standard deviation  

        \)  # Tuple end

        (?P<tag> , [^,()]+ )?  # Optional tag

        \)  # ufloat() closing parenthesis
        ''',
        ufloat_tuple_fix,
        line,
        re.VERBOSE)

    # Update of ufloat(string,...):
    line = re.sub(
        r'''
        ufloat \(

        (?P<string>
            (?P<quote>['"])  # String start  # "'  # For emacs fontify
            [^'"]+  # String contents
            (?P=quote)  # String end
        )  # End of string
        
        (<?P<tag> , [^,()]+ )?  # Optional tag

        \)  # ufloat() closing parenthesis
        ''',  # "  # For emacs fontify
        ufloat_str_fix,
        line,
        re.VERBOSE)
    
    return line

def update_file(file_path, backup=True):
    '''
    Updates code so that the changes introduced in versions 2.0+ are
    taken into account.

    Updated parts:
    - calls to std_dev(),
    - calls to ufloat().

    The update is done through update_line().
    
    file_path -- path to the file to be updated.
    
    backup -- True if backup files should be created.
    '''

    # This function is only exceptionally used: imports are only done
    # here, not globally:
    import tempfile
    import re
    
    # The new code goes to a temporary file: the input file is not
    # modified even in case of failure:
    tmp_file = tempfile.NamedTemporaryFile(delete=False)
    tmp_file_path = tmp_file.name
    tmp_file.close()  # The file will be reused
    
    with open(file_path, 'rbU') as input_file, \
        open(tmp_file_path, 'wb') as out_file:

        for line in input_file:

            out_file.write(update_line(line))

    # Backup:
    if backup:
        shutil.copy2(file_path, file_path + '.bak')

    # Replacement with the updated code:
    shutil.move(tmp_file_path, file_path)
        
            
def test_update_line():
    '''
    Test of update_line().
    '''

    tests = {
        # Tuples:
        'ufloat((3, 0.14))': 'ufloat(3, 0.14)',
        'ufloat((3, 0.14), "pi")': 'ufloat(3, 0.14, "pi")',
        "ufloat((3, 0.14), 'pi')": "ufloat(3, 0.14, 'pi')",
        "x = ufloat((3, 0.14), tag='pi')": "x = ufloat(3, 0.14, tag='pi')",

        # Strings:
        'ufloat("-1.23(3.4)")': 'ufloat_from_str("-1.23(3.4)")',
        "ufloat('-1.23(3.4)')": "ufloat_from_str('-1.23(3.4)')",
        'ufloat("-1.23(3.4)", "var")': 'ufloat_from_str("-1.23(3.4)", "var")',
        'ufloat("-1.23(3.4)", tag="var")':
            'ufloat_from_str("-1.23(3.4)", tag="var")',

        # Simple expressions that can be transformed:
        'ufloat((n, s), tag="var")': 'ufloat(n, s, tag="var")'
        
        # Simple expressions that cannot be transformed:
        'ufloat(str_repr, tag="var")': 'ufloat(str_repr, tag="var")',
        'ufloat(*tuple_repr, tag="var")': 'ufloat(*tuple_repr, tag="var")',
        'ufloat(*t[0, 0])': 'ufloat(*t[0, 0])',        
        }

    for (input_str, out_str) in tests.items():
        assert update_line(input_str) == out_str
        

    
if __name__ == '__main__':
    pass #!!!!!!! take file names
