#!/usr/bin/env python

'''
Unit tests for the uncertainties.lib1to2 code update package.

Meant to be run through nosetests.

(c) 2013-2020 by Eric O. LEBIGOT (EOL).
'''

# Code inspired by:
#
# - lib2to3.tests.test_fixers.py

from builtins import str
import sys
import os

# !! Would it be possible to use an import hook so as to stop the
# import if the Python version is not high enough, instead of having
# like here a whole indented block?


if sys.version_info < (2, 7) or "TRAVIS" in os.environ or "APPVEYOR" in os.environ:
    
    # This package uses lib2to3, which requires Python 2.6+.
    
    # lib2to3.tests.support is missing from 2.7.3 Travis Python packages.

    # !!  Nosetests for Python 2.6 also fails (it looks like it tries
    # to run tests via lib2to3/tests/test_refactor.py):
    
    pass

else:

    import os
    try:
        # lib2to3 test support seems to have moved to a new place in 2013:
        import test.test_lib2to3.support as support
    except ImportError:
        # Pre-~2013 path for lib2to3 test support
        import lib2to3.tests.support as support

    # The lib1to2.fixes package given to lib2to3 is the *local* package
    # (not to another installed module). This is important for the
    # __import__() used via support.get_refactorer().
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

    def check_refactor(refactorer, source, expected):
        """
        Raises an AssertionError if the given
        lib2to3.refactor.RefactoringTool does not refactor 'source' into
        'expected'.

        source, expected -- strings (typically with Python code).
        """

        # !! str() is from future's builtins and is only needed for Python 2,
        # where it is mostly equivalent to unicode():
        new = str(
            refactorer.refactor_string(support.reformat(source), '<string>'))

        assert support.reformat(expected) == new, (
            "Refactoring failed: '{}' => '{}' instead of '{}'".format(
            source, new.strip(), expected))

        # print 'Checked:', source, '=>', expected
        
    def check_all(fixer, tests):
        '''
        Takes a fixer name (module from fixes) and a mapping that maps
        code using the obsolete syntax into updated code, and checks
        whether the code is correctly updated.
        '''    

        refactorer = support.get_refactorer(
            fixer_pkg='lib1to2', fixers=[fixer])

        for (input_str, out_str) in tests.items():
            check_refactor(refactorer, input_str, out_str)

    def test_fix_std_dev():
        'Tests the transformation of std_dev() into std_dev.'


        tests = {
            'x.std_dev()': 'x.std_dev',
            'y.std_dev();  unc.std_dev(z)': 'y.std_dev;  unc.std_dev(z)',
            'uncertainties.std_dev(x)': 'uncertainties.std_dev(x)',
            'std_dev(x)': 'std_dev(x)',
            'obj.x.std_dev()': 'obj.x.std_dev',

            """
            long_name.std_dev(
            # No argument!
            )""":
            """
            long_name.std_dev""",

            # set_std_dev => .std_dev:
            'x.set_std_dev(3)': 'x.std_dev = 3',
            'y = set_std_dev(3)': 'y = set_std_dev(3)',  # None
            'func = x.set_std_dev': 'func = x.set_std_dev',
            'obj.x.set_std_dev(sin(y))': 'obj.x.std_dev = sin(y)'
        }

        check_all('std_dev', tests)

    def test_ufloat():
        '''
        Test of the transformation of ufloat(tuple,...) and
        ufloat(string,...) into ufloat(nominal_value, std_dev, tag=...).
        '''

        tests = {
            # Tuples:
            'ufloat((3, 0.14))': 'ufloat(3, 0.14)',
            'ufloat((3, 0.14), "pi")': 'ufloat(3, 0.14, "pi")',
            "ufloat((3, 0.14), 'pi')": "ufloat(3, 0.14, 'pi')",
            "x = ufloat((3, 0.14), tag='pi')": "x = ufloat(3, 0.14, tag='pi')",

            # Simple expressions that can be transformed:
            'ufloat((n, s), tag="var")': 'ufloat(n, s, tag="var")',

            # Simple expressions that cannot be transformed automatically:
            'ufloat(str_repr, tag="var")': 'ufloat(str_repr, tag="var")',
            'ufloat(*tuple_repr, tag="var")': 'ufloat(*tuple_repr, tag="var")',
            'ufloat(*t[0, 0])': 'ufloat(*t[0, 0])',        

            # Strings:
            'ufloat("-1.23(3.4)")': 'ufloat_fromstr("-1.23(3.4)")',
            "ufloat('-1.23(3.4)')": "ufloat_fromstr('-1.23(3.4)')",
            'ufloat("-1.23(3.4)", "var")':
            'ufloat_fromstr("-1.23(3.4)", "var")',
            'ufloat("-1.23(3.4)", tag="var")':
                'ufloat_fromstr("-1.23(3.4)", tag="var")'

        }

        # Automatic addition of a dotted access:
        tests.update(dict(
            # !! Dictionary comprehension usable with Python 2.7+
            (orig.replace('ufloat', 'unc.ufloat'),
             new.replace('ufloat', 'unc.ufloat'))
            for (orig, new) in tests.items()))

        # Test for space consistency:
        tests[' t  =  u.ufloat("3")'] = ' t  =  u.ufloat_fromstr("3")'

        # Exponentiation test:
        tests.update(dict(
            # !! Dictionary comprehension usable with Python 2.7+
            (orig+'**2', new+'**2')
            for (orig, new) in tests.items()))

        # Exponent test:
        tests['2**ufloat("3")'] = '2**ufloat_fromstr("3")'

        # Opposite test:
        tests['-ufloat("3")'] = '-ufloat_fromstr("3")'

        check_all('ufloat', tests)

    def test_uarray_umatrix():
        '''
        Test of the transformation of uarray(tuple,...) into
        uarray(nominal_values, std_devs). Also performs the same tests
        on umatrix().
        '''
        
        tests = {
            'uarray((arange(3), std_devs))': 'uarray(arange(3), std_devs)',
            'uarray(tuple_arg)': 'uarray(*tuple_arg)',
            # Unmodified, correct code:
            'uarray(values, std_devs)': 'uarray(values, std_devs)',
            # Spaces tests:
            'uarray( ( arange(3),  std_devs ) ) ':
            'uarray( arange(3),  std_devs) ',
            'uarray(  tuple_arg )': 'uarray(*  tuple_arg)'

        }

        # Automatic addition of a dotted access:
        tests.update(dict(
            # !! Dictionary comprehension usable with Python 2.7+
            (orig.replace('uarray', 'un.uarray'),
             new.replace('uarray', 'un.uarray'))
            for (orig, new) in tests.items()))
                             
        # Exponentiation test:
        tests.update(dict(
            # !! Dictionary comprehension usable with Python 2.7+
            (orig+'**2', new+'**2')
            for (orig, new) in tests.items()))

        # Test for space consistency:
        tests[' t  =  u.uarray(args)'] = ' t  =  u.uarray(*args)'

        # Same tests, but for umatrix:
        tests.update(dict(
            (orig.replace('uarray', 'umatrix'),
             new.replace('uarray', 'umatrix'))
            for (orig, new) in tests.items()))
        
        check_all('uarray_umatrix', tests)

