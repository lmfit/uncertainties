#!/usr/bin/env python
# -*- coding: utf-8 -*-

# !! This program must run with all version of Python since 2.3 included.

import sys
import os

min_version = (2, 3)
error_msg = ("Sorry, this package is for Python %d.%d and higher only."
             % min_version)

try:
    if sys.version_info < min_version:
        sys.exit(error_msg)
except AttributeError:  # sys.version_info was introduced in Python 2.0
    sys.exit(error_msg)

# Determination of the directory that contains the source code:
if os.path.exists('uncertainties'):
    # Case of a direct download of a Python-version-specific Git
    # branch:
    package_dir = 'uncertainties'
else:
    # Case of a PyPI package download:
    if sys.version_info >= (2, 6):
        package_dir = 'uncertainties-py26'
    else:
        package_dir = 'uncertainties-py23'

#! The following code was intended to automatically fetch the version
# number; however, it fails when run from Python3 if the downloaded
# code is not the Python 3 version.  An alternative approach would be
# to run 2to3 just before, instead of using build_py_2to3 (which does
# not modify the source), but care should be taken so that users can
# run the setup.py script many times anyway.
## Access to the local uncertainties package (and not to an already
## installed uncertainties package):
# sys.path.insert(0, package_dir)
# uncertainties = __import__(package_dir)

# Common options for distutils/setuptools's setup():
setup_options = dict(
    name='uncertainties',
    version='2.4.6',
    author='Eric O. LEBIGOT (EOL)',
    author_email='eric.lebigot@normalesup.org',
    url='http://pythonhosted.org/uncertainties/',
      
    license='''\
This software can be used under one of the following two licenses: \
(1) The Revised BSD License. \
(2) Any other license, as long as it is obtained from the original \
author.''',
      
    description=('Transparent calculations with uncertainties on the'
                 ' quantities involved (aka "error propagation");'
                 ' fast calculation of derivatives'),
    
    long_description='''\
Overview
========

``uncertainties`` allows **calculations** such as (2 +/- 0.1)*2 = 4 +/- 
0.2 to be **performed transparently**.  Much more complex mathematical 
expressions involving numbers with uncertainties can also be evaluated 
directly.

The ``uncertainties`` package **takes the pain and complexity out** 
of uncertainty calculations.

**Detailed information** about this package can be found on its `main
website`_.

Basic examples
==============

::

    >>> from uncertainties import ufloat
    
    >>> x = ufloat(2, 0.25)
    >>> x
    2.0+/-0.25
    
    >>> square = x**2  # Transparent calculations
    >>> square
    4.0+/-1.0
    >>> square.nominal_value
    4.0
    >>> square.std_dev  # Standard deviation
    1.0

    >>> square - x*x
    0.0  # Exactly 0: correlations taken into account

    >>> from uncertainties.umath import *  # sin(), etc.
    >>> sin(1+x**2)
    -0.95892427466313845+/-0.2836621854632263
    
    >>> print (2*x+1000).derivatives[x]  # Automatic calculation of derivatives
    2.0
    
    >>> from uncertainties import unumpy  # Array manipulation
    >>> random_vars = unumpy.uarray([1, 2], [0.1, 0.2])
    >>> print random_vars
    [1.0+/-0.1 2.0+/-0.2]
    >>> print random_vars.mean()
    1.50+/-0.11
    >>> print unumpy.cos(random_vars)
    [0.540302305868+/-0.0841470984808 -0.416146836547+/-0.181859485365]

Main features
=============

- **Transparent calculations with uncertainties**: **no or little
  modification of existing code** is needed.  Similarly, the Python_ (or
  IPython_) shell can be used as **a powerful calculator** that
  handles quantities with uncertainties (``print`` statements are
  optional, which is convenient).

- **Correlations** between expressions are correctly taken into
  account.  Thus, ``x-x`` is exactly zero, for instance (most
  implementations found on the web yield a non-zero uncertainty for
  ``x-x``, which is incorrect).

- **Almost all mathematical operations** are supported, including most
  functions from the standard math_ module (sin,...).  Comparison
  operators (``>``, ``==``, etc.) are supported too.
  
- Many **fast operations on arrays and matrices** of numbers with
  uncertainties are supported.

- **Extensive support for printing** numbers with uncertainties
  (including LaTeX support and pretty-printing).
  
- Most uncertainty calculations are performed **analytically**.
  
- This module also gives access to the **derivatives** of any 
  mathematical expression (they are used by error
  propagation theory, and are thus automatically calculated by this
  module).


Installation or upgrade
=======================

Installation instructions are available on the `main web site
<http://pythonhosted.org/uncertainties/#installation-and-download>`_
for this package.

Contact
=======

Please send **feature requests, bug reports, or feedback** to
`Eric O. LEBIGOT (EOL)`_.

Version history
===============

Main changes:

- 2.4.5: String formatting now works for ``NaN+/-...`` numbers.
- 2.4.4: The documentation license now allows its commercial use.
- 2.4.2: `NumPy 1.8 compatibility <https://github.com/numpy/numpy/issues/4063>`_.
- 2.4.1: In ``uncertainties.umath``, functions ``ceil()``, ``floor()``, \
         ``isinf()``, ``isnan()`` and ``trunc()`` now return values of \
         the same type as the corresponding ``math`` module function \
         (instead of generally returning a value with a zero uncertainty \
         ``...+/-0``).
- 2.4: Extensive support for the formatting_ of numbers with uncertainties. \
       A zero uncertainty is now explicitly displayed as the integer 0. \
       The new formats are generally understood by ``ufloat_fromstr()``. \
       Abbreviations for the nominal value (``n``) and the standard \
       deviation (``s``) are now available.
- 2.3.6:  Full support for limit cases of the power operator \
          ``umath.pow()``.
- 2.3.5: Uncertainties and derivatives can now be NaN (not-a-number). \
         Full support for numbers with a zero uncertainty \
         (``sqrt(ufloat(0, 0))`` now works). \
         Full support for limit cases of the power operator (``x**y``).
- 2.3: Functions wrapped \
       so that they accept numbers with uncertainties instead of floats \
       now have full keyword arguments support \
       (improved ``wrap()`` function). Incompatible change: \
       ``wrap(..., None)`` should be replaced by ``wrap(...)`` or \
       ``wrap(..., [])``.
- 2.2: Creating arrays and matrices of numbers with uncertainties \
       with ``uarray()`` and ``umatrix()`` now requires two simple arguments \
       (nominal values and standard deviations) instead of a tuple argument. \
       This is consistent with the new, simpler ``ufloat()`` interface. \
       The previous \
       usage will be supported for some time. Users are encouraged to update \
       their code, for instance through the newly provided `code updater`_, \
       which in addition now automatically converts ``.set_std_dev(v)`` to \
       ``.std_dev = v``.
- 2.1: Numbers with uncertainties are now created more directly like \
       ``ufloat(3, 0.1)``, ``ufloat(3, 0.1, "pi")``, \
       ``ufloat_fromstr("3.0(1)")``, or ``ufloat_fromstr("3.0(1)", "pi")``. \
       The previous ``ufloat((3, 0.1))`` and ``ufloat("3.0(1)")`` forms \
       will be supported for some time. Users are encouraged to update \
       their code, for instance through the newly provided `code updater`_.
- 2.0: The standard deviation is now obtained more directly without an \
       explicit \
       call (``x.std_dev`` instead of ``x.std_dev()``). ``x.std_dev()`` \
       will be supported for some time. Users are encouraged to update \
       their code. The standard deviation of a variable can now be \
       directly updated with ``x.std_dev = 0.1``. As a consequence, \
       ``x.set_std_dev()`` is deprecated.
- 1.9.1: Support added for pickling subclasses of ``UFloat`` (= ``Variable``).
- 1.9: Added functions for handling correlation matrices: \
       ``correlation_matrix()`` and \
       ``correlated_values_norm()``. (These new functions mirror the \
       covariance-matrix based ``covariance_matrix()`` and \
       ``correlated_values()``.) ``UFloat.position_in_sigmas()`` is \
       now named ``UFloat.std_score()``, so as to follow the common \
       naming convention (`standard score \
       <http://en.wikipedia.org/wiki/Standard_score>`_).  Obsolete \
       functions were removed (from the main module: \
       ``NumberWithUncert``, ``num_with_uncert``, ``array_u``, \
       ``nominal_values``, ``std_devs``).
- 1.8: Compatibility with Python 3.2 added.
- 1.7.2: Compatibility with Python 2.3, Python 2.4, Jython 2.5.1 and \
         Jython 2.5.2 added.
- 1.7.1: New semantics: ``ufloat("12.3(78)")`` now represents 12.3+/-7.8 \
         instead of 12.3+/-78.
- 1.7: ``ufloat()`` now raises ValueError instead of a generic Exception, \
       when given an incorrect \
       string representation, like ``float()`` does.
- 1.6: Testing whether an object is a number with uncertainty should now \
       be done with ``isinstance(..., UFloat)``. \
       ``AffineScalarFunc`` is not imported by ``from uncertainties import *`` \
       anymore, but its new alias ``UFloat`` is.
- 1.5.5: The first possible license is now the Revised BSD License \
         instead of GPLv2, which \
         makes it easier to include this package in other projects.
- 1.5.4.2: Added ``umath.modf()`` and ``umath.frexp()``.
- 1.5.4: ``ufloat`` does not accept a single number (nominal value) anymore. \
       This removes some potential confusion about \
       ``ufloat(1.1)`` (zero uncertainty) being different from \
       ``ufloat("1.1")`` (uncertainty of 1 on the last digit).
- 1.5.2: ``float_u``, ``array_u`` and ``matrix_u`` renamed ``ufloat``, \
       ``uarray`` and ``umatrix``, for ease of typing.
- 1.5:  Added functions ``nominal_value`` and ``std_dev``, and \
       modules ``unumpy`` (additional support for NumPy_ arrays and \
       matrices) and ``unumpy.ulinalg`` (generalization of some \
       functions from ``numpy.linalg``). \
       Memory footprint of arrays of numbers with uncertainties \
       divided by 3. \
       Function ``array_u`` is 5 times faster. \
       Main function ``num_with_uncert`` renamed \
       ``float_u``, for consistency with ``unumpy.array_u`` and \
       ``unumpy.matrix_u``, with the added benefit of a shorter name.
- 1.4.5: Added support for the standard ``pickle`` module.
- 1.4.2: Added support for the standard ``copy`` module.
- 1.4: Added utilities for manipulating NumPy_ arrays of numbers with\
       uncertainties (``array_u``, ``nominal_values`` and ``std_devs``).
- 1.3: Numbers with uncertainties are now constructed with \
  ``num_with_uncert()``, which replaces ``NumberWithUncert()``.  This \
  simplifies the class hierarchy by removing the ``NumberWithUncert`` class.
- 1.2.5: Numbers with uncertainties can now be entered as \
         ``NumberWithUncert("1.23+/-0.45")`` too.
- 1.2.3: ``log(x, base)`` is now supported by ``umath.log()``, in addition \
         to ``log(x)``.
- 1.2.2: Values with uncertainties are now output like 3+/-1, in order \
         to avoid confusing 3+-1 with 3+(-1).
- 1.2: A new function, ``wrap()``, is exposed, which allows non-Python \
       functions (e.g. Fortran or C used through a module such as SciPy) to \
       handle numbers with uncertainties.
- 1.1: Mathematical functions (such as cosine, etc.) are in a new \
       uncertainties.umath module; \
       they do not override functions from the ``math`` module anymore.
- 1.0.12: Main class (``Number_with_uncert``) renamed ``NumberWithUncert`` \
          so as to follow `PEP 8`_.
- 1.0.11: ``origin_value`` renamed more appropriately as \
          ``nominal_value``.
- 1.0.9: ``correlations()`` renamed more appropriately as \
         ``covariance_matrix()``.

.. _Python: http://docs.python.org/tutorial/interpreter.html
.. _IPython: http://ipython.scipy.org/
.. _NumPy: http://numpy.scipy.org/
.. _math: http://docs.python.org/library/math.html
.. _PEP 8: http://www.python.org/dev/peps/pep-0008/
.. _error propagation theory: http://en.wikipedia.org/wiki/Propagation\
_of_uncertainty
.. _Eric O. LEBIGOT (EOL): mailto:eric.lebigot@normalesup.org
.. _PayPal: https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=4TK7KNDTEDT4S
.. _main website: http://pythonhosted.org/uncertainties/
.. _code updater: http://pythonhosted.org/uncertainties/#migration-from-version-1-to-version-2
.. _formatting: http://pythonhosted.org/uncertainties/user_guide.html#printing
''',
      
    keywords=['error propagation', 'uncertainties',
              'uncertainty calculations',
              'standard deviation',
              'derivatives', 'partial derivatives', 'differentiation'],
    
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Other Audience',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.3',
        'Programming Language :: Python :: 2.4',
        'Programming Language :: Python :: 2.5',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        # Python 3.1 failed because of a problem with NumPy 1.6.1 (whereas
        # everything was fine with Python 3.2 and Python 2.7).
        'Programming Language :: Python :: 3.1',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: Implementation :: Jython',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Topic :: Education',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Software Development',
        'Topic :: Software Development :: Libraries',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Utilities'
    ],

    # Where to find the source code:
    package_dir={'uncertainties': package_dir},

    # Files are defined in MANIFEST (which is automatically created by
    # python setup.py sdist):
    packages=['uncertainties', 'uncertainties.unumpy',
              'uncertainties.lib1to2', 'uncertainties.lib1to2.fixes']
    )

# The best available setup() is used (some users do not have
# setuptools):
try:
    from setuptools import setup

    # Some setuptools-specific options can be added:
    
    addtl_setup_options = {
        
        # Allows python setup.py nosetests to do the right thing:
        'use_2to3': True,

        'tests_require': ['nose', 'numpy'],
        # Optional dependencies install using:
        # `easy_install uncertainties[optional]`
        'extras_require': {
            'optional': ['numpy'],
            'docs': ['sphinx'],
        }
    }

    # easy_install uncertainties[tests] option:
    addtl_setup_options['extras_require']['tests'] = (
        addtl_setup_options['tests_require'])
    
    # easy_install uncertainties[all] option: all dependencies are
    # gathered
    addtl_setup_options['extras_require']['all'] = set(
        sum(addtl_setup_options['extras_require'].values(), []))

    setup_options.update(addtl_setup_options)
    
except ImportError:
    from distutils.core import setup

    # distutils.core.setup is not like setuptools: it does not have an
    # option for automatically calling 2to3 if Python 3 is used, so the
    # conversion is done here:
    if sys.version_info >= (3,):
        import subprocess
        subprocess.check_call(["2to3", "-w", "."])

# End of setup definition

setup(**setup_options)
