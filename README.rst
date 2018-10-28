uncertainties
=============


.. image:: https://travis-ci.org/lebigot/uncertainties.svg?branch=master
   :target: https://travis-ci.org/lebigot/uncertainties
.. image:: https://ci.appveyor.com/api/projects/status/j5238244myqx0a0r?svg=true
   :target: https://ci.appveyor.com/project/lebigot/uncertainties
.. image:: https://codecov.io/gh/lebigot/uncertainties/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/lebigot/uncertainties/
.. image:: https://readthedocs.org/projects/uncertainties-python-package/badge/?version=latest
   :target: http://uncertainties-python-package.readthedocs.io/en/latest/?badge=latest
.. image:: https://img.shields.io/pypi/v/uncertainties.svg
   :target: https://pypi.org/project/uncertainties/

   
This is the ``uncertainties`` Python package, which performs **transparent
calculations with uncertainties** (aka "error propagation"):

    >>> from uncertainties import ufloat
    >>> from uncertainties.umath import *  # sin(), etc.
    >>> x = ufloat(1, 0.1)  # x = 1+/-0.1
    >>> print 2*x
    2.00+/-0.20
    >>> sin(2*x)  # In a Python shell, "print" is optional
    0.9092974268256817+/-0.08322936730942848

This package also automatically calculates derivatives:

    >>> (2*x+1000).derivatives[x]
    2.0

Some useful links:

* Documentation: http://uncertainties-python-package.readthedocs.io/
* Issues: https://github.com/lebigot/uncertainties/issues/
* Python Package Index entry: http://pypi.python.org/pypi/uncertainties/
* Code: https://github.com/lebigot/uncertainties/

GitHub
------

The ``release`` branch is the latest stable release for Python 2.7+ (including Python 3+ through
``2to3``), while the ``release_python2.3`` branch is the same but for Python 2.3 to
2.6 (with unit tests only run with Python 2.6). They should pass the tests.


``master*`` branches in the Github repository are bleeding-edge, and do not necessarily pass the tests. The ``master`` and ``master_python2.3`` are the latest, relatively stable versions (while other ``master*`` branches are more experimental).

Other branches might be present in the GitHub repository, but they are
also temporary and represent work in progress that does not necessarily run
properly yet.

License
-------

This package and its documentation are released under the `Revised BSD
License <LICENSE.txt>`_.
