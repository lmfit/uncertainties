uncertainties
=============

.. image:: https://readthedocs.org/projects/uncertainties/badge/?version=latest
   :target: https://uncertainties.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. image:: https://img.shields.io/pypi/v/uncertainties.svg
   :target: https://pypi.org/project/uncertainties/
.. image:: https://pepy.tech/badge/uncertainties/week
   :target: https://pepy.tech/project/uncertainties
.. image:: https://codecov.io/gh/lmfit/uncertainties/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/lmfit/uncertainties/
.. image:: https://img.shields.io/github/actions/workflow/status/lmfit/uncertainties/python-package.yml?logo=github%20actions
   :target: https://github.com/lmfit/uncertainties/actions/workflows/python-package.yml

This is the ``uncertainties`` Python package, which performs **transparent
calculations with uncertainties** (aka "error propagation"):

    >>> from uncertainties import ufloat
    >>> from uncertainties.umath import *  # sin(), etc.
    >>> x = ufloat(1, 0.1)  # x = 1+/-0.1
    >>> print(2*x)
    2.00+/-0.20
    >>> sin(2*x)  # In a Python shell, "print" is optional
    0.9092974268256817+/-0.08322936730942848

This package also **automatically calculates derivatives of arbitrary functions**:

    >>> (2*x+1000).derivatives[x]
    2.0

The main documentation is available at
https://uncertainties.readthedocs.io/.

Git branches
------------

The ``release`` branch is the latest stable release. It should pass the tests.


``master*`` branches in the Github repository are bleeding-edge, and do not
necessarily pass the tests. The ``master`` branch is the latest, relatively
stable versions (while other ``master*`` branches are more experimental).

Other branches might be present in the GitHub repository, but they are
typically temporary and represent work in progress that does not necessarily run
properly yet.

License
-------

This package and its documentation are released under the `Revised BSD
License <LICENSE.txt>`_.

Voluntary donations
-------------------
If you find this open-source software useful (e.g. in saving you time or helping you produce
something valuable), please consider `donating $10 or more <https://www.paypal.com/donate/?cmd=_s-xclick&hosted_button_id=4TK7KNDTEDT4S>`_.

History
-------

..
   Note from Eric Lebigot: I would like the origin of the package to
   remain documented for its whole life. Thanks!

This package was created back around 2009 by `Eric O. LEBIGOT <https://github.com/lebigot>`_.

Ownership of the package was taken over by the `lmfit GitHub organization <https://github.com/lmfit>`_ in 2024.
