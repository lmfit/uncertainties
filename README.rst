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

``uncertainties`` allows **calculations** such as (2 +/- 0.1)*2 = 4 +/-
0.2 to be **performed transparently**.  Much more complex mathematical
expressions involving numbers with uncertainties can also be evaluated
directly.

The ``uncertainties`` package **takes the pain and complexity out**
of uncertainty calculations.

**Detailed information** about this package can be found on its `main
website`_.

Basic examples
--------------

.. code-block:: python

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
-------------

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
  mathematical expression (they are used by `error
  propagation theory`_, and are thus automatically calculated by this
  module).


Installation or upgrade
-----------------------

Installation instructions are available on the `main web site
<http://uncertainties-python-package.readthedocs.io/en/latest/index.html#installation-and-download>`_
for this package.



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

.. _Python: http://docs.python.org/tutorial/interpreter.html
.. _IPython: http://ipython.readthedocs.io/en/stable/
.. _math: http://docs.python.org/library/math.html
.. _error propagation theory: http://en.wikipedia.org/wiki/Propagation_of_uncertainty
.. _main website: http://uncertainties-python-package.readthedocs.io/
