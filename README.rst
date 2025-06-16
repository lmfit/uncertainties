uncertainties
=============

.. image:: https://readthedocs.org/projects/uncertainties/badge/?version=latest
   :target: https://uncertainties.readthedocs.io/en/latest/?badge=latest
.. image:: https://img.shields.io/pypi/v/uncertainties.svg
   :target: https://pypi.org/project/uncertainties/
.. image:: https://pepy.tech/badge/uncertainties/week
   :target: https://pepy.tech/project/uncertainties
.. image:: https://codecov.io/gh/lmfit/uncertainties/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/lmfit/uncertainties/
.. image:: https://img.shields.io/github/actions/workflow/status/lmfit/uncertainties/python-package.yml?logo=github%20actions
   :target: https://github.com/lmfit/uncertainties/actions/workflows/python-package.yml

The ``uncertainties`` package allows calculations with values that have
uncertaintes, such as (2 +/- 0.1)*2 = 4 +/- 0.2.  ``uncertainties`` takes the
pain and complexity out of error propagation and calculations of values with
uncertainties.  For more information, see https://uncertainties.readthedocs.io/

Basic examples
--------------

.. code-block:: python

    >>> from uncertainties import UFloat
    >>> x = UFloat(2, 0.25)
    >>> print(x)
    2.0+/-0.25

    >>> print(x**2)
    >>> square = x**2
    >>> print(square)
    4.0+/-1.0
    >>> print(square.nominal_value)
    4.0
    >>> print(square.std_dev)  # Standard deviation
    1.0

    >>> print(square - x*x)
    0.0  # Exactly 0: `uncertainties` is aware of correlations

    >>> from uncertainties.umath import sin, cos  # and many more.
    >>> print(sin(1+x**2))
    -0.95892427466313845+/-0.2836621854632263

    >>> print(2*x+1000).derivatives[x]  # Automatic calculation of derivatives
    2.0

    >>> from uncertainties import unumpy  # Array manipulation
    >>> varr = unumpy.uarray([1, 2], [0.1, 0.2])
    >>> print(varr)
    [1.0+/-0.1 2.0+/-0.2]
    >>> print(varr.mean())
    1.50+/-0.11
    >>> print(unumpy.cos(varr))
    [0.540302305868+/-0.0841470984808 -0.416146836547+/-0.181859485365]

Main features
-------------

- **Transparent calculations with uncertainties**: Little or
  no modification of existing code is needed to convert calculations of floats
  to calculations of values with uncertainties.

- **Correlations** between expressions are correctly taken into
  account.  Thus, ``x-x`` is exactly zero.

- **Many  mathematical operations** are supported, including many
  functions from the standard math_ module (sin,...).

- Many **fast operations on arrays and matrices** of numbers with
  uncertainties are supported.

- **Extensive support for formatting** numbers with uncertainties
  (including LaTeX support and pretty-printing).

- Most uncertainty calculations are performed **analytically**.

- This module gives access to the error components and weights which contributed to the
   total uncertainty on any :class:`UFloat` object.

Installation or upgrade
-----------------------

To install `uncertainties`, use::

     pip install uncertainties


Further details are in the `on-line documentation
<https://uncertainties.readthedocs.io/en/latest/install.html>`_.


Git branches
------------

The GitHub ``main`` branch is the latest development version, and is intended to be a
stable pre-release version. It will be experimental, but should pass all tests.  Tagged
releases will be available on GitHub, and correspond to the releases to PyPI.  The
GitHub ``gh-pages`` branch will contain a stable test version of the documentation that
can be viewed at `<https://lmfit.github.io/uncertainties/>`_.  Other Github branches
should be treated as unstable and in-progress development branches.


License
-------

This package and its documentation are released under the `Revised BSD
License <LICENSE.txt>`_.


History
-------

..
   Note from Eric Lebigot: I would like the origin of the package to
   remain documented for its whole life. Thanks!

This package was created back around 2009 by `Eric O. LEBIGOT <https://github.com/lebigot>`_.

Ownership of the package was taken over by the `lmfit GitHub organization <https://github.com/lmfit>`_ in 2024.

.. _IPython: https://ipython.readthedocs.io/en/stable/
.. _math: https://docs.python.org/library/math.html
.. _error propagation theory: https://en.wikipedia.org/wiki/Propagation_of_uncertainty
.. _main website: https://uncertainties.readthedocs.io/
