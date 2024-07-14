.. meta::
   :description: The uncertainties Python package
   :keywords: error propagation, uncertainties, error calculations, Python,
              calculator, library, package

Uncertainties
=================

The `uncertainties package`_ is an open source Python library for doing
calculations on numbers that have uncertainties (like 3.14±0.01) that are
common in many scientific fields.  The calculations done with this package will
propagate the uncertainties to the result of mathematical calculations.
The :mod:`uncertainties` package takes the pain and
complexity out of uncertainty calculations and error propagation.  Here is a
quick taste of how to use :mod:`uncertainties`:

>>> from uncertainties import ufloat
>>> x = ufloat(2, 0.1)   # x = 2+/-0.1
>>> y = ufloat(3, 0.2)   # y = 3+/-0.2
>>> print(2*x)
4.00+/-0.20
>>> print(x+y)
5.00+/-0.22
>>> print(x*y)
6.0+/-0.5

The :mod:`uncertainties` library calculates uncertainties using linear `error
propagation theory`_ by automatically :ref:`calculating derivatives
<derivatives>` and analytically propagating these to the results.  Correlations
between variables are automatically handled.  This library can also yield the
derivatives of any expression with respect to the variables that have uncertain
values.  For other approaches, see soerp_ (using higher-order terms) and mcerp_
(using a Monte-Carlo approach).

The `source code`_ for the uncertainties package is licensed under the `Revised
BSD License`_.  This documentation is licensed under the `CC-SA-3 License`_.

.. _uncertainties package: https://pypi.python.org/pypi/uncertainties/
.. _error propagation theory: https://en.wikipedia.org/wiki/Propagation_of_uncertainty
.. _soerp: https://pypi.python.org/pypi/soerp
.. _mcerp: https://pypi.python.org/pypi/mcer
.. _Revised BSD License: https://opensource.org/licenses/BSD-3-Clause
.. _CC-SA-3 License: https://creativecommons.org/licenses/by-sa/3.0
.. _source code:   https://github.com/lmfit/uncertainties/
.. _version history: https://pypi.python.org/pypi/uncertainties#version-history

.. _Pint: https://pypi.python.org/pypi/Pint/
.. _future: https://pypi.org/project/future/


Table of Contents
=================

.. toctree::
   :maxdepth: 2

   install
   user_guide
   numpy_guide
   formatting
   tech_guide
