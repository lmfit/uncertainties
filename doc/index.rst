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

Version 4.0 Coming Soon
=======================

The :mod:`uncertainties` team is working on developing a new major release.
This release modernizes the :mod:`uncertainties` codebase and simplifies the underlying
architecture used for uncertainty propagation. These changes greatly improve the
maintainability of this code base which helps guarantee its functionality and support
into the future. At the heart of the architecture change is the idea that
:class:`UFloat` objects should closely resemble mathematical random variables. Here are
some specific benefits from the architecture update.

- It will now be possible to hash :class:`UFloat` objects in a well-controlled way.
- :class:`UFloat` objects will exhibit more intuitive behavior under copying and
   pickling.
- There is a clear path towards implementing correlation-preserving serialization and
   deserialization of :class:`UFloat` objects.
- Previously when users generated a number with uncertainties, the resulting object may
   have been a :class:`Variable` or :class:`AffineScalarFunc` instance depending on how
   the object was created. Now all numbers with uncertainty are uniformly represented by
   :class:`UFloat` instances. :class:`AffineScalarFunc` remains as a legacy name for
   :class:`UFloat`.

In addition to reformating the underlying architecture, this update makes a number of
changes and simplifications to the API to improve consistency. Some of these API change
break backwards compatibility.
Notable changes are

- Elimination of the :class:`Variable` class.
- Elimination of the :attr:`AffineScalarFunc.derivatives` property.
- Replacing the :func:`AffineScalarFunc.error_components()` method with the
   :attr:`UFloat.error_components` property
- Elimination of a number of :class:`UFloat` instance methods and :mod:`umath` module
   functions which do not make sense for random variables such as comparison methods
   and modulo methods and functions.

There will continue to be support for the latest version 3 release as the community
transitions to the new version.


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
