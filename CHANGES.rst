Change Log
===================

3.2.2   2024-July-08
-----------------------

Fixes:

 - fix support for Numpy 2.0 (#245).  Note: `uncertainties.unumpy` still
    provides `umatrix` based on `numpy.matrix`.  With `numpy.matrix`
    discouraged, `umatrix` is too, and will be dropped in a  future release.
 - fix automated running and reporting of code coverage with tests (#246)
 - use `setuptools-scm` for setting version number from git tag  (#247)

 3.2.1   2024-June-08
-----------------------

Fixes for build, deployment, and docs

 - Use explicit package list to make sure unumpy is included (#232)
 - Use setuptools-scm to make sure all files are in the source distribution (#235)
 - updates to configuration for and links to readthedocs documentation. (#239)
 - use double backticks more uniformly in docs. (#240)
 - fixes to README.rst to allow it to render (needed for PyPI upload) (#243)

3.2.0   2024-June-02
-----------------------

Version 3.2.0 is the first release of Uncertainties in nearly two years and the
first minor release in over five years. It marks the beginning of an effort to
refresh and update the project with a new and expanded team of maintainers.

* Main Changes

  - Moved code development to lmfit organization, with 4 maintainers.
  - Update documentation.
  - Drop future dependency. Uncertainties now has no external dependencies when
     not using Numpy integration (Drop official support for Python versions before 3.8 #200).
  - Drop support for Python versions before 3.8, including Python 2 (Drop official support for Python versions before 3.8 #200)
  - remove 1to2 and deprecations (remove 1to2 and depreciations #214)

* Developer related changes

  - Moved from setup.py to pyproject.toml (Transition from setup.py to pyproject.toml #199)
  - Move tests to tests folder (Move tests to tests folder #216)
  - Update unumpy test to be compatible with numpy 2
  - Mark docstrings with backslashes as raw strings in tests (Mark docstrings with backslashes as raw strings #226)



Older Version history
------------------------

Main changes:
- 3.1.6: The pretty-print and LaTeX format can now be customized.
- 3.1.5: Added a "p" formatting option, that makes sure that there are always
  parentheses around the … ± … part of printed numbers.
- 3.1.4: Python 2.7+ is now required.
- 3.1.2: Fix for NumPy 1.17 and ``unumpy.ulinalg.pinv()``.
- 3.1: Variables built through a correlation or covariance matrix, and that
  have uncertainties that span many orders of magnitude are now
  calculated more accurately (improved ``correlated_values()`` and
  ``correlated_values_norm()`` functions).
- 3.0: Massive speedup for some operations involving large numbers of numbers with uncertainty, like ``sum(ufloat(1, 1) for _ in xrange(100000))`` (this is about 5,000 times faster than before).
- 2.4.8: Friendlier completions in Python shells, etc.: internal functions should not appear anymore (for the user modules: ``uncertainties``, ``uncertainties.umath`` and  ``uncertainties.unumpy``). Parsing the shorthand notation (e.g. ``3.1(2)``) now works with infinite values (e.g. ``-inf(inf)``); this mirrors the ability to print such numbers with uncertainty. The Particle Data Group rounding rule is applied in more cases (e.g. printing 724.2±26.2 now gives ``724±26``). The shorthand+LaTeX formatting of numbers with an infinite nominal value is fixed. ``uncertainties.unumpy.matrix`` now uses ``.std_devs`` instead of ``.std_devs()``, for consistency with floats with uncertainty (automatic conversion of code added to ``uncertainties.1to2``).
- 2.4.7: String formatting now works for ``(-)inf+/-...`` numbers.
- 2.4.5: String formatting now works for ``NaN+/-...`` numbers.
- 2.4.4: The documentation license now allows its commercial use.
- 2.4.2: `NumPy 1.8 compatibility <https://github.com/numpy/numpy/issues/4063>`_.
- 2.4.1: In ``uncertainties.umath``, functions ``ceil()``, ``floor()``,
  ``isinf()``, ``isnan()`` and ``trunc()`` now return values of
  the same type as the corresponding ``math`` module function
  (instead of generally returning a value with a zero uncertainty
  ``...+/-0``).
- 2.4: Extensive support for the formatting_ of numbers with uncertainties.
  A zero uncertainty is now explicitly displayed as the integer 0.
  The new formats are generally understood by ``ufloat_fromstr()``.
  Abbreviations for the nominal value (``n``) and the standard
  deviation (``s``) are now available.
- 2.3.6:  Full support for limit cases of the power operator
  ``umath.pow()``.
- 2.3.5: Uncertainties and derivatives can now be NaN (not-a-number).
  Full support for numbers with a zero uncertainty
  (``sqrt(ufloat(0, 0))`` now works).
  Full support for limit cases of the power operator (``x**y``).
- 2.3: Functions wrapped
  so that they accept numbers with uncertainties instead of floats
  now have full keyword arguments support
  (improved ``wrap()`` function). Incompatible change:
  ``wrap(..., None)`` should be replaced by ``wrap(...)`` or
  ``wrap(..., [])``.
- 2.2: Creating arrays and matrices of numbers with uncertainties
  with ``uarray()`` and ``umatrix()`` now requires two simple arguments
  (nominal values and standard deviations) instead of a tuple argument.
  This is consistent with the new, simpler ``ufloat()`` interface.
  The previous
  usage will be supported for some time. Users are encouraged to update
  their code, for instance through the newly provided `code updater`_,
  which in addition now automatically converts ``.set_std_dev(v)`` to
  ``.std_dev = v``.
- 2.1: Numbers with uncertainties are now created more directly like
  ``ufloat(3, 0.1)``, ``ufloat(3, 0.1, "pi")``,
  ``ufloat_fromstr("3.0(1)")``, or ``ufloat_fromstr("3.0(1)", "pi")``.
  The previous ``ufloat((3, 0.1))`` and ``ufloat("3.0(1)")`` forms
  will be supported for some time. Users are encouraged to update
  their code, for instance through the newly provided `code updater`_.
- 2.0: The standard deviation is now obtained more directly without an
  explicit
  call (``x.std_dev`` instead of ``x.std_dev()``). ``x.std_dev()``
  will be supported for some time. Users are encouraged to update
  their code. The standard deviation of a variable can now be
  directly updated with ``x.std_dev = 0.1``. As a consequence,
  ``x.set_std_dev()`` is deprecated.
- 1.9.1: Support added for pickling subclasses of ``UFloat`` (= ``Variable``).
- 1.9: Added functions for handling correlation matrices:
  ``correlation_matrix()`` and
  ``correlated_values_norm()``. (These new functions mirror the
  covariance-matrix based ``covariance_matrix()`` and
  ``correlated_values()``.) ``UFloat.position_in_sigmas()`` is
  now named ``UFloat.std_score()``, so as to follow the common
  naming convention (`standard score
  <http://en.wikipedia.org/wiki/Standard_score>`_).  Obsolete
  functions were removed (from the main module:
  ``NumberWithUncert``, ``num_with_uncert``, ``array_u``,
  ``nominal_values``, ``std_devs``).
- 1.8: Compatibility with Python 3.2 added.
- 1.7.2: Compatibility with Python 2.3, Python 2.4, Jython 2.5.1 and
  Jython 2.5.2 added.
- 1.7.1: New semantics: ``ufloat("12.3(78)")`` now represents 12.3+/-7.8
  instead of 12.3+/-78.
- 1.7: ``ufloat()`` now raises ValueError instead of a generic Exception,
  when given an incorrect
  string representation, like ``float()`` does.
- 1.6: Testing whether an object is a number with uncertainty should now
  be done with ``isinstance(..., UFloat)``.
  ``AffineScalarFunc`` is not imported by ``from uncertainties import *``
  anymore, but its new alias ``UFloat`` is.
- 1.5.5: The first possible license is now the Revised BSD License
  instead of GPLv2, which
  makes it easier to include this package in other projects.
- 1.5.4.2: Added ``umath.modf()`` and ``umath.frexp()``.
- 1.5.4: ``ufloat`` does not accept a single number (nominal value) anymore.
  This removes some potential confusion about
  ``ufloat(1.1)`` (zero uncertainty) being different from
  ``ufloat("1.1")`` (uncertainty of 1 on the last digit).
- 1.5.2: ``float_u``, ``array_u`` and ``matrix_u`` renamed ``ufloat``,
  ``uarray`` and ``umatrix``, for ease of typing.
- 1.5:  Added functions ``nominal_value`` and ``std_dev``, and
  modules ``unumpy`` (additional support for NumPy arrays and
  matrices) and ``unumpy.ulinalg`` (generalization of some
  functions from ``numpy.linalg``).
  Memory footprint of arrays of numbers with uncertainties
  divided by 3.
  Function ``array_u`` is 5 times faster.
  Main function ``num_with_uncert`` renamed
  ``float_u``, for consistency with ``unumpy.array_u`` and
  ``unumpy.matrix_u``, with the added benefit of a shorter name.
- 1.4.5: Added support for the standard ``pickle`` module.
- 1.4.2: Added support for the standard ``copy`` module.
- 1.4: Added utilities for manipulating NumPy arrays of numbers with
  uncertainties (``array_u``, ``nominal_values`` and ``std_devs``).
- 1.3: Numbers with uncertainties are now constructed with
  ``num_with_uncert()``, which replaces ``NumberWithUncert()``.  This
  simplifies the class hierarchy by removing the ``NumberWithUncert`` class.
- 1.2.5: Numbers with uncertainties can now be entered as
  ``NumberWithUncert("1.23+/-0.45")`` too.
- 1.2.3: ``log(x, base)`` is now supported by ``umath.log()``, in addition
  to ``log(x)``.
- 1.2.2: Values with uncertainties are now output like 3+/-1, in order
  to avoid confusing 3+-1 with 3+(-1).
- 1.2: A new function, ``wrap()``, is exposed, which allows non-Python
  functions (e.g. Fortran or C used through a module such as SciPy) to
  handle numbers with uncertainties.
- 1.1: Mathematical functions (such as cosine, etc.) are in a new
  uncertainties.umath module;
  they do not override functions from the ``math`` module anymore.
- 1.0.12: Main class (``Number_with_uncert``) renamed ``NumberWithUncert``
  so as to follow `PEP 8`_.
- 1.0.11: ``origin_value`` renamed more appropriately as
  ``nominal_value``.
- 1.0.9: ``correlations()`` renamed more appropriately as
  ``covariance_matrix()``.

.. _math: http://docs.python.org/library/math.html
.. _PEP 8: http://www.python.org/dev/peps/pep-0008/
.. _code updater: http://uncertainties-python-package.readthedocs.io/en/latest/index.html#migration-from-version-1-to-version-2
.. _formatting: http://uncertainties-python-package.readthedocs.io/en/latest/user_guide.html#printing
