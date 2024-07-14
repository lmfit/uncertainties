.. index:: formatting Variables
.. _formatting guide:

========================================
Formatting Variables with uncertainties
========================================

.. index::
   printing
   formatting

Printing
========

.. Overview:

Numbers with uncertainties can be printed conveniently:

>>> print(x)
0.200+/-0.010

The resulting form can generally be parsed back with
:func:`ufloat_fromstr` (except for the LaTeX form).

.. Precision matching:

The nominal value and the uncertainty always have the **same
precision**: this makes it easier to compare them.

Standard formats
----------------

.. Formatting method:

More **control over the format** can be obtained (in Python 2.6+)
through the usual :func:`format` method of strings:

>>> print('Result = {:10.2f}'.format(x))
Result =       0.20+/-      0.01


.. Legacy formats and base syntax of the format specification:

**All the float format specifications** are accepted, except those
with the ``n`` format type. In particular, a fill character, an
alignment option, a sign or zero option, a width, or the ``%`` format
type are all supported.

The usual **float formats with a precision** retain their original
meaning (e.g. ``.2e`` uses two digits after the decimal point): code
that works with floats produces similar results when running with
numbers with uncertainties.

Precision control
-----------------

.. Precision control:

It is possible to **control the number of significant digits of the
uncertainty** by adding the precision modifier ``u`` after the
precision (and before any valid float format type like ``f``, ``e``,
the empty format type, etc.):

>>> print('1 significant digit on the uncertainty: {:.1u}'.format(x))
1 significant digit on the uncertainty: 0.20+/-0.01
>>> print('3 significant digits on the uncertainty: {:.3u}'.format(x))
3 significant digits on the uncertainty: 0.2000+/-0.0100
>>> print('1 significant digit, exponent notation: {:.1ue}'.format(x))
1 significant digit, exponent notation: (2.0+/-0.1)e-01
>>> print('1 significant digit, percentage: {:.1u%}'.format(x))
1 significant digit, percentage: (20+/-1)%

When :mod:`uncertainties` must **choose the number of significant
digits on the uncertainty**, it uses the `Particle
Data Group
<http://PDG.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf>`_ rounding
rules (these rules keep the number of digits small, which is
convenient for reading numbers with uncertainties, and at the same
time prevent the uncertainty from being displayed with too few
digits):

>>> print('Automatic number of digits on the uncertainty: {}'.format(x))
Automatic number of digits on the uncertainty: 0.200+/-0.010
>>> print(x)
0.200+/-0.010

Custom options
--------------

.. Options:

:mod:`uncertainties` provides even more flexibility through custom
formatting options. They can be added at the end of the format string:

- ``P`` for **pretty-printing**:

  >>> print('{:.2e}'.format(x))
  (2.00+/-0.10)e-01
  >>> print(u'{:.2eP}'.format(x))
  (2.00±0.10)×10⁻¹

  The pretty-printing mode thus uses "±", "×" and superscript
  exponents.

- ``S`` for the **shorthand notation**:

  >>> print('{:+.1uS}'.format(x))  # Sign, 1 digit for the uncertainty, shorthand
  +0.20(1)

  In this notation, the digits in parentheses represent the
  uncertainty on the last digits of the nominal value.

- ``L`` for a **LaTeX** output:

  >>> print(x*1e7)
  (2.00+/-0.10)e+06
  >>> print('{:L}'.format(x*1e7))  # Automatic exponent form, LaTeX
  \left(2.00 \pm 0.10\right) \times 10^{6}

- ``p`` is for requiring that parentheses be always printed around the …±… part
  (without enclosing any exponent or trailing "%", etc.). This can for instance
  be useful so as to explicitly factor physical units:

    >>> print('{:p} kg'.format(x))  # Adds parentheses
    (0.200+/-0.010) kg
    >>> print("{:p} kg".format(x*1e7))  # No parentheses added (exponent)
    (2.00+/-0.10)e+06 kg

These custom formatting options **can be combined** (when meaningful).

Details
-------

.. Common exponent:

A **common exponent** is automatically calculated if an exponent is
needed for the larger of the nominal value (in absolute value) and the
uncertainty (the rule is the same as for floats). The exponent is
generally **factored**, for increased legibility:

>>> print(x*1e7)
(2.00+/-0.10)e+06

When a *format width* is used, the common exponent is not factored:

>>> print('Result = {:10.1e}'.format(x*1e-10))
Result =    2.0e-11+/-   0.1e-11

(Using a (minimal) width of 1 is thus a way of forcing exponents to not
be factored.) Thanks to this feature, each part (nominal value and
standard deviation) is correctly aligned across multiple lines, while the
relative magnitude of the error can still be readily estimated thanks to
the common exponent.

.. Special cases:

An uncertainty which is *exactly* **zero** is always formatted as an
integer:

>>> print(ufloat(3.1415, 0))
3.1415+/-0
>>> print(ufloat(3.1415e10, 0))
(3.1415+/-0)e+10
>>> print(ufloat(3.1415, 0.0005))
3.1415+/-0.0005
>>> print('{:.2f}'.format(ufloat(3.14, 0.001)))
3.14+/-0.00
>>> print('{:.2f}'.format(ufloat(3.14, 0.00)))
3.14+/-0

**All the digits** of a number with uncertainty are given in its
representation:

>>> y = ufloat(1.23456789012345, 0.123456789)
>>> print(y)
1.23+/-0.12
>>> print(repr(y))
1.23456789012345+/-0.123456789
>>> y
1.23456789012345+/-0.123456789


Global formatting
-----------------

It is sometimes useful to have a **consistent formatting** across
multiple parts of a program. Python's `string.Formatter class
<https://docs.python.org/3/library/string.html#custom-string-formatting>`_
allows one to do just that. Here is how it can be used to consistently
use the shorthand notation for numbers with uncertainties:

.. code-block:: python

   class ShorthandFormatter(string.Formatter):

       def format_field(self, value, format_spec):
           if isinstance(value, uncertainties.UFloat):
               return value.format(format_spec+'S')  # Shorthand option added
           # Special formatting for other types can be added here (floats, etc.)
           else:
               # Usual formatting:
               return super(ShorthandFormatter, self).format_field(
                   value, format_spec)

   frmtr = ShorthandFormatter()

   print(frmtr.format("Result = {0:.1u}", x))  # 1-digit uncertainty

prints with the shorthand notation: ``Result = 0.20(1)``.


Customizing the pretty-print and LaTeX outputs
----------------------------------------------

The pretty print and LaTeX outputs themselves can be customized.

For example, the pretty-print representation of numbers with uncertainty can
display multiplication with a centered dot (⋅) instead of the default symbol
(×), like in ``(2.00±0.10)⋅10⁻¹``; this is easily done through the global
setting ``uncertainties.core.MULT_SYMBOLS["pretty-print"] = "⋅"``.

Beyond this multiplication symbol, the "±" symbol, the parentheses and the
exponent representations can also be customized globally. The details can be
found in the documentation of :func:`uncertainties.core.format_num`.
