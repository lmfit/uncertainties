"""
Mathematical operations that generalize many operations from the
standard math module so that they also work on numbers with
uncertainties.

Examples:

  from umath import sin

  # Manipulation of numbers with uncertainties:
  x = uncertainties.ufloat(3, 0.1)
  print sin(x)  # prints 0.141120008...+/-0.098999...

  # The umath functions also work on regular Python floats:
  print sin(3)  # prints 0.141120008...  This is a Python float.

Importing all the functions from this module into the global namespace
is possible.  This is encouraged when using a Python shell as a
calculator.  Example:

  import uncertainties
  from uncertainties.umath import *  # Imports tan(), etc.

  x = uncertainties.ufloat(3, 0.1)
  print tan(x)  # tan() is the uncertainties.umath.tan function

The numbers with uncertainties handled by this module are objects from
the uncertainties module, from either the Variable or the
AffineScalarFunc class.

(c) 2009-2016 by Eric O. LEBIGOT (EOL) <eric.lebigot@normalesup.org>.
Please send feature requests, bug reports, or feedback to this address.

This software is released under a dual license.  (1) The BSD license.
(2) Any other license, as long as it is obtained from the original
author."""

from .umath_core import *  # noqa
