Welcome to the uncertainties package
====================================

The `uncertainties package`_ handles calculations that involve
**numbers with uncertainties** (like 3.14±0.01).  It also
transparently yields the :ref:`derivatives <derivatives>` of any
expression (these derivatives are :ref:`used <linear_method>` for
calculating uncertainties).

Whatever the complexity of the calculation, this package returns the
result with its uncertainty as predicted by linear `error propagation
theory`_.  In particular, it handles **correlations** between
variables, which sets it apart from many existing error propagation
codes.

Calculations involving numbers with uncertainties are made **very
simple** thanks to this package.  In fact, it transparently calculates
the `numerous derivatives`_ required by linear error propagation
theory.

Calculations of results with uncertainties, or of derivatives, can
either be performed in an **interactive session**, or in programs
written in the Python_ programming language.  Existing calculation
code can **run with no or little change**.

Let's now see how to use these unique features!

.. index:: calculator

An easy-to-use calculator
=========================

Calculations involving **numbers with uncertainties** can be performed
even without knowing anything about the Python_ programming language.
After `installing this package`_ and `invoking the Python
interpreter`_, calculations with automatic error propagation can be
performed directly and transparently:

  >>> from uncertainties import ufloat
  >>> from uncertainties.umath import *  # sin(), etc.
  >>> x = ufloat((1, 0.1))  # x = 1+/-0.1
  >>> print 2*x
  2.0+/-0.2
  >>> sin(2*x)  # In a Python shell, "print" is optional
  0.90929742682568171+/-0.083229367309428481

Thus, existing calculation code designed for floats can be run with
numbers with uncertainties with no or little modification.

.. index:: correlations; simple example

Another strength of this package is its correct handling of
**correlations**.  For instance, the following quantity is exactly
zero even though ``x`` has an uncertainty:

  >>> x-x
  0.0

Many other error propagation codes return the incorrect value
0±0.1414… because they assume that the two subtracted quantities are
*independent* random variables.

**Derivatives** are similarly very easy to obtain:

  >>> (2*x+1000).derivatives[x]
  2.0

They are calculated with a :ref:`fast method <differentiation method>`.

.. index:: installation

.. _installing this package:

Installation and download
=========================

Automatic install
-----------------

The :mod:`uncertainties` package may be automatically downloaded and
installed (or updated) with

.. code-block:: sh

   easy_install -U uncertainties

Whether this works depends on your system (this does not require any
manual download, but does require setuptools_).  Under Unix, it may be
necessary to prefix this command with ``sudo``, so that the
installation program has sufficient access rights to the system.

Manual download and install
---------------------------

Alternatively, you can simply download_ the program from the
Python Package Index (PyPI) and, after unpacking it, install it
with:

.. code-block:: sh

   python setup.py install

or, for an installation in a custom directory my_directory:

.. code-block:: sh

   python setup.py install --prefix my_directory

or, if additional access rights are needed (Unix):

.. code-block:: sh

   sudo python setup.py install

You can also simply **copy** the :file:`uncertainties/` directory to a
location that Python can import from (directory in which scripts using
uncertainties are run, etc.).

The **code** and the **documentation source** are available on GitHub_.

Available documentation
=======================

In addition to this introduction, details on the features of the
:mod:`uncertainties` package can be found in the :doc:`user_guide`.

Support for arrays of numbers with uncertainties is described in
:doc:`numpy_guide`.

Technical details about the mathematics behind this package are given
in the :doc:`tech_guide`.

In addition, the documentation for any of the modules defined in this
package is available through the pydoc_ command.

.. index:: license

License
=======

This software is released under a **dual license**; one of the
following options can be chosen:

1. The `GNU General Public License version 2`_.
2. Any other license, as long as it is obtained from the original author.

.. index:: support

Contact
=======

Please send feature requests, bug reports, or feedback to `Eric
O. LEBIGOT (EOL)`_.

.. figure:: _static/eol.*
   :height: 64
   :width:  64

   `Visit EOL's page`_

Please `support the development`_ of this program by donating $5 or
more through PayPal!



.. toctree::
   :hidden:
   :maxdepth: 1

   Overview <self>
   user_guide
   numpy_guide
   tech_guide

.. _Python: http://python.org/
.. _error propagation theory: http://en.wikipedia.org/wiki/Propagation_of_uncertainty
.. _invoking the Python interpreter: http://docs.python.org/tutorial/interpreter.html
.. _setuptools: http://pypi.python.org/pypi/setuptools
.. _download: http://pypi.python.org/pypi/uncertainties/
.. _numerous derivatives: http://en.wikipedia.org/wiki/Propagation_of_uncertainty#Non-linear_combinations

.. _Eric O. LEBIGOT (EOL): mailto:eric.lebigot@normalesup.org
.. _support the development: https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=4TK7KNDTEDT4S
.. _GNU General Public License version 2: http://creativecommons.org/licenses/GPL/2.0/
.. _uncertainties package: http://pypi.python.org/pypi/uncertainties/
.. _Visit EOL's page: http://lebigot.pip.verisignlabs.com/
.. _GitHub: http://github.com/lebigot/uncertainties
.. _pydoc: http://docs.python.org/library/pydoc.html
