====================================
Welcome to the uncertainties package
====================================

The `uncertainties package`_ handles calculations that involve
**numbers with uncertainties** (like 3.14±0.01).  It also
transparently yields the **derivatives** of any expression (these
derivatives are used for calculating uncertainties).

Whatever the complexity of the calculation, this package returns the
result with its uncertainty as predicted by linear `error propagation
theory`_.  In particular, it handles **correlations** between
variables, which sets it apart from many existing error propagation
codes.

Calculations involving numbers with uncertainties are made **very
simple** thanks to this package.  In fact, it :ref:`transparently
<derivatives>` calculates the `numerous derivatives`_ required by
linear error propagation theory.

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
numbers with uncertainties with :ref:`no or little modification <user guide>`.

.. index:: correlations; simple example

Another strength of this package is its correct handling of
**correlations**.  For instance, the following quantity is exactly
zero even though ``x`` has an uncertainty:

  >>> x-x
  0.0

Many other error propagation codes return the incorrect value
0±0.1414… because they assume that the two subtracted quantities are
*independent* random variables.

Operations on **arrays** of numbers with uncertainties are
:ref:`transparently handled <simple_array_use>` too.


**Derivatives** are similarly very :ref:`easy to obtain <derivatives>`:

  >>> (2*x+1000).derivatives[x]
  2.0

They are calculated with a :ref:`fast method <differentiation method>`.

.. index:: installation

.. _installing this package:

Installation and download
=========================

Automatic install
-----------------

One of the automatic installation procedures below might work on your
system, if you have a Python package installer.

Under **Unix**, it may be necessary to prefix the installation command
with ``sudo``, so that the installation program has sufficient access
rights to the system.

If you have setuptools_, you can try to automatically install or
upgrade this package with

.. code-block:: sh

   easy_install --upgrade uncertainties

If you have `pip <http://pip.openplans.org/>`_, you can try to
do

.. code-block:: sh

   pip install --upgrade uncertainties



Manual download and install
---------------------------

Alternatively, you can simply download_ the package archive from the
Python Package Index (PyPI) and unpack it.  The package can then be
installed by going into the unpacked :file:`uncertainties` directory,
and running the provided :file:`setup.py` program with

.. code-block:: sh

   python setup.py install

or, for an installation in a custom directory :file:`my_directory`:

.. code-block:: sh

   python setup.py install --prefix my_directory

or, if additional access rights are needed (Unix):

.. code-block:: sh

   sudo python setup.py install

You can also simply **copy** the :file:`uncertainties` directory to a
location that Python can import from (directory in which scripts using
uncertainties are run, etc.).

Source code
-----------

The **code** and the **documentation source** are available on
GitHub_.  The :mod:`uncertainties` package is written in pure Python,
and contains about 4000 lines of code.  75 % of those lines are
documentation strings and comments.  The remaining 25 % are equally
split between unit tests and the calculation code proper.
:mod:`uncertainties` is thus a **lightweight package** with abundant
documentation and tests.

Available documentation
=======================

The :doc:`user_guide` details many of the features of this package.

The part :doc:`numpy_guide` describes how arrays of numbers with
uncertainties can be created and used.

The :doc:`tech_guide` gives advanced technical details.

In addition to this web documentation, the pydoc_ gives access to many
of the documentation strings included in the code.

.. index:: license

What others say
===============

- "*An awesome python package*" (`Jason Moore
  <http://biosport.ucdavis.edu/blog/2010/05/07/uncertainty-analysis>`_)
- "*Your package is brilliant and I love it.*" (e-mail)
- "*A gift of the gods for the work I\'m doing*" (e-mail)
- "*PyPI\'s uncertainties rocks!*" (`identi.ca
  <http://identi.ca/notice/23330742>`_)
- "*Holy f\*\*\* this would have saved me so much f\*\*\*ing time last
  semester*." (`reddit
  <http://www.reddit.com/r/Python/comments/am84v/now_you_can_do_calculations_with_uncertainties_5/>`_)

License
=======

This software is released under a **dual license**; one of the
following options can be chosen:

1. The `BSD license`_.
2. Any other license, as long as it is obtained from the creator of
   this package.

.. index:: support

Contact
=======

Please send feature requests, bug reports, or feedback to the creator
of :mod:`uncertainties`, `Eric O. LEBIGOT (EOL)`_.

.. figure:: _static/eol.*
   :height: 64
   :width:  64
   :target: http://lebigot.pip.verisignlabs.com/
   :align: center
   :alt: Eric O. LEBIGOT (EOL)

Please `support the development`_ of this program by donating $5 or
more through PayPal!

Acknowledgments
===============

The author wishes to thank Arnaud Delobelle, Pierre Cladé, and
Sebastian Walter for very useful technical input.  I would like to
thank Joaquim Abian, Jason Moore, and the many other users for their
feedback and suggestions, which greatly helped improve this program.


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
.. _BSD license: http://creativecommons.org/licenses/BSD/
.. _uncertainties package: http://pypi.python.org/pypi/uncertainties/
.. _GitHub: http://github.com/lebigot/uncertainties
.. _pydoc: http://docs.python.org/library/pydoc.html
