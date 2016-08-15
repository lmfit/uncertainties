.. meta::
   :description: The uncertainties Python package
   :keywords: error propagation, uncertainties, error calculations, Python,
              calculator, library, package

====================================
Welcome to the uncertainties package
====================================

The `uncertainties package`_ is a free, cross-platform program that
**transparently** handles calculations with **numbers with uncertainties**
(like 3.14±0.01).  It can also yield the **derivatives** of any
expression.

The :mod:`uncertainties` package **takes the pain and complexity out**
of uncertainty calculations. Error propagation is not to be feared
anymore!

Calculations of results with uncertainties, or of derivatives, can be
performed either in an **interactive session** (as with a calculator),
or in **programs** written in the Python_ programming language.
Existing calculation code can **run with little or no change**.

Whatever the complexity of a calculation, this package returns its
result with an uncertainty as predicted by linear `error propagation
theory`_. It automatically :ref:`calculates derivatives <derivatives>`
and uses them for calculating uncertainties. Almost all uncertainty
calculations are performed **analytically**.

**Correlations** between variables are automatically handled, which
sets this module apart from many existing error propagation codes.

You may want to check the following related uncertainty calculation
Python packages to see if they better suit your needs: soerp_
(higher-order approximations) and mcerp_ (Monte-Carlo approach).

.. index:: calculator

An easy-to-use calculator
=========================

Calculations involving **numbers with uncertainties** can be performed
even without knowing anything about the Python_ programming language.
After `installing this package`_ and `invoking the Python interpreter`_,
calculations with **automatic error propagation** can be performed
**transparently** (i.e., through the usual syntax for mathematical
formulas):

>>> from uncertainties import ufloat
>>> from uncertainties.umath import *  # sin(), etc.
>>> x = ufloat(1, 0.1)  # x = 1+/-0.1
>>> print 2*x
2.00+/-0.20
>>> sin(2*x)  # In a Python shell, "print" is optional
0.9092974268256817+/-0.08322936730942848

Thus, existing calculation code designed for regular numbers can run
with numbers with uncertainties with :ref:`no or little modification
<user guide>`.

.. index:: correlations; simple example

Another strength of this package is its correct handling of
**correlations**.  For instance, the following quantity is exactly
zero even though :data:`x` has an uncertainty:

>>> x-x
0.0+/-0

Many other error propagation codes return the incorrect value 0±0.1414…
because they wrongly assume that the two subtracted quantities are
*independent* random variables.

**Arrays** of numbers with uncertainties are :ref:`transparently
handled <simple_array_use>` too.


**Derivatives** are similarly very :ref:`easy to obtain <derivatives>`:

>>> (2*x+1000).derivatives[x]
2.0

They are calculated with a :ref:`fast method <differentiation method>`.

Available documentation
=======================

The :doc:`user_guide` details many of the features of this package.

The part :doc:`numpy_guide` describes how arrays of numbers with
uncertainties can be created and used.

The :doc:`tech_guide` gives advanced technical details.

.. only:: html

   A :download:`PDF version <_build/latex/uncertaintiesPythonPackage.pdf>`
   of the documentation is also available.

Additional information is available through the pydoc_ command, which
gives access to many of the documentation strings included in the code.

.. index:: installation

.. _installing this package:

Installation and download
=========================

Important note
--------------

The installation commands below should be **run in a DOS or Unix
command shell** (*not* in a Python shell).

Under Windows (version 7 and earlier), a command shell can be obtained
by running ``cmd.exe`` (through the Run… menu item from the Start
menu). Under Unix (Linux, Mac OS X,…), a Unix shell is available when
opening a terminal (in Mac OS X, the Terminal program is found in the
Utilities folder, which can be accessed through the Go menu in the
Finder).

Automatic install or upgrade
----------------------------

One of the automatic installation or upgrade procedures below might work
on your system, if you have a Python package installer or use certain
Linux distributions.

Under Unix, it may be necessary to prefix the commands below with
``sudo``, so that the installation program has **sufficient access
rights to the system**.

If you have `pip <http://pip.openplans.org/>`_, you can try to install
the latest version with

.. code-block:: sh

   pip install --upgrade uncertainties

If you have setuptools_, you can try to automatically install or
upgrade this package with

.. code-block:: sh

   easy_install --upgrade uncertainties

The :mod:`uncertainties` package is also available for **Windows**
through the `Python(x,y)`_ distribution. It may also be included in
Christoph Gohlke's Base distribution of `scientific Python packages`_.

**Mac OS X** users who use the `MacPorts package manager
<http://www.macports.org/>`_ can install :mod:`uncertainties` with
``sudo port install py**-uncertainties``, and upgrade it with ``sudo
port upgrade py**-uncertainties`` where ``**`` represents the desired
Python version (``27``, ``33``, etc.).

The :mod:`uncertainties` package is also available through the
following **Linux** distributions and software platforms: `Ubuntu
<https://launchpad.net/ubuntu/+source/uncertainties>`_, `Fedora
<http://pkgs.org/fedora-18/rpm-sphere-i586/python-uncertainties-1.8.dev418-4.1.noarch.rpm.html>`_,
`openSUSE
<https://build.opensuse.org/package/show?package=python-uncertainties&project=home%3Aocefpaf>`_,
`Debian
<http://packages.debian.org/search?keywords=python-uncertainties>`_
and `Maemo <http://maemo.org/packages/view/python-uncertainties/>`_.


Manual download and install
---------------------------

Alternatively, you can simply download_ the package archive from the
Python Package Index (PyPI) and unpack it.  The package can then be
installed by **going into the unpacked directory**
(:file:`uncertainties-…`), and running the provided :file:`setup.py`
program with

.. code-block:: sh

   python setup.py install

(where the default ``python`` interpreter must generally be replaced
by the version of Python for which the package should be installed:
``python3``, ``python3.3``, etc.).

For an installation with Python 2.6+ in the *user* Python library
(no additional access rights needed):

.. code-block:: sh

   python setup.py install --user

For an installation in a custom directory :file:`my_directory`:

.. code-block:: sh

   python setup.py install --install-lib my_directory

If additional access rights are needed (Unix):

.. code-block:: sh

   sudo python setup.py install

You can also simply **move** the :file:`uncertainties-py*` directory
that corresponds best to your version of Python to a location that
Python can import from (directory in which scripts using
:mod:`uncertainties` are run, etc.); the chosen
:file:`uncertainties-py*` directory should then be renamed
:file:`uncertainties`. Python 3 users should then run ``2to3 -w .``
from inside this directory so as to automatically adapt the code to
Python 3.

Source code
-----------

The latest, bleeding-edge but working `code
<https://github.com/lebigot/uncertainties/tree/master/uncertainties>`_
and `documentation source
<https://github.com/lebigot/uncertainties/tree/master/doc/>`_ are
available `on GitHub <https://github.com/lebigot/uncertainties/>`_.
The :mod:`uncertainties` package is written in pure Python and has no
external dependency (the `NumPy`_ package is optional).  It contains
about 7000 lines of code.  75 % of these lines are documentation
strings and comments.  The remaining 25 % are split between unit tests
(15 % of the total) and the calculation code proper (10 % of the
total).  :mod:`uncertainties` is thus a **lightweight, portable
package** with abundant documentation and tests.


Migration from version 1 to version 2
=====================================

Some **incompatible changes** were introduced in version 2 of
:mod:`uncertainties` (see the `version history`_). While the version 2
line will support the version 1 syntax for some time, it is
recommended to **update existing programs** as soon as possible. This
can be made easier through the provided **automatic updater**.

The automatic updater works like Python's `2to3
<http://docs.python.org/2/library/2to3.html>`_ updater. It can be run
(in a Unix or DOS shell) with:

.. code-block:: sh

   python -m uncertainties.1to2

For example, updating a single Python program can be done with

.. code-block:: sh

   python -m uncertainties.1to2 -w example.py

All the Python programs contained under a directory ``Programs``
(including in nested sub-directories) can be automatically updated
with

.. code-block:: sh

   python -m uncertainties.1to2 -w Programs

Backups are automatically created, unless the ``-n`` option is given.

Some **manual adjustments** might be necessary after running the
updater (incorrectly modified lines, untouched obsolete syntax).

While the updater creates backup copies by default, it is generally
useful to **first create a backup** of the modified directory, or
alternatively to use some `version control
<http://en.wikipedia.org/wiki/Version_control_system>`_
system. Reviewing the modifications with a `file comparison tool
<http://en.wikipedia.org/wiki/File_comparison>`_ might also be useful.

What others say
===============

- "*Superb,*" "*wonderful,*" "*It's like magic.*" (`Joaquin Abian
  <http://blog.garlicsim.org/post/1266209646/cool-python-module-uncertainties#comment-85154147>`_)
- "*pretty amazing*" (`John Kitchin <http://kitchingroup.cheme.cmu.edu/blog/2013/03/07/Another-approach-to-error-propagation/>`_)
- "*An awesome python package*" (`Jason Moore
  <http://biosport.ucdavis.edu/blog/2010/05/07/uncertainty-analysis>`_)
- "*Utterly brilliant*" (`Jeffrey Simpson
  <http://twitter.com/#!/GeekyJeffrey>`_)
- "*An amazing time saver*" (`Paul Nakroshis
  <http://scipyscriptrepo.com/wp/?p=41>`_)
- "*Seems to be the gold standard for this kind of thing*" (`Peter Williams
  <http://newton.cx/~peter/work/?p=660>`_)
- "*This package has a great interface and makes error propagation
  something to stop fearing.*" (`Dr Dawes
  <http://dawes.wordpress.com/2011/01/02/scientific-python/>`_)
- "*uncertainties makes error propagation dead simple.*" (`enrico
  documentation <http://readthedocs.org/docs/enrico/en/latest/setup.html>`_)
- "*many inspiring ideas*" (`Abraham Lee
  <https://pypi.python.org/pypi/soerp#acknowledgements>`_)
- "*Those of us working with experimental data or simulation results
  will appreciate this.*" (`Konrad Hinsen
  <http://khinsen.wordpress.com/2010/07/12/euroscipy-2010/>`_)
- "*PyPI\'s uncertainties rocks!*" (`Siegfried Gevatter
  <http://identi.ca/notice/23330742>`_)
- "*A very cool Python module*" (`Ram Rachum
  <http://blog.garlicsim.org/post/1266209646/cool-python-module-uncertainties>`_)
- "*Holy f\*\*\* this would have saved me so much f\*\*\*ing time last
  semester*." (`reddit
  <http://www.reddit.com/r/Python/comments/am84v/now_you_can_do_calculations_with_uncertainties_5/>`_)


Future developments
===================

Planned future developments include (starting from the most requested
ones):

- handling of complex numbers with uncertainties;

- increased support for `NumPy`_: Fourier Transform with
  uncertainties, automatic wrapping of functions that accept or
  produce arrays, standard deviation of arrays, more convenient matrix
  creation, new linear algebra methods (eigenvalue and QR
  decompositions, determinant,…), input of arrays with uncertainties
  as strings (like in NumPy),…;

- `JSON <http://docs.python.org/library/json.html>`_ support;

- addition of :attr:`real` and :attr:`imag` attributes, for increased
  compatibility with existing code (Python numbers have these attributes);

- addition of new functions from the :mod:`math` module;

- fitting routines that conveniently handle data with uncertainties;

- a re-correlate function that puts correlations back between data
  that was saved in separate files;

- support for multi-precision numbers with uncertainties.

**Call for contributions**: I got multiple requests for complex
numbers with uncertainties, Fourier Transform support, and the
automatic wrapping of functions that accept or produce arrays. Please
contact me if you are interested in contributing. Patches are
welcome. They must have a high standard of legibility and quality in
order to be accepted (otherwise it is always possible to create a new
Python package by branching off this one, and I would still be happy
to help with the effort).

**Please support the continued development of this program** by
`donating $10`_ or more through PayPal (no PayPal account
necessary). I love modern board games, so this will go towards giving
my friends and I some special gaming time!

.. index:: support

Contact
=======

**Feature requests, bug reports, or feedback are much welcome.** They
can be sent_ to the creator of :mod:`uncertainties`, `Eric O. LEBIGOT
(EOL)`_.

.. figure:: _static/eol.*
   :height: 64
   :width:  64
   :target: http://linkedin.com/pub/eric-lebigot/22/293/277
   :align: center
   :alt: Eric O. LEBIGOT (EOL)

How to cite this package
========================

If you use this package for a publication (in a journal, on the web,
etc.), please cite it by including as much information as possible
from the following: *Uncertainties: a Python package for calculations
with uncertainties*, Eric O. LEBIGOT,
`<http://pythonhosted.org/uncertainties/>`_.  Adding the version
number is optional.


Acknowledgments
===============

The author wishes to thank all the people who made generous
`donations`_: they help keep this project alive by providing positive
feedback.

I greatly appreciated getting key technical input from Arnaud
Delobelle, Pierre Cladé, and Sebastian Walter.  Patches by Pierre
Cladé, Tim Head, José Sabater Montes, Martijn Pieters, Ram Rachum,
Christoph Deil, and Gabi Davar are gratefully acknowledged.

I would also like to thank users who contributed with feedback and
suggestions, which greatly helped improve this program: Joaquin Abian,
Jason Moore, Martin Lutz, Víctor Terrón, Matt Newville, Matthew Peel,
Don Peterson, Mika Pflueger, Albert Puig, Abraham Lee, Arian Sanusi,
Martin Laloux, Jonathan Whitmore, Federico Vaggi, Marco A. Ferra,
Hernan Grecco, David Zwicker, and many others.

I am also grateful to Gabi Davar and Pierre Raybaut for including it
in `Python(x,y)`_, to Christoph Gohlke for including it in his Base
distribution of `scientific Python packages`_ for Windows, and to the
Mac OS X and Linux distribution maintainers of this package (Jonathan
Stickel, David Paleino, Federico Ceratto, Roberto Colistete Jr, and
Filipe Pires Alvarenga Fernandes).

.. index:: license

License
=======

This software is released under a **dual license**; one of the
following options can be chosen:

1. The `Revised BSD License`_ (© 2010–2016, Eric O. LEBIGOT [EOL]).
2. Any other license, as long as it is obtained from the creator of
   this package.

.. _Python: http://python.org/
.. _Python(x,y): https://code.google.com/p/pythonxy/
.. _scientific Python packages: http://www.lfd.uci.edu/~gohlke/pythonlibs/
.. _error propagation theory: http://en.wikipedia.org/wiki/Propagation_of_uncertainty
.. _invoking the Python interpreter: http://docs.python.org/tutorial/interpreter.html
.. _setuptools: http://pypi.python.org/pypi/setuptools
.. _download: http://pypi.python.org/pypi/uncertainties/#downloads
.. _donations: https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=4TK7KNDTEDT4S
.. _Eric O. LEBIGOT (EOL): http://linkedin.com/pub/eric-lebigot/22/293/277
.. _sent: mailto:eric.lebigot@normalesup.org
.. _Revised BSD License: http://opensource.org/licenses/BSD-3-Clause
.. _uncertainties package: http://pypi.python.org/pypi/uncertainties/
.. _pydoc: http://docs.python.org/library/pydoc.html
.. _NumPy: http://numpy.scipy.org/
.. _donating $10: donations_
.. _version history: https://pypi.python.org/pypi/uncertainties#version-history
.. _soerp: https://pypi.python.org/pypi/soerp
.. _mcerp: https://pypi.python.org/pypi/mcerp
.. _Pint: https://pypi.python.org/pypi/Pint/
