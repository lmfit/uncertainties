.. index:: installation
.. index:: credits
.. _installation:

====================================
Installation and Credits
====================================

Download and Installation
=========================

The :mod:`uncertainties` package supports Python versions 3.8 and higher.
Earlier versions of Python are not tested, but may still work.  Development
version of Python (currently, 3.13) are likely to work, but are not regularly
tested.

To install :mod:`uncertainties`, use:

.. code-block:: sh

    pip install uncertainties

You can upgrade from an older version of :mod:`uncertainties` with:

.. code-block:: sh

  pip install --upgrade uncertainties


Other packaging systems such as `Anaconda <https://www.anaconda.com>`_,
`MacPorts <http://www.macports.org/>`_, or Linux package manager may also
maintain packages for :mod:`uncertainties`, so that you may also be able to
install using something like

.. code-block:: sh

   conda install -c conda-forge uncertainties

.. code-block:: sh

  sudo port install py**-uncertainties

or

.. code-block:: sh

  sudo apt get python-uncertainties

depending on your platform and installation of Python.  For all installations
of Python, using `pip` should work and is therefore recommended.


Source code and Development Version
=====================================

.. _download:  https://pypi.python.org/pypi/uncertainties/#files
.. _GitHub releases: https://github.com/lmfit/uncertainties/releases
.. _NumPy: http://numpy.scipy.org/

You can `download`_ the latest source package archive from the Python Package
Index (PyPI) and unpack it, or from the `GitHub releases`_ page.  This package
can be unpacked using `unzip`, `tar xf` , or other similar utilities, and then
installed with

.. code-block:: sh

   python -m pip install .

To work with the development version, use `git` to fork or clone the code:

.. code-block:: sh

   git clone git@github.com:lmfit/uncertainties.git

The :mod:`uncertainties` package is written in pure Python and has no external
dependencies.  If available (and recommended), the `NumPy`_ package can be
used.  Running the test suite requires `pytest` and `pytest_cov`, and building
these docs requires `sphinx`.  To install these optional packages, use one of:

.. code-block:: sh

    pip install ".[arrays]"    # to install numpy
    pip install ".[test]"      # to enable running the tests
    pip install ".[doc]"       # to enable building the docs
    pip install ".[all]"       # to enable all of these options

Getting Help
=================

.. _GitHub Discussions: https://github.com/lmfit/uncertainties/discussions
.. _GitHub Issues: https://github.com/lmfit/uncertainties/issues
.. _lmfit GitHub organization: https://github.com/lmfit/

If you have questions about :mod:`uncertainties` or run into trouble, use the
`GitHub Discussions`_ page.  For bug reports, use the `GitHub Issues`_ pages.


Credits
================

.. _Eric O. LEBIGOT (EOL): http://linkedin.com/pub/eric-lebigot/22/293/277

The :mod:`uncertainties` package was written and developed by `Eric O. LEBIGOT
(EOL)`_.  EOL also maintained the package until 2024, when the GitHub project
was moved to the `lmfit GitHub organization`_ to allow more sustainable
development and maintenance.  Current members of the devlopment and
maintenance team include `Andrew G Savage <https://github.com/andrewgsavage>`_,
`Justin Gerber <https://github.com/jagerber48>`_,
`Eric O Legibot <https://github.com/lebigot>`_,
`Matt Newville <https://github.com/newville>`_,
and `Will Shanks <https://github.com/wshanks>`_.  Contributions and suggestions
for development are welcome.


How to cite this package
========================

If you use this package for a publication, please cite it as *Uncertainties: a
Python package for calculations with uncertainties*, Eric O. LEBIGOT.  A
version number can be added, but is optional.


Acknowledgments
===============

.. _Python(x,y): https://python-xy.github.io/
.. _scientific Python packages: http://www.lfd.uci.edu/~gohlke/pythonlibs/

Eric O. LEBIGOT (EOL) thanks all the people who made generous donations: that
help to keep this project alive by providing positive feedback.

EOL greatly appreciates having gotten key technical input from Arnaud Delobelle,
Pierre Cladé, and Sebastian Walter.  Patches by Pierre Cladé, Tim Head, José
Sabater Montes, Martijn Pieters, Ram Rachum, Christoph Deil, Gabi Davar, Roman
Yurchak and Paul Romano are gratefully acknowledged.

EOL also thanks users who contributed with feedback and
suggestions, which greatly helped improve this program: Joaquin Abian,
Jason Moore, Martin Lutz, Víctor Terrón, Matt Newville, Matthew Peel,
Don Peterson, Mika Pflueger, Albert Puig, Abraham Lee, Arian Sanusi,
Martin Laloux, Jonathan Whitmore, Federico Vaggi, Marco A. Ferra,
Hernan Grecco, David Zwicker, James Hester, Andrew Nelson, and many others.

EOL is grateful to the Anaconda, macOS and Linux distribution maintainers
of this package (Jonathan Stickel, David Paleino, Federico Ceratto,
Roberto Colistete Jr, Filipe Pires Alvarenga Fernandes, and Felix Yan)
and also to Gabi Davar and Pierre Raybaut for including it in
`Python(x,y)`_ and to Christoph Gohlke for including it in his Base
distribution of `scientific Python packages`_ for Windows.

.. index:: license

License
=======

.. _Revised BSD License: http://opensource.org/licenses/BSD-3-Clause

This software is released under the  `Revised BSD License`_ (© 2010–2024,
Eric O. LEBIGOT [EOL]).
