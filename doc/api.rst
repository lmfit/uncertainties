.. index:: technical details

=========================
Advanced Topics
=========================


This page gives more in-depth technical description of the
:mod:`uncertainties` package.

 .. index:: api

.. _api_funcs:


API: Application Programming Interface
==============================================

.. module:: uncertainties

The most important object is the :class:`UFloat` object.
The :class:`UFloat` object can also be constructed using the :class:`ufloat` alias.
The :func:`ufloat_from_str` function is an alternate constructor.
The :func:`wrap` function wraps user-provided functions to support :class:`UFloat`
objects with uncertainty propagation.

.. autoclass:: UFloat

.. autofunction:: ufloat_fromstr

.. autofunction:: wrap

Special Technical Topics
============================================================
