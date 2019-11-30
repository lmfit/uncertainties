#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
from setuptools import setup

PY_VER = sys.version_info
if PY_VER.major < 3:
    min_version = (2, 7)
else:
    min_version = (3, 5)


error_msg = ("Sorry, this package is for Python %d.%d and higher only." %
             min_version)

try:
    if PY_VER < min_version:
        sys.exit(error_msg)
except AttributeError:  # sys.version_info was introduced in Python 2.0
    sys.exit(error_msg)

# Common options for distutils/setuptools's setup():
setup_options = dict(
    name='uncertainties',
    version='3.1.2',
    author='Eric O. LEBIGOT (EOL)',
    author_email='eric.lebigot@normalesup.org',
    url='http://uncertainties-python-package.readthedocs.io/',
    license='''\
This software can be used under one of the following two licenses: \
(1) The Revised BSD License. \
(2) Any other license, as long as it is obtained from the original \
author.''',
    description=('Transparent calculations with uncertainties on the'
                 ' quantities involved (aka error propagation);'
                 ' fast calculation of derivatives'),
    long_description=open('README.rst').read(),
    keywords=[
        'error propagation', 'uncertainties', 'uncertainty calculations',
        'standard deviation', 'derivatives', 'partial derivatives',
        'differentiation'
    ],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Other Audience',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: Implementation :: Jython',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Topic :: Education',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Software Development',
        'Topic :: Software Development :: Libraries',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Utilities'
    ],

    # Files are defined in MANIFEST (which is automatically created by
    # python setup.py sdist bdist_wheel):
    packages=[
        'uncertainties', 'uncertainties.unumpy', 'uncertainties.lib1to2',
        'uncertainties.lib1to2.fixes'])

if PY_VER.major < 3:
    setup_options["install_requires"] = ['future']

# Some setuptools-specific options can be added:
addtl_setup_options = {
    # Allows python setup.py test to do the right thing:
    'tests_require': ['pytest', 'numpy'],
    # Optional dependencies install using:
    # `easy_install uncertainties[optional]`
    'extras_require': {
        'optional': ['numpy'],
        'docs': ['sphinx'],
    }
}

# easy_install uncertainties[tests] option:
addtl_setup_options['extras_require']['tests'] = (
    addtl_setup_options['tests_require'])

# easy_install uncertainties[all] option: all dependencies are
# gathered
addtl_setup_options['extras_require']['all'] = set(
    sum(list(addtl_setup_options['extras_require'].values()), []))

setup_options.update(addtl_setup_options)


# End of setup definition

setup(**setup_options)
