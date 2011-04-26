#!/bin/sh

# This script must be run before packaging (python setup.py sdist upload).

# WARNING: this script erases any uncertainties-py23 and uncertainties-py25
# found in the current directory.

## Only committed versions are packaged, to help with debugging published code:
git commit -a

## Getting the Python 2.5 version:

rm -rf uncertainties-py25 && \
git archive master uncertainties > /tmp/u.tar && \
tar -C /tmp -xf /tmp/u.tar && \
mv /tmp/uncertainties uncertainties-py25 && \
echo "Python 2.5 version imported"

## Getting the Python 2.3 version:

rm -rf uncertainties-py23 && \
git archive python2.3 uncertainties > /tmp/u.tar && \
tar -C /tmp -xf /tmp/u.tar && \
mv /tmp/uncertainties uncertainties-py23 && \
echo "Python 2.3 version imported"

# Packaging:
python setup.py sdist && \
echo "Package created.  The package can be uploaded with setup.py sdist upload."
