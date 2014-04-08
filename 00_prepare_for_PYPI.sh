#!/bin/sh

# This script must be run before packaging (python setup.py sdist upload).

# WARNING: this script erases any uncertainties-py23 and uncertainties-py26
# found in the current directory.

## Only committed versions are packaged, to help with debugging published code:
git checkout release
git commit -a

# The Python 2.3-2.5 version should always be up to date:
git checkout release_python2.3
git merge release

# Default branch for working on the code:
git checkout release

## Getting the Python 2.6+ version:

rm -rf uncertainties-py26 && \
git archive --output /tmp/u.tar master uncertainties && \
tar -C /tmp -xf /tmp/u.tar && \
mv /tmp/uncertainties uncertainties-py26 && \
echo "Python 2.6+ version imported"

## Getting the Python 2.3 version:

rm -rf uncertainties-py23 && \
git archive --output /tmp/u.tar release_python2.3 uncertainties && \
tar -C /tmp -xf /tmp/u.tar && \
mv /tmp/uncertainties uncertainties-py23 && \
echo "Python 2.3 version imported"

# Packaging:
python setup.py sdist && \
echo "Package created.  The package can be uploaded with setup.py sdist upload."
echo "WARNING: current git branch is:"
git branch | grep '^\*'