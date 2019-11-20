#!/bin/sh

# This script must be run before packaging (twine upload dist/*).

# WARNING: this script erases any uncertainties-py23 and uncertainties-py27
# found in the current directory.

## Only committed versions are packaged, to help with debugging published code:
git checkout release
git commit -a

# The Python 2.3-2.5 version should always be up to date:
git checkout release_python2.3
git merge release

# Default branch for working on the code:
git checkout release

## Getting the Python 2.7+ version:

rm -rf uncertainties-py27 && \
git archive --output /tmp/u.tar master uncertainties && \
tar -C /tmp -xf /tmp/u.tar && \
mv /tmp/uncertainties uncertainties-py27 && \
echo "Python 2.7+ version imported"

## Getting the Python 2.3 version:

rm -rf uncertainties-py23 && \
git archive --output /tmp/u.tar release_python2.3 uncertainties && \
tar -C /tmp -xf /tmp/u.tar && \
mv /tmp/uncertainties uncertainties-py23 && \
echo "Python 2.3 version imported"

# Packaging:
python setup.py sdist bdist_wheel && \
echo "Package created.  The package can be uploaded with twine upload dist/*."
echo "WARNING: current git branch is:"
git branch | grep '^\*'
