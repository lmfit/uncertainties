#!/bin/sh

# This script prepares the package for PyPI. It must be run
# before uploading it on PyPI.

# This script must be run from its directory.

# Fail the script at the first failed command (HOWEVER, maybe when there are
# no commits to be done during the merges, the commands fail?):
#set -e 

echo "****************************************************************"
echo "WARNING: if any commit fails, RESOLVE IT before running this"
echo "script again. Otherwise conflict marks will be committed by the"
echo "second run!"
echo "****************************************************************"

## Only committed versions are packaged, to help with debugging published code:
git commit -a

# We make sure that the release and master branches are merged (changes
# may have been made on both sides):
git checkout master
git merge release

git checkout release
git merge master

# Default branch for working on the code:
git checkout release

# Packaging. We include wheels because it makes it easier to install,
# in some cases (https://github.com/lebigot/uncertainties/pull/108,
# https://discourse.slicer.org/t/problems-installing-lmfit-python-package/9210/6):
python setup.py sdist bdist_wheel
echo "Package created.  The package can be uploaded with twine upload dist/...*"
echo "where ...* is the new versions."
echo "WARNING: current git branch is:"
git branch | grep '^\*'
