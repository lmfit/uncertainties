# This script must be run before packaging (python setup.py sdist upload).

# WARNING: this script erases any uncertainties-py23 and uncertainties-py25
# found in the current directory.

rm -rf uncertainties-py25

# Only committed version are packaged, to help with debugging:
git commit -a

# Let's move the original files out of the way:
mv uncertainties uncertainties-py25

# Getting the Python <2.5 version:

# git read-tree -u --prefix=uncertainties-py23/ python-pre-2.5:uncertainties/ && \
rm -rf uncertainties-py23
git checkout python2.3 -- uncertainties
mv uncertainties uncertainties-py23
echo "Python 2.3 version imported"

# Getting the Python 3+ version:
#git co python-3 uncertainties/
#mv uncertainties uncertainties-py30

# Removing files added to the index by the checkout:
git reset

# The original directory is put back:
cp -r uncertainties-py25 uncertainties

# Packaging:
python setup.py sdist
