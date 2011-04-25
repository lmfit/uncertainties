# This script must be run before packaging (python setup.py sdist upload).

# WARNING: this script erases any uncertainties-py23 and uncertainties-py25
# found in the current directory.

## Only committed versions are packaged, to help with debugging published code:
git commit -a

## Let's move the original files out of the way:
mv uncertainties uncertainties-orig

## Getting the Python 2.5 version:

rm -rf uncertainties-py25
git checkout master -- uncertainties
mv uncertainties uncertainties-py25
echo "Python 2.5 version imported"

## Getting the Python 2.3 version:

rm -rf uncertainties-py23
git checkout python2.3 -- uncertainties
mv uncertainties uncertainties-py23
echo "Python 2.3 version imported"

# The original directory is put back:
mv uncertainties-orig uncertainties

# Removing files added to the index by the checkouts:
git reset

# Packaging:
python setup.py sdist
echo "Package created.  The package can be uploaded with setup.py upload."
