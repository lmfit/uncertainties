* Some installation methods:

  python setup.py install

or, for an installation in the user Python library (no additional access 
rights needed):

  python setup.py install --user

or, for an installation in a custom directory my_directory:

  python setup.py install --install-lib my_directory

or, if additional access rights are needed (Unix):

  sudo python setup.py install

You can also simply move the uncertainties/ or uncertainties-py*/ 
directory to a location that Python can import from (directory in which 
scripts using uncertainties are run, etc.), and then rename it 
uncertainties/ if necessary (with no Python version number).

* The tests programs (test_*.py) are meant to be run through the Nose 
testing framework.  This can be achieved for instance with a command 
like

  nosetests -sv uncertainties/

or simply

  nosetests uncertainties/

(for a less verbose output).
