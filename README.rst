==========================================
bioy: a collection of bioinformatics tools
==========================================

Bio-y
    (pronounced "bio-ee") The adjective form of the noun "Bio"

.. contents:: Table of Contents

authors
=======

* Noah Hoffman
* Chris Rosenthal
* Tyler Land

dependencies
============

* Tested primarily on Ubuntu 12.04
* Python 2.7.x
* setuptools

Some functions require

* numpy
* pandas
* biopython

installation
============

Using setup.py or pip from the project directory::

  % cd bioy
  % python setup.py install
  # or
  % pip install -U .

A virtualenv containing a complete execution environment can be
created using `dev/bootstrap.sh`::

  % dev/bootstrap.sh -h
  Create a virtualenv and install all pipeline dependencies
  Options:
  --venv            - path of virtualenv [bioy-env]
  --python          - path to the python interpreter [/usr/local/bin/python]
  --wheelstreet     - path to directory containing python wheels; wheel files will be
  in a subdirectory named according to the python interpreter version;
  uses WHEELSTREET if defined.
  (a suggested location is ~/wheelstreet) []
  --requirements    - a file listing python packages to install [requirements.txt]


execution
=========

The ``bioy`` script provides the user interface, and uses standard
UNIX command line syntax. Note that for development, it is convenient
to run ``bioy`` from within the project directory by specifying the
relative path to the script::

  % ./bioy

Commands are constructed as follows. Every command starts with the
name of the script, followed by an "action" followed by a series of
required or optional "arguments". The name of the script, the action,
and options and their arguments are entered on the command line
separated by spaces. Help text is available for both the ``bioy``
script and individual actions using the ``-h`` or ``--help`` options.

versions
========

We use abbrevited git sha hashes to identify the software version::

  % ./bioy --version
  0128.9790c13

The version information is saved in ``bioy_pkg/data`` when ``setup.py``
is run (on installation, or even by executing ``python setup.py
-h``).

unit tests
==========

Unit tests are implemented using the ``unittest`` module in the Python
standard library. The ``tests`` subdirectory is itself a Python
package that imports the local version (ie, the version in the project
directory, not the version installed to the system) of the
package. All unit tests can be run like this::

    % ./testall
    ...........
    ----------------------------------------------------------------------
    Ran 11 tests in 0.059s

    OK

A single unit test can be run by referring to a specific module,
class, or method within the ``tests`` package using dot notation::

    % ./testone -v tests.test_utils

license
=======

Copyright (c) 2012 Noah Hoffman

Released under the GPLv3 License:

TODO: include license text
