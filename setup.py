"""
Setup.py template. Try this:

    sed 's/bioy_pkg/newpackagename/g;s/bioy/newscriptname/g' setup.py
"""

import subprocess
from setuptools import setup, find_packages
from os.path import join

subprocess.call('git log --pretty=format:%h -n 1 > bioy_pkg/data/sha', shell = True)
subprocess.call('git shortlog --format="XXYYXX%h" | grep -c XXYYXX > bioy_pkg/data/ver', shell = True)

from bioy_pkg import __version__

params = {'author': 'Noah Hoffman',
          'author_email': 'ngh2@uw.edu',
          'description': 'A collection of bioinformatics tools',
          'name': 'bioy_pkg',
          'packages': find_packages(),
          'package_dir': {'bioy_pkg': 'bioy_pkg'},
          'scripts': ['bioy'],
          'version': __version__,
          'package_data': {'bioy_pkg': [join('data',f) for f in ['sha','ver']]}
          }

setup(**params)

