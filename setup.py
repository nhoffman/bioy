import subprocess

from setuptools import setup, find_packages

subprocess.call('git log --pretty=format:%h -n 1 > bioy_pkg/data/sha', shell = True)
subprocess.call('git shortlog --format="XXYYXX%h" | grep -c XXYYXX > bioy_pkg/data/ver', shell = True)

version = subprocess.check_output('git describe --tags', shell = True).strip()

setup(author = 'Noah Hoffman',
      author_email = 'ngh2@uw.edu',
      description = 'A collection of bioinformatics tools',
      name = 'bioy',
      packages = find_packages(),
      scripts = ['bioy'],
      version = version,
      url = 'https://bitbucket.org/uwlabmed/bioy/',
      requires = ['python (>= 2.7.5)'],
      install_requires = [
          'numpy==1.8.1',
          'pandas==0.14.0',
          'biopython>=1.6.3'
      ])
