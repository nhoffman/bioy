import subprocess

from setuptools import setup, find_packages

try:
    subprocess.check_call(
        'git describe --tags --long > bioy_pkg/data/version',
        shell = True)
except subprocess.CalledProcessError:
    version = ''
else:
    with open('bioy_pkg/data/version') as f:
        version = f.read().strip()

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
          'pandas==0.13.1',
          'biopython>=1.6.3'
      ])
