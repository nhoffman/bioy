import subprocess

from setuptools import setup, find_packages

try:
    subprocess.check_call(
        'git describe --tags --long > bioy_pkg/data/version',
        shell=True)
except subprocess.CalledProcessError:
    version = ''
else:
    with open('bioy_pkg/data/version') as f:
        version = f.read().strip().split('-')[0]

setup(author='Noah Hoffman',
      author_email='noah.hoffman@gmail.com',
      description='A collection of bioinformatics tools',
      name='bioy',
      packages=find_packages(),
      scripts=['bioy'],
      version=version,
      url='https://github.com/nhoffman/bioy',
      package_data={'bioy_pkg': ['data/*']},
      requires=['python (>= 2.7.5)'],
      install_requires=[
          'numpy>=1.8.1',
          'pandas>=0.17.1',
          'biopython>=1.6.3',
          'matplotlib'
      ])
