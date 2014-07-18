"""
Tools for microbial sequence analysis and classification.
"""

from os.path import join, dirname

_data = join(dirname(__file__), 'data')

try:
    with open(join(_data, 'version')) as v:
        __version__ = v.read().strip()
except Exception, e:
    __version__ = ''
