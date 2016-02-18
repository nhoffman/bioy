# This file is part of Bioy
#
#    Bioy is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Bioy is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Bioy.  If not, see <http://www.gnu.org/licenses/>.

import os
import bz2
import gzip
import logging
import pandas
import re
import shutil
import sys
import contextlib
import tempfile

from itertools import takewhile, izip_longest, groupby
from csv import DictReader
from collections import Iterable, OrderedDict
from os import path

log = logging.getLogger(__name__)


def apply_df_status(func, df, msg=''):
    """
    """
    tmp_column = 'index_number'
    row_count = float(len(df))
    df[tmp_column] = xrange(int(row_count))
    msg += ' {:.0%}\r'

    def apply_func(item, msg):
        sys.stderr.write(msg.format(item[tmp_column] / row_count))
        return func(item)

    df = df.apply(apply_func, args=[msg], axis=1)
    return df.drop(tmp_column, axis=1)


def flattener(iterable):
    """
    Flatten nested iterables (not strings or dict-like objects).

    Poached from http://stackoverflow.com/questions/2158395
                       /flatten-an-irregular-list-of-lists-in-python
    """
    for el in iterable:
        if isinstance(el, Iterable) and \
                not (isinstance(el, basestring) or hasattr(el, 'get')):
            for sub in flattener(el):
                yield sub
        else:
            yield el


def chunker(seq, size, combine_last=None):
    """
    Break sequence seq into lists of length `size`. If the length of
    the final list is < 'combine_last', it is appended to the end of
    the penultimate element.
    """

    chunks = [seq[pos:pos + size] for pos in xrange(0, len(seq), size)]
    if combine_last and len(chunks[-1]) < combine_last:
        chunks[-2].extend(chunks.pop(-1))

    return iter(chunks)


def grouper(n, iterable, pad=True):
    """
    Return sequence of n-tuples composed of successive elements of
    iterable; last tuple is padded with None if necessary. Not safe
    for iterables with None elements.
    """

    args = [iter(iterable)] * n
    iterout = izip_longest(fillvalue=None, *args)

    if pad:
        return iterout
    else:
        return (takewhile(lambda x: x is not None, c) for c in iterout)


def cast(val):
    for func in [int, float, lambda x: x.strip()]:
        try:
            return func(val)
        except ValueError:
            pass


def mkdir(dirpath, clobber=False):
    """
    Create a (potentially existing) directory without errors. Raise
    OSError if directory can't be created. If clobber is True, remove
    dirpath if it exists.
    """

    if clobber:
        shutil.rmtree(dirpath, ignore_errors=True)

    try:
        os.mkdir(dirpath)
    except OSError:
        pass

    if not path.exists(dirpath):
        raise OSError('Failed to create %s' % dirpath)

    return dirpath


def parse_extras(s, numeric=True):
    """
    Return an OrderedDict parsed from a string in the format
    "key1:val1,key2:val2"
    """

    # allow for escaped quoted text
    commas = re.compile(r"""((?:[^,"']|"[^"]*"|'[^']*')+)""")
    colons = re.compile(r"""((?:[^:"']|"[^"]*"|'[^']*')+)""")

    extras = commas.split(s)[1::2]
    extras = (colons.split(e)[1::2] for e in extras)
    extras = ((k, cast(v) if numeric else v) for k, v in extras)
    extras = OrderedDict(extras)

    return extras


class Opener(object):

    """Factory for creating file objects

    Keyword Arguments:
        - mode -- A string indicating how the file is to be opened. Accepts the
            same values as the builtin open() function.
        - bufsize -- The file's desired buffer size. Accepts the same values as
            the builtin open() function.
    """

    def __init__(self, mode='r', bufsize=-1):
        self._mode = mode
        self._bufsize = bufsize

    def __call__(self, string):
        if string is sys.stdout or string is sys.stdin:
            return string
        elif string == '-':
            return sys.stdin if 'r' in self._mode else sys.stdout
        elif string.endswith('.bz2'):
            return bz2.BZ2File(string, self._mode, self._bufsize)
        elif string.endswith('.gz'):
            return gzip.open(string, self._mode, self._bufsize)
        else:
            return open(string, self._mode, self._bufsize)

    def __repr__(self):
        args = self._mode, self._bufsize
        args_str = ', '.join(repr(arg) for arg in args if arg != -1)
        return '{}({})'.format(type(self).__name__, args_str)


def opener(pth, mode='r', bufsize=-1):
    return Opener(mode, bufsize)(pth)


class Csv2Dict(object):

    """Easy way to convert a csv file into a
    dictionary using the argparse type function

    If no arguments the first column of the csv
    will be the key and every column
    will be the value in an OrderedDict.

    Keyword Arguments:
        - index -- csv column to key index the dictionary
        - value -- csv column to value the dictionary
        - fieldnames -- csv column names
    """

    def __init__(self, index=None, value=None, *args, **kwds):
        self.index = index
        self.value = value
        self.args = args
        self.kwds = kwds

    def __call__(self, pth):
        reader = DictReader(opener(pth), *self.args, **self.kwds)

        if not self.index:
            self.index = reader.fieldnames[0]

        results = {}
        for r in reader:
            key = r[self.index]
            if self.value:
                results[key] = r[self.value]
            else:
                fields = lambda k: reader.fieldnames.index(k[0])
                results[key] = OrderedDict(sorted(r.items(), key=fields))

        return results


def csv2dict(pth, index=None, value=None, *args, **kwds):
    return Csv2Dict(index, value, args, kwds)(pth)


def groupbyl(li, key=None, as_dict=False):
    groups = sorted(li, key=key)
    groups = groupby(groups, key=key)
    groups = ((g, list(l)) for g, l in groups)
    if as_dict:
        return(dict(groups))
    else:
        return groups


def read_csv(filename, compression=None, limit=None, **kwargs):
    """Read a csv file using pandas.read_csv with compression defined by
    the file suffix unless provided.
    """

    suffixes = {'.bz2': 'bz2', '.gz': 'gzip'}
    compression = compression or suffixes.get(path.splitext(filename)[-1])
    kwargs['compression'] = compression

    return pandas.read_csv(filename, **kwargs)


@contextlib.contextmanager
def named_tempfile(*args, **kwargs):
    """Near-clone of tempfile.NamedTemporaryFile, but the file is deleted
    when the context manager exits, rather than when it's closed.

    """

    kwargs['delete'] = False
    tf = tempfile.NamedTemporaryFile(*args, **kwargs)
    try:
        with tf:
            yield tf
    finally:
        os.unlink(tf.name)
