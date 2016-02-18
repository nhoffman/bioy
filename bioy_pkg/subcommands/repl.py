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

"""
Replace strings in one or more files.

Yes, all of this to replace a sed one-liner. Couldn't figure out how
to make a command the worked for both GNU and BSD sed.
"""

import fileinput
import sys
import re
import logging

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('infiles',
            nargs='+', help='Input files')
    parser.add_argument('-r','--replace',
            action = 'append',
            help="""String of the form 'pattern:repl' -
            each instance of pattern will be replaced with repl.
            -r can be used multiple times.""")
    parser.add_argument('-n','--dry-run',
            action = 'store_false', default = True, dest = 'inplace',
            help = 'Print replaced lines to stdout without modifying files.')

def action(args):
    replacements = [val.split(':') for val in args.replace]
    for line in fileinput.input(args.infiles, inplace = args.inplace):
        for pattern, repl in replacements:
            line = line.replace(pattern, repl)
        sys.stdout.write(line)
