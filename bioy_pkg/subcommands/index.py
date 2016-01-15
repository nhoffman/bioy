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
Add simple indices to an sqlite database
"""

import logging
import sys
import sqlite3

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('db',
                        help = 'An sqlite database')
    parser.add_argument('colnames',
                        help = 'Comma-delimited list of column names')
    parser.add_argument('-t', '--tables',
                        help = ('One or more table names (comma-delimited); if not provided, '
                                'attempts to add an index to any '
                                'table with a column in --colnames'))
    parser.add_argument('-c', '--clobber', action='store_true', default=False,
                        help='drop specified indices before attempted creation')

def action(args):

    con = sqlite3.connect(args.db)
    cur = con.cursor()

    colnames = args.colnames.split(',')

    if args.tables:
        tables = args.tables.split(',')
    else:
        cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = [x[0] for x in cur.fetchall()]

    for table in tables:
        for column in colnames:
            index = '{}_{}'.format(table, column)

            if args.clobber:
                cmd = 'DROP INDEX IF EXISTS "{}"'.format(index)
                log.info(cmd)
                cur.execute(cmd)

            cmd = 'CREATE INDEX IF NOT EXISTS {index} ON {table} ({column})'.format(
                index=index, table=table, column=column)
            log.info(cmd)

            try:
                cur.execute(cmd)
            except sqlite3.OperationalError, e:
                log.info(e)
