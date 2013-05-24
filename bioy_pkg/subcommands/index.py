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
