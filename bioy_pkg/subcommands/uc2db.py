"""
(deprecated) Parses usearch_global blast results and inserts into specified database

Use csvsql --db sqlite:///path/to/database.db --insert /path/to/csv/file.csv
http://csvkit.readthedocs.org/en/latest/
"""

import sys
import sqlite3
import argparse
import csv

from romperroom import uclust

def build_parser(parser):
    parser.add_argument('infile', 
            nargs = '?',
            default = sys.stdin,
            help = "A uclust formatted input file")
    parser.add_argument('-d', '--db', 
            help = "sqlite3 database",
            default = 'species.db')
    parser.add_argument('-i', '--taxinfo', 
            help = 'A csv format of taxonomy info to retrieve species names',
            type = argparse.FileType('r'),
            required = True)
    parser.add_argument('-t', '--threshold', 
            type = int)
    parser.add_argument('-n', '--table-name',
            help = 'A database table name',
            default = 'uclust')

def action(args):
    #"seqname","tax_id","accession","description","length","ambig_count","is_type","rdp_lineage"
    #"S000438419","53635","U75647","Acidimicrobium ferrooxidans (T); ICP (DSM 10331).","1465","0",
    #"FALSE","Root;Bacteria;Actinobacteria;Actinobacteria;Acidimicrobidae;
    #         Acidimicrobiales;Acidimicrobineae;Acidimicrobiaceae;Acidimicrobium"
    taxinfo = {row['seqname']:row for row in csv.DictReader(args.taxinfo)}

    con = sqlite3.connect(args.db)
    cur = con.cursor()

    cur.execute("""
            create table if not exists %s (
            q_name text, 
            cluster_size integer, 
            t_name text, 
            tax_id text, 
            species_name text, 
            pct_id real, 
            threshold integer)""" % args.table_name)

    for col in ['q_name', 'threshold']:
        cur.execute('create index if not exists %s_%s on %s(%s)' % (args.table_name, col, args.table_name, col))

    cmd = """insert into %s values (?, ?, ?, ?, ?, ?, ?)""" % args.table_name

    for a in uclust.parse_uclust_out(args.infile):
        if a.type == 'H':
            cur.execute(cmd, ([a.query_label,
                a.query_label.split('|')[1],
                a.target_label,
                taxinfo[a.target_label]['tax_id'],
                taxinfo[a.target_label]['description'],
                a.pct_id,
                args.threshold,
                ]))
            con.commit()
    con.close()

