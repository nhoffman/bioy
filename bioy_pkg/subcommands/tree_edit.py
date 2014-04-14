"""
Tree leaf name editor that wraps BioPython.
"""

import logging
import re
import sys

from Bio import Phylo
from csv import DictReader

from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('tree',
            default = sys.stdin,
            nargs = '?',
            type = Opener('r'),
            help='sequences in fasta format')
    parser.add_argument('--info',
            type = Opener('r'),
            metavar = 'CSV',
            help = """selectively replace leaves from a two column
                     csv [seqname,*newname*] *optional column*""")
    parser.add_argument('-o', '--out', default='/dev/stdout',
            help='output file [%(default)s]')
    parser.add_argument('--add-prefix', default = '',
            metavar = 'PRE',
            help = 'append a prefix string to all names')
    parser.add_argument('--add-suffix', default = '',
            metavar = 'SUF',
            help = 'append a suffix string to all names')
    parser.add_argument('--tree-type', default = 'newick',
            help = 'tree type to parse')
    parser.add_argument('--remove-word', metavar = 'REGEX',
            help = 'remove a word from a ')

def action(args):
    def newname(leaf, newname):
        leaf.name = newname
        return leaf

    tree = Phylo.parse(args.tree, args.tree_type).next()
    leafs = (leaf for leaf in tree.get_terminals())

    if args.info:
        info = DictReader(args.info, fieldnames = ['seqname','newname'])
        info = {i['seqname']:i['newname'] for i in info}

        # for newick trees :s will be replaced by |s
        if args.tree_type == 'newick':
            info = {s.replace(':', '|'):n for s,n in info.items()}

        leafs = (l for l in leafs if l.name in info)
        leafs = (newname(l, info[l.name]) for l in leafs)

    if args.remove_word:
        leafs = (newname(l, re.sub(args.remove_word, '', l.name)) for l in leafs)
        leafs = (newname(l, l.name.strip()) for l in leafs)

    leafs = (newname(l, args.add_prefix + l.name) for l in leafs)
    leafs = (newname(l, l.name + args.add_suffix) for l in leafs)

    # do this last
    if args.tree_type == 'newick':
        leafs = (newname(l, l.name.replace(' ', '_')) for l in leafs)

    # execute changes and write tree
    list(leafs)
    Phylo.write(tree, args.out, args.tree_type)

