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

"""Classify sequences by grouping blast output by matching taxonomic names

Optional grouping by specimen and query sequences

Running the program
-------------------

::

    positional arguments:
      blast_file            CSV tabular blast file of query and subject hits.
      seq_info              File mapping reference seq name to tax_id
      taxonomy              Table defining the taxonomy for each tax_id

    optional arguments:
      -h, --help            show this help message and exit
      --threads NUM         Number of threads (CPUs). Can also specify with
                            environment variable THREADS_ALLOC. [32]
      --copy-numbers CSV    Estimated 16s rRNA gene copy number for each
                            tax_ids (CSV file with columns: tax_id, median)
      --rank-thresholds CSV
                            Columns [tax_id,ranks...]
      --specimen-map CSV    CSV file with columns (name, specimen) assigning
                            sequences to groups. The default behavior is to
                            treat all query sequences as
                            belonging to one specimen.
      -w CSV, --weights CSV
                            Optional headless csv file with columns 'seqname',
                            'count' providing weights for each query sequence
                            described in the blast input (used, for example, to
                            describe cluster sizes for corresponding cluster
                            centroids).
      -o FILE, --out FILE   Classification results.
      -O FILE, --details-out FILE
                            Optional details of taxonomic assignments.
      --details-full        do not limit out_details to only larget cluster per
                            assignment
      --group-def INTEGER   define a group threshold for a particular rank
                            overriding --target-max-group-size. example:
                            genus:2 (NOT IMPLEMENTED)
      --has-header          specify this if blast data has a header
      --min-identity PERCENT
                            minimum identity threshold for accepting matches
                            [>= 0.0]
      --max-identity PERCENT
                            maximum identity threshold for accepting matches
                            [<= 100.0]
      --min-cluster-size INTEGER
                            minimum cluster size to include in classification
                            output [1]
      --min-coverage PERCENT
                            percent of alignment coverage of blast result [0.0]
      --specimen LABEL      Single group label for reads
      --starred PERCENT     Names of organisms for which at least one reference
                            sequence has pairwise identity with a query
                            sequence of at least PERCENT will be marked with an
                            asterisk[100.0]
      --target-max-group-size INTEGER
                            group multiple target-rank assignments that excede
                            a threshold to a higher rank [3]
      --target-rank TARGET_RANK
                            Rank at which to classify. Default: "species"

Positional arguments
++++++++++++++++++++

blast_file
==========

A csv file with columns **qseqid**, **sseqid**, **pident**,
**qstart**, **qend** and **qlen**.

.. note:: The actual header is optional and if
          present make sure to use the --has-header switch

seq_info
========

A csv file with minimum columns **seqname** and **tax_id**.  Additional
columns will be included in the details output.

taxonomy
========

A csv file with columns **tax_id**, **rank** and **tax_name**, plus at least
one additional rank column(s) creating a taxonomic tree such as **species**,
**genus**, **family**, **class**, **pylum**, **kingdom** and/or **root**.
The rank columns also give an order of specificity from right to left,
least specific to most specific respectively.

Optional input
++++++++++++++

rank-thresholds
===============

TODO

copy-numbers
============

Below is an *example* copy numbers csv with the required columns:

    ====== ==================== ======
    tax_id tax_name             median
    ====== ==================== ======
    155977 Acaryochloris        2.00
    155978 Acaryochloris marina 2.00
    434    Acetobacter          5.00
    433    Acetobacteraceae     3.60
    ====== ==================== ======

weights
=======

Headerless file containing two columns specifying the seqname (clustername) and
weight (or number of sequences in the cluster).

Output
++++++

out
===

A csv with columns and headers as in the example below:

    =========== =============== ======================================
     specimen    assignment_id   assignment
    =========== =============== ======================================
      039_3      0               Pseudomonas mendocina;Pseudonocardia
      039_3      1               Rhizobiales
      039_3      2               Alcaligenes faecalis*
      039_3      3               [no blast result]
    =========== =============== ======================================

    ======= ============= =============
     low     max_percent   min_percent
    ======= ============= =============
     95.00   99.02         95.74
     95.00   98.91         95.31
     99.00   100.00        99.00

    ======= ============= =============

    ============= ======= =========== ===========
     target_rank   reads   pct_reads   clusters
    ============= ======= =========== ===========
     species       6       35.29       1
     genus         5       29.41       1
     species       5       29.41       1
                   1       5.88        1
    ============= ======= =========== ===========

details-out
===========

A csv that is basically a blast results breakdown of the `out`_ output.

Internal functions
------------------

"""

import sys
import logging

from os import path

import pandas as pd
import math

from bioy_pkg import sequtils, _data as datadir
from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)


def freq(df):
    """Calculate pct_reads column.
    """
    df['pct_reads'] = df['reads'] / df['reads'].sum() * 100
    return df


def read_csv(filename, compression=None, **kwargs):
    """Read a csv file using pandas.read_csv with compression defined by
    the file suffix unless provided.
    """

    suffixes = {'.bz2': 'bz2', '.gz': 'gzip'}
    compression = compression or suffixes.get(path.splitext(filename)[-1])
    kwargs['compression'] = compression

    return pd.read_csv(filename, **kwargs)


def star(df, starred):
    """Assign boolean if any items in the
    dataframe are above the star threshold.
    """

    df['starred'] = df.pident.apply(lambda x: x >= starred).any()
    return df


def condense_ids(df, tax_dict, ranks, max_group_size):
    """Create mapping from tax_id to its
    condensed id and set assignment hash.
    """

    condensed = sequtils.condense_ids(
        df.tax_id.unique(),
        tax_dict,
        ranks=ranks,
        max_size=max_group_size)
    condensed = pd.DataFrame(
        condensed.items(),
        columns=['tax_id', 'condensed_id']).set_index('tax_id')
    assignment_hash = hash(frozenset(condensed.condensed_id.unique()))
    condensed['assignment_hash'] = assignment_hash
    return df.join(condensed, on='tax_id')


def assign(df, tax_dict):
    """Create str assignment based on tax_ids str and starred boolean.
    """
    ids_stars = df.groupby(by=['condensed_id', 'starred']).groups.keys()
    df['assignment'] = sequtils.compound_assignment(ids_stars, tax_dict)
    return df


def agg_columns(df):
    """Create aggregate columns for assignments.
    """
    agg = {}
    agg['max_percent'] = df['pident'].max()
    agg['min_percent'] = df['pident'].min()
    return pd.Series(agg)


def pct_corrected(df):
    """Calculate pct_corrected values.
    """
    df['pct_corrected'] = df['corrected'] / df['corrected'].sum() * 100
    return df


def assignment_id(df):
    """Resets and drops the current dataframe's
    index and sets it to the assignment_hash
    """
    df = df.reset_index(drop=True)  # specimen is retained in the group key
    df.index.name = 'assignment_id'
    return df


def load_rank_thresholds(
        path=path.join(datadir, 'rank_thresholds.csv'), usecols=None):
    """Load a rank-thresholds file.  If no argument is specified the default
    rank_threshold_defaults.csv file will be loaded.
    """
    return read_csv(
        path,
        comment='#',
        usecols=['tax_id'] + usecols,
        dtype=dict(tax_id=str)).set_index('tax_id')


def find_tax_id(series, ranks, rank):
    """Return the most taxonomic specific tax_id available for the given
    Series.
    """

    index = ranks.index(rank)
    series = series[ranks[index:]]
    series = series[~series.isnull()]
    return series.iloc[0]


def assign_tax_ids(df, ranks):
    for r in ranks:
        valid = df[df['{}_threshold'.format(r)] < df['pident']]
        if not valid.empty:
            tax_ids = df[r]
            na_ids = tax_ids.isnull()
            if na_ids.any():
                # bump up missing taxids
                found_ids = df[na_ids].apply(
                    find_tax_id, args=(ranks, r), axis=1)
                have_ids = df[tax_ids.notnull()][r]
                tax_ids = have_ids.append(found_ids)
            valid['tax_id'] = tax_ids
            return valid


def build_parser(parser):
    # required inputs
    parser.add_argument(
        'blast_file',
        help='CSV tabular blast file of query and subject hits.')
    parser.add_argument(
        'seq_info',
        help='File mapping reference seq name to tax_id')
    parser.add_argument(
        'taxonomy',
        help="""Table defining the taxonomy for each tax_id""")

    # optional inputs
    parser.add_argument(
        '--copy-numbers', metavar='CSV',
        help="""Estimated 16s rRNA gene copy number for each tax_ids
        (CSV file with columns: tax_id, median)""")
    parser.add_argument(
        '--rank-thresholds', metavar='CSV',
        help="""Columns [tax_id,ranks...]""")
    parser.add_argument(
        '--specimen-map', metavar='CSV',
        help="""CSV file with columns (name, specimen) assigning sequences to
        groups. The default behavior is to treat all query sequences
        as belonging to one specimen.""")
    parser.add_argument(
        '-w', '--weights', metavar='CSV',
        help="""Optional headless csv file with columns 'seqname',
        'count' providing weights for each query sequence described in
        the blast input (used, for example, to describe cluster sizes
        for corresponding cluster centroids).""")

    # common outputs
    parser.add_argument(
        '-o', '--out', default=sys.stdout, type=Opener('w'),
        metavar='FILE',
        help="Classification results.")
    parser.add_argument(
        '-O', '--details-out', type=Opener('w'), metavar='FILE',
        help="""Optional details of taxonomic assignments.""")

    # switches and options
    parser.add_argument(
        '--details-full', action='store_true',
        help='do not limit out_details to only larget cluster per assignment')
    parser.add_argument(
        '--group-def', metavar='INTEGER', action='append',
        default=[], help="""define a group threshold for a
        particular rank overriding --target-max-group-size. example:
        genus:2 (NOT IMPLEMENTED)""")
    parser.add_argument(
        '--has-header', action='store_true',
        help='specify this if blast data has a header')
    parser.add_argument(
        '--min-identity', default=0.0, metavar='PERCENT', type=float,
        help="""minimum identity threshold
        for accepting matches [>= %(default)s]""")
    parser.add_argument(
        '--max-identity', default=100.0, metavar='PERCENT', type=float,
        help="""maximum identity threshold for
        accepting matches [<= %(default)s]""")
    parser.add_argument(
        '--min-cluster-size', default=1, metavar='INTEGER', type=int,
        help="""minimum cluster size to include in
        classification output [%(default)s]""")
    parser.add_argument(
        '--min-coverage', default=0.0, type=float, metavar='PERCENT',
        help='percent of alignment coverage of blast result [%(default)s]')
    parser.add_argument(
        '--specimen', metavar='LABEL',
        help="""Single group label for reads""")
    parser.add_argument(
        '--starred', default=100.0, metavar='PERCENT', type=float,
        help="""Names of organisms for which at least one reference
        sequence has pairwise identity with a query sequence of at
        least PERCENT will be marked with an asterisk [%(default)s]""")
    parser.add_argument(
        '--target-max-group-size', metavar='INTEGER', default=3, type=int,
        help="""group multiple target-rank assignments that excede a
        threshold to a higher rank [%(default)s]""")
    parser.add_argument(
        '--target-rank', default='species',
        help='Rank at which to classify. Default: "%(default)s"')


def action(args):
    # for debugging:
    pd.set_option('display.max_columns', None)
    # pd.set_option('display.max_rows', None)

    # format blast data and add additional available information
    names = None if args.has_header else sequtils.BLAST_HEADER
    header = 0 if args.has_header else None
    usecols = ['qseqid', 'sseqid', 'pident', 'coverage']
    blast_results = read_csv(
        args.blast_file,
        dtype=dict(qseqid=str, sseqid=str, pident=float, coverage=float),
        names=names,
        na_filter=True,  # False is faster
        header=header,
        usecols=usecols)

    # get a set of qseqids for identifying [no blast hits] after filtering
    qseqids = blast_results[['qseqid']].drop_duplicates().set_index('qseqid')

    # run raw hi, low and coverage filters
    blast_results = blast_results[
        blast_results['coverage'] >= args.min_coverage]

    blast_results = blast_results[
        blast_results['pident'] <= args.max_identity]

    blast_results = blast_results[
        blast_results['pident'] >= args.min_identity]

    # load seq_info as a bridge to the sequence taxonomy.  Additional
    # columns can be specified to be included in the details-out file
    # such as accession number
    seq_info = read_csv(
        args.seq_info,
        usecols=['seqname', 'tax_id', 'accession'],
        dtype=dict(seqname=str, tax_id=str, accession=str),
        index_col='seqname')
    # rename index to match blast results column name
    seq_info.index.name = 'sseqid'

    # merge blast results with seq_info - do this early so that
    # refseqs not represented in the blast results are discarded in
    # the merge.
    blast_results = blast_results.join(seq_info, on='sseqid', how='inner')

    # load the full taxonomy table.  Rank specificity as ordered from
    # left (less specific) to right (more specific)
    taxonomy = read_csv(args.taxonomy, dtype=str)
    # set index after assigning dtype
    taxonomy = taxonomy.set_index('tax_id')

    # get the a list of rank columns ordered by specificity (see above)
    # NOTE: we are assuming the rank columns
    #       are last N columns staring with 'root'
    rank_cols = taxonomy.columns.tolist()
    rank_cols = rank_cols[rank_cols.index('root'):]

    # now combine just the rank columns to the blast results
    blast_results = blast_results.join(taxonomy[rank_cols], on='tax_id')

    # load the default rank thresholds
    rank_thresholds = load_rank_thresholds(usecols=rank_cols)

    # and any additional thresholds specified by the user
    if args.rank_thresholds:
        rank_thresholds = rank_thresholds.append(
            load_rank_thresholds(path=args.rank_thresholds, usecols=rank_cols))

    blast_results = blast_results.join(
        rank_thresholds, on='tax_id', rsuffix='_threshold')

    #  tax_id will be later based on pident and thresholds
    blast_results = blast_results.drop('tax_id', axis=1)

    # Keep most specific hits.  This step is required so taxonomic
    # assignments are not masked in the results by hits to their
    # more generic taxonomic parents
    blast_results = blast_results.groupby(
        by=['qseqid'], sort=False, group_keys=False).apply(
            assign_tax_ids, list(reversed(rank_cols)))

    # join with taxonomy for tax_name and rank
    blast_results = blast_results.join(
        taxonomy[['tax_name', 'rank']], on='tax_id')

    # merge qseqids that have no hits back into blast_results
    blast_results = blast_results.join(qseqids, on='qseqid', how='outer')

    # assign specimen groups
    specimens = blast_results[['qseqid']].drop_duplicates().set_index('qseqid')

    # load specimen-map and assign specimen names
    if args.specimen_map:
        # if a specimen_map is defined and a qseqid is not included in the map
        # hits to that qseqid will be dropped
        spec_map = read_csv(
            args.specimen_map,
            names=['qseqid', 'specimen'],
            index_col='qseqid')
        specimens = specimens.join(spec_map)
    elif args.specimen:
        specimens['specimen'] = args.specimen
    else:
        specimens['specimen'] = specimens.index  # qseqid

    # join specimen labels onto blast_results
    blast_results = blast_results.join(specimens, on='qseqid')

    # assign seqs that had no results to [no blast_result]
    no_hits = blast_results[blast_results.sseqid.isnull()]
    no_hits['assignment'] = '[no blast result]'
    no_hits['assignment_hash'] = 0

    # move on to seqs that have blast hits
    blast_results = blast_results[blast_results.sseqid.notnull()]

    # TODO: this is relatively slow, need to integrate
    # pandas into sequtils.condense_ids
    tax_dict = {i: t.to_dict() for i, t in taxonomy.fillna('').iterrows()}

    # create condensed assignment hashes by qseqid
    blast_results = blast_results.groupby(
        by=['qseqid'], sort=False).apply(
            condense_ids, tax_dict, rank_cols, args.target_max_group_size)

    # star condensed ids if one hit meets star threshold
    blast_results = blast_results.groupby(
        by=['assignment_hash', 'condensed_id'],
        sort=False).apply(star, args.starred)

    blast_results = blast_results.groupby(
        by=['assignment_hash'], sort=False).apply(assign, tax_dict)

    # put assignments and no assignments back together
    blast_results = pd.concat([blast_results, no_hits])

    # concludes our blast details, on to output summary

    output = blast_results.groupby(
        by=['specimen', 'assignment_hash', 'assignment'])
    output = output.apply(agg_columns)
    output = output.reset_index(level='assignment')

    # qseqid cluster stats
    clusters = blast_results[['qseqid', 'specimen', 'assignment_hash']]
    clusters = clusters.drop_duplicates().set_index('qseqid')

    if args.weights:
        weights = read_csv(
            args.weights,
            names=['qseqid', 'weight'],
            index_col='qseqid')
        clusters = clusters.join(weights)
        # enforce weight dtype as int and unlisted qseq's to weight of 1
        clusters['weight'] = clusters['weight'].fillna(1).astype(int)
    else:
        clusters['weight'] = 1

    clusters = clusters.groupby(by=['specimen', 'assignment_hash'], sort=False)

    output['reads'] = clusters['weight'].sum()

    output = output.groupby(level='specimen').apply(freq)

    output['clusters'] = clusters.size()

    # copy number corrections
    if args.copy_numbers:
        copy_numbers = read_csv(
            args.copy_numbers,
            dtype=dict(tax_id=str),
            usecols=['tax_id', 'median']).set_index('tax_id')

        # get root out (taxid: 1) and set it as the
        # default correction value with index nan
        default = copy_numbers.get_value('1', 'median')
        default = pd.DataFrame(default, index=[None], columns=['median'])
        copy_numbers = copy_numbers.append(default)

        # do our copy number correction math
        corrections = blast_results[['tax_id', 'specimen', 'assignment_hash']]
        corrections = corrections.drop_duplicates().set_index('tax_id')
        corrections = corrections.join(copy_numbers)
        corrections = corrections.groupby(
            by=['specimen', 'assignment_hash'], sort=False)
        corrections = corrections['median'].mean()
        output['corrected'] = output['reads'] / corrections

        # reset corrected counts to int before calculating pct_corrected
        output['corrected'] = output['corrected'].apply(math.ceil)
        output['corrected'] = output['corrected'].fillna(1).astype(int)

        output = output.groupby(level='specimen').apply(pct_corrected)

    # round up any pct < 0.01 for before sorting
    round_up = lambda x: max(0.01, x)
    output['pct_reads'] = output['pct_reads'].map(round_up)
    if args.copy_numbers:
        output['pct_corrected'] = output['pct_corrected'].map(round_up)

    # sort by:
    # 1) specimen
    # 2) read/corrected count
    # 3) cluster count
    # 4) alpha assignment
    columns = ['corrected'] if args.copy_numbers else ['reads']
    columns += ['clusters', 'assignment']
    output = output.sort(columns=columns, ascending=False)
    # nobody cares about assignment_hash
    output = output.reset_index(level='assignment_hash', drop=True)
    output = output.sort_index()

    # create assignment ids by specimen
    output = output.groupby(level="specimen", sort=False).apply(assignment_id)

    # output results
    with args.out as out:
        output.to_csv(out, index=True, float_format='%.2f')

    # output to details.csv.bz2
    if args.details_out:
        blast_results = blast_results.merge(output.reset_index(), how='left')

        if not args.details_full:
            largest = clusters.apply(lambda x: x['weight'].nlargest(1))
            blast_results = blast_results.merge(largest.reset_index())

        # nobody cares about assignment_hash
        blast_results = blast_results.drop(labels='assignment_hash', axis=1)

        with args.details_out as out_details:
            blast_results.to_csv(
                out_details, header=True, index=False, float_format='%.2f')
