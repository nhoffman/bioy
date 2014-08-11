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

"""Create a full taxonomic table of thresholds based on a csv of defaults.

Child tax_ids not specified in the defaults receive the parent threshold
"""

import sys
import logging

from os import path

import pandas as pd

from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)


def read_csv(filename, compression=None, **kwargs):
    """Read a csv file using pandas.read_csv with compression defined by
    the file suffix unless provided.
    """

    suffixes = {'.bz2': 'bz2', '.gz': 'gzip'}
    compression = compression or suffixes.get(path.splitext(filename)[-1])
    kwargs['compression'] = compression

    return pd.read_csv(filename, **kwargs)


def build_parser(parser):
    # required inputs
    parser.add_argument(
        'taxonomy',
        help='must have tax_id column and rank columns')
    parser.add_argument(
        'thresholds',
        help='with required columns tax_id, low, target_rank')

    # common outputs
    parser.add_argument(
        '-o', '--out',
        default=sys.stdout, type=Opener('w'), metavar='FILE',
        help="Classification results.")


def action(args):
    # for debugging:
    # pd.set_option('display.max_columns', None)
    # pd.set_option('display.max_rows', None)

    # load default_tresholds
    taxonomy = read_csv(
        args.taxonomy,
        comment='#',
        dtype=str,
        na_filter=True,  # False is faster
        ).set_index('tax_id')

    # load default_tresholds
    default_thresholds = read_csv(
        args.thresholds,
        comment='#',
        dtype=dict(tax_id=str, low=float, target_rank=str),
        na_filter=True,  # False is faster
        usecols=['tax_id', 'low', 'target_rank']
        )

    tax_cols = taxonomy.columns.tolist()
    tax_cols = tax_cols[tax_cols.index('root'):]

    # out output data structure
    full_tree = pd.DataFrame(index=taxonomy.index)

    # start with the root column and move right
    for index, rank in enumerate(tax_cols):
        defaults = default_thresholds[
            default_thresholds['target_rank'] == rank]
        defaults = defaults.set_index('tax_id')

        # iterate taxonomy most to least specificity and join with defaults
        target_thresholds = []
        for specificity in reversed(tax_cols):
            thresholds = taxonomy[[specificity]].join(
                defaults[['low']], on=specificity, how='inner')
            # append just the low column
            target_thresholds.append(thresholds['low'])

        # concat and take just the first(), most specific threshold
        target_thresholds = pd.concat(target_thresholds)
        target_thresholds = target_thresholds.groupby(
            target_thresholds.index, sort=False).first()

        full_tree[rank] = target_thresholds

        # fill in tax holes
        index = max(index-1, 0)
        blanks = full_tree[full_tree[rank].isnull()]
        full_tree[rank] = full_tree[rank].fillna(
            blanks[tax_cols[index]])

    full_tree.to_csv(args.out)
