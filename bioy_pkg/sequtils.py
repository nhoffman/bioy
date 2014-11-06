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

import csv
import contextlib
import tempfile
import logging
import numpy
import re
import subprocess
import utils

from cStringIO import StringIO
from itertools import tee, izip_longest, groupby, takewhile, izip
from collections import Counter, defaultdict, namedtuple
from operator import itemgetter
from subprocess import Popen, PIPE

log = logging.getLogger(__name__)

BLAST_HEADERS = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

# use BLAST_FORMAT as input to blastn -outfmt
BLAST_FORMAT = "6 qseqid sseqid pident qstart qend qlen"
BLAST_HEADER = BLAST_FORMAT.split()[1:] + ['coverage']

ERRORS = ['snp', 'indel', 'homoindel', 'compound']

UCLUST_HEADERS = ['type', 'cluster_number', 'size', 'pct_id', 'strand',
                  'query_start', 'seed_start', 'alignment', 'query_label',
                  'target_label']

FETCH_HEADERS = ['seqid', 'gi', 'gb', 'taxid', 'fullname']
RANKS = ['root', 'superkingdom', 'phylum', 'class',
         'order', 'family', 'genus', 'species']

gap = '-'

homogap = '='

IUPAC = {('A',): 'A',
         ('C',): 'C',
         ('G',): 'G',
         ('T',): 'T',
         (gap,): gap,
         (homogap,): homogap,
         ('A', 'C'): 'M',
         ('A', 'G'): 'R',
         ('A', 'T'): 'W',
         ('C', 'G'): 'S',
         ('C', 'T'): 'Y',
         ('G', 'T'): 'K',
         ('A', 'C', 'G'): 'V',
         ('A', 'C', 'T'): 'H',
         ('A', 'G', 'T'): 'D',
         ('C', 'G', 'T'): 'B',
         ('A', 'C', 'G', 'T'): 'N', }

CAPUI = {v: set(k) for k, v in IUPAC.items()}

# provides criteria for defining matching tax_ids as "unclassified"
UNCLASSIFIED_REGEX = re.compile(
    r'' + r'|'.join(frozenset(['-like',
                               'Taxon'
                               '\d\d',
                               'acidophile',
                               'actinobacterium',
                               'aerobic',
                               r'\b[Al]g(um|a)\b',
                               r'\b[Bb]acteri(um|a)',
                               'Barophile',
                               'cyanobacterium',
                               'Chloroplast',
                               '[Cc]lone',
                               '[Cc]loning',
                               'cluster',
                               '-containing',
                               'environmental',
                               'epibiont',
                               # 'et al',
                               'eubacterium',
                               r'\b[Gg]roup\b',
                               'halophilic',
                               r'hydrothermal\b',
                               'isolate',
                               'marine',
                               'methanotroph',
                               'microorganism',
                               'mollicute',
                               'pathogen',
                               '[Pp]hytoplasma',
                               'proteobacterium',
                               'putative',
                               r'\bsp\.',
                               'species',
                               'spirochete',
                               r'str\.'
                               'strain',
                               'symbiont',
                               'taxon',
                               'unicellular',
                               'uncultured',
                               'unclassified',
                               'unidentified',
                               'unknown',
                               'vector\b',
                               r'vent\b',
                               ])))


def homoencode(seq):
    """Run length encode a string

    Returns a 2-tuple (rle_seq, counts), where rle_seq is the run
    length encoded sequence with no counts, and counts is a list of
    integers giving the lengths of the homopolymers.

    For example, homoencode('AATGGGC') ==> ('ATGC', [2,1,3,1])
    """

    assert gap not in seq

    seq = seq.upper()
    seq = groupby(seq)
    seq = ((c, len(list(g))) for c, g in seq)

    chars, counts = izip(*seq)

    return ''.join(chars), list(counts)


def homodecodealignment(seq1, counts1, seq2, counts2, insertion=homogap):
    """Decode a pair of aligned, run length encoded sequences.

    This maps an alignment in run length encoding to an alignment in
    normal encoding. *seq1* and *seq2* should be aligned, run length
    encoded sequences. *counts1* and *counts2* are their corresponding
    counts.
    """

    assert len(seq1.replace(gap, '').replace(insertion, '')) == len(counts1)
    assert len(seq2.replace(gap, '').replace(insertion, '')) == len(counts2)

    if not seq1:
        return seq1, homodecode(seq2, counts2)
    if not seq2:
        return homodecode(seq1, counts1), seq2

    assert not seq1[0] == seq2[0] == gap

    if seq1[0] == gap:
        suff1, suff2 = homodecodealignment(
            seq1[1:], counts1, seq2[1:], counts2[1:])
        deseq1 = gap * counts2[0] + suff1
        deseq2 = seq2[0] * counts2[0] + suff2
    elif seq2[0] == gap:
        suff1, suff2 = homodecodealignment(
            seq1[1:], counts1[1:], seq2[1:], counts2)
        deseq1 = seq1[0] * counts1[0] + suff1
        deseq2 = gap * counts1[0] + suff2
    else:
        m = max(counts1[0], counts2[0])
        suff1, suff2 = homodecodealignment(
            seq1[1:], counts1[1:], seq2[1:], counts2[1:])
        deseq1 = seq1[0] * counts1[0] + insertion * (m - counts1[0]) + suff1
        deseq2 = seq2[0] * counts2[0] + insertion * (m - counts2[0]) + suff2

    return deseq1, deseq2


def homodecode(seq, counts, insertion=homogap):
    """Expand a run length encoded sequence.

    *seq* should be a run length encoded string, minus homopolymer
    lengths, and *counts* should be a list of integers giving the
    length. homodecode is the inverse of homoencode, so::

        homodecode(*homoencode(seq)) == seq

    *counts* should have the same length as the number of bases in
    *seq* that are not - or = (both representations for gaps.
    """

    if not seq or not counts:
        return seq
    elif seq[0] == gap or seq[0] == insertion:
        return gap + homodecode(seq[1:], counts)
    else:
        return seq[0] * counts[0] + homodecode(seq[1:], counts[1:])


def to_ascii(nums):
    """
    Encode a list of integers as an ascii string. Max allowed value is
    126-48 = 78 (avoids many special characters that might cause
    trouble). The purpose of this encoding is to store run-length
    encodings, so we're asuming that values > 78 are not plausible.
    """

    if max(nums) > 78:
        raise ValueError('values over 78 are not allowed')

    return ''.join([chr(i + 48) for i in nums])


def from_ascii(chars):
    """
    Decode an ascii-encoded list of integers.
    """
    return [ord(c) - 48 for c in chars]


def cons_rle(c, ceiling=False):
    """
    Choose a consensus run length given counts of run lengths in
    Counter object `c`. Choose the most common run length count. In
    the case of ties, choose the smaller value.
    """

    # TODO: return counts for only chosen cons base
    most_common = c.most_common()
    _, top_freq = most_common[0]
    most_common = [cnt for cnt, freq in most_common if freq == top_freq]
    most_common = sorted(most_common, reverse=ceiling)

    return most_common[0]


def cons_char(c, gap_ratio=0.5):
    """
    Choose a consensus character given counts of characters in Counter
    object `c`.
    """

    # if gap char freq more than half then return gaps, else remove gap char
    if float(c[gap]) / sum(c.values()) > gap_ratio:
        return gap

    del c[gap]

    most_common = c.most_common()
    _, top_freq = most_common[0]
    most_common = [ch.upper() for ch, freq in most_common if freq == top_freq]
    most_common = sorted(most_common)

    return IUPAC[tuple(most_common)]


def get_char_counts(seqs):
    """
    Return a list of Counter objects for each position of a sequence
    of iterables in `seqs`.
    """

    counter = defaultdict(Counter)
    for seq in seqs:
        for i, c in enumerate(seq.seq):
            for b in CAPUI[c]:
                counter[i][b] += 1

    # counters, sorted by position
    return [count for position, count in sorted(counter.items())]


def get_rle_counts(seqs, rlelist):
    """
    Return a list of two-tuples: (char_count, rle_count) where
    char_count and rle_count are Counter objects tallying characters
    and run length counts, respectively for each position of sequences
    in in `seqs`.

     * seqs - sequence of SeqRecord objects
     * rlelist - sequence of lists of run length
                 counts correponding to seach sequence
    """

    char_counter = defaultdict(Counter)
    rle_counter = defaultdict(Counter)
    for seq, rle_counts in izip_longest(seqs, rlelist):
        assert len(seq.seq.replace(gap, '')) == len(rle_counts)

        rle_counts = iter(rle_counts)

        # for each non-gap character in seq, tally the run length; if
        # c is a gap, save a tally of 1
        for i, c in enumerate(seq.seq):
            for b in CAPUI[c]:
                char_counter[i][b] += 1
            rle_counter[i][rle_counts.next() if c != gap else 1] += 1

    # tuples of counters, sorted by position
    return [(char_counter[i], rle_counter[i])
            for i in xrange(len(char_counter))]


def consensus(seqs, rlelist=None, degap=True):
    """
    Calculate a consensus for an iterable of SeqRecord objects. seqs
    are decoded using corresponding lists of read length counts in
    `rlelist` if provided. Gaps are removed if degap is True.
    """

    if rlelist:
        cons = [cons_char(c) * cons_rle(n)
                for c, n in get_rle_counts(seqs, rlelist)]
    else:
        cons = [cons_char(c) for c in get_char_counts(seqs)]

    return ''.join(cons).replace(gap, '') if degap else ''.join(cons)


@contextlib.contextmanager
def fasta_tempfile(seqs, dir=None):
    """Creates a temporary FASTA file representing an iterable of
    objects *seqs* containing attributes `description` and `seq` (eg,
    SeqRecord or fastalite objects). Returns the name of a temorary
    file and, then deletes it at the end of the with block.
    """

    handle = tempfile.NamedTemporaryFile(
        mode='w', suffix='.fasta', dir=dir)

    handle.write(''.join('>{s.id}\n{s.seq}\n'.format(s=s) for s in seqs))
    handle.flush()

    try:
        yield handle.name
    finally:
        handle.close()


@contextlib.contextmanager
def run_ssearch(query, target, outfile=None, cleanup=True,
                ssearch='ssearch36', max_hits=None,
                full_length=True, m10=True, dna=True,
                forward_only=True, args=None):
    """Align sequences in fasta-format files ``query`` and ``target``
    using ssearch36. Returns a file-like object open for
    reading. This is meant to be run in a with block.
    Other options are follows:

    * max_hits      if provided an integer value, specify option '-d <value>'
    * full_length   if True, specify '-a'
    * m10           if True, specify '-m 10'
    * dna           if True, specify '-n'
    * forward_only  if True, specify '-3'

    `args` is a list of options whicn overrides all of the above if
    provided.

    Example:

        with run_ssearch(query, target) as f:
            aligns = parse_ssearch36(f)

    """

    if outfile:
        filename = outfile
    else:
        handle = tempfile.NamedTemporaryFile(
            mode='rw', suffix='.ssearch')
        filename = handle.name

    cmd = [ssearch]
    if args:
        cmd.extend(args)
    else:
        if max_hits:
            cmd.extend(['-d', str(max_hits)])
        if full_length:
            cmd.append('-a')
        if m10:
            cmd.extend(['-m', '10'])
        if dna:
            cmd.append('-n')
        if forward_only:
            cmd.append('-3')

    cmd.extend([query, target, '>', filename])
    command = ' '.join(cmd)

    log.info(command)

    try:
        subprocess.check_call(command, shell=True)
        if outfile:
            handle = open(filename, 'rU')
        else:
            handle.seek(0)
        yield handle
    finally:
        if not outfile and cleanup:
            handle.close()


def all_pairwise(seqs):
    """
    Perform all pairwise alignments among
    sequences in list of SeqRecords `seqs`
    """

    if not hasattr(seqs, 'len'):
        seqs = list(seqs)

    for i in xrange(len(seqs) - 1):
        with fasta_tempfile([seqs[i]]) as target:
            with fasta_tempfile(seqs[i + 1:]) as query:
                with run_ssearch(query, target, max_hits=1) as aligned:
                    for d in parse_ssearch36(aligned):
                        yield (d['t_name'], d['q_name'], float(d['sw_ident']))


def names_from_pairs(pairs):

    try:
        first, second, _ = pairs.pop(0)
    except AttributeError:
        first, second, _ = pairs.next()

    yield first
    yield second
    for q, t, _ in pairs:
        if q == first:
            yield t


def run_muscle(seqs, tmpdir=None, keep_order=True):
    """
    Write an iterable of SeqRecord objects to a temporary file, align
    with muscle, and return an iterable of aligned SeqRecords. Output
    order is same as input if `keep_order` is True.
    """

    # note that muscle has limited support for ambuguity codes:
    # http://www.drive5.com/muscle/muscle_userguide3.8.html
    # 3.1.2 Nucleotide sequences - The usual letters A, G, C, T and U
    # stand for nucleotides. The letters T and U are equivalent as far
    # as MUSCLE is concerned. N is the wildcard meaning "unknown
    # nucleotide". R means A or G, Y means C or T/U. Other wildcards,
    # such as those used by RFAM, are not understood in this version
    # and will be replaced by Ns. If you would like support for other
    # DNA / RNA alphabets, please let me know.

    if keep_order:
        seqs, seqs1 = tee(seqs)
        sortdict = {s.id: i for i, s in enumerate(seqs1)}

    with fasta_tempfile(seqs, tmpdir) as f:
        command = ['muscle', '-quiet', '-seqtype', 'dna', '-in', f]
        pipe = Popen(command, stdout=PIPE)
        (seqstr, _) = pipe.communicate()

    aligned = fastalite(StringIO(seqstr))

    if keep_order:
        aligned = iter(sorted(aligned, key=lambda s: sortdict[s.id]))

    return aligned


def parse_uc(infile):
    """
    Return two dicts: {read_name: cluster} and {cluster: cluster_size}
    """

    cluster_ids = {}
    cluster_sizes = {}

    rows = csv.DictReader(infile, delimiter='\t', fieldnames=UCLUST_HEADERS)
    for row in rows:
        cluster = int(row['cluster_number'])
        if row['type'] == 'C':
            cluster_sizes[cluster] = int(row['size'])
        else:
            cluster_ids[row['query_label']] = cluster

    return cluster_ids, cluster_sizes


def parse_blast(blast, extras=[]):
    for line in blast:
        if line and not line[0] == '#':
            yield dict(zip(BLAST_HEADERS + extras, line.split()))


def itemize_errors(ref, query):
    """
    Itemize differences between aligned strings `ref` and `query`

     * ref - reference sequence
     * query - sequence on which error calculation is performed

     Return an iterable of dicts with keys

      * i - position relative to ref
      * ref - base or bases in ref
      * query - base or bases in query

    """

    # shift query to preserve alignment
    query = query[_find_homochar_length(ref, gap):]
    ref = ref.strip(gap)  # strip terminal gaps

    t = h = 0
    m = min(len(query), len(ref))
    results = []

    while t < m:
        h += max(_find_homochar_length(ref[t:]),
                 _find_homochar_length(query[t:]))

        assert (t < h)  # strange chars found or miss-alignment

        results.append({
            'i': len(ref[:t].replace(gap, '').replace(homogap, '')),
            'ref': ref[t:h],
            'query': query[t:h]
        })

        t = h

    return results


def _find_homochar_length(s, char='', ignore=[gap, homogap]):
    '''
    Returns the length of repeating chars from the left side of a string
    '''
    i = 0
    if s:
        char = char or ('' if s[0] in ignore else s[0])
        while char == s[i:i + 1]:
            i += 1
    return i


def error_category(e, gap=gap, homogap=homogap, errors=ERRORS):
    r, q = e['ref'], e['query']
    if r == q:
        return 'equal'
    elif gap in r or gap in q:
        return errors[1]
    elif len(r) == 1:
        return errors[0]
    elif set(c for c in r if c != homogap) == \
            set(c for c in q if c != homogap):
        return errors[2]
    else:
        return errors[3]


def error_count(errors):
    cnt = Counter()
    for e in errors:
        cnt[error_category(e)] += 1
        cnt['length'] += len(e['ref'].strip('-='))
    return cnt


def parse_ssearch36(lines, numeric=False):
    """
    Parse output of 'ssearch36 -m 10 query.fasta library.fasta'

    Return an iterator of dicts containing alignment data read from
    file-like object `lines`. Coerce strings to ints or floats if
    numeric is True.

    Each alignment of a single query sequence to multiple targets (ie,
    different target sequences or multiple regions within the same
    target) is represented by an element in the output; use
    'groupby(results, key = lambda hit: hit['q_name'])' to group by
    query.

    Note: use 'ssearch36 -a' to retain full sequence.
    """

    query_count = 0
    hit_count = 0
    keeplines = False
    prefix = None
    hit = defaultdict()

    for line in lines:
        line = line.rstrip('\n')

        if line.startswith('>>><<<'):
            # query end
            keeplines = False
        elif line.startswith('>>>'):
            # start of a new hit
            if not line.startswith('>>>///'):
                query_count += 1
            q_name = line.lstrip('>').split(',')[0]
        elif line.startswith('>>') or line.startswith('>--'):
            # hit-specific results; keep results starting here
            if prefix:
                yield hit
            hit_count += 1
            if line.startswith('>>'):
                t_description = line[2:]
                t_name = t_description.split()[0]
            prefix = ''
            keeplines = True
            hit = {'q_name': q_name,
                   't_name': t_name,
                   't_description': t_description,
                   'q_seq': '',
                   't_seq': ''}
        elif line.startswith('>'):
            prefix = 't_' if prefix else 'q_'
        elif line.startswith(';') and keeplines:
            k, v = line.lstrip('; ').split(':', 1)
            k = k.replace(gap, '').replace(' ', '_').lower()
            if k == 'al_cons':
                hit[k] = ''
            else:
                hit[prefix + k] = utils.cast(v) if numeric else v.strip()
        elif prefix and keeplines:
            if 'al_cons' in hit:
                hit['al_cons'] += line
            else:
                hit[prefix + 'seq'] += line.strip()

    if hit:
        yield hit

    log.info('%s queries, %s hits' % (query_count, hit_count))


def wrap(text, width=60):
    """
    Wraps input string [text] to [width] characters. Return a list of
    substrings.
    """

    r1 = range(0, len(text) + 1, width)
    r2 = range(width, len(text) + width + 1, width)
    return [text[f:t] for f, t in zip(r1, r2)]


def format_alignment(seq1, seq2, name1='', name2='',
                     seqwidth=80, namewidth=15):

    izl = lambda s: izip_longest(*s, fillvalue='')

    seqs = [seq1,
            ''.join((':' if c1 == c2 else ' ')
                    for c1, c2 in izl([seq1, seq2])),
            seq2]

    names = [name1, '', name2]

    fstr = '%%%ss %%s' % namewidth
    ss = []
    for strands in izl([wrap(s, seqwidth) for s in seqs]):
        for name, strand in zip(names, strands):
            ss.append(fstr % (name[:namewidth], strand))
        ss.append('')

    return '\n'.join(ss)


def grouper(n, iterable, pad=False):
    """
    Return sequence of n-tuples composed of successive elements of
    iterable; If pad is True, the last tuple is padded with None to
    make it the same length as the others. Not safe for iterables with
    None elements.
    """

    args = [iter(iterable)] * n
    iterout = izip_longest(fillvalue=None, *args)

    if pad:
        return iterout
    else:
        return (takewhile(lambda x: x is not None, c) for c in iterout)


def show_errors(errorlist, width=40):
    errors = sorted(errorlist, key=lambda e: e['i'])
    out = ''
    for group in grouper(width, errors):
        ii, rr, qq, cc = [], [], [], []
        for d in group:
            i, r, q = str(d['i']), d['ref'], d['query']
            fstr = '%-' + str(max([len(r), len(q), len(i)])) + 's'
            ii.append(fstr % i)
            rr.append(fstr % r)
            qq.append(fstr % q)
            c = error_category(d)[0]
            cc.append(fstr % ('' if c == 'e' else c))

        out += '\n  i:   ' + '|'.join(ii)
        out += '\nref:   ' + '|'.join(rr)
        out += '\nquery: ' + '|'.join(qq)
        out += '\n     : ' + ' '.join(cc)
        out += '\n'

    return out


def format_taxonomy(names, selectors, asterisk='*'):
    """
    Create a friendly formatted string of taxonomy names. Names will
    have an asterisk value appended *only* if the cooresponding
    element in the selectors evaluates to True.
    """

    names = izip_longest(names, selectors)
    names = ((n, asterisk if s else '')
             for n, s in names)  # add asterisk to selected names
    names = set(names)
    names = sorted(names)  # sort by the name plus asterisk
    names = groupby(names, key=itemgetter(0))  # group by just the names
    # prefer asterisk names which will be at the bottom
    names = (list(g)[-1] for _, g in names)
    names = (n + a for n, a in names)  # combine names with asterisks
    # assume species names have exactly two words
    species = lambda n: len(n.split()) == 2
    names = sorted(names, key=species)
    names = groupby(names, key=species)

    tax = []

    for species, assigns in names:
        if species:
            # take the species word and combine them with a '/'
            assigns = (a.split() for a in assigns)
            # group by genus name
            assigns = groupby(assigns, key=itemgetter(0))
            assigns = ((k, map(itemgetter(1), g))
                       for k, g in assigns)  # get a list of the species names
            assigns = ('{} {}'.format(k, '/'.join(g))
                       for k, g in assigns)  # combine species names with '/'

        tax.extend(assigns)

    return ';'.join(sorted(tax))


def compound_assignment(assignments, taxonomy):
    """
    Create taxonomic names based on 'assignmnets', which are a set of
    two-tuples: {(tax_id, is_starred), ...} where each tax_id is a key
    into taxdict, and is_starred is a boolean indicating whether at
    least one reference sequence had a parirwise alignment identity
    score meeting some thresholed. 'taxdict' is a dictionary keyed by
    tax_id and returning a dict of taxonomic data corresponding to a
    row from the taxonomy file. If 'include_stars' is False, ignore
    the second element of each tuple in 'assignments' and do not
    include asterisks in the output names.

    assignments = [(tax_id, is_starred),...]
    taxonomy = {taxid:taxonomy}

    Functionality: see format_taxonomy
    """

    if not taxonomy:
        raise TypeError('taxonomy must not be empty or NoneType')

    assignments = ((taxonomy[i]['tax_name'], a) for i, a in assignments)
    assignments = zip(*assignments)

    return format_taxonomy(*assignments, asterisk='*')


def condense_ids(assignments,
                 taxonomy,
                 ranks=RANKS,
                 floor_rank=None,
                 ceiling_rank=None,
                 max_size=3,
                 rank_thresholds={}):
    """
    assignments = [tax_ids...]
    taxonomy = {taxid:taxonomy}

    Functionality: Group items into taxonomic groups given max rank sizes.
    """

    floor_rank = floor_rank or ranks[-1]
    ceiling_rank = ceiling_rank or ranks[0]

    if not taxonomy:
        raise TypeError('taxonomy must not be empty or NoneType')

    if floor_rank not in ranks:
        msg = '{} not in ranks: {}'.format(floor_rank, ranks)
        raise TypeError(msg)

    if ceiling_rank not in ranks:
        msg = '{} not in ranks: {}'.format(ceiling_rank, ranks)
        raise TypeError(msg)

    if ranks.index(floor_rank) < ranks.index(ceiling_rank):
        msg = '{} cannot be lower rank than {}'.format(
            ceiling_rank, floor_rank)
        raise TypeError(msg)

    # set rank to ceiling
    try:
        assignments = {a: taxonomy[a][ceiling_rank] for a in assignments}
    except KeyError:
        print assignments
        error = ('Assignment id not found in taxonomy.')
        raise KeyError(error)

    def condense(groups, ceiling_rank=ceiling_rank, max_size=max_size):
        new_groups = {}

        for a, r in groups.items():
            new_groups[a] = taxonomy[a][ceiling_rank] or r

        num_groups = len(set(new_groups.values()))

        if rank_thresholds.get(ceiling_rank, max_size) < num_groups:
            return groups

        groups = new_groups

        # return if we hit the floor
        if ceiling_rank == floor_rank:
            return groups

        # else move down a rank
        ceiling_rank = ranks[ranks.index(ceiling_rank) + 1]

        # recurse each branch down the tax tree
        for _, g in utils.groupbyl(groups.items(), itemgetter(1)):
            g = condense(dict(g),
                         ceiling_rank,
                         max_size - num_groups + 1)
            groups.update(g)

        return groups

    return condense(assignments)


def correct_copy_numbers(assignments, copy_numbers):
    """Return a float representing a copy number adjustment factor given
    'assignments' (a set of two-tuples: {(tax_id, is_starred), ...})
    and .

    """

    if not assignments:
        raise TypeError('assignments must not be empty or NoneType')

    if copy_numbers:
        assignments = map(itemgetter(0), assignments)
        corrections = [float(copy_numbers[i]) for i in assignments]
        return numpy.mean(corrections)
    else:
        return 1

SeqLite = namedtuple('SeqLite', 'id, description, seq')


def fastalite(handle, limit=None):
    """
    Return a sequence of namedtupe objects given fasta format open
    file-like object `handle`. Sequence is a list if `readfile` is
    True, an iterator otherwise.
    """
    limit = limit or -1

    name, seq = '', ''
    for line in handle:
        if line.startswith('>'):
            if limit != 0:
                limit -= 1
            else:
                break

            if name:
                yield SeqLite(name.split()[0], name, seq)

            name, seq = line[1:].strip(), ''
        else:
            seq += line.strip()

    if name and seq:
        yield SeqLite(name.split()[0], name, seq)


# Taken from Connor McCoy's Deenurp


def tax_of_genbank(gb):
    """
    Get the tax id from a genbank record, returning None if no taxonomy is
    available.
    """
    # Check for bad name
    try:
        source = next(i for i in gb.features if i.type == 'source')
        taxon = source.qualifiers.get('db_xref', [])
        taxon = next(i[6:] for i in taxon if i.startswith('taxon:'))
        taxon = re.findall('\d+', taxon)[0]
        return taxon
    except StopIteration:
        return None


def count_ambiguous(seq):
    s = frozenset('ACGT')
    return sum(i not in s for i in seq)


def is_type(gb):
    """
    Returns a boolean indicating whether a
    sequence is a member of a type strain,
    as indicated by the presence of the string
    '(T)' within the record description.
    """
    return '(T)' in gb.description
###
