import csv
import contextlib
import tempfile
import os
import logging
import re

from cStringIO import StringIO
from itertools import tee, izip_longest, groupby, takewhile
from collections import Counter, defaultdict, namedtuple
from subprocess import Popen, PIPE
from utils import cast

log = logging.getLogger(__name__)

BLAST_HEADERS = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

# use BLAST_FORMAT as input to blastn -outfmt
BLAST_FORMAT = "6 qseqid sseqid pident qstart qend qlen"
BLAST_HEADER = BLAST_FORMAT.split()[1:]

ERRORS = ['snp', 'indel', 'homoindel', 'compound']

UCLUST_HEADERS = ['type', 'cluster_number', 'size', 'pct_id', 'strand',
        'query_start', 'seed_start', 'alignment', 'query_label', 'target_label']

IUPAC = {('A',): 'A',
        ('A', 'C'): 'M',
        ('A', 'C', 'G'): 'V',
        ('A', 'C', 'T'): 'H',
        ('A', 'G'): 'R',
        ('A', 'G', 'T'): 'D',
        ('A', 'T'): 'W',
        ('C',): 'C',
        ('C', 'G'): 'S',
        ('C', 'G', 'T'): 'B',
        ('C', 'T'): 'Y',
        ('G',): 'G',
        ('A', 'C', 'G', 'T'): 'N',
        ('G', 'T'): 'K',
        ('T',): 'T'}

gap = '-'

homogap = '='

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

    chars, counts = zip(*[(c, len(list(g))) for c, g in groupby(seq.upper())])
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
        suff1, suff2 = homodecodealignment(seq1[1:], counts1, seq2[1:], counts2[1:])
        deseq1 = gap * counts2[0] + suff1
        deseq2 = seq2[0] * counts2[0] + suff2
    elif seq2[0] == gap:
        suff1, suff2 = homodecodealignment(seq1[1:], counts1[1:], seq2[1:], counts2)
        deseq1 = seq1[0] * counts1[0] + suff1
        deseq2 = gap * counts1[0] + suff2
    else:
        m = max(counts1[0], counts2[0])
        suff1, suff2 = homodecodealignment(seq1[1:], counts1[1:], seq2[1:], counts2[1:])
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

def cons_rle(c, ceiling = False):
    """
    Choose a consensus run length given counts of run lengths in
    Counter object `c`. Choose the most common run length count. In
    the case of ties, choose the smaller value.
    """

    # TODO: return counts for only chosen cons base
    most_common = c.most_common()
    _,top_freq = most_common[0]
    most_common = [cnt for cnt,freq in most_common if freq == top_freq]
    most_common = sorted(most_common, reverse = ceiling)

    return most_common[0]

def cons_char(c, gap_ratio = 0.5):
    """
    Choose a consensus character given counts of characters in Counter
    object `c`.
    """

    # if gap char freq more than half then return gaps, else remove gap char
    if float(c[gap]) / sum(c.values()) > gap_ratio:
        return gap

    del c[gap]

    most_common = c.most_common()
    _,top_freq = most_common[0]
    most_common = [ch.upper() for ch,freq in most_common if freq == top_freq]
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
            counter[i][c] += 1

    # counters, sorted by position
    return [count for position, count in sorted(counter.items())]

def get_rle_counts(seqs, rlelist):
    """
    Return a list of two-tuples: (char_count, rle_count) where
    char_count and rle_count are Counter objects tallying characters
    and run length counts, respectively for each position of sequences
    in in `seqs`.

     * seqs - sequence of SeqRecord objects
     * rlelist - sequence of lists of run length counts correponding to seach sequence
    """

    char_counter = defaultdict(Counter)
    rle_counter = defaultdict(Counter)
    for seq, rle_counts in izip_longest(seqs, rlelist):
        assert len(seq.seq.replace(gap, '')) == len(rle_counts)

        rle_counts = iter(rle_counts)

        # for each non-gap character in seq, tally the run length; if
        # c is a gap, save a tally of 1
        for i, c in enumerate(seq.seq):
            char_counter[i][c] += 1
            rle_counter[i][rle_counts.next() if c != gap else 1] += 1

    # tuples of counters, sorted by position
    return [(char_counter[i], rle_counter[i]) for i in xrange(len(char_counter))]

def consensus(seqs, rlelist = None, degap = True):
    """
    Calculate a consensus for an iterable of SeqRecord objects. seqs
    are decoded using corresponding lists of read length counts in
    `rlelist` if provided. Gaps are removed if degap is True.
    """

    if rlelist:
        cons = [cons_char(c) * cons_rle(n) for c, n in get_rle_counts(seqs, rlelist)]
    else:
        cons = [cons_char(c) for c in get_char_counts(seqs)]

    return ''.join(cons).replace(gap,'') if degap else ''.join(cons)

@contextlib.contextmanager
def fasta_tempfile(seqs, tmpdir = None, cleanup = True):
    """Creates a temporary FASTA file representing an iterable of
    objects *seqs* containing attributes `description` and `seq` (eg,
    SeqRecord or fastalite objects). Returns the name of a temorary
    file and, then deletes it at the end of the with block.
    """

    (db_fd, db_name) = tempfile.mkstemp(text=True, dir=tmpdir)
    db_handle = os.fdopen(db_fd, 'w')
    db_handle.write(''.join('>{s.description}\n{s.seq}\n'.format(s=s) for s in seqs))
    db_handle.close()

    try:
        yield db_name
    finally:
        if cleanup:
            os.unlink(db_name)

def run_muscle(seqs, tmpdir = None, keep_order = True):
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
        sortdict = {s.id: i for i,s in enumerate(seqs1)}

    with fasta_tempfile(seqs, tmpdir) as f:
        command = ['muscle', '-quiet', '-seqtype', 'dna', '-in', f]
        pipe = Popen(command, stdout=PIPE)
        (seqstr, _) = pipe.communicate()

    aligned = fastalite(StringIO(seqstr))

    if keep_order:
        aligned = iter(sorted(aligned, key = lambda s: sortdict[s.id]))

    return aligned

def parse_uc(infile):
    """
    Return two dicts: {read_name: cluster} and {cluster: cluster_size}
    """

    cluster_ids = {}
    cluster_sizes = {}

    rows = csv.DictReader(infile, delimiter = '\t', fieldnames = UCLUST_HEADERS)
    for row in rows:
        cluster = int(row['cluster_number'])
        if row['type'] == 'C':
            cluster_sizes[cluster] = int(row['size'])
        else:
            cluster_ids[row['query_label']] = cluster

    return cluster_ids, cluster_sizes

def parse_blast(blast, extras = []):
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

    query = query[_find_homochar_length(ref, gap):] # shift query to preserve alignment
    ref = ref.strip(gap) # strip terminal gaps

    t = h = 0
    m = min(len(query), len(ref))
    results = []

    while t < m:
        h += max(_find_homochar_length(ref[t:]), _find_homochar_length(query[t:]))

        assert (t < h) # strange chars found or miss-alignment

        results.append({
            'i':len(ref[:t].replace(gap, '').replace(homogap, '')),
            'ref':ref[t:h],
            'query':query[t:h]
            })

        t = h

    return results

def _find_homochar_length(s, char = '', ignore = [gap, homogap]):
    '''
    Returns the length of repeating chars from the left side of a string
    '''
    i = 0
    if s:
        char = char or ('' if s[0] in ignore else s[0])
        while char == s[i:i+1]: i += 1
    return i

def error_category(e, gap = gap, homogap = homogap, errors = ERRORS):
    i, r, q = e['i'], e['ref'], e['query']
    if r == q:
        return 'equal'
    elif gap in r or gap in q:
        return errors[1]
    elif len(r) == 1:
        return errors[0]
    elif set(c for c in r if c != homogap) == set(c for c in q if c != homogap):
        return errors[2]
    else:
        return errors[3]

def error_count(errors):
    cnt = Counter()
    for e in errors:
        cnt[error_category(e)] += 1
        cnt['length'] += len(e['ref'].strip('-='))

        ### homoindel counts..
        #lr = len(d['ref'].strip(homogap))
        #lq = len(d['query'].strip(homogap))
        #cnt[(lr if lr <= args.max else gtceil, lq if lq <= args.max else gtceil)] += 1

    return cnt

def parse_ssearch36(lines, numeric = False):
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
        if line.startswith('>>><<<'):
            # query end
            keeplines = False
        elif line.startswith('>>>'):
            # start of a new hit
            if not line.startswith('>>>///'):
                query_count +=1
            q_name = line.lstrip('>').split(',')[0]
        elif line.startswith('>>') or line.startswith('>--'):
            # hit-specific results; keep results starting here
            if prefix:
                yield hit
            hit_count += 1
            if line.startswith('>>'):
                t_description = line[2:].rstrip('\n')
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
            k = k.replace(gap,'').replace(' ','_').lower()
            if k == 'al_cons':
                hit[k] = ''
            else:
                hit[prefix + k] = cast(v) if numeric else v.strip()
        elif prefix and keeplines:
            if 'al_cons' in hit:
                hit['al_cons'] += line.strip()
            else:
                hit[prefix + 'seq'] += line.strip()

    yield hit
    log.info('%s queries, %s hits' % (query_count, hit_count))

def parse_primer_alignments(hits, lprimer = 'lprimer', rprimer = 'rprimer'):
    """
    Hits is a sequence of dictionaries containing alignment results
    for a single query. 'lprimer' and 'rprimer' are left and right
    primer names, respectively. It is assumed that 'q_*' identifies
    attributes of the query sequence, and 't_*' identifies attributes
    of the primer.
    """

    d = {'l':{}, 'r':{}}
    # save the first (shoud be highest-scoring) alignment for each
    # primer
    for hit in hits:
        for label, primer in [('l', lprimer), ('r', rprimer)]:
            if not d[label] and hit['t_name'] == primer:
                d[label]['start']= int(hit['q_al_start']) - 1
                d[label]['stop'] = int(hit['q_al_stop'])
                d[label]['align'] = (hit['q_seq'], hit['t_seq'])
                d[label]['sw_expect'] = float(hit['sw_expect'])
                d[label]['sw_ident'] = float(hit['sw_ident'])
                d[label]['sw_zscore'] = float(hit['sw_zscore'])
                d[label]['sw_overlap'] = int(hit['sw_overlap'])

                # for key in ['sw_expect','sw_zscore','sw_ident','sw_overlap']:
                #     d[label][key] = cast(hit[key])

    d['name'] = hit['q_name']
    return d

def wrap(text, width=60):
    """
    Wraps input string [text] to [width] characters. Return a list of
    substrings.
    """

    r1 = range(0, len(text)+1, width)
    r2 = range(width, len(text) + width +1, width)
    return [text[f:t] for f,t in zip(r1,r2)]


def format_alignment(seq1, seq2, name1 = '', name2 = '', seqwidth = 80, namewidth = 15):

    izl = lambda s: izip_longest(*s, fillvalue = '')

    seqs = [seq1,
            ''.join((':' if c1 == c2 else ' ') for c1,c2 in izl([seq1, seq2])),
            seq2]

    names = [name1, '', name2]

    fstr = '%%%ss %%s' % namewidth
    ss = []
    for strands in izl([wrap(s, seqwidth) for s in seqs]):
        for name, strand in zip(names, strands):
            ss.append(fstr % (name[:namewidth], strand))
        ss.append('')

    return '\n'.join(ss)

def grouper(n, iterable, pad = False):

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

def show_errors(errorlist, width = 40):
    errors = sorted(errorlist, key = lambda e: e['i'])
    out = ''
    for group in grouper(width, errors):
        ii, rr, qq, cc = [],[],[],[]
        for d in group:
            i, r, q = str(d['i']), d['ref'], d['query']
            fstr = '%-' + str(max([len(r), len(q), len(i)])) + 's'
            ii.append(fstr % i)
            rr.append(fstr % r)
            qq.append(fstr % q)
            c = error_category(d)[0]
            cc.append(fstr % ('' if c =='e' else c))

        out += '\n  i:   ' + '|'.join(ii)
        out += '\nref:   ' + '|'.join(rr)
        out += '\nquery: ' + '|'.join(qq)
        out += '\n     : ' + ' '.join(cc)
        out += '\n'

    return out

def format_taxonomy(names, selectors, asterisk = '*'):
    """
    Create a friendly formatted string of taxonomy names. Names will
    have an asterisk value appended *only* if the cooresponding
    element in the selectors evaluates to True.
    """

    taxons = defaultdict(lambda: defaultdict(bool))

    for i,name in enumerate(names):
        name = name.split(None, 1)
        subject = name[0]
        predicate = name[1] if name[1:] else None
        selector = selectors[i] if selectors[i:] else False
        taxons[subject][predicate] |= selector

    taxonomy = []
    for tax, preds in sorted(taxons.items()):
        onomy = []
        for pred, has_asterisk in sorted(preds.items()):
            mark = asterisk if has_asterisk else ''
            if pred:
                onomy.append(pred + mark)
            else:
                tax += mark
        taxonomy.append('{} {}'.format(tax, '/'.join(onomy)).strip())

    return ';'.join(taxonomy)

def _readfasta(handle):
    """
    Lightweight fasta parser. Returns iterator of namedtuple instances
    with fields (id, description, seq) given file-like object `handle`.
    """

    seqlite = namedtuple('SeqLite', 'id, description, seq')
    for h in handle.read().lstrip('> ').split('\n>'):
        part = h.find('\n')
        yield seqlite(h[:part].split()[0], h[:part], ''.join(h[part:].split()))

def _iterfasta(handle):
    """
    Lightweight fasta parser. Returns iterator of namedtuple instances
    with fields (id, description, seq) given file-like object `handle`.
    """

    seqlite = namedtuple('SeqLite', 'id, description, seq')
    name, seq = '', ''
    for line in handle:
        if line.startswith('>'):
            if name:
                yield seqlite(name.split()[0], name, seq)
            name, seq = line[1:].strip(), ''
        else:
            seq += line.strip()

    yield seqlite(name.split()[0], name, seq)

def fastalite(handle, readfile = True):
    return _readfasta(handle) if readfile else _iterfasta(handle)

### Taken from Connor McCoy's Deenurp
def tax_of_genbank(gb):
    """
    Get the tax id from a genbank record, returning None if no taxonomy is
    available.
    """
    # Check for bad name
    try:
        source = next(i for i in gb.features if i.type == 'source')
        return re.findall('\d+', next(i[6:] for i in source.qualifiers.get('db_xref', [])
                     if i.startswith('taxon:')))[0]
    except StopIteration:
        return None

def count_ambiguous(seq):
    s = frozenset('ACGT')
    return sum(i not in s for i in seq)

def is_type(gb):
    """
    Returns a boolean indicating whether a sequence is a member of a type strain,
    as indicated by the presence of the string '(T)' within the record description.
    """
    return '(T)' in gb.description
###
