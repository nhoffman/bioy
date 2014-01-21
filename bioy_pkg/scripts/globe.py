import os
from multiprocessing import cpu_count

def parse_args(parser):
    parser.add_argument('--threads',
            metavar = 'NUM',
            default = int(os.environ.get('THREADS_ALLOC') or cpu_count()),
            type = int,
            help = """Number of threads (CPUs). Can also specify
                      with environment variable THREADS_ALLOC.
                      (Default = %(default)s)""")
    return parser
