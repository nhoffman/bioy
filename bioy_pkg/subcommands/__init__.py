import os
from multiprocessing import cpu_count


# global actions
def action(args):
    return args


# global args
def parse_args(parser):
    default_threads = int(os.environ.get('THREADS_ALLOC') or cpu_count())

    parser.add_argument('--threads',
                        metavar='NUM',
                        default=default_threads,
                        type=int,
                        help="""Number of threads (CPUs). Can also specify
                                with environment variable THREADS_ALLOC.
                                [%(default)s]""")
    return parser
