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
Tools for microbial sequence analysis and classification.

Assembles subcommands and provides top-level script.
"""

import argparse
import importlib
import logging
import os
import pkgutil
import subcommands
import sys

log = logging.getLogger(__name__)

_data = os.path.join(os.path.dirname(__file__), 'data')

try:
    with open(os.path.join(_data, 'version')) as v:
        __version__ = v.read().strip()
except Exception, e:
    log.warn(e)
    __version__ = ''


def main(argv):
    action, arguments = parse_arguments(argv)

    loglevel = {
        0: logging.ERROR,
        1: logging.WARNING,
        2: logging.INFO,
        3: logging.DEBUG,
    }.get(arguments.verbosity, logging.DEBUG)

    if arguments.verbosity > 1:
        logformat = '%(levelname)s %(module)s %(lineno)s %(message)s'
    else:
        logformat = '%(message)s'

    # set up logging
    logging.basicConfig(file=sys.stdout, format=logformat, level=loglevel)

    return action(arguments)


def parse_arguments(argv):
    """
    Create the argument parser
    """

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('-V', '--version', action='version',
                        version=__version__,
                        help='Print the version number and exit')

    parser.add_argument('-v', '--verbose',
                        action='count', dest='verbosity', default=1,
                        help='Increase verbosity of screen output '
                             '(eg, -v is verbose, -vv more so)')
    parser.add_argument('-q', '--quiet',
                        action='store_const', dest='verbosity', const=0,
                        help='Suppress output')

    ##########################
    # Setup all sub-commands #
    ##########################

    subparsers = parser.add_subparsers(dest='subparser_name')

    # Begin help sub-command
    parser_help = subparsers.add_parser(
        'help', help='Detailed help for actions using `help <action>`')
    parser_help.add_argument('action', nargs=1)
    # End help sub-command

    # Organize submodules by argv
    modules = [
        name for _, name, _ in pkgutil.iter_modules(subcommands.__path__)]
    commands = [m for m in modules if m in argv]

    actions = {}

    # `commands` will contain the module corresponding to a single
    # subcommand if provided; otherwise, generate top-level help
    # message from all submodules in `modules`.
    for name in commands or modules:
        # set up subcommand help text. The first line of the dosctring
        # in the module is displayed as the help text in the
        # script-level help message (`script -h`). The entire
        # docstring is displayed in the help message for the
        # individual subcommand ((`script action -h`))
        # if no individual subcommand is specified (run_action[False]),
        # a full list of docstrings is displayed
        try:
            imp = '{}.{}'.format(subcommands.__name__, name)
            mod = importlib.import_module(imp)
        except Exception, e:
            log.error('error importing module {}'.format(imp))
            log.error(e)
            continue

        subparser = subparsers.add_parser(
            name,
            help=mod.__doc__.lstrip().split('\n', 1)[0],
            description=mod.__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
        subparser = subcommands.parse_args(subparser)
        mod.build_parser(subparser)
        actions[name] = mod.action
    # Determine we have called ourself (e.g. "help <action>")
    # Set arguments to display help if parameter is set
    #           *or*
    # Set arguments to perform an action with any specified options.
    arguments = parser.parse_args(argv)
    # Determine which action is in play.
    action = arguments.subparser_name

    # Support help <action> by simply having this function call itself and
    # translate the arguments into something that argparse can work with.
    if action == 'help':
        return parse_arguments([str(arguments.action[0]), '-h'])  # loop back
    else:
        return actions[action], arguments
