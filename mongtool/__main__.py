import sys
import os

from mongtool import __author__, __copyright__, __version__
from mongtool.cli import get_main_parser
from mongtool.main import OptionsParser

def print_help():
    print('''

                    ...::: Mongtool v%s :::...
Author(s): %s

Description:
    This software is a mongodb tool that fetches, inserts and 
    removes specific sample data. Furthermore, it validates new
    pipeline data against that of the old data stored in mongodb 
    database.

Usage: mongtool <method> <options> 

Information:
    -V, --version       Display the version of Mongtool and exit.
    -h,  --help         Print help.

Methods:
    insert            Compare new pipeline data with old data.
                Compare new pipeline data with old data.
    validate            Compare new pipeline data with old data.
    
''' % (__version__, __author__))

def main():
    args = None
    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {'-V', '--version'}:
        print(f"Mongtool version {__version__} {__copyright__} {__author__}")
        sys.exit(0)
    elif sys.argv[1] in {'-h', '--h', '-help', '--help'}:
        print_help()
        sys.exit(0)
    else:
        args = get_main_parser().parse_args()
        ts_parser = OptionsParser(__version__)
        ts_parser.parse_options(args)
    try:
        print("Done")
    except SystemExit:
        print('Controlled exit resulting from early termination.')
        sys.exit(1)
    except KeyboardInterrupt:
        print('Controlled exit resulting from interrupt signal.')
        sys.exit(1)
    except Exception as e:
        error_message = 'Uncontrolled exit resulting from an unexpected error.\n\n'
        error_message += '-' * 80 + '\n'
        error_message += 'EXCEPTION: {}\n'.format(type(e).__name__)
        error_message += 'MESSAGE: {}\n'.format(e)
        error_message += '-' * 80 + '\n\n'
        print(error_message)
        sys.exit(1)

if __name__ == "__main__":
    main()