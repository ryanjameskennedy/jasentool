"""__main__ file that handles help and cli execution"""

import sys

from jasentool import __author__, __copyright__, __version__
from jasentool.cli import get_main_parser
from jasentool.main import OptionsParser

def print_help():
    """Print help string for jasentool software"""
    print(f'''

                    ...::: Jasentool v{__version__} :::...
Author(s): {__author__}

Description:
    This software is a mongodb tool that fetches, inserts and 
    removes specific sample data. Furthermore, it validates new
    pipeline data against that of the old data stored in mongodb 
    database.

Usage: jasentool <method> <options> 

Information:
    -V, --version       Display the version of Jasentool and exit.
    -h,  --help         Print help.

Methods:
    find                Find samples in given db.
    insert              Compare new pipeline data with old data.
    missing             Find data missing from Bonsai compared to cgviz.
    validate            Compare new pipeline data with old data.
    convert             Convert cgmlst.org target files to bed files.
    fix                 Fix output files from bjorn.
    converge            Converge tuberculosis mutation catlogues.
    qc                  Extract QC values after alignment.
''')

def main():
    """Main function that handles cli"""
    args = None
    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {'-V', '--version'}:
        print(f"Jasentool version {__version__} {__copyright__} {__author__}")
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
    except Exception as error_code:
        error_message = 'Uncontrolled exit resulting from an unexpected error.\n\n'
        error_message += '-' * 80 + '\n'
        error_message += f'EXCEPTION: {type(error_code).__name__}\n'
        error_message += f'MESSAGE: {error_code}\n'
        error_message += '-' * 80 + '\n\n'
        print(error_message)
        sys.exit(1)

if __name__ == "__main__":
    main()
