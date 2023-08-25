import argparse
from contextlib import contextmanager

@contextmanager
def subparser(parser, name, desc):
    yield parser.add_parser(name, conflict_handler='resolve', help=desc, formatter_class=argparse.RawDescriptionHelpFormatter)

@contextmanager
def mutex_group(parser, required):
    group = parser.add_argument_group(f'mutually exclusive {"required" if required else "optional"} arguments')
    yield group.add_mutually_exclusive_group(required=required)

@contextmanager
def arg_group(parser, name):
    yield parser.add_argument_group(name)

def __query(group, required):
    group.add_argument('-q', '--query', required=required, nargs='+', help='sample query')

def __input_dir(group, required, help):
    group.add_argument('--input_dir', required=required, help=help)

def __input_file(group, required):
    group.add_argument('-i', '--input_file', nargs='+', help='input filepath(s)')

def __output_file(group, required, help):
    group.add_argument('-o', '--output_file', required=required, type=str, help=help)

def __output_dir(group, required):
    group.add_argument('--output_dir', required=required, type=str, help='directory to output files')

def __uri(group):
    group.add_argument('--address', '--uri', default='mongodb://localhost:27017/', help='Mongodb host address. Use: `sudo lsof -iTCP -sTCP:LISTEN | grep mongo` to get address')

def __db_name(group, required):
    group.add_argument('--db_name', required=required, help='Mongodb database name address. Use: `show dbs` to get db name')

def __db_collection(group, required):
    group.add_argument('--db_collection', required=required, help='Mongodb collection name. Use: `show collections` to get db collection')

def __prefix(group):
    group.add_argument('--prefix', type=str, default='mongtool_results_', help='prefix for all output files')

def __combined_output(group):
    group.add_argument('--combined_output', dest='combined_output', action='store_true', help='combine all of the outputs into one output')

def __help(group):
    group.add_argument('-h', '--help', action='help', help='show help message')

def get_main_parser():
    main_parser = argparse.ArgumentParser(prog='mongtool', conflict_handler='resolve')
    sub_parsers = main_parser.add_subparsers(help='--', dest='subparser_name')
    with subparser(sub_parsers, 'find', 'Find sample from given mongo db') as parser:
        with mutex_group(parser, required=True) as group:
            __output_dir(group, required=False)
            __output_file(group, required=False, help='path to mongo db output file')
        with arg_group(parser, 'required named arguments') as group:
            __query(group, required=True)
            __db_name(group, required=True)
            __db_collection(group, required=True)
        with arg_group(parser, 'optional arguments') as group:
            __uri(group)
            __prefix(group)
            __help(group)

    with subparser(sub_parsers, 'insert', 'Insert sample(s) into db') as parser:
        with mutex_group(parser, required=True) as group:
            __input_file(group, required=False)
            __input_dir(group, required=False, help='path to directory containing sample files')
        with arg_group(parser, 'required named arguments') as group:
            __db_name(group, required=True)
            __db_collection(group, required=True)
        with arg_group(parser, 'optional arguments') as group:
            __uri(group)
            __help(group)

    with subparser(sub_parsers, 'validate', 'Compare results from new pipeline to old results') as parser:
        with mutex_group(parser, required=True) as group:
            __input_file(group, required=False)
            __input_dir(group, required=False, help='path to directory containing sample files')
        with mutex_group(parser, required=True) as group:
            __output_dir(group, required=False)
            __output_file(group, required=False, help='path to mongo db output file')
        with arg_group(parser, 'required named arguments') as group:
            __db_name(group, required=True)
            __db_collection(group, required=True)
        with arg_group(parser, 'optional arguments') as group:
            __combined_output(group)
            __uri(group)
            __prefix(group)
            __help(group)

    return main_parser