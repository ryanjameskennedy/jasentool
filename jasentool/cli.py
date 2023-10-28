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

def __input_file(group, required, help):
    group.add_argument('-i', '--input_file', nargs='+', help=help)

def __output_file(group, required, help):
    group.add_argument('-o', '--output_file', required=required, type=str, help=help)

def __output_dir(group, required):
    group.add_argument('--output_dir', required=required, type=str, help='directory to output files')

def __analysis_dir(group, required):
    group.add_argument('--analysis_dir', required=required, type=str, help='analysis results dir containing jasen results')

def __restore_dir(group, required):
    group.add_argument('--restore_dir', required=required, type=str, default='/fs2/seqdata/restored', help='directory user wishes spring files to be restored to')

def __restore_file(group, required):
    group.add_argument('--restore_file', required=required, type=str, help='filepath bash shell script (.sh) to be output')

def __missing_log(group, required):
    group.add_argument('--missing_log', required=required, type=str, default='missing_samples.log', help='file containing missing files')

def __assay(group, required):
    group.add_argument('--assay', required=required, type=str, default='jasen-saureus-dev', help='assay for jasen to run')

def __platform(group, required):
    group.add_argument('--platform', required=required, type=str, default='illumina', help='sequencing platform for jasen to run')

def __uri(group):
    group.add_argument('--address', '--uri', default='mongodb://localhost:27017/', help='Mongodb host address. Use: `sudo lsof -iTCP -sTCP:LISTEN | grep mongo` to get address')

def __db_name(group, required):
    group.add_argument('--db_name', required=required, help='Mongodb database name address. Use: `show dbs` to get db name')

def __db_collection(group, required):
    group.add_argument('--db_collection', required=required, help='Mongodb collection name. Use: `show collections` to get db collection')

def __out_format(group, required):
    group.add_argument('-f', '--out_format', required=required, type=str, default="bed", help='output format')

def __accession(group, required):
    group.add_argument('-a', '--accession', required=required, type=str, help='accession number')

def __prefix(group):
    group.add_argument('--prefix', type=str, default='jasentool_results_', help='prefix for all output files')

def __combined_output(group):
    group.add_argument('--combined_output', dest='combined_output', action='store_true', help='combine all of the outputs into one output')

def __sample_sheet(group, required):
    group.add_argument('--sample_sheet', required=required, dest='sample_sheet', action='store_true', help='sample sheet input')

def __help(group):
    group.add_argument('-h', '--help', action='help', help='show help message')

def get_main_parser():
    main_parser = argparse.ArgumentParser(prog='jasentool', conflict_handler='resolve')
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
            __combined_output(group)
            __uri(group)
            __prefix(group)
            __help(group)

    with subparser(sub_parsers, 'insert', 'Insert sample(s) into db') as parser:
        with mutex_group(parser, required=True) as group:
            __input_file(group, required=False, help='path to json file to be inserted into db')
            __input_dir(group, required=False, help='path to directory containing sample files')
        with arg_group(parser, 'required named arguments') as group:
            __db_name(group, required=True)
            __db_collection(group, required=True)
        with arg_group(parser, 'optional arguments') as group:
            __combined_output(group)
            __uri(group)
            __help(group)

    with subparser(sub_parsers, 'validate', 'Compare results from new pipeline to old results') as parser:
        with mutex_group(parser, required=True) as group:
            __input_file(group, required=False, help='input filepath(s)')
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

    with subparser(sub_parsers, 'missing', 'Find missing sample data from old runs') as parser:
        with arg_group(parser, 'required named arguments') as group:
            __input_file(group, required=False, help='''path to cgviz meta csv file, created using: mongoexport --quiet --db=cgviz --collection=sample --type=csv --fields=id,mlst.sequence_type,aribavir.lukF_PV.present,aribavir.lukS_PV.present,missing,metadata.QC,metadata.Comment,run --query='{"metadata.QC":"OK"}' | grep -v FOHM | sed "1s/id,mlst.sequence_type,aribavir.lukF_PV.present,aribavir.lukS_PV.present,missing,metadata.QC,metadata.Comment,run/id,mlst,lukF_PV,lukS_PV,missing,QC,Comment,run/" > cgviz_meta.csv''')
            __output_file(group, required=False, help='path to mongo db output file')
        with arg_group(parser, 'optional arguments') as group:
            __analysis_dir(group, required=False)
            __restore_dir(group, required=False)
            __restore_file(group, required=False)
            __missing_log(group, required=False)
            __assay(group, required=False)
            __platform(group, required=False)
            __sample_sheet(group, required=False)
            __help(group)

    with subparser(sub_parsers, 'convert', 'Convert file format') as parser:
        with arg_group(parser, 'required named arguments') as group:
            __input_file(group, required=False, help='path to targets tsv file')
            __output_file(group, required=False, help='path to mongo db output file')
        with arg_group(parser, 'optional arguments') as group:
            __out_format(group, required=False)
            __accession(group, required=False)
            __help(group)

    return main_parser