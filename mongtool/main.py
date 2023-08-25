import os
import sys
import json

from mongtool.database import Database
from mongtool.validate import Validate

class OptionsParser(object):
    def __init__(self, version):
        self.version = version
        self._check_python()

    def _check_python(self):
        if sys.version_info.major < 3:
            print('Python 2 is no longer supported.')
            sys.exit(1)

    def _traverse_input_dir(self, input_dir):
        return [os.path.join(input_dir, filename) for filename in os.listdir(input_dir) if filename.endswith("result.json")]

    def _input_to_process(self, input_file, input_dir):
        if input_dir:
            valid_inputs = self._traverse_input_dir(input_dir)
        elif input_file:
            valid_inputs = input_file
        return valid_inputs

    def _get_output_fpaths(self, input_files, output_dir, output_file, prefix, combined_output):
        output_fpaths = []
        if output_dir:
            output_dir = os.path.expanduser(output_dir)
            if combined_output: 
                output_fpaths = [os.path.join(output_dir, prefix + "combined_outputs")]
            else:
                output_fpaths = [os.path.join(output_dir, prefix + os.path.basename(os.path.splitext(input_fpath)[0])) for input_fpath in input_files]
        elif output_file:
            if len(input_files) > 1:
                print(f'ERROR: You have input {len(input_files)} input_files and provided an outfile instead of an out directory. Use --outdir instead.')
                sys.exit(1)
            output_fpaths = [os.path.splitext(output_file)[0]]
        return output_fpaths

    def find(self, options):
        Database.initialize(options.db_name)
        output_fpaths = self._get_output_fpaths(options.query, options.output_dir, options.output_file, options.prefix)
        for query_idx, query in enumerate(options.query):
            find = json.dumps(Database.find(options.db_collection, query))
            print(find)
            with open(output_fpaths[query_idx], 'w+') as fout:
                json.dump(find, fout)

    def insert(self, options):
        Database.initialize(options.db_name)
        input_files = self._input_to_process(options.input_file, options.input_dir)
        for input_file in input_files:
            with open(input_file, 'r') as fin:
                input_sample = json.load(fin)
                Database.insert(options.db_collection, input_sample)

    def validate(self, options):
        Database.initialize(options.db_name)
        input_files = self._input_to_process(options.input_file, options.input_dir)
        output_fpaths = self._get_output_fpaths(input_files, options.output_dir, options.output_file, options.prefix, options.combined_output)
        validate = Validate()
        validate.run(input_files, output_fpaths, options.db_collection, options.combined_output)

    def parse_options(self, options):
        if options.subparser_name == 'find':
            self.find(options)

        elif options.subparser_name == 'insert':
            self.insert(options)

        elif options.subparser_name == 'validate':
            self.validate(options)