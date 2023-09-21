import os
import sys
import json
import pprint

from mongtool.database import Database
from mongtool.validate import Validate
from mongtool.utils import Utils
from mongtool.missing import Missing

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
        output_fpaths = self._get_output_fpaths(options.query, options.output_dir, options.output_file, options.prefix, options.combined_output)
        for query_idx, query in enumerate(options.query):
            find = list(Database.find(options.db_collection, {"id": query}))
            if not find:
                find = list(Database.find(options.db_collection, {"sample_id": query}))
            pp = pprint.PrettyPrinter(indent=4)
            pp.pprint(find)
            #with open(output_fpaths[query_idx], 'w+') as fout:
                #json.dump(find, fout)

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

    def missing(self, options):
        utils = Utils()
        missing = Missing()
        if options.sample_sheet:
            csv_dict = missing.parse_sample_sheet(options.input_file[0], options.restore_dir)
            utils.write_out_csv(csv_dict, options.assay, options.platform, options.output_file)
        if options.analysis_dir:
            log_fpath = os.path.splitext(options.missing_log)[0] + ".log"
            empty_fpath = os.path.splitext(options.output_file)[0] + "_empty.csv"
            meta_dict = missing.parse_mongodb_csv(options.input_file[0])
            analysis_dir_fnames = missing.parse_dir(options.analysis_dir)
            csv_dict, missing_samples_txt = missing.find_missing(meta_dict, analysis_dir_fnames, options.restore_dir)
            empty_files_dict, csv_dict = missing.remove_empty_files(csv_dict)
            utils.write_out_csv(csv_dict, options.assay, options.platform, options.output_file)
            utils.write_out_csv(empty_files_dict, options.assay, options.platform, empty_fpath)
            utils.write_out_txt(missing_samples_txt, log_fpath)
        if options.restore_file:
            bash_fpath = os.path.splitext(options.restore_file)[0] + ".sh"
            bash_script = missing.create_bash_script(csv_dict, options.restore_dir)
            utils.write_out_txt(bash_script, bash_fpath)

    def parse_options(self, options):
        if options.subparser_name == 'find':
            self.find(options)

        elif options.subparser_name == 'insert':
            self.insert(options)

        elif options.subparser_name == 'validate':
            self.validate(options)

        elif options.subparser_name == 'missing':
            self.missing(options)