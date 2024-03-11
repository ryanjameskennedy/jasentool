"""Module for executing each module/class"""

import os
import sys
import json
import pprint

from jasentool.database import Database
from jasentool.validate import Validate
from jasentool.utils import Utils
from jasentool.missing import Missing
from jasentool.convert import Convert
from jasentool.fix import Fix
from jasentool.converge import Converge
from jasentool.qc import QC

class OptionsParser:
    """Class that parses through cli arguments and executes respective modules"""
    def __init__(self, version):
        """Initiate OptionsParser class"""
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
        """Find entry in mongodb"""
        Database.initialize(options.db_name)
        output_fpaths = self._get_output_fpaths(options.query, options.output_dir,
                                                options.output_file, options.prefix,
                                                options.combined_output)
        for query_idx, query in enumerate(options.query):
            find = list(Database.find(options.db_collection, {"id": query}, {}))
            if not find:
                find = list(Database.find(options.db_collection, {"sample_id": query}, {}))
            sample_pp = pprint.PrettyPrinter(indent=4)
            sample_pp.pprint(find)
            with open(output_fpaths[query_idx], 'w+', encoding="utf-8") as fout:
                json.dump(find, fout)

    def insert(self, options):
        """Insert entry in mongodb"""
        Database.initialize(options.db_name)
        input_files = self._input_to_process(options.input_file, options.input_dir)
        for input_file in input_files:
            with open(input_file, 'r', encoding="utf-8") as fin:
                input_sample = json.load(fin)
                Database.insert(options.db_collection, input_sample)

    def validate(self, options):
        """Execute validation of old vs new pipeline results"""
        Database.initialize(options.db_name)
        input_files = self._input_to_process(options.input_file, options.input_dir)
        output_fpaths = self._get_output_fpaths(input_files, options.output_dir,
                                                options.output_file, options.prefix,
                                                options.combined_output)
        validate = Validate()
        validate.run(input_files, output_fpaths, options.db_collection, options.combined_output)

    def missing(self, options):
        """Execute search for missing samples from new pipeline results"""
        utils = Utils()
        missing = Missing()
        db = Database()
        db.initialize(options.db_name)
        if options.sample_sheet:
            csv_dict = missing.parse_sample_sheet(options.input_file[0], options.restore_dir)
            utils.write_out_csv(csv_dict, options.assay, options.platform, options.output_file)
        if options.analysis_dir:
            log_fpath = os.path.splitext(options.missing_log)[0] + ".log"
            empty_fpath = os.path.splitext(options.output_file)[0] + "_empty.csv"
            meta_dict = db.find(options.db_collection, {"metadata.QC": "OK"}, db.get_meta_fields())
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

    def convert(self, options):
        """Execute conversion of file formats"""
        utils = Utils()
        convert = Convert()
        input_file = options.input_file[0]
        output_fpath = os.path.splitext(options.output_file)[0] + "." + options.out_format
        in_format = os.path.splitext(input_file)[1].lstrip(".")
        if in_format == "tsv" and options.out_format == "bed":
            output_txt = convert.targets2bed(input_file, options.accession)
            utils.write_out_txt(output_txt, output_fpath)

    def fix(self, options):
        """Execute fixing of file to desired format(s)"""
        utils = Utils()
        fix = Fix()
        csv_files, assays = fix.fix_csv(options.csv_file, options.output_file)
        batch_files = fix.fix_sh(options.sh_file, options.output_file, assays)
        if (options.remote or options.auto_start) and batch_files:
            utils.copy_batch_and_csv_files(batch_files, csv_files, options.remote_dir, options.remote_hostname, options.auto_start or options.remote)
            if options.auto_start:
                utils.start_remote_pipelines(batch_files, options.remote_dir)

    def converge(self, options):
        """Execute convergence of mutation catalogues"""
        converge = Converge(options.output_dir)
        converge.run(options.save_dbs)

    def qc(self, options):
        """Execute retrieval of qc results"""
        qc = QC(options)
        json_result = qc.run()
        qc.write_json_result(json_result, options.output_file)

    def parse_options(self, options):
        """Options parser"""
        if options.subparser_name == 'find':
            self.find(options)

        elif options.subparser_name == 'insert':
            self.insert(options)

        elif options.subparser_name == 'validate':
            self.validate(options)

        elif options.subparser_name == 'missing':
            self.missing(options)

        elif options.subparser_name == 'convert':
            self.convert(options)

        elif options.subparser_name == 'fix':
            self.fix(options)

        elif options.subparser_name == 'converge':
            self.converge(options)

        elif options.subparser_name == 'qc':
            self.qc(options)
