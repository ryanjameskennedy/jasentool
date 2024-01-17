"""Module that fixes csv and shell scripts"""

import os
import pandas as pd
from jasentool.utils import Utils

class Fix:
    """Class that fixes csvs for start_nextflow_analysis.pl"""
    @staticmethod
    def fix_csv(input_file, output_fpath):
        """Convert the provided bjorn csvs into new jasen-compatible csvs"""
        assays = []
        out_fpaths = []
        with open(input_file, 'r', encoding="utf-8") as csvfile:
            samples = pd.read_csv(csvfile)
            samples['assay'] = samples['species']
            for assay, df_assay in samples.groupby('assay'):
                out_fpath = f'{os.path.splitext(output_fpath)[0]}_{assay}.csv'
                df_assay.to_csv(out_fpath, encoding='utf-8', index=False)
                out_fpaths.append(out_fpath)
                assays.append(assay)
        return out_fpaths, assays

    @staticmethod
    def fix_sh(input_file, output_fpath, assays):
        """Fix the shell scripts"""
        utils = Utils()
        output_content = ""
        out_fpaths = []
        with open(input_file, 'r', encoding="utf-8") as shfile:
            for line in shfile:
                line = line.rstrip()
                if line.startswith('/fs2/sw/bnf-scripts/start_nextflow_analysis.pl'):
                    for assay in assays:
                        output_txt = ""
                        line = '/fs2/sw/bnf-scripts/start_nextflow_analysis.pl ' + \
                            f'$SCRIPTPATH/{os.path.splitext(output_fpath)[0]}_{assay}.csv'
                        out_fpath = f'{os.path.splitext(output_fpath)[0]}_{assay}.sh'
                        output_txt += output_content+line+'\n'
                        utils.write_out_txt(output_txt, out_fpath)
                        out_fpaths.append(out_fpath)
                    break
                output_content += line+'\n'

        return out_fpaths
