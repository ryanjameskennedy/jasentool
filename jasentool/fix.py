import os
import pandas as pd
from jasentool.utils import Utils

class Fix(object):
    @staticmethod
    def fix_csv(input_file, output_fpath):
        assays = []
        out_fpaths = []
        with open(input_file, 'r') as csvfile:
            df = pd.read_csv(csvfile)
            df['assay'] = df['species']
            for assay, df_assay in df.groupby('assay'):
                out_fpath = f'{os.path.splitext(output_fpath)[0]}_{assay}.csv'
                df_assay.to_csv(out_fpath, encoding='utf-8')
                out_fpaths.append(out_fpath)
                assays.append(assay)
        return out_fpaths, assays

    @staticmethod
    def fix_sh(input_file, output_fpath, assays):
        utils = Utils()
        output_txt = ""
        out_fpaths = []
        with open(input_file, 'r') as shfile:
            for line in shfile:
                line = line.rstrip()
                if line.startswith('/fs2/sw/bnf-scripts/start_nextflow_analysis.pl'):
                    for assay in assays:
                        line = f'/fs2/sw/bnf-scripts/start_nextflow_analysis.pl  $SCRIPTPATH/{os.path.splitext(output_fpath)[0]}_{assay}.csv'
                        out_fpath = f'{os.path.splitext(output_fpath)[0]}_{assay}.sh'
                        output_txt += line+'\n'
                        utils.write_out_txt(output_txt, out_fpath)
                        out_fpaths.append(out_fpath)
                    break
                output_txt += line

        return out_fpaths