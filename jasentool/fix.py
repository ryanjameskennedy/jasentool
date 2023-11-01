import os
import pandas as pd

class Fix(object):
    @staticmethod
    def fix_csv(input_file, output_fpath):
        with open(input_file, 'rb') as csvfile:
            df = pd.read_csv(csvfile)
            df['assay'] = df['species']
            for assay, df_assay in df.groupby('assay'):
                out_fpath = f'{os.path.splitext(output_fpath)[0]}_{assay}.csv'
                df_assay.to_csv(out_fpath, encoding='utf-8')
