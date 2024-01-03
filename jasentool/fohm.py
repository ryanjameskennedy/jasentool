import os
import pandas as pd

class Fohm:
    def __init__(self, download_dir):
        self.download_dir = download_dir
        self.fohm_filepath = os.path.join(download_dir, "fohm.csv")

    def convert_colour():
        from openpyxl import load_workbook
        excel_file = 'color_codes.xlsx' 
        wb = load_workbook(excel_file, data_only = True)
        sh = wb['Sheet1']
        color_in_hex = sh['A2'].fill.start_color.index # this gives you Hexadecimal value of the color
        print ('HEX =',color_in_hex) 
        print('RGB =', tuple(int(color_in_hex[i:i+2], 16) for i in (0, 2, 4))) # Color in RGB
    
    def read_file(self, csv_filepath, xlsx_filepath):
        catalogue = pd.read_csv(csv_filepath, header=True)
        catalogue = pd.read_excel(xlsx_filepath, sheet_name='Mutation_catalogue', header=[0,1]).set_index([('variant (common_name)', 'Unnamed: 2_level_1')])
        return catalogue

    def convert2hgvs(self, mutation):
            if mutation[:3].isalpha() and mutation[0].isupper():
                return f'p.{mutation}'
            elif mutation[0].isalpha() and mutation[0].islower() and not mutation[1].isalpha():
                if 'Stop' in mutation:
                    mutation.replace('Stop', '*')
                ref = mutation[0].upper()
                alt = mutation[-1].upper()
                pos = mutation[1:-1]
                return f'c.{pos}{ref}>{alt}'
            else:
                return mutation
    
    def _parse(self):
        catalogue = pd.read_csv(self.fohm_filepath, header=True)
        catalogue['Mutation'] = catalogue.Mutation.apply(self.convert2hgvs)
