import os
import re
import pandas as pd
from tqdm import tqdm
from jasentool.utils import Utils

class WHO(object):
    def __init__(self):
        self.aa_dict_1 = self.get_aa_dict()
        self.aa_dict_2 = self.inv_dict()
        self.nucleotide_complements = self.get_nt_complements()
        self.drug_dict = self.get_drug_dict()
        self.re_c, self.re_p, self.re_d, self.re_i = self.setup_re()
        self.re_attr = re.compile('Name=([^;]+).*locus_tag=([^;|\n]+)')

    def get_aa_dict(self):
        return {
            'Ala': 'A',
            'Arg': 'R',
            'Asn': 'N',
            'Asp': 'D',
            'Asx': 'B',
            'Cys': 'C',
            'Glu': 'E',
            'Gln': 'Q',
            'Glx': 'Z',
            'Gly': 'G',
            'His': 'H',
            'Ile': 'I',
            'Leu': 'L',
            'Lys': 'K',
            'Met': 'M',
            'Phe': 'F',
            'Pro': 'P',
            'Ser': 'S',
            'Thr': 'T',
            'Trp': 'W',
            'Tyr': 'Y',
            'Val': 'V',
            '*': '!',
        }
    
    def inv_dict(self):
        return {v: k for k, v in self.aa_dict_1.items()}

    def get_nt_complements(self):
        return {
            'C': 'G',
            'G': 'C',
            'T': 'A',
            'A': 'T',
        }

    def get_drug_dict(self):
        return {
            'RIF': 'rifampicin',
            'INH': 'isoniazid',
            'EMB': 'ethambutol',
            'PZA': 'pyrazinamide',
            'LEV': 'levofloxacin',
            'MXF': 'moxifloxacin',
            'BDQ': 'bedaquiline',
            'LZD': 'linezolid',
            'CFZ': 'clofazimine',
            'DLM': 'delamanid',
            'AMI': 'amikacin',
            'STM': 'streptomycin',
            'ETH': 'ethionamide',
            'KAN': 'kanamycin',
            'CAP': 'capreomycin'
        }

    def setup_re(self):
        # Setup the regular expressions
        re_c = re.compile('^(\w+)_([actg])(-*\d+)([actg])$') #regex pattern for 
        re_p = re.compile('^(\w+)_([A-Z])(\d+)([A-Z!])$') #regex pattern for protein
        re_d = re.compile('^(\w+)_(-*\d+)_del_(\d+)_([actg]+)_([actg]+)$') #regex pattern for deletions
        re_i = re.compile('^(\w+)_(-*\d+)_ins_(\d+)_([actg]+)_([actg]+)$') #regex pattern for insertions
        return re_c, re_p, re_d, re_i

    def lower_row(self, row):
        return row.str.lower()

    def read_files(self, gff_filepath, xlsx_filepath, h37rv_filepath):
        # Load the reference GFF file
        gff = pd.read_csv(gff_filepath, names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'], sep='\t', header=None)
        # Load the WHO catalogue
        catalogue = pd.read_excel(xlsx_filepath, sheet_name='Catalogue_master_file', header=2)
        # Load the reference genome to impute missing data from deletions
        h37rv = ''
        with open(h37rv_filepath, 'r') as fin:
            for line in fin:
                h37rv += line.replace('\n', '')
        return gff, catalogue, h37rv

    def process_variant(self, variant, gff_dict):
        '''Translates variants in the WHO catalogue format to HGVS'''

        c_match = self.re_c.match(variant)
        if c_match:
            if gff_dict[c_match[1]]['type'] == 'rRNA':
                v_type = 'n'
                ref = c_match[2].upper()
                alt = c_match[4].upper()
            else:
                v_type = 'c'
                ref = c_match[2].upper()
                alt = c_match[4].upper()
            return (c_match[1], v_type, '{}.{}{}>{}'.format(v_type, c_match[3], ref, alt), False, None)

        p_match = self.re_p.match(variant)
        if p_match:
            return (p_match[1], 'p', 'p.{}{}{}'.format(self.aa_dict_2[p_match[2].upper()], p_match[3], self.aa_dict_2[p_match[4].upper()]), False, None)

        d_match = self.re_d.match(variant)
        if d_match:
            if int(d_match[3]) != len(d_match[4]) - len(d_match[5]):
                return (None, None, None, True, 'length mismatch')

            starts = [pos for pos in range(1, len(d_match[4]) + 1 - int(d_match[3])) if d_match[4][:pos]+d_match[4][pos+int(d_match[3]):] == d_match[5]]
            if not starts:
                return (None, None, None, True, 'invalid indel')
            if not gff_dict[d_match[1]]['strand']:
                hgvs = []
                for start in starts:
                    if int(d_match[3]) == 1:
                        hgvs.append('c.{}del'.format(int(d_match[2])+start))
                    else:
                        hgvs.append('c.{}_{}del'.format(int(d_match[2])+start, int(d_match[2])+start-1+int(d_match[3])))
                return (d_match[1], 'c', '|'.join(hgvs), False, None)
            else:
                hgvs = []
                for start in starts:
                    if int(d_match[3]) == 1:
                        hgvs.append('c.{}del'.format(int(d_match[2]) - start - int(d_match[3]) + 1))
                    else:
                        v = 'c.{}_{}del'.format(int(d_match[2]) - start - int(d_match[3]) + 1, int(d_match[2]) - start)
                        hgvs.append(v)
                return (d_match[1], 'c', '|'.join(hgvs), False, None)

        i_match = self.re_i.match(variant)
        if i_match:
            if int(i_match[3]) != len(i_match[5]) - len(i_match[4]):
                return (None, None, None, True, 'length mismatch')
            starts = [pos for pos in range(1, len(i_match[4]) + 1) if i_match[4][:pos]+i_match[5][pos:pos+int(i_match[3])]+i_match[4][pos:] == i_match[5]]
            if not starts:
                return (None, None, None, True, 'invalid indel')
            if not gff_dict[i_match[1]]['strand']:
                hgvs = []
                for start in starts:
                    hgvs.append('c.{}_{}ins{}'.format(int(i_match[2])+start-1, int(i_match[2])+start, ''.join([i.upper() for i in i_match[5][start:start+int(i_match[3])]])))
                return (i_match[1], 'c', '|'.join(hgvs), False, None)
            else:
                hgvs = []
                for start in starts:
                    v = 'c.{}_{}ins{}'.format(int(i_match[2])-start, int(i_match[2]) - start+1, ''.join([self.nucleotide_complements[i.upper()] for i in i_match[5][start:start+int(i_match[3])][::-1]]))
                    hgvs.append(v)
                return (i_match[1], 'c', '|'.join(hgvs), False, None)

        return (None, None, None, True, 'does not match indel or variant')

    def extract_info(self, info_string):
        if pd.notna(info_string):
            match = self.re_attr.search(info_string)
            if match:
                if len(match.groups([1, 2])) != 2:
                    print(match.groups([1, 2]))
                return pd.Series(match.groups([1, 2]))
        return pd.Series([None, None])

    def get_gene_info(self, gff):
        # Get the gene information from the GFF file
        # Apply the function to the 'attributes' column
        gff[['locus_tag', 'name']] = gff.attributes.apply(self.extract_info)

        gff_dict = {}
        for _, row in gff.iterrows():
            gene = {}
            gene['seqid'] = row.seqid
            gene['source'] = row.source
            gene['type'] = row.type
            gene['start'] = row.start
            gene['end'] = row.end
            gene['score'] = row.score
            gene['strand'] = 0 if row.strand == '+' else 1
            gene['attributes'] = row.attributes
            gene['locus_tag'] = row.locus_tag
            gene['name'] = row['name']
            gff_dict[row['locus_tag']] = gene
            gff_dict[row['name']] = gene
        return gff_dict

    def prep_catalogue(self, catalogue):
        # Prepare the WHO catalogue dataframe
        classified = []
        v = re.compile('^(.*) \((.*)\)')
        for var, row in catalogue[catalogue[('FINAL CONFIDENCE GRADING', 'Unnamed: 51_level_1')].apply(lambda conf: conf != 'combo')].iterrows():
            drug_key = row[('drug', 'Unnamed: 0_level_1')]
            drug = self.drug_dict[drug_key]
            v_match = v.match(var)
            if v_match:
                # Include all variants listed
                variants = [v_match[1]] + [i.strip() for i in v_match[2].split(',')]
                # Include only first variant (the one on which the analysis was performed)
                variants = [v_match[1]]
            else:
                variants = [var]
            category = ' '.join(row[('FINAL CONFIDENCE GRADING', 'Unnamed: 51_level_1')].split(' ')[1:])
            genome_pos = '{:.0f}'.format(row[('Genome position', 'Unnamed: 3_level_1')])
            for variant in variants:
                classified.append([variant, drug, 'resistance', '', 'https://www.who.int/publications/i/item/9789240028173', category])
        classified = pd.DataFrame(classified, columns=['variant', 'Drug', 'Confers', 'Interaction', 'Literature', 'WHO Confidence'])
        return classified

    def var2hgvs(self, classified, gff_dict):
        # Convert the variants to HGVS format
        for idx, row in tqdm(classified.iterrows(), total=classified.shape[0]):
            gene, var_type, variant, fail, fail_reason = self.process_variant(row.variant, gff_dict)
            classified.loc[idx, 'gene'] = gene
            classified.loc[idx, 'hgvs'] = variant
            classified.loc[idx, 'type'] = var_type
            classified.loc[idx, 'fail'] = fail
            classified.loc[idx, 'fail_reason'] = fail_reason
        return classified

    def impute_del(self, classified, gff_dict, h37rv):
        # Impute missing data for deletions
        length_mismatch = classified[classified.fail_reason == 'length mismatch'].sort_values(by='variant', key=self.lower_row)

        for idx, row in tqdm(length_mismatch.iterrows(), total=length_mismatch.shape[0]):
            d_match = self.re_d.match(row.variant)
            if d_match:
                if not gff_dict[d_match[1]]['strand']:
                    indexing_correction = -1 if int(d_match[2]) < 0 else -2 # correct for 0 based python indexing (-1 if promotor, -2 if within gene)
                    start = int(gff_dict[d_match[1]]['start']) + int(d_match[2]) + int(indexing_correction)
                    end = start + int(d_match[3]) + len(d_match[5]) # add the lenght of the alt allele to account for the bases not part of the indel
                    try:
                        complete_variant = '{}_{}_del_{}_{}_{}'.format(d_match[1], d_match[2], d_match[3], h37rv[start:end].lower(), d_match[5])
                    except TypeError:
                        print(f"{start}: {type(start)}\n{end}: {type(end)}")
                    classified.loc[idx, 'complete_variant'] = complete_variant
                    classified.loc[idx, 'complete_variant_fail'] = False
                else:
                    indexing_correction = -1 if int(d_match[2]) < 0 else 0 # correct for 0 based python indexing (-1 if promotor, 0 if within gene)
                    start = int(gff_dict[d_match[1]]['end']) - int(d_match[2]) + int(indexing_correction) # subtract d_match[2] instead of adding as this is the opposite strand
                    end = start + int(d_match[3]) + len(d_match[5]) # add the lenght of the alt allele to account for the bases not part of the indel
                    complete_variant = '{}_{}_del_{}_{}_{}'.format(d_match[1], d_match[2], d_match[3], h37rv[start:end].lower(), d_match[5])
                    classified.loc[idx, 'complete_variant'] = complete_variant
                    classified.loc[idx, 'complete_variant_fail'] = False
                continue

            i_match = self.re_i.match(row.variant)
            if i_match:
                classified.loc[idx, 'complete_variant_fail'] = True
                classified.loc[idx, 'complete_variant_fail_reason'] = 'Not assuming for insertions'
                if not gff_dict[i_match[1]]['strand']:
                    pass
                else:
                    pass
                continue
        return classified
    
    def imp2hgvs(self, classified, gff_dict):
        # Convert imputed deletions to HGVS format
        for idx, row in tqdm(classified[classified.complete_variant_fail == False].iterrows(), total=classified[classified.complete_variant_fail == False].shape[0]):
            if row.complete_variant_fail:
                continue
            gene, var_type, variant, fail, fail_reason = self.process_variant(row.complete_variant, gff_dict)
            classified.loc[idx, 'gene'] = gene
            classified.loc[idx, 'hgvs'] = variant
            classified.loc[idx, 'type'] = var_type
            classified.loc[idx, 'fail'] = fail
            classified.loc[idx, 'fail_reason'] = fail_reason
        return classified

    def write_out_csv(self, classified, csv_outpath):
        # Write results to csv file
        classified.to_csv(csv_outpath, index=False)

    def _parse(self, fasta_filepath, gff_filepath, download_dir):
        utils = Utils()
        #who_url = "https://apps.who.int/iris/bitstream/handle/10665/341906/WHO-UCN-GTB-PCI-2021.7-eng.xlsx"
        who_url = "https://raw.githubusercontent.com/GTB-tbsequencing/mutation-catalogue-2023/main/Final%20Result%20Files/WHO-UCN-TB-2023.5-eng.xlsx"
        who_filepath = os.path.join(download_dir, "who.xlsx")
        utils.download_and_save_file(who_url, who_filepath)
        gff, catalogue, h37rv = self.read_files(gff_filepath, who_filepath, fasta_filepath)
        gff_dict = self.get_gene_info(gff)
        catalogue.columns = catalogue.columns.str.title()
        catalogue.rename(columns={'Final Confidence Grading': 'WHO Confidence'}, inplace=True)
        catalogue['Confers'] = 'resistance'
        catalogue['Interaction'] = ''
        catalogue['Literature'] = 'https://www.who.int/publications/i/item/9789240082410'
        catalogue['WHO Confidence'] = catalogue['WHO Confidence'].apply(lambda x: ' '.join(x.split(' ')[1:]))
        catalogue['Drug'] = catalogue['Drug'].apply(lambda x: x.lower())
        catalogue = catalogue.loc[:, ["Drug","Confers","Interaction","Literature","WHO Confidence","Gene","Mutation"]]
        csv_outpath = os.path.join(download_dir, "who.csv")
        self.write_out_csv(catalogue, csv_outpath)
        return catalogue
