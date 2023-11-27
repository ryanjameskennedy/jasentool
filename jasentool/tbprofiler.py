import os
import re
import sys
import pandas as pd
from tqdm import tqdm
from jasentool.utils import Utils

class Tbprofiler(object):
    def __init__(self, tbdb_dir):
        self.tbdb_filepath = os.path.join(tbdb_dir, "tbdb.csv")
        self.chr_name = "Chromosome"
        self.aa_long2short = self.get_aa_dict()

    def fasta2dict(self, filepath):
        fa_dict = {}
        seq_name = ""
        with open(filepath, 'r') as fin:
            for line in fin:
                line = line.rstrip()
                if line.startswith(">"):
                    seq_name = line[1:].split()[0]
                    fa_dict[seq_name] = []
                else:
                    fa_dict[seq_name].append(line)
        return {seq: "".join(fa_dict[seq]) for seq in fa_dict}

    def reverse_complement(self, seq):
        """Return reverse complement of a sequence"""
        def complement(seq):
            basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
            letters = [basecomplement[base] for base in seq]
            return ''.join(letters)
        return complement(seq[::-1])

    def write_gene_pos(self, infile, genes, outfile):
        output_txt = ""
        with open(infile, "r") as fin:
            for line in fin:
                row = line.strip().split()
                rv, gene, chr_start, chr_end, gene_start, gene_end = [row[0], row[1]]+[int(row[i]) for i in range(2,6)]
                if rv in genes:
                    y = 0
                    for i, chr_pos in enumerate(range(chr_start, chr_end+1)):
                        x = 1 if gene_start< gene_end else -1
                        if gene_start+(x*i) == 0:
                            y = 1 if gene_start< gene_end else -1
                        output_txt += "%s\t%s\t%s\t%s\n" % (self.chr_name, chr_pos, rv, gene_start+(x*i)+y)
        with open(outfile, "w") as fout:
            fout.write(output_txt)

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
            "Stop":"*",
            "-":"-"
        }
    
    def parse_mutation(self, mut, gene, fasta_dict, gene_info):
        # AA change
        re_obj = re.search("p.([A-Z][a-z][a-z])([0-9]+)([A-Z][a-z][a-z])", mut)
        if re_obj:
            ref_aa = self.aa_long2short[re_obj.group(1)]
            alt_aa = self.aa_long2short[re_obj.group(3)]
            codon_num = re_obj.group(2)
            return ["%s%s>%s%s" % (codon_num, ref_aa, codon_num, alt_aa)]
        # Stop codon
        re_obj = re.search("p.([A-Z][a-z][a-z])([0-9]+)(\*)", mut)
        if re_obj:
            ref_aa = self.aa_long2short[re_obj.group(1)]
            alt_aa = re_obj.group(3)
            codon_num = re_obj.group(2)
            return ["%s%s>%s%s" % (codon_num, ref_aa, codon_num, alt_aa)]
        # Deletion single base
        re_obj = re.search("c.([\-0-9]+)del", mut)
        if re_obj:
            gene_start_nt = int(re_obj.group(1))
            strand = gene_info[gene]["strand"]
            if strand == "-":
                chr_start_nt = gene_info[gene]["end"] + gene_info[gene]["gene_end"] - gene_start_nt + (1 if gene_info[gene]["gene_end"]<0 else 0)
            else:
                chr_start_nt = gene_info[gene]["start"] - gene_info[gene]["gene_start"] + gene_start_nt - (0 if gene_start_nt<0 else 1)
            seq = fasta_dict["Chromosome"][chr_start_nt-2:chr_start_nt]
            return ["%s%s>%s" % (chr_start_nt-1,seq,seq[0])]
        # Deletion multi base
        re_obj = re.search("c.([\-0-9]+)_([\-0-9]+)del", mut)
        if re_obj:
            gene_start_nt = int(re_obj.group(1))
            gene_end_nt = int(re_obj.group(2))
            del_len = gene_end_nt-gene_start_nt+1
            strand = gene_info[gene]["strand"]
            if strand == "-":
                chr_start_nt = gene_info[gene]["end"] + gene_info[gene]["gene_end"] - gene_start_nt - (del_len-1) + (1 if gene_info[gene]["gene_end"]<0 else 0)
            else:
                chr_start_nt = gene_info[gene]["start"] - gene_info[gene]["gene_start"] + gene_start_nt - (0 if gene_start_nt<0 else 1)
            chr_end_nt = chr_start_nt+del_len-1
            seq = fasta_dict["Chromosome"][chr_start_nt-2:chr_end_nt]
            return ["%s%s>%s" % (chr_start_nt-1, seq, seq[0])]
        # Insertion
        re_obj = re.search("c.([0-9]+)_([0-9]+)ins([A-Z]+)", mut)
        if re_obj:
            gene_start_nt = int(re_obj.group(1))
            seq_ins = re_obj.group(3)
            strand = gene_info[gene]["strand"]
            if strand == "-":
                chr_start_nt = gene_info[gene]["end"] + gene_info[gene]["gene_end"] - gene_start_nt + 1
                seq_ins = self.reverse_complement(seq_ins)
            else:
                chr_start_nt = gene_info[gene]["start"] - gene_info[gene]["gene_start"] + gene_start_nt - 1
            seq_start = fasta_dict["Chromosome"][chr_start_nt-1]
            return ["%s%s>%s" % (chr_start_nt,seq_start,seq_start+seq_ins)]
        # Promoter Mutation
        ## c.-16G>C
        re_obj = re.search("c.(\-[0-9]+)([A-Z])>([A-Z])",mut)
        if re_obj:
            nt_pos = int(re_obj.group(1))
            ref_nt = re_obj.group(2)
            alt_nt = re_obj.group(3)
            strand = gene_info[gene]["strand"]

            if strand == "+":
                chr_pos = gene_info[gene]["start"] - (gene_info[gene]["gene_start"] - nt_pos)
                return ["%s%s>%s" % (nt_pos,ref_nt,alt_nt)]
            else:
                chr_pos = gene_info[gene]["end"] + (gene_info[gene]["gene_end"] - nt_pos)
                return ["%s%s>%s" % (nt_pos, self.reverse_complement(ref_nt), self.reverse_complement(alt_nt))]
        # ncRNA Mutation
        ## r.514a>c
        re_obj = re.search("r.([0-9]+)([a-z]+)>([a-z]+)",mut)
        if re_obj:
            nt_pos = re_obj.group(1)
            ref_nt = re_obj.group(2)
            alt_nt = re_obj.group(3)
            return ["%s%s>%s" % (nt_pos,ref_nt.upper(),alt_nt.upper())]
        # frameshift
        re_obj = re.search("frameshift",mut)
        if re_obj:
            return ["frameshift"]
        # premature_stop
        re_obj = re.search("premature_stop",mut)
        if re_obj:
            return ["premature_stop"]
        # Codon range
        ## any_missense_codon_425_452
        re_obj = re.search("any_missense_codon_([0-9]+)_([0-9]+)",mut)
        if re_obj:
            start = int(re_obj.group(1))
            end = int(re_obj.group(2))
            return ["any_missense_codon_%s" % i for i in range(start,end+1)]
        # Codon single
        ## any_missense_codon_425
        re_obj = re.search("any_missense_codon_([0-9]+)",mut)
        if re_obj:
            start = int(re_obj.group(1))
            return ["any_missense_codon_%s" % start]
        # Indel range
        re_obj = re.search("any_indel_nucleotide_([0-9]+)_([0-9]+)",mut)
        if re_obj:
            start = int(re_obj.group(1))
            end = int(re_obj.group(2))
            return ["any_indel_nucleotide_%s" % i for i in range(start,end+1)]
        # large_deletion
        re_obj = re.search("large_deletion",mut)
        if re_obj:
            return ["large_deletion"]
        sys.exit("%s is not a valid formatted mutation... Exiting!" % mut)
    
    def _parse(self, fasta_filepath, gff_filepath, download_dir):
        utils = Utils()
        tbdb_url = "https://raw.githubusercontent.com/jodyphelan/tbdb/master/tbdb.csv"
        tbdb_filepath = os.path.join(download_dir, "tbdb.csv")
        utils.download_and_save_file(tbdb_url, tbdb_filepath)
        tbdb_df = pd.read_csv(tbdb_filepath, header=0)
        return tbdb_df
