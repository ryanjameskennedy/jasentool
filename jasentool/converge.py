import os
import pandas as pd
from jasentool.who import WHO
from jasentool.genome import Genome
from jasentool.tbprofiler import Tbprofiler

class Converge(object):
    @staticmethod
    def merge_dfs(df1, df2):
        pd.merge(df1,df2, on='name', how='inner')

    @staticmethod
    def run(tbdb_dir, download_dir):
        # Download the genome
        mycobacterium_genome = Genome("NC_000962.3", "AL123456.3", download_dir, "h37rv")
        fasta_filepath = mycobacterium_genome.download_fasta()
        gff_filepath = mycobacterium_genome.download_gff()
        who = WHO()
        tbprofiler = Tbprofiler(tbdb_dir)
        #h37rv_gb_filepath = mycobacterium_genome.download_genbank()
        who_df = who._parse(fasta_filepath, gff_filepath, download_dir)
        tbdb_df = tbprofiler._parse(fasta_filepath, gff_filepath, download_dir)

