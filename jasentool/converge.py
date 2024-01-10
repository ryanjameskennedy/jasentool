import os
import pandas as pd
from jasentool.who import WHO
from jasentool.genome import Genome
from jasentool.tbprofiler import Tbprofiler

class Converge(object):
    def __init__(self, tbdb_dir, download_dir):
        self.tbdb_dir = tbdb_dir
        self.download_dir = download_dir
        self.intersection_outfpath = os.path.join(download_dir, "intersection.csv")
        self.unique_tbdb_outfpath = os.path.join(download_dir, "unique_tbdb.csv")
        self.unique_who_outfpath = os.path.join(download_dir, "unique_who.csv")

    def merge_dfs(self, df1, df2):
        pd.merge(df1,df2, on='name', how='inner')

    def compare_columns(self, tbdb_df, who_df, column_names):
        """Return a list of all of the unique and common variants in each dataframe"""

        # Merge DataFrames on both columns with indicator set to True
        merged = tbdb_df.merge(who_df, on=column_names, how='outer', indicator=True)

        # Filter rows that are unique to each DataFrame
        unique_tbdb_df = merged[merged['_merge'] == 'left_only']
        unique_who_df = merged[merged['_merge'] == 'right_only']

        # Filter rows that intersect
        intersecting_rows = merged[merged['_merge'] == 'both']

        drop_columns = ['Confers_y', 'Interaction_y', 'Literature_y', 'WHO Confidence_y', '_merge']
        column_mapping = {
            'Confers_x': 'Confers',
            'Interaction_x': 'Interaction',
            'Literature_x': 'Literature',
            'WHO Confidence_x': 'WHO Confidence',
        }

        # Drop the indicator column
        unique_tbdb_df = unique_tbdb_df.drop(columns=drop_columns).rename(columns=column_mapping)
        unique_who_df = unique_who_df.drop(columns=drop_columns).rename(columns=column_mapping)
        intersection_df = intersecting_rows.drop(columns=drop_columns).rename(columns=column_mapping)
        return intersection_df, unique_tbdb_df, unique_who_df

    def run(self):
        # Download the genome
        mycobacterium_genome = Genome("NC_000962.3", "AL123456.3", self.download_dir, "h37rv")
        fasta_filepath = mycobacterium_genome.download_fasta()
        gff_filepath = mycobacterium_genome.download_gff()
        who = WHO()
        tbprofiler = Tbprofiler(self.tbdb_dir)
        #h37rv_gb_filepath = mycobacterium_genome.download_genbank()
        who_df = who._parse(fasta_filepath, gff_filepath, self.download_dir)
        tbdb_df = tbprofiler._parse(fasta_filepath, gff_filepath, self.download_dir)
        #tbdb_df, who_df = pd.read_csv("/data/bnf/dev/ryan/pipelines/jasen/converge/tbdb.csv"), pd.read_csv("/data/bnf/dev/ryan/pipelines/jasen/converge/who.csv")
        intersection_df, unique_tbdb_df, unique_who_df = self.compare_columns(tbdb_df, who_df, ['Drug', 'Gene', 'Mutation'])
        intersection_df.to_csv(self.intersection_outfpath, index=False)
        unique_tbdb_df.to_csv(self.unique_tbdb_outfpath, index=False)
        unique_who_df.to_csv(self.unique_who_outfpath, index=False)
