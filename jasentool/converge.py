"""Module to converge mutation catalogues"""

import os
import shutil
import pandas as pd
from jasentool.who import WHO
from jasentool.genome import Genome
from jasentool.tbprofiler import Tbprofiler
from jasentool.utils import Utils

class Converge:
    """Class that converges mutation catalogues"""
    def __init__(self, download_dir):
        self.download_dir = download_dir
        self.fohm_fpath = os.path.join(os.path.dirname(__file__), "data/dbs/fohm.csv")
        self.intersection_outfpath = os.path.join(download_dir, "intersection.csv")
        self.unique_tbdb_outfpath = os.path.join(download_dir, "unique_tbdb.csv")
        self.unique_who_outfpath = os.path.join(download_dir, "unique_who.csv")
        self.fohm_tbdb_outfpath = os.path.join(download_dir, "fohm_tbdb.csv")
        self.convereged_outfpath = os.path.join(download_dir, "converged_who_fohm_tbdb.csv")
        self.tbdb_filepath = os.path.join(download_dir, "tbdb.csv")
        self.tbdb_url = "https://raw.githubusercontent.com/jodyphelan/tbdb/master/tbdb.csv"

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

    def rm_intermediary_files(self):
        files = os.listdir(self.download_dir)
        files.remove('converged_who_fohm_tbdb.csv')
        files.remove('unique_tbdb.csv')
        files.remove('unique_who.csv')
        files.remove('fohm.csv')
        for filename in files:
            filepath = os.path.join(self.download_dir, filename)
            if os.path.isfile(filepath):
                os.remove(filepath)
            elif os.path.isdir(filepath):
                shutil.rmtree(filepath)

    def run(self):
        """Run the retrieval and convergance of mutation catalogues"""
        utils = Utils()
        # Download the genome
        mycobacterium_genome = Genome("NC_000962.3", "AL123456.3", self.download_dir, "h37rv")
        fasta_filepath = mycobacterium_genome.download_fasta()
        gff_filepath = mycobacterium_genome.download_gff()
        utils.download_and_save_file(self.tbdb_url, self.tbdb_filepath)
        who = WHO()
        tbprofiler = Tbprofiler(self.tbdb_filepath)
        #h37rv_gb_filepath = mycobacterium_genome.download_genbank()
        who_df = who._parse(fasta_filepath, gff_filepath, self.download_dir)
        tbdb_df = tbprofiler._parse(self.download_dir)
        #tbdb_df = pd.read_csv("/data/bnf/dev/ryan/pipelines/jasen/converge/tbdb.csv")
        #who_df = pd.read_csv("/data/bnf/dev/ryan/pipelines/jasen/converge/who.csv")
        fohm_df = pd.read_csv(self.fohm_fpath)
        column_names = ['Drug', 'Gene', 'Mutation']
        intersection_df, unique_tbdb_df, unique_who_df = self.compare_columns(tbdb_df, who_df, column_names)
        dfs_to_concat = [intersection_df, unique_tbdb_df, fohm_df]
        fohm_tbdb_df = pd.concat(dfs_to_concat, ignore_index=True).drop_duplicates()
        intersection_df.to_csv(self.intersection_outfpath, index=False)
        unique_tbdb_df.to_csv(self.unique_tbdb_outfpath, index=False)
        unique_who_df.to_csv(self.unique_who_outfpath, index=False)
        fohm_tbdb_df.to_csv(self.fohm_tbdb_outfpath, index=False)
        dfs_to_converge = [intersection_df, unique_tbdb_df, unique_who_df, fohm_df]
        converged_df = pd.concat(dfs_to_converge, ignore_index=True).drop_duplicates()
        converged_df.to_csv(self.convereged_outfpath, index=False)
        self.rm_intermediary_files()
