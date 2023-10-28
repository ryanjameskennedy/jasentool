#!/usr/bin/env python3

import csv
import pandas as pd

class Utils(object):
    @staticmethod
    def write_out_csv(csv_dict, assay, platform, out_fpath):
        with open(out_fpath, 'w+') as csvfile:
            fieldnames = ["id", "group", "species", "assay", "platform", "read1", "read2"] #header
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for sample in csv_dict:
                row_dict = {"id":sample, "group": csv_dict[sample][0], "species": csv_dict[sample][1], "assay": assay, "platform": platform, "read1": csv_dict[sample][2][0], "read2": csv_dict[sample][2][1]} #write rows to CSV
                writer.writerow(row_dict)

    @staticmethod
    def write_out_txt(output_txt, out_fpath):
        with open(out_fpath, 'w+') as fout:
            fout.write(output_txt)
