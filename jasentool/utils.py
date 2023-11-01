#!/usr/bin/env python3

import os
import csv
import shutil
import pathlib
import subprocess
import pandas as pd
from time import sleep

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

    @staticmethod
    def copy_batch_and_csv_files(batch_files, csv_files, output_dir, remote_hostname, remote=False):
        if remote:
            # Copy files to remote server using ssh/scp
            process = subprocess.run(
                f'ssh {remote_hostname} mkdir -p {output_dir}',
                shell=True
            )
            process = subprocess.run(
                f'scp {" ".join(batch_files)} {" ".join(csv_files)} {remote_hostname}:{output_dir}',
                shell=True,
                stdout=subprocess.PIPE,
                universal_newlines=True
            )
        else:
            # Copy files to a local directory
            pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
            for fn in batch_files + csv_files:
                shutil.copy(fn, output_dir)

    def start_remote_pipelines(batch_files, remote_hostname, remote_dir):
        for batch_file in batch_files:
            sleep(10.0) # Avoid maxing SSH auth connections
            process = subprocess.Popen(
                ["ssh", remote_hostname, "bash", f"{remote_dir}/{os.path.basename(batch_file)}"],
                close_fds=True
            )
