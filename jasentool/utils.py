#!/usr/bin/env python3

import os
import csv
import shutil
import pathlib
import requests
import subprocess
import pandas as pd
from time import sleep
from zipfile import ZipFile 

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
    def pipeline_ready(batch_file):
        assays = ['saureus']
        for assay in assays:
            if assay in batch_file:
                return True
        return False

    @staticmethod
    def copy_batch_and_csv_files(batch_files, csv_files, remote_dir, remote_hostname, remote=False):
        if remote:
            # Copy files to remote server using ssh/scp
            process = subprocess.run(
                f'ssh {remote_hostname} mkdir -p {remote_dir}',
                shell=True
            )
            process = subprocess.run(
                f'scp {" ".join(batch_files)} {" ".join(csv_files)} {remote_hostname}:{remote_dir}',
                shell=True,
                stdout=subprocess.PIPE,
                universal_newlines=True
            )
        else:
            # Copy files to a local directory
            pathlib.Path(remote_dir).mkdir(parents=True, exist_ok=True)
            for fn in batch_files + csv_files:
                shutil.copy(fn, remote_dir)

    @staticmethod
    def start_remote_pipelines(batch_files, remote_hostname, remote_dir):
        for batch_file in batch_files:
            if Utils.pipeline_ready(batch_file):
                sleep(10.0) # Avoid maxing SSH auth connections
                process = subprocess.Popen(
                    ["ssh", remote_hostname, "bash", f"{remote_dir}/{os.path.basename(batch_file)}"],
                    close_fds=True
                )

    @staticmethod
    def download_and_save_file(url, output_filepath):
        try:
            # Make a request to the URL
            response = requests.get(url, stream=True)
            response.raise_for_status()  # Raise an error for bad responses

            # Open the output file in binary write mode
            with open(output_filepath, 'wb') as output_file:
                # Iterate over the content in chunks and write to the file
                for chunk in response.iter_content(chunk_size=8192):
                    output_file.write(chunk)

            print(f"File downloaded and saved to: {output_filepath}")

        except requests.exceptions.RequestException as e:
            print(f"Error downloading the file: {e}")

    @staticmethod
    def unzip(zip_file, outdir):
        with ZipFile(zip_file, 'r') as zip_object:
            zip_object.extractall(path=outdir)
