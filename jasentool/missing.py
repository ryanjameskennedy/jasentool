import os
import re
import csv
import pandas as pd

class Missing(object):
    @staticmethod
    def rm_double_dmltplx(read_files):
        first_reads = read_files[0]
        for read_file in read_files[1:]:
            errors = 0
            for i in range(len(first_reads)):
                if first_reads[i] != read_file[i]:
                    errors += 1
            if errors == 1:
                return [first_reads, read_file]
        return read_files

    @staticmethod
    def find_files(search_term, parent_dir):
        try:
            search_files = os.listdir(parent_dir)
        except FileNotFoundError:
            print(f"WARN: {parent_dir} does not exist! Trying to fix.")
        finally:
            search_files = os.listdir(parent_dir)
            found_files = sorted([os.path.join(parent_dir, search_file) for search_file in search_files if re.search(search_term, search_file) and not search_file.endswith("~")])
            return found_files

    @staticmethod
    def edit_read_paths(reads, restore_dir):
        restore_dirs = set([restore_dir.rstrip("/"), "/fs2/seqdata/restored"])
        filename = os.path.join(restore_dir, reads.split("BaseCalls/")[1])
        read1, read2 = [filename.rstrip(".spring") + f"_R{i}_001.fastq.gz" for i in [1, 2]]
        return os.path.join(restore_dir, reads.split("BaseCalls/")[1]), [read1, read2]
    
    @staticmethod
    def check_file_cp(reads, restore_dir):
        checked_reads = []
        restore_dirs = set([restore_dir.rstrip("/"), "/fs2/seqdata/restored"])
        for filepath in reads:
            filename = os.path.basename(filepath)
            if filepath.startswith("/fs") and os.path.exists(filepath):
                checked_reads.append(filepath)
            else:
                for dir in restore_dirs:
                    read_fpath = os.path.join(dir, filename)
                    if os.path.exists(read_fpath) and not os.path.isdir(read_fpath) and len(checked_reads) != 2:
                        checked_reads.append(read_fpath)
        if len(checked_reads) == 0:
            checked_reads = [os.path.join(restore_dir, os.path.basename(read_filepath)) for read_filepath in reads]
        return checked_reads

    @staticmethod
    def parse_sample_sheet(sample_sheet, restore_dir):
        csv_dict = {}
        with open(sample_sheet, "r") as fin:
            for line in fin:
                if line.endswith("saureus\n"):
                    line = line.rstrip()
                    sample_id = line.split(",")[-1].split("_")[1]
                    species = line.split(",")[-1].split("_")[2]
                    try:
                        clarity_id = line.split(",")[0].split(":")[1]
                    except IndexError:
                        clarity_id = line.split(",")[0]
                    try:
                        clarity_group_id = clarity_id.split("_")[1]
                    except IndexError:
                        clarity_group_id = clarity_id
                    if ":" in line:
                        parent_dir = os.path.join(line.split(":")[0].rstrip("SampleSheet.csv"), "Data/Intensities/BaseCalls/")
                    else:
                        parent_dir = os.path.join(os.path.dirname(sample_sheet), "Data/Intensities/BaseCalls/")
                    try:
                        paired_reads = Missing.find_files(r'^' + clarity_id, parent_dir)
                        if len(paired_reads) == 2 and paired_reads[0].endswith(".gz"):
                            restored_reads_fpaths = Missing.check_file_cp(paired_reads, restore_dir)
                            csv_dict[sample_id] = [clarity_group_id, species, restored_reads_fpaths, None, paired_reads]
                        elif len(paired_reads) == 1 and paired_reads[0].endswith(".spring"):
                            spring_fpaths = paired_reads
                            (restored_spring_fpaths, paired_reads) = list(map(Missing.edit_read_paths, spring_fpaths, [restore_dir]*len(spring_fpaths)))[0]
                            csv_dict[sample_id] = [clarity_group_id, species, paired_reads, spring_fpaths, restored_spring_fpaths]
                        elif len(paired_reads) == 4 and paired_reads[0].endswith(".gz"):
                            paired_reads = Missing.rm_double_dmltplx(paired_reads)
                            if len(paired_reads) == 2:
                                restored_reads_fpaths = Missing.check_file_cp(paired_reads, restore_dir)
                                csv_dict[sample_id] = [clarity_group_id, species, restored_reads_fpaths, None, paired_reads]
                            elif len(paired_reads) == 4:
                                paired_reads_string = '\n-'.join(paired_reads)
                                print(f"There are 4 sets of reads related to sample {sample_id} from the {parent_dir}: \n-{paired_reads_string}\n")
                        elif len(paired_reads) == 3:
                            paired_reads = [paired_read for paired_read in paired_reads if paired_read.endswith(".fastq.gz")]
                            restored_reads_fpaths = Missing.check_file_cp(paired_reads, restore_dir)
                            csv_dict[sample_id] = [clarity_group_id, species, restored_reads_fpaths, None, paired_reads]
                        elif len(paired_reads) == 6:
                            paired_reads = [paired_read for paired_read in paired_reads if paired_read.endswith(".fastq.gz")]
                            restored_reads_fpaths = Missing.check_file_cp(paired_reads, restore_dir)
                            csv_dict[sample_id] = [clarity_group_id, species, restored_reads_fpaths, None, paired_reads]
                        #elif len(paired_reads) == 0:
                            #print(f"The sample {sample_id} doesn't have read/spring files in the {parent_dir} ({paired_reads}).")
                        #else:
                            #print(len(paired_reads))
                    except FileNotFoundError:
                        print(f"WARNING: {parent_dir} does not exist regarding {sample_id}.")
                        print(sample_sheet)

        return csv_dict

    @staticmethod
    def check_format(fpath):
        if fpath.startswith("/fs1") and not os.path.exists(os.path.join(fpath, "Data/Intensities/BaseCalls")):
            print(f"WARN: {fpath} does not exist! Fixing by removing '/fs1' prefix.")
            fpath = fpath.replace("/fs1", "")
        if fpath.startswith("NovaSeq"):
            fpath = "/seqdata/" + fpath
            print(f"WARN: {fpath} does not exist! Fixing by adding '/seqdata/' as a prefix.")
        if not fpath.startswith("/data"):
            fs2_fpath = "/fs2" + fpath
            isilon_fpath = "/media/isilon/backup_hopper" + fpath
            data_fpath = "/data" + fpath
            if os.path.exists(os.path.join(fs2_fpath, "Data/Intensities/BaseCalls")):
                return fs2_fpath
            elif os.path.exists(os.path.join(isilon_fpath, "Data/Intensities/BaseCalls")):
                return isilon_fpath
            elif os.path.exists(os.path.join(data_fpath, "Data/Intensities/BaseCalls")):
                return data_fpath
            elif os.path.exists(fpath):
                return fpath.rstrip("Data/Intensities/BaseCalls/")#.replace("Data/Intensities/BaseCalls", "")
            else:
                print(f"WARN: Base calls for {fpath} cannot be found.")
        return fpath

    @staticmethod
    def parse_dir(dir_fpath):
        return [filename.split("_")[0] for filename in os.listdir(dir_fpath)]
    
    @staticmethod
    def filter_csv_dict(csv_dict, missing_samples):
        filtered_csv_dict = {}
        not_found = []
        for missing_sample in missing_samples:
            try:
                filtered_csv_dict[missing_sample] = csv_dict[missing_sample]
            except KeyError:
                not_found.append(missing_sample)
                #print(f"{missing_sample} could not be found")
        print(f"{len(not_found)} samples could not be found")
        print(f"{len(filtered_csv_dict.keys())} samples remain after filtering")
        return filtered_csv_dict, not_found
    
    @staticmethod
    def find_missing(meta_dict, analysis_dir_fnames, restore_dir):
        sample_run = ""
        missing_samples = []
        csv_dict = {}
        #print(f"{len(list(meta_dict))} samples found in the meta dictionary")
        for sample in meta_dict:
            if sample["id"] not in analysis_dir_fnames:
                missing_samples.append(sample["id"])
                if sample_run != sample["run"]: #if sample run changes based on 
                    ss_dict = {}
                    sample_run_dir = Missing.check_format(sample["run"])
                    sample_sheets = Missing.find_files(r'.csv$', sample_run_dir)
                    if sample_sheets:
                        for sample_sheet in sample_sheets:
                            ss_dict |= Missing.parse_sample_sheet(sample_sheet, restore_dir)
                        if not sample_sheet:
                            print(f"sample sheets yieded nothing from {sample['run']}")
                        csv_dict |= ss_dict
                    else:
                        print(f"WARN: No sample sheets exist in the following path path {sample['run']}!")
                    sample_run = sample["run"]

        print(f"{len(csv_dict.keys())} samples found")
        print(f"{len(missing_samples)} samples missing")
        print(f"{len(missing_samples)-len(set(missing_samples))} duplicate sample ids")
        filtered_csv_dict, not_found = Missing.filter_csv_dict(csv_dict, missing_samples)
        return filtered_csv_dict, "\n".join(not_found)

    @staticmethod
    def create_bash_script(csv_dict, restore_dir):
        spring_command = ""
        shell_script_path = 'SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"\n'
        shell_fail_count = "FAIL=0\n"
        shell_for_loop = 'for job in $PIDS; do wait $job || let "FAIL+=1"; done \nif [ "$FAIL" != "0" ]; \nthen \n\techo Failed to restore from backup \n\texit 2 \nfi \n'
        for sample in csv_dict:
            jcp_command = ""
            unspring_command = ""
            try:
                spring_fpaths, restored_fpaths = csv_dict[sample][3][0], csv_dict[sample][4]
                read1, read2 = csv_dict[sample][2]
                if not os.path.exists(restored_fpaths) and not os.path.exists(read1):
                    jcp_command = f'/fs2/sw/bnf-scripts/jcp {spring_fpaths} {restore_dir}/ && '
                    unspring_command = f'/fs2/sw/bnf-scripts/unspring_file.pl {restored_fpaths} {restore_dir}/ WAIT &\nPIDS="$PIDS $!"\n'
                spring_command = spring_command + jcp_command + unspring_command
            except TypeError:
                for read_fpath in csv_dict[sample][4]:
                    jcp_command = f'/fs2/sw/bnf-scripts/jcp {read_fpath} {restore_dir}/ WAIT &\nPIDS="$PIDS $!"\n'
                    spring_command = spring_command + jcp_command
        bash_script = shell_script_path + shell_fail_count + spring_command + shell_for_loop
        return bash_script

    @staticmethod
    def remove_empty_files(csv_dict):
        empty_files_dict = {}
        for sample in csv_dict:
            try:
                file_size_r1 = os.path.getsize(csv_dict[sample][2][0]) / (1024 * 1024)
                file_size_r2 = os.path.getsize(csv_dict[sample][2][1]) / (1024 * 1024)
                if file_size_r1 < 10 or file_size_r2 < 10:
                    empty_files_dict[sample] = csv_dict[sample]
            except FileNotFoundError:
                print(f"WARN: {sample} read files ({csv_dict[sample][2][0]} and/or {csv_dict[sample][2][1]}) could not be found!")
            except IndexError:
                print(csv_dict[sample])
        for empty_file in list(empty_files_dict.keys()):
            csv_dict.pop(empty_file, None)
        return empty_files_dict, csv_dict
