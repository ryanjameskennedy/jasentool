import os
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
        search_files = os.listdir(parent_dir)
        return sorted([os.path.join(parent_dir, search_file) for search_file in search_files if search_file.startswith(search_term)])

    @staticmethod
    def edit_read_paths(reads, restore_dir):
        filename = os.path.join(restore_dir, reads.split("BaseCalls/")[1])
        read1 = filename.rstrip(".spring") + "_R1_001.fastq.gz"
        read2 = filename.rstrip(".spring") + "_R2_001.fastq.gz"
        return os.path.join(restore_dir, reads.split("BaseCalls/")[1]), [read1, read2]

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
                        paired_reads = Missing.find_files(clarity_id, parent_dir)
                        if len(paired_reads) == 2 and paired_reads[0].endswith(".gz"):
                            csv_dict[sample_id] = [clarity_group_id, species, paired_reads]
                        elif len(paired_reads) == 1 and paired_reads[0].endswith(".spring"):
                            spring_fpaths = paired_reads
                            (restored_spring_fpaths, paired_reads) = list(map(Missing.edit_read_paths, spring_fpaths, [restore_dir]*len(spring_fpaths)))[0]
                            csv_dict[sample_id] = [clarity_group_id, species, paired_reads, spring_fpaths, restored_spring_fpaths]
                        elif len(paired_reads) == 4 and paired_reads[0].endswith(".gz"):
                            paired_reads = Missing.rm_double_dmltplx(paired_reads)
                            if len(paired_reads) == 2:
                                csv_dict[sample_id] = [clarity_group_id, species, paired_reads]
                            elif len(paired_reads) == 4:
                                paired_reads_string = '\n-'.join(paired_reads)
                                print(sample_id)
                                print(f"There are 4 sets of reads related to sample {sample_id} from the {parent_dir}: \n-{paired_reads_string}\n")
                        #else:
                            #print(f"The sample {sample_id} doesn't have read/spring files in the {parent_dir} ({paired_reads}).")
                    except FileNotFoundError:
                        print(f"WARNING: {parent_dir} does not exist regarding {sample_id}.")
                        print(sample_sheet)

        return csv_dict

    @staticmethod
    def check_format(fpath):
        if not fpath.startswith("/data"):
            isilon_fpath = "/media/isilon/backup_hopper" + fpath
            fs2_fpath = "/fs2" + fpath
            if os.path.exists(os.path.join(isilon_fpath, "Data/Intensities/BaseCalls")):
                return isilon_fpath
            elif os.path.exists(os.path.join(fs2_fpath, "Data/Intensities/BaseCalls")):
                return fs2_fpath
        return fpath

    @staticmethod
    def parse_mongodb_csv(input_fpath):
        with open(input_fpath, "rb") as csvfile:
            meta = pd.read_csv(csvfile)
            meta = meta.drop(columns=["mlst", "lukF_PV", "lukS_PV", "missing", "QC", "Comment"])
            meta["run"] = meta["run"].apply(lambda x: Missing.check_format(x))
            meta_dict = meta.to_dict(orient="records")
        return meta_dict
    
    @staticmethod
    def parse_dir(dir_fpath):
        return [filename.split("_")[0] for filename in os.listdir(dir_fpath)]
    
    @staticmethod
    def find_missing(meta_dict, analysis_dir_fnames, restore_dir):
        sample_run = ""
        missing_samples = []
        csv_dict = {}
        for sample in meta_dict:
            if sample["id"] not in analysis_dir_fnames:
                if sample_run != sample["run"]: #if sample run changes based on 
                    ss_dict = {}
                    sample_sheets = Missing.find_files("SampleSheet", sample["run"])
                    if sample_sheets:
                        for sample_sheet in sample_sheets:
                            ss_dict |= Missing.parse_sample_sheet(sample_sheet, restore_dir)
                        if not sample_sheet:
                            print(f"sample sheets yieded nothing from {sample['run']}")
                        csv_dict |= ss_dict
                    else:
                        print(f"Note: No sample sheets exist in the following path path {sample['run']}!")
                    sample_run = sample["run"]

                missing_samples.append(sample["id"])
        print(f"{len(csv_dict.keys())} samples found")
        print(f"{len(missing_samples)} samples missing")
        return csv_dict, "\n".join(missing_samples)

    @staticmethod
    def create_bash_script(csv_dict, restore_dir):
        spring_command = ""
        shell_script_path = 'SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"\n'
        shell_fail_count = "FAIL=0\n"
        shell_for_loop = 'for job in $PIDS; do wait $job || let "FAIL+=1"; done \nif [ "$FAIL" != "0" ]; \nthen \n\techo Failed to restore from backup \n\texit 2 \nfi \n'
        for sample in csv_dict:
            try:
                spring_command = spring_command + f'/fs2/sw/bnf-scripts/jcp {csv_dict[sample][3][0]} {restore_dir}/ && /fs2/sw/bnf-scripts/unspring_file.pl {csv_dict[sample][4]} {restore_dir}/ WAIT &\nPIDS="$PIDS $!"\n'
            except IndexError:
                continue
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
        for empty_file in list(empty_files_dict.keys()):
            csv_dict.pop(empty_file, None)
        return empty_files_dict, csv_dict