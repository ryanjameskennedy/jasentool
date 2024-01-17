"""Module for retrieving qc results"""

import os
import json
import subprocess

class QC:
    """Class for retrieving qc results"""
    def __init__(self, args):
        self.results = {}
        self.bam = args.bam
        self.bed = args.bed
        self.sample_id = args.sample_id
        self.cpus = args.cpus
        self.baits = args.baits
        self.reference = args.reference
        self.paired = self.is_paired()

    def write_json_result(self, json_result, output_filepath):
        """Write out json file"""
        with open(output_filepath, 'w', encoding="utf-8") as json_file:
            json_file.write(json_result)

    def parse_basecov_bed(self, basecov_fpath, thresholds):
        """Parse base coverage bed file"""
        with open(basecov_fpath, "r", encoding="utf-8") as cov_fh:
            head_str = cov_fh.readline().strip().lstrip("#")
            head = head_str.split("\t")
            cov_field = head.index("COV")

            tot_bases = 0
            above_cnt = {min_val: 0 for min_val in thresholds}

            tot, cnt = 0, 0
            levels = {}
            for line in cov_fh:
                line = line.strip().split("\t")
                tot += int(line[2])
                cnt += 1
                tot_bases += 1
                for min_val in thresholds:
                    if int(line[cov_field]) >= min_val:
                        above_cnt[min_val] += 1

            above_pct = {min_val: 100 * (above_cnt[min_val] / tot_bases) for min_val in thresholds}

            mean_cov = tot / cnt

            # Calculate the inter-quartile range / median (IQR/median)
            q1_num = cnt / 4
            q3_num = 3 * cnt / 4
            median_num = cnt / 2
            sum_val = 0
            quartile1, quartile3, median = None, None, None
            iqr_median = "9999"
            for level in sorted(levels):
                sum_val += levels[level]
                if sum_val >= q1_num and not quartile1:
                    quartile1 = level
                if sum_val >= median_num and not median:
                    median = level
                if sum_val >= q3_num and not quartile3:
                    quartile3 = level

            if quartile1 and quartile3 and median:
                iqr_median = (quartile3 - quartile1) / median

            return above_pct, mean_cov, iqr_median

    def is_paired(self):
        """Check if reads are paired"""
        line = subprocess.check_output(f"samtools view {self.bam} | head -n 1| awk '{{print $2}}'", shell=True, text=True)
        remainder = int(line) % 2
        is_paired = 1 if remainder else 0
        return is_paired

    def system_p(self, *cmd):
        """Execute subproces"""
        print(f"RUNNING: {' '.join(cmd)}")
        print()
        subprocess.run(cmd, check=True)

    def run(self):
        """Run QC info extraction"""
        if self.baits and self.reference:
            print("Calculating HS-metrics...")
            dict_file = self.reference
            if not dict_file.endswith(".dict"):
                dict_file += ".dict"
            if not os.path.isfile(f"{self.bed}.interval_list"):
                self.system_p(f"picard BedToIntervalList -I {self.bed} -O {self.bed}.interval_list -SD {dict_file}")
            if not os.path.isfile(f"{self.baits}.interval_list"):
                self.system_p(f"picard BedToIntervalList -I {self.baits} -O {self.baits}.interval_list -SD {dict_file}")
            self.system_p(f"picard CollectHsMetrics -I {self.bam} -O {self.bam}.hsmetrics -R {self.reference} -BAIT_INTERVALS {self.baits}.interval_list -TARGET_INTERVALS {self.bed}.interval_list")

            with open(f"{self.bam}.hsmetrics", "r", encoding="utf-8") as fin:
                for line in fin:
                    if line.startswith("## METRICS CLASS"):
                        next(fin)
                        vals = next(fin).split("\t")
                        self.results['pct_on_target'] = vals[18]
                        self.results['fold_enrichment'] = vals[26]
                        self.results['median_coverage'] = vals[23]
                        self.results['fold_80'] = vals[33]

        print("Collecting basic stats...")
        flagstat = subprocess.check_output(f"sambamba flagstat {'-t '+str(self.cpus) if self.cpus else ''} {self.bam}", shell=True, text=True).splitlines()
        num_reads = int(flagstat[0].split()[0])
        dup_reads = int(flagstat[3].split()[0])
        mapped_reads = int(flagstat[4].split()[0])

        if self.paired:
            print("Collect insert sizes...")
            self.system_p(f"picard CollectInsertSizeMetrics -I {self.bam} -O {self.bam}.inssize -H {self.bam}.ins.pdf -STOP_AFTER 1000000")
            with open(f"{self.bam}.inssize", "r", encoding="utf-8") as ins:
                for line in ins:
                    if line.startswith("## METRICS CLASS"):
                        next(ins)
                        vals = next(ins).split("\t")
                        self.results['ins_size'] = vals[0]
                        self.results['ins_size_dev'] = vals[1]

            os.remove(f"{self.bam}.inssize")
            os.remove(f"{self.bam}.ins.pdf")

        out_prefix = f"{self.bam}_postalnQC"
        thresholds = [1, 10, 30, 100, 250, 500, 1000]

        print("Collecting depth stats...")
        self.system_p(f"sambamba depth base -c 0 {'-t '+str(self.cpus) if self.cpus else ''} -L {self.bed} {self.bam} > {out_prefix}.basecov.bed")
        pct_above, mean_cov, iqr_median = self.parse_basecov_bed(f"{out_prefix}.basecov.bed", thresholds)
        os.remove(f"{out_prefix}.basecov.bed")

        self.results['pct_above_x'] = pct_above
        self.results['tot_reads'] = num_reads
        self.results['mapped_reads'] = mapped_reads
        self.results['dup_reads'] = dup_reads
        self.results['dup_pct'] = dup_reads / mapped_reads
        self.results['sample_id'] = self.sample_id
        self.results['mean_cov'] = mean_cov
        self.results['iqr_median'] = iqr_median

        json_result = json.dumps(self.results, indent=4)
        return json_result
