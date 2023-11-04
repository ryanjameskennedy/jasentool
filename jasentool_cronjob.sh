#!/bin/bash
seqrunid=$(head -2 /data/tmp/multi_microbiology.csv | tail -1 | cut -d',' -f7 | cut -d'/' -f5)

conda-exec -n jasen jasentool fix --csv_file /data/tmp/multi_microbiology.csv --sh_file /data/tmp/multi_microbiology.sh -o ${seqrunid}_jasen.csv --remote_dir /fs1/ryan/pipelines/jasen/bjorn/ --remote --auto-start
