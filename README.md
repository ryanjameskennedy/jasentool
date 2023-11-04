# Jasentool v0.0.1: A mongodb validation tool for comaparing pipeline outputs
## Dependencies (latest)
* python=3.11
* pymongo

## Using Jasentool
### Use the help argument for information regarding the Jasentool's methods
```
jasentool -h
```

### Use the method help argument for information regarding the input for each of Jasentool's methods (`find`, `insert`, `remove` and `validate`)
```
jasentool <method> -h
```

### Validate pipeline data
```
jasentool validate (-i INPUT_FILE [INPUT_FILE ...] | --input_dir INPUT_DIR) --db_name DB_NAME --db_collection DB_COLLECTION -o OUTPUT_FILE [--address ADDRESS] [-h]
```

### Get samplesheet inputs
```
grep saureus /fs2/seqdata/*/*/SampleSheet.csv > saureus_fs2_sample_sheet.csv
grep saureus /data/*/*/SampleSheet.csv > saureus_data_sample_sheet.csv
grep saureus /media/isilon/backup_hopper/seqdata/*/*/SampleSheet.csv > saureus_isilon_sample_sheet.csv
```

### Find missing samples
```
jasentool missing --db_name <db_name> --db_collection <db_collection> --analysis_dir <jasen_analysis_results_dir> --restore_dir <restore_dir> --restore_file <restore_file.sh> -o <output_file.csv>
```

### Fix bjorn csv
```
jasentool fix --csv_file /data/tmp/multi_microbiology.csv --sh_file /data/tmp/multi_microbiology.sh -o <flow_cell_id>_jasen.csv --remote_dir /fs1/ryan/pipelines/jasen/bjorn/ --remote
```
