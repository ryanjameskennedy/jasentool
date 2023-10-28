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

### Get cgviz meta csv
```
mongoexport --quiet --db=cgviz --collection=sample --type=csv --fields=id,mlst.sequence_type,aribavir.lukF_PV.present,aribavir.lukS_PV.present,missing,metadata.QC,metadata.Comment,run --query='{"metadata.QC":"OK"}' | grep -v FOHM | sed "1s/id,mlst.sequence_type,aribavir.lukF_PV.present,aribavir.lukS_PV.present,missing,metadata.QC,metadata.Comment,run/id,mlst,lukF_PV,lukS_PV,missing,QC,Comment,run/" > cgviz_meta.csvmongoexport --quiet --db=cgviz --collection=sample --type=csv --fields=id,mlst.sequence_type,aribavir.lukF_PV.present,aribavir.lukS_PV.present,missing,metadata.QC,metadata.Comment,run --query='{"metadata.QC":"OK"}' | grep -v FOHM | sed "1s/id,mlst.sequence_type,aribavir.lukF_PV.present,aribavir.lukS_PV.present,missing,metadata.QC,metadata.Comment,run/id,mlst,lukF_PV,lukS_PV,missing,QC,Comment,run/" > cgviz_meta.csv
```

### Find missing samples
```
jasentool missing -i <cgviz_meta.csv> --analysis_dir <jasen_analysis_results_dir> --restore_dir <restore_dir> --restore_file <restore_file.sh> -o <output_file.csv>
```