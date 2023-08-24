# Mongtool v0.0.1: A mongodb validation tool for comaparing pipeline outputs
## Dependencies (latest)
* python=3.11
* pymongo

## Using Mongtool
### Use the help argument for information regarding the Mongtool's methods
```
mongtool -h
```

### Use the method help argument for information regarding the input for each of Mongtool's methods (`find`, `insert`, `remove` and `validate`)
```
mongtool <method> -h
```

### Validate pipeline data
```
mongtool validate (-i INPUT_FILE [INPUT_FILE ...] | --input_dir INPUT_DIR) --db_name DB_NAME --db_collection DB_COLLECTION -o OUTPUT_FILE [--address ADDRESS] [-h]
```
