def write_out(output_fpath, output):
    with open(f"{output_fpath}.csv", 'w+') as fout:
        fout.write(output)