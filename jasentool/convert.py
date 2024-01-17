"""Module that converts file type"""

class Convert:
    """Convert class for converting files into desired format"""
    @staticmethod
    def targets2bed(target_file, accn):
        """Convert cgmlst locus targets to bed file format"""
        bed_output = ""
        with open(target_file, 'r', encoding="utf-8") as fin:
            for line in fin:
                if line.startswith("Locus"):
                    continue
                line_split = line.split("\t")
                start = int(line_split[3]) - 1
                length = int(line_split[4])
                end = start + length
                bed_output += f"{accn}\t{start}\t{end}\n"
        return bed_output
