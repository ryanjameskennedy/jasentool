import os
import sys
import json
from jasentool.database import Database
from jasentool.utils import Utils

class Validate(object):
    def get_sample_id(self, results):
        return results["sample_id"]

    def get_species_name(self, results):
        return results["species_prediction"][0]["scientific_name"]

    def _check_exists(self, db_collection, sample_id):
        return (True if list(Database.find(db_collection, {"id": sample_id}, {})) else False)

    def search(self, search_query, search_kw, search_list):
        return [element for element in search_list if element[search_kw] == search_query]

    def get_virulence_results(self, results):
        return self.search("VIRULENCE", "type", results["element_type_result"])

    def get_pvl(self, results):
        virulence_results = self.get_virulence_results(results)
        return (True if self.search("lukS-PV", "gene_symbol", virulence_results[0]["result"]["genes"]) else False)

    def get_mlst(self, results):
        return self.search("mlst", "type", results["typing_result"])

    def get_cgmlst(self, results):
        return self.search("cgmlst", "type", results["typing_result"])

    def get_mdb_cgv_data(self, db_collection, sample_id):
        mdb_pvl = list(Database.get_pvl(db_collection, {"id": sample_id, "metadata.QC": "OK"}))
        mdb_mlst = list(Database.get_mlst(db_collection, {"id": sample_id, "metadata.QC": "OK"}))
        mdb_cgmlst = list(Database.get_cgmlst(db_collection, {"id": sample_id, "metadata.QC": "OK"}))
        try:
            mdb_pvl_present = int(mdb_pvl[0]["aribavir"]["lukS_PV"]["present"])
            mdb_mlst_seqtype = str(mdb_mlst[0]["mlst"]["sequence_type"]) if mdb_mlst[0]["mlst"]["sequence_type"] != "-" else str(None)
            mdb_mlst_alleles = mdb_mlst[0]["mlst"]["alleles"]
            mdb_cgmlst_alleles = mdb_cgmlst[0]["alleles"]
            return {"pvl": mdb_pvl_present, "mlst_seqtype": mdb_mlst_seqtype, "mlst_alleles": mdb_mlst_alleles, "cgmlst_alleles": mdb_cgmlst_alleles}
        except IndexError:
            return False

    def get_fin_data(self, sample_json):
        fin_pvl_present = self.get_pvl(sample_json)
        fin_mlst = self.get_mlst(sample_json)
        fin_cgmlst = self.get_cgmlst(sample_json)
        fin_mlst_seqtype = str(fin_mlst[0]["result"]["sequence_type"])
        fin_mlst_alleles = fin_mlst[0]["result"]["alleles"]
        fin_cgmlst_alleles = list(fin_cgmlst[0]["result"]["alleles"].values())
        return {"pvl": fin_pvl_present, "mlst_seqtype": fin_mlst_seqtype, "mlst_alleles": fin_mlst_alleles, "cgmlst_alleles": fin_cgmlst_alleles}

    def compare_mlst_alleles(self, old_mlst_alleles, new_mlst_alleles):
        match_count, total_count = 0, 0
        for allele in old_mlst_alleles:
            if str(old_mlst_alleles[allele]) == str(new_mlst_alleles[allele]):
                match_count += 1
            total_count += 1
        return 100*(match_count/total_count)

    def compare_cgmlst_alleles(self, old_cgmlst_alleles, new_cgmlst_alleles):
        match_count, total_count = 0, 0
        for allele in range(0, len(old_cgmlst_alleles)):
            if str(old_cgmlst_alleles[allele]) == str(new_cgmlst_alleles[allele]):
                match_count += 1
            total_count += 1
        return 100*(match_count/total_count)

    def compare_data(self, sample_id, old_data, new_data):
        pvl_comp = int(old_data["pvl"] == new_data["pvl"])
        mlst_seqtype_comp = int(old_data["mlst_seqtype"] == new_data["mlst_seqtype"])
        if mlst_seqtype_comp == 0:
            mlst_at_list = [f'{old_data["mlst_alleles"][gene]},{new_data["mlst_alleles"][gene]}' for gene in sorted(old_data["mlst_alleles"].keys())]
            mlst_at_str = ",".join(mlst_at_list)
            print(f'{sample_id},{old_data["mlst_seqtype"]},{new_data["mlst_seqtype"]},{mlst_at_str}')
        mlst_alleles = self.compare_mlst_alleles(old_data["mlst_alleles"], new_data["mlst_alleles"])
        cgmlst_alleles = self.compare_cgmlst_alleles(old_data["cgmlst_alleles"], new_data["cgmlst_alleles"])
        return f"{sample_id},{pvl_comp},{mlst_seqtype_comp},{mlst_alleles},{cgmlst_alleles}"

    def run(self, input_files, output_fpaths, db_collection, combined_output):
        utils = Utils()
        csv_output = "sample_id,pvl,mlst_seqtype,mlst_allele_matches(%),cgmlst_allele_matches(%)"
        for input_idx, input_file in enumerate(input_files):
            with open(input_file, 'r') as fin:
                sample_json = json.load(fin)
                sample_id = self.get_sample_id(sample_json)
                if not self._check_exists(db_collection, sample_id):
                    print(f"The sample provided ({sample_id}) does not exist in the provided database ({Database.db_name}) or collection ({db_collection}).")
                    continue
                mdb_data_dict = self.get_mdb_cgv_data(db_collection, sample_id)
                if mdb_data_dict:
                    species_name = self.get_species_name(sample_json)
                    fin_data_dict = self.get_fin_data(sample_json)
                    compared_data_output = self.compare_data(sample_id, mdb_data_dict, fin_data_dict)
                    csv_output += "\n" + compared_data_output
            if not combined_output:
                utils.write_out_txt(csv_output, f"{output_fpaths[input_idx]}.csv")
                csv_output = "pvl,mlst_seqtype,mlst_allele_matches(%),cgmlst_allele_matches(%)\n"

        if combined_output:
            utils.write_out_txt(csv_output, f"{output_fpaths[0]}.csv", )
