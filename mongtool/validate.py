import json
from mongtool.database import Database

class Validate(object):
    def get_sample_id(self, results):
        return results["sample_id"]

    def get_species_name(self, results):
        return results["species_prediction"][0]["scientific_name"]

    def search(self, search_query, search_kw, search_list):
        return [element for element in search_list if element[search_kw] == search_query]

    def get_virulence_results(self, results):
        return self.search("VIRULENCE", "type", results["element_type_result"])

    def get_pvl(self, results):
        virulence_results = self.get_virulence_results(results)
        return (True if self.search("lukS-PV", "gene_symbol", virulence_results[0]["result"]["genes"]) else False)

    def get_mlst(self, results):
        return self.search("mlst", "type",results["typing_result"])

    def get_cgmlst(self, results):
        return self.search("cgmlst", "type", results["typing_result"])

    def get_mdb_data(self, db_collection, sample_id):
        mdb_pvl = list(Database.get_pvl(db_collection, {"id": sample_id}))
        mdb_mlst = list(Database.get_mlst(db_collection, {"id": sample_id}))
        mdb_cgmlst = list(Database.get_cgmlst(db_collection, {"id": sample_id}))
        mdb_pvl_present = bool(mdb_pvl[0]["aribavir"]["lukS_PV"]["present"])
        mdb_mlst_seqtype = int(mdb_mlst[0]["mlst"]["sequence_type"])
        mdb_mlst_alleles = mdb_mlst[0]["mlst"]["alleles"]
        mdb_cgmlst_alleles = mdb_cgmlst[0]["alleles"]
        return {"pvl": mdb_pvl_present, "mlst_seqtype": mdb_mlst_seqtype, "mlst_alleles": mdb_mlst_alleles, "cgmlst_alleles": mdb_cgmlst_alleles}
    
    def get_fin_data(self, sample_json):
        fin_pvl_present = self.get_pvl(sample_json)
        fin_mlst = self.get_mlst(sample_json)
        fin_cgmlst = self.get_cgmlst(sample_json)
        fin_mlst_seqtype = fin_mlst[0]["result"]["sequence_type"]
        fin_mlst_alleles = fin_mlst[0]["result"]["alleles"]
        fin_cgmlst_alleles = list(fin_cgmlst[0]["result"]["alleles"].values())
        return {"pvl": fin_pvl_present, "mlst_seqtype": fin_mlst_seqtype, "mlst_alleles": fin_mlst_alleles, "cgmlst_alleles": fin_cgmlst_alleles}

    def run(self, input_files, output_fpaths, db_collection):
        for input_idx, input_file in enumerate(input_files):
            with open(input_file, 'r') as fin:
                sample_json = json.load(fin)
                sample_id = self.get_sample_id(sample_json)
                mdb_data_dict = self.get_mdb_data(db_collection, sample_id)
                species_name = self.get_species_name(sample_json)
                fin_data_dict = self.get_fin_data(sample_json)
