import os
from Bio import Entrez, SeqIO
from jasentool.utils import Utils

class Genome:
    def __init__(self, refseq_accn, genbank_accn, download_dir, prefix, email="rjkennedyy@gmail.com"):
        Entrez.email = email
        self.refseq_accn = refseq_accn
        self.genbank_accn = genbank_accn
        self.download_dir = download_dir
        self.zip_filepath = os.path.join(download_dir, "GCF_000195955.2.zip")
        self.fasta_filepath = os.path.join(download_dir, f"{prefix}.fasta")
        self.genbank_filepath = os.path.join(download_dir, f"{prefix}.gb")
        self.gff_filepath = os.path.join(download_dir, f"{prefix}.gff")

    def download_fasta(self):
        try:
            # Fetch the fasta record from NCBI
            fasta_handle = Entrez.efetch(db="nucleotide", id=self.refseq_accn, rettype="fasta", retmode="text")
            fasta_record = SeqIO.read(fasta_handle, "fasta")
            fasta_handle.close()

            # Save the fasta record to a file
            SeqIO.write(fasta_record, self.fasta_filepath, "fasta")

            print(f"Fasta downloaded and saved to {self.fasta_filepath}")

        except Exception as e:
            print(f"Error downloading the genome: {e}")
        return self.fasta_filepath

    def download_genbank(self):
        try:
            # Fetch the GenBank record from NCBI
            genbank_handle = Entrez.efetch(db="nucleotide", id=self.genbank_accn, rettype="gb", retmode="text")
            genbank_record = SeqIO.read(genbank_handle, "genbank")
            genbank_handle.close()

            # Save the GenBank record to a file
            SeqIO.write(genbank_record, self.genbank_filepath, "genbank")

            print(f"Genbank file downloaded and saved to {self.genbank_filepath}")

        except Exception as e:
            print(f"Error downloading the genbank file: {e}")
        return self.genbank_filepath
    
    def download_gff(self):
        utils = Utils()
        h37rv_url = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000195955.2/download?include_annotation_type=GENOME_GFF&filename=GCF_000195955.2.zip"
        try:
            utils.download_and_save_file(h37rv_url, self.zip_filepath)
            utils.unzip(self.zip_filepath, self.download_dir)
            source = os.path.join(self.download_dir, "ncbi_dataset/data/GCF_000195955.2/genomic.gff")
            destination = os.path.join(self.download_dir, "h37rv.gff")
            utils.copy_file(source, destination)
        except Exception as e:
            print(f"Error downloading the gff file: {e}")
        return self.gff_filepath
