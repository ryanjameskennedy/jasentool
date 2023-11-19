import os
from Bio import Entrez, SeqIO

class Genome:
    def __init__(self, refseq_accn, genbank_accn, download_dir, prefix, email="rjkennedyy@gmail.com"):
        Entrez.email = email
        self.refseq_accn = refseq_accn
        self.genbank_accn = genbank_accn
        self.fasta_filepath = os.path.join(download_dir, f"{prefix}.fasta")
        self.genbank_filepath = os.path.join(download_dir, f"{prefix}.gb")

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