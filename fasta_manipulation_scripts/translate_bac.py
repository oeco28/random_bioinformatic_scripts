from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Define bacterial codon table
# Table 11 from https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
bacterial_codon_table = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S",
    "TCC": "S", "TCA": "S", "TCG": "S", "TAT": "Y", "TAC": "Y",
    "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C", "TGA": "*",
    "TGG": "W", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAT": "H",
    "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "ATT": "I", "ATC": "I", "ATA": "I",
    "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S",
    "AGC": "S", "AGA": "R", "AGG": "R", "GTT": "V", "GTC": "V",
    "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A",
    "GCG": "A", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}

def translate_dna_to_protein(dna_seq):
    protein = ""
    unrecognized_codons = set()  # To collect codons not found in the table

    for i in range(0, len(dna_seq), 3):
        codon = dna_seq[i:i+3].upper()  # Ensure the codon is in uppercase

        # Check if codon exists in our table
        if codon in bacterial_codon_table:
            protein += bacterial_codon_table[codon]
        else:
            protein += 'X'
            unrecognized_codons.add(codon)

    if unrecognized_codons:
        print(f"Unrecognized codons found: {', '.join(unrecognized_codons)}")

    return protein

# Rest of the script remains the same...


def main(input_file, output_file):
    translated_sequences = []

    with open(input_file, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            dna_seq = str(record.seq)
            protein_seq = translate_dna_to_protein(dna_seq)
            protein_record = SeqRecord(Seq(protein_seq), id=record.id, description="")
            translated_sequences.append(protein_record)

    with open(output_file, "w") as output:
        SeqIO.write(translated_sequences, output, "fasta")

if __name__ == "__main__":
    input_file = "virulent_genes.fas"  # Modify this with the path to your input FASTA file
    output_file = "virulent_genes.aa.fas"  # Modify this with the path where you want the output FASTA file
    main(input_file, output_file)

