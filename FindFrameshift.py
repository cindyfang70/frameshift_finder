from Bio.Blast import NCBIWWW
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
import re


# fetches genomes based on the genome refseq id and writes them to a gb file
def get_genomes(_id: str) -> None:
    Entrez.email = "cindyfang70@gmail.com"
    with open("my_genbank.gb", "w") as out_handle:
        result_handle = Entrez.efetch(db="nucleotide", id=_id, rettype="gb",
                                      retmode="text")
        out_handle.write(result_handle.read())
        result_handle.close()


# extracts gene ids from the genbank file of the genomes
def locate_gene_id(filename: str) -> list:
    gene_ids = []
    protein_ids = []
    for seq_record in SeqIO.parse(filename, "genbank"):
        for seq_feature in seq_record.features:
            if seq_feature.type == "CDS":
                if protein_seq in seq_feature.qualifiers["protein_id"]:
                    protein_ids.extend(seq_feature.qualifiers["protein_id"])
                    gene_ids.extend(seq_feature.qualifiers["db_xref"])
        # print(seq_record.annotations)
    return gene_ids


# formats the gene ids into an appropriate format
def format_gene_id(gene_ids: list) -> list:
    formatted_genes = []
    for gene in gene_ids:
        formatted_genes.append(gene.split(":")[1])
    return formatted_genes


# fetches the region of the genome corresponding to the gene based on the gene id
def find_genome_region(geneid: str) -> list:
    Entrez.email = "cindyfang70@gmail.com"
    with open("genome.gb", "w") as out_handle:
        request = Entrez.efetch(db="gene", id=geneid, rettype="fasta",
                                retmode="text")
        out_handle.write(request.read())
        request.close()
    file = open("genome.gb", "r")
    region = ""
    for line in file:
        if "Annotation" in line:
            line = line.split("(")
            region = line[1]
            region = region.split(")")[0]
            region = region.split("..")
    for reg in region:
        if "complement" in reg:
            region.append("c")
    for i in range(len(region)):
        n = len(region[i])
        for j in range(len(region[i])):
            if not region[i][j].isdigit():
                n = j
                break
        region[i] = region[i][:n]
    return region


# returns the desired portion of the DNA sequence
def find_extended_DNA_seq(genomeid: str, region: list) -> Seq:
    Entrez.email = "cindyfang70@gmail.com"
    with open("seq.txt", "w") as out_handle:
        request = Entrez.efetch(db="nuccore", id=genomeid, rettype="fasta",
                                retmode="text", seq_start=str(int(region[0])-308),
                                seq_stop=str(int(region[1]) + 272))
        out_handle.write(request.read())
        request.close()
    file = open("seq.txt", "r")
    seq = ""
    for line in file:
        if ">" not in line:
            seq = seq + line
            seq = seq.replace("\n", "")
    sequence = Seq(seq)
    if len(region) == 3:
        sequence = sequence.reverse_complement()
    return sequence


def find_actual_DNA_seq(genomeid: str, region: list) -> Seq:
    Entrez.email = "cindyfang70@gmail.com"
    with open("seq.txt", "w") as out_handle:
        request = Entrez.efetch(db="nuccore", id=genomeid, rettype="fasta",
                                retmode="text",
                                seq_start=region[0],
                                seq_stop=region[1])
        out_handle.write(request.read())
        request.close()
    file = open("seq.txt", "r")
    seq = ""
    for line in file:
        if ">" not in line:
            seq = seq + line
            seq = seq.replace("\n", "")
    sequence = Seq(seq)
    if len(region) == 3:
        sequence = sequence.reverse_complement()
    return sequence


def DNA_to_protein(coding_dna: Seq) -> tuple:
    mrna = coding_dna.transcribe()
    protein = mrna.translate()
    return protein


def find_orf(dna: Seq, actual_length: int) -> list:
    """CUT THEM ALL OFF AT THE FIRST STOP CODON THEN FIND THE LONGEST ONE"""
    dna = str(dna[actual_length - 30:])
    seqs = []
    for i in range(29):
        rna = Seq(dna[i:]).transcribe()
        protein = rna.translate()
        protein = protein.split("*")[0]
        print(protein)
        seqs.append(protein)
    return seqs

def find_longest(proteins: list, g_sequence: Seq) -> list:
    concatenated_proteins = []
    for i in range(len(proteins)):
        if proteins[i] not in g_sequence.translate():
            concatenated_proteins.append(g_sequence.translate() + proteins[i])
    max_length = 0
    max_orf = ""
    for i in range(len(concatenated_proteins)):
        print(max_length)
        if len(concatenated_proteins[i]) >= max_length:
            max_length = len(concatenated_proteins[i])
            max_orf = concatenated_proteins[i]
            print((len(max_orf), len(concatenated_proteins[i])))
    print(max_orf)
    print(len(max_orf))
    return max_orf


def find_frameshift(desired_protein: Seq, orfs: list, g_sequence: Seq) -> tuple:
    start_index = 0
    for orf in orfs:
        if orf == desired_protein:
            start_index = g_sequence.find(orf)
            print(start_index)
            break
    return g_sequence[start_index-60:]


def find_fs_protein(fusion_prot: Seq) -> Seq:
    start_index = fusion_prot.find("*")
    return fusion_prot[start_index + 1:]


if __name__ == "__main__":
    # protein_seq = "NP_040593.1"
    # genome_refseq = "NC_001416"
    # get_genomes(genome_refseq)
    # gene_ids = locate_gene_id("my_genbank.gb")
    # formatted_genes = format_gene_id(gene_ids)
    # print(formatted_genes)
    g_region = find_genome_region("26634148")
    extended_g_sequence = find_extended_DNA_seq("NC_028896.1", g_region)
    g_sequence = find_actual_DNA_seq("NC_028896.1", g_region)
    print(g_sequence.transcribe())
    print(extended_g_sequence.transcribe())
    print(g_sequence[len(g_sequence)-30:].transcribe())
    orfs = find_orf(extended_g_sequence, len(g_sequence))
    proteins = []
    for orf in orfs:
        proteins.append(DNA_to_protein(orf))
    print("****LONGEST****")
    fusion_prot = find_longest(proteins, g_sequence)
    print(fusion_prot)
    print("******FRAMESHIFT AREA******")
    desired_prot = find_fs_protein(fusion_prot)
    print(desired_prot)
    frameshift_region = find_frameshift(desired_prot, orfs, g_sequence)[0]
    print(frameshift_region)


