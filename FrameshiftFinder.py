from Bio.Blast import NCBIWWW
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
import re

protein_seq = "NP_040593.1"
genome_refseq = "NC_001416"


def blast(protein_sequence) -> str:
    result_handle = NCBIWWW.qblast("blastp", "nr", protein_sequence)
    xml = result_handle.read()
    return xml


def get_genomes(_id: str) -> None:
    Entrez.email = "cindyfang70@gmail.com"
    with open("my_genbank.gb", "w") as out_handle:
        result_handle = Entrez.efetch(db="nucleotide", id=_id, rettype="gb",
                                      retmode="text")
        out_handle.write(result_handle.read())
        result_handle.close()


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


def format_gene_id(gene_ids: list) -> list:
    formatted_genes = []
    for gene in gene_ids:
        formatted_genes.append(gene.split(":")[1])
    return formatted_genes


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
    return region


def find_DNA_seq(genomeid: str, region: list) -> Seq:
    Entrez.email = "cindyfang70@gmail.com"
    with open("seq.txt", "w") as out_handle:
        request = Entrez.efetch(db="nuccore", id=genomeid, rettype="fasta",
                                retmode="text", seq_start=region[0],
                                seq_stop=str(int(region[1])))
        out_handle.write(request.read())
        request.close()
    file = open("seq.txt", "r")
    seq = ""
    for line in file:
        if ">" not in line:
            seq = seq + line
            seq = seq.replace("\n", "")
    sequence = Seq(seq)
    return sequence


def DNA_to_protein(coding_dna: Seq) -> tuple:
    mrna = coding_dna.transcribe()
    protein = mrna.translate()
    return protein


def find_orf(dna: Seq):
    dna = str(dna)
    pattern = re.compile(r'(?=((ATG|GTG)(?:...)*?)(?=TAG|TGA|TAA))')
    orfs = set(pattern.findall(dna))
    seqs = []
    for orf in orfs:
        seqs.append(DNA_to_protein(Seq(orf[0])))
    return seqs



if __name__ == "__main__":
    get_genomes(genome_refseq)
    gene_ids = locate_gene_id("my_genbank.gb")
    formatted_genes = format_gene_id(gene_ids)
    v_region = find_genome_region("2703487")
    g_region = find_genome_region("2703488")
    h_region = find_genome_region("2703511")
    t_region = find_genome_region("2703489")
    # v_sequence = find_DNA_seq("NC_001416", v_region)
    g_sequence = find_DNA_seq("NC_001416", g_region)
    # h_sequence = find_DNA_seq("NC_001416", h_region)
    # t_sequence = find_DNA_seq("NC_001416", t_region)

    # get_genomes("NC_028896")
    # gene_ids = locate_gene_id("my_genbank.gb")
    # formatted_genes = format_gene_id(gene_ids)
    # g_region = find_genome_region("26634148")
    # h_region = find_genome_region("26634146")
    # t_region = find_genome_region("2703489")
    # g_sequence = find_DNA_seq("NC_028896", g_region)
    # h_sequence = find_DNA_seq("NC_028896", h_region)
    # t_sequence = find_DNA_seq("NC_001416", t_region)

    # print(t_sequence in (g_sequence + h_sequence))
    # max_orf = find_orf(v_sequence + g_sequence + h_sequence[:15])
    # print(v_sequence + g_sequence + h_sequence[:15])
    # print(max_orf)
    # full_seq = v_sequence + g_sequence + h_sequence
    #
    # print(str(t_sequence).replace("G", "A", 1) in full_seq)
    # print(str(t_sequence).replace("G", "A", 1))
    # print(full_seq)
    # print(full_seq.translate())
    print(g_sequence)
    print(g_sequence.translate())

    # print(blast(max_orf))

"""
stops as soon as it sees a stop codon
frameshift mutation 
need to try starting at +1 or -1? see if you can get it to skip over the stop 
because regex doesn't detect all possible reading frames?? 


how to find the whole sequence if you don't know where the frameshift is? 
could just keep going back to try and make the orf longer and longer? 
take gpG sequence and concatenate it with the rest of the bases until the next
orf 
could blast it to see which one it matches up with
then see how many AAs are missing
then turn that into how many bases are missing
and find the whole sequence that way 
"""
