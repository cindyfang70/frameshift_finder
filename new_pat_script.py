import xlrd
import pandas as pd
from FindFrameshift import *
import urllib.error
from Bio.SeqRecord import SeqRecord


class Organism:

    def __init__(self, species_name, genome_refseq, gpg_refseq=None,
                 gpt_refseq=None, gph_refseq=None):
        self.name = species_name
        self.genome_refseq = genome_refseq
        # tail assembly chaperone
        self.tg_refseq = gpg_refseq
        # tail assembly chaperone c-terminal frameshift
        self.th_refseq = gpt_refseq
        # tape measure
        self.tm_refseq = gph_refseq
        self.tg_geneid = None
        self.tm_geneid = None
        self.th_geneid = None
        self.extended_tg_sequence = None
        self.tg_length = None
        self.tg_sequence = None
        self.th_sequence = None
        self.fusion_protein = None
        self.iscDNA = False

    def print_organism(self):
        print(self.name, self.genome_refseq, self.tg_refseq, self.th_refseq,
              self.tm_refseq, self.tg_geneid, self.tm_geneid,
              self.th_geneid)


th_df = pd.read_excel("gpt_dataa.xlsx", sheet_name="gpT")
tg_df = pd.read_excel("gpt_dataa.xlsx", sheet_name="gpG")
tm_df = pd.read_excel("gpt_dataa.xlsx", sheet_name="gpH")

th_df.drop_duplicates(subset="Species", keep="first", inplace=True)
tg_df.drop_duplicates(subset="Species", keep="first", inplace=True)
tm_df.drop_duplicates(subset="Species", keep="first", inplace=True)

th_df.reset_index(inplace=True)
tg_df.reset_index(inplace=True)


 tm_df.reset_index(inplace=True)
th_len = len(th_df)
tg_len = len(tg_df)
tm_len = len(tm_df)

th_organisms = []
for i in range(th_len):
    species_name = th_df["Species"][i]
    genome_refseq = th_df["Genome RefSeq"][i]
    protein_refseq = th_df["Protein RefSeq"][i]
    organism = Organism(species_name, genome_refseq, None, protein_refseq)
    th_organisms.append(organism)

tg_organisms = []
for i in range(tg_len):
    species_name = tg_df["Species"][i]
    genome_refseq = tg_df["Genome RefSeq"][i]
    protein_refseq = tg_df["Protein RefSeq"][i]
    tg_organism = Organism(species_name, genome_refseq)
    tg_organism.tg_refseq = protein_refseq
    tg_organisms.append(tg_organism)

tm_organisms = []
for i in range(tm_len):
    species_name = tm_df["Species"][i]
    genome_refseq = tm_df["Genome RefSeq"][i]
    protein_refseq = tm_df["Protein RefSeq"][i]
    tm_organism = Organism(species_name, genome_refseq, None, None)
    tm_organism.tm_refseq = protein_refseq
    tm_organisms.append(tm_organism)

shared_by_all = []
for tg_organism in tg_organisms:
    for tm_organism in tm_organisms:
        if tg_organism.genome_refseq == tm_organism.genome_refseq:
            tg_organism.tm_refseq = tm_organism.tm_refseq
            shared_by_all.append(tg_organism)

gene_ids = []
for tg_organism in shared_by_all:
    Entrez.email = "cindyfang70@gmail.com"
    try:
        with open("pat_organism_genome.gb", "a") as out_handle:
            result_handle = Entrez.efetch(db="nucleotide",
                                          id=tg_organism.genome_refseq,
                                          rettype="gb", retmode="text")
            out_handle.write(result_handle.read())
            result_handle.close()
    except urllib.error.HTTPError:
        print(tg_organism.genome_refseq)
gene_ids = []
translations = []
for seq_record in SeqIO.parse("pat_organism_genome.gb", "genbank"):
    for seq_feature in seq_record.features:
        if seq_feature.type == "CDS" and "db_xref" in seq_feature.qualifiers:
            gene_ids.append((seq_feature.qualifiers["db_xref"],
                             seq_feature.qualifiers["protein_id"]))

        if "translation" in seq_feature.qualifiers:
            translations.append((seq_feature.qualifiers["translation"],
                                 seq_feature.qualifiers["protein_id"]))

for gene_id in gene_ids:
    for tg_organism in shared_by_all:
        if tg_organism.tg_refseq in gene_id[1]:
            gene = gene_id[0][0]
            gene = gene.split(":")[1]
            tg_organism.tg_geneid = gene
        if tg_organism.tm_refseq in gene_id[1]:
            gene = gene_id[0][0]
            gene = gene.split(":")[1]
            tg_organism.tm_geneid = gene
        if tg_organism.th_refseq in gene_id[1]:
            gene = gene_id[0][0]
            gene = gene.split(":")[1]
            tg_organism.th_geneid = gene
print("finished checking geneids")
for translation in translations:
    for tg_organism in shared_by_all:
        if tg_organism.th_refseq in translation[1]:
            tg_organism.th_sequence = translation[0]

num_match = 0
num_mismatch = 0
frameshift_regions = []
TG_proteins = []
desired_prots = []
for tg_organism in shared_by_all:
    print("***********************")
    try:
        g_region = find_genome_region(tg_organism.tg_geneid)
        h_region = find_genome_region(tg_organism.tm_geneid)
    except AttributeError:
        print(tg_organism.tg_geneid)
        print(tg_organism.tm_geneid)
    if len(g_region) == 3 or len(h_region) == 3:
        tg_organism.iscDNA = True
    print(g_region, h_region)
    try:
        h_sequence = find_actual_DNA_seq(tg_organism.genome_refseq, h_region)
    except urllib.error.HTTPError:
        pass
    except IndexError:
        tg_organism.print_organism()
    try:
        g_sequence = find_actual_DNA_seq(tg_organism.genome_refseq, g_region)
        tg_organism.tg_sequence = g_sequence
        tg_protein = g_sequence.transcribe().translate()
        TG_proteins.append(SeqRecord(tg_protein, tg_organism.name+ "|" + tg_organism.genome_refseq))
        tg_organism.tg_length = len(g_sequence)
    except (RuntimeError, urllib.error.HTTPError, IndexError):
        pass
    # if len(g_region) == 3:
    #     extended_g_region = (int(g_region[0]) - len(h_sequence), g_region[1])
    # else:
    try:
        extended_g_region = (g_region[0], int(g_region[1]) + len(h_sequence))
    except IndexError:
        pass
    try:
        extended_g_sequence = find_actual_DNA_seq(tg_organism.genome_refseq,
                                                  extended_g_region)
        tg_organism.extended_tg_sequence = extended_g_sequence
    except urllib.error.HTTPError:
        pass
    if tg_organism.tg_length != None:
        orfs = find_orf(tg_organism.extended_tg_sequence, tg_organism.tg_length)
    proteins = []
    for orf in orfs:
        proteins.append(orf)
    if tg_organism.tg_sequence != None:
        fusion_prot = find_longest(proteins, tg_organism.tg_sequence)
        tg_organism.fusion_protein = fusion_prot
        print("fusion: ", fusion_prot)
        desired_prot = find_fs_protein(fusion_prot)
        desired_prots.append(SeqRecord(desired_prot,
                                       tg_organism.name + "|" + tg_organism.genome_refseq))
        print("desired prot:", desired_prot)
        frameshift_region = find_frameshift(desired_prot, orfs,
                                            tg_organism.extended_tg_sequence)

        fs_start = tg_organism.extended_tg_sequence.find(frameshift_region[:60])
        print("fs_start: ", fs_start)
        try:
            g_region_seq = find_actual_DNA_seq(tg_organism.genome_refseq, g_region)
            if tg_organism.iscDNA:
                g_region_seq = g_region_seq.reverse_complement()
        except urllib.error.HTTPError:
            pass

        fs_region = (int(g_region[0]) + fs_start,
                     int(g_region[0]) + fs_start + len(
                         frameshift_region[:60]) - 1)
        frameshift_regions.append(SeqRecord(frameshift_region[:60],
                                            tg_organism.name + "|" + tg_organism.genome_refseq + " " + "(" + str(fs_region[0]) + ", " + str(fs_region[1]) + ")",
                                            description="Tail Assembly Chaperone-C-Terminal Frameshift"))
        try:
            found_region = find_actual_DNA_seq(tg_organism.genome_refseq,
                                               fs_region)
        except urllib.error.HTTPError:
            pass

# alignment = MultipleSeqAlignment(frameshift_regions)
# print(alignment)

with open("frameshift.fasta", "w") as output_handle:
    SeqIO.write(frameshift_regions, output_handle, "fasta")

with open("frameshifted_proteins.fasta", "w") as output_handle:
    SeqIO.write(desired_prots, output_handle, "fasta")

with open("TG_proteins.fasta", "w") as output_handle:
    SeqIO.write(TG_proteins, output_handle, "fasta")
