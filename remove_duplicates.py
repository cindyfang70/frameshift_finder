from Bio import SeqIO
from Bio.Blast import NCBIWWW

ids = []
records = []
seqs = []
for record in SeqIO.parse("frameshift.fasta", "fasta"):
    if record.id not in ids:
        records.append(record)
        ids.append(record.id)
        seqs.append(record.seq)

with open("frameshift_no_duplicates.fasta", "w") as output_handle:
    SeqIO.write(records, output_handle, "fasta")
print(len(ids))
# count = 0
# for record in SeqIO.parse("frameshift_proteins_no_duplicates.fasta", "fasta"):
#     count += 1
#     protein = record.seq
#     result_handle = NCBIWWW.qblast("tblastn", "nr", record.seq)
#     with open("blast_results.xml", "w") as out_handle:
#         out_handle.write(result_handle.read())
#     result_handle.close()
#     print(count)
