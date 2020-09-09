from Bio import SeqIO
import glob

def keep_only_longest_sequence(all_fasta, longest_seq_fasta):
    myList = []
    for seq_record in SeqIO.parse(all_fasta, "fasta"):
        myList.append([seq_record.id, str(seq_record.seq), len(seq_record)])

    myList.sort(key=lambda x: x[2])

    with open(longest_seq_fasta, "w") as fout:
        print(">", myList[-1][0], sep='', file=fout)
        print(myList[-1][1], file=fout)


dir="/home/leandro/git/pandora1_paper/snippy_debug/refs"
files = glob.glob(f"{dir}/*_all.fa")

for file in files:
    longest_seq_file = file.replace("_all.fa", "_chrom.fa")
    keep_only_longest_sequence(file, longest_seq_file)