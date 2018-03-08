import sys
from Bio import SeqIO

with open(sys.argv[2], "r") as nucmer_file:
    lines = nucmer_file.readlines()
    record_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[1], "fasta"))

    with open(sys.argv[3], "w") as output_handle:
        if len(lines) == 1 or  lines[1].startswith("CONTIG"):
            pass
        else:
            SeqIO.write([record_dict[lines[4].split("\t")[0].strip()]], output_handle, "fasta")
