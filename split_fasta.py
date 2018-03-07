import sys
from Bio import SeqIO

import os

directory = sys.argv[2]
try:
    os.stat(directory)
except:
    os.mkdir(directory)


current_index = 1
for record in SeqIO.parse(sys.argv[1], "fasta"):
    with open(directory + os.sep + str(current_index) + ".fasta", "w") as output_handle:

        SeqIO.write([record], output_handle, "fasta")
        current_index += 1