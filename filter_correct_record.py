import sys
from Bio import SeqIO

with open(sys.argv[2], "r") as nucmer_file:
    lines = nucmer_file.readlines()
    record_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[1], "fasta"))
    if len(lines) == 1:
        exit()

    records_to_output = []
    with open(sys.argv[3], "w") as output_handle:
        for line in lines:
            split = line.split("\t")
            if len(split) < 5 or split[0] == "S1":
                pass
            else:
                if split[4].strip() not in records_to_output:
                    records_to_output.append(split[4].strip())
        for record in records_to_output:
            SeqIO.write([record_dict[record]], output_handle, "fasta")
