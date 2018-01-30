#!/usr/bin/env python
import sys
import re
from Bio import SeqIO

def parse_coverage(ctg_name):
   m = re.search("cov_([\d\.]+)", ctg_name) 
   return float(m.group(1))
    

if len(sys.argv) < 4: 
    print("Usage: %s <coverage> <contigs_file> <output>" % sys.argv[0])
    sys.exit(1) 

threshold = float(sys.argv[1])
f_n = sys.argv[2]
print("Filtering %s for contigs with coverage higher than %.1f" % (f_n, threshold))
input_seq_iterator = SeqIO.parse(open(f_n, "r"), "fasta")
filtered_iterator = (record for record in input_seq_iterator \
                      if parse_coverage(record.name) > threshold)
 
output_handle = open(sys.argv[3], "w")
SeqIO.write(filtered_iterator, output_handle, "fasta")
output_handle.close()
