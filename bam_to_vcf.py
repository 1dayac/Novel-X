import pysam
#import vcf
import sys
import os
from Bio.Seq import Seq

def Usage():
    print("bam_to_vcf.py iterates through all sam-files in a given directory and constructs VCF files with information abou insertions")
    print("Usage: bam_to_vcf.py <directory with sam-files> <path_to_vcf_file>")


if sys.argc != 2:
    Usage()
    exit(0)

onlyfiles = [os.join(sys.argv[1], f) for f in os.listdir(sys.argv[1]) if f.endswith(".sam")]
print(onlyfiles)
with open(sys.argv[2], "w") as vcf:
    for file in onlyfiles:
        with pysam.AlignmentFile(sys.argv[1], "rs") as samfile_in:
            for aligned_segment in samfile_in:
                ref_id = aligned_segment.reference_id
                cigartuples = aligned_segment.cigartuples
                ref_start = aligned_segment.reference_start
                read_start = 0
                seq = Seq(aligned_segment.query_sequence)
                if aligned_segment.is_reverse:
                    seq = seq.reverse_complement()
                for tuple in cigartuples:
                    if tuple[0] == 0:
                        ref_start += tuple[1]
                        read_start += tuple[1]
                    if tuple[0] == 1:
                        vcf.write(ref_id + "\t" + ref_start + "\t" + "." + "\t" + seq[read_start] + "\t" + seq[read_start:read_start + tuple[1]] + "\t" + "." + "\t" + "PASS" + "\t" + "DP=100" + "\n")
                        read_start += tuple[1]

                    if tuple[0] == 2:
                        ref_start += tuple[1]

                    if tuple[0] == 4:
                        read_start += tuple[1]