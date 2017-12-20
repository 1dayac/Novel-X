import pysam
import sys

samfile_in = pysam.AlignmentFile(sys.argv[1], "rb")
samfile_out = pysam.AlignmentFile(sys.argv[2], "wb", header=samfile_in.header)

current_name = ""

for r in samfile_in:
    if current_name == "":
        current_name = r.query_name
        prev_record = r
        continue
    if current_name == r.query_name:
        samfile_out.write(prev_record)
        samfile_out.write(r)
        current_name = ""
        continue
    if current_name != r.query_name:
        current_name = r.query_name
        prev_record = r
        continue

samfile_in.close()
samfile_out.close()
