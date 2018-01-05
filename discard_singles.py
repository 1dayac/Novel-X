import pysam
import sys

samfile_in = pysam.AlignmentFile(sys.argv[1], "rb")
samfile_out = pysam.AlignmentFile(sys.argv[2], "wb", header=samfile_in.header)

current_name = ""
records = []
for r in samfile_in:
    if current_name == "":
        current_name = r.query_name
        records.append(r)
        continue
    if current_name == r.query_name:
        records.append(r)
        continue
    if current_name != r.query_name:
        if len(records) != 1:
            for record in records:
                samfile_out.write(record)
        records[:] = []

        current_name = r.query_name
        records.append(r)
        continue

if len(records) != 1:
    for r in records:
        samfile_out.write(r)
records[:] = []

samfile_in.close()
samfile_out.close()
