import pysam
import sys

samfile_in = pysam.AlignmentFile(sys.argv[1], "rb")
samfile_out = pysam.AlignmentFile(sys.argv[2], "wb", header=samfile_in.header)

def HasBothReads(records):
    has_first = False
    has_second = False
    for record in records:
        if record.is_read1:
            has_first = True
        if record.is_read2:
            has_second = True

    return has_first and has_second


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
            if HasBothReads(records):
                for record in records:
                    samfile_out.write(record)
        records[:] = []

        current_name = r.query_name
        records.append(r)
        continue

if len(records) != 1:
    if HasBothReads(records):
        for r in records:
            samfile_out.write(r)

records[:] = []

samfile_in.close()
samfile_out.close()
