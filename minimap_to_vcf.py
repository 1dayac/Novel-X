import sys
import re
from Bio import SeqIO



current_set = []
current_node = ""
current_chrom = ""


unique_insertions = []
non_unique_insertions = []
unique_deletions = []
non_unique_deletions = []


class Deletion(object):
    def __init__(self, node, contig_pos, rc, ref_id, ref_pos1, ref_pos2, anchor1, anchor2):
        self.node = node
        self.contig_pos = contig_pos
        self.rc = rc
        self.ref_id = ref_id
        self.ref_pos1 = ref_pos1
        self.ref_pos2 = ref_pos2
        self.anchor1 = anchor1
        self.anchor2 = anchor2

    def __eq__(self, other):
        return (abs(self.ref_pos1 - other.ref_pos1) < 100 or abs(self.ref_pos2 - other.ref_pos2) < 100 or abs(self.ref_pos2 - other.ref_pos1) < 100 or abs(self.ref_pos1 - other.ref_pos2) < 100)

    def size(self):
        return  abs(self.ref_pos2 - self.ref_pos1) - 1

class Insertion(object):
    def __init__(self, node, pos1, pos2, rc, ref_id, ref_pos, anchor1, anchor2):
        self.node = node
        self.pos1 = pos1
        self.pos2 = pos2
        self.rc = rc
        self.ref_id = ref_id
        self.ref_pos = ref_pos
        self.anchor1 = anchor1
        self.anchor2 = anchor2

    def __eq__(self, other):
        return abs(self.ref_pos - other.ref_pos) < 100

    def size(self):
        return abs(self.pos2 - self.pos1) - 1

def Overlaps(al1, al2):
    cont11 = al1.cont1
    cont12 = al1.cont2

    cont21 = al2.cont1
    cont22 = al2.cont2

    if cont11 > cont12:
        cont11, cont12 = cont12, cont11

    if cont21 > cont22:
        cont21, cont22 = cont22, cont21

    return (cont11 >= cont21 and cont11 <= cont22) or (cont12 <= cont21 and cont12 >= cont22)


def Near(sv1, sv2):
    return sv2 - sv1 < 50 and sv2 - sv1 > -10

def BigAnchors(al1, al2):
    return abs(al1.cont1 - al1.cont2) > 300 or abs(al2.cont1 - al2.cont2) > 300



def IsDeletion(al1, al2):
    ref11 = al1.ref1
    ref12 = al1.ref2

    is_rc = ref12 - ref11 < 0
    ref21 = al2.ref1
    ref22 = al2.ref2
    is_rc2 = ref22 - ref21 < 0

    if is_rc != is_rc2:
        return False
    if ref11 > ref12:
        ref11, ref12 = ref12, ref11

    if ref21 > ref22:
        ref21, ref22 = ref22, ref21

    if is_rc:
        contig_breakpoint = al1.cont1
    else:
        contig_breakpoint = al1.cont2

    if cont21 - cont12 < 5:
        if ref21 - ref12 > 5:
            deletion = Deletion(al1.node, contig_breakpoint, is_rc, al1.chrom, ref12, ref21, abs(cont22 - cont21), abs(cont12 - cont11))
            if deletion.size() > 1000:
                print(deletion.size())
                print("Deletion: " + str(al1) + " \t " + str(al2))
            if deletion not in unique_deletions:
                unique_deletions.append(deletion)
            else:
                non_unique_deletions.append(deletion)
            return True

    if cont22 < cont11:
        deletion = Deletion(al1.node, contig_breakpoint, is_rc, al1.chrom, ref12, ref21, abs(cont22 - cont21), abs(cont12 - cont11))
        if deletion.size() > 1000:
            print(deletion.size())
            print("Deletion: " + str(al1) + " \t " + str(al2))

        if deletion in non_unique_deletions:
            return False
        if cont11 - cont22 > 5:
            if deletion not in unique_deletions:
                unique_deletions.append(deletion)
            else:
                non_unique_deletions.append(deletion)
            return True
    return False

def IsInsertion(al1, al2):
    cont11 = al1.cont1
    cont12 = al1.cont2

    is_rc = cont12 - cont11 < 0
    cont21 = al2.cont1
    cont22 = al2.cont2
    is_rc2 = cont22 - cont21 < 0

    if is_rc != is_rc2:
        return False
    if cont11 > cont12:
        cont11, cont12 = cont12, cont11

    if cont21 > cont22:
        cont21, cont22 = cont22, cont21
    if is_rc:
        ref_breakpoint = al1.ref1
    else:
        ref_breakpoint = al1.ref2
    if cont12 < cont21:
        if cont21 - cont12 > 5:
            ins = Insertion(al1.node, cont12, cont21, is_rc, al1.chrom, ref_breakpoint, abs(cont22 - cont21), abs(cont12 - cont11))
            if ins.size() > 500:
                print(ins.size())
                print("Insertion: " + str(al1) + " \t " + str(al2))
            if ins not in unique_insertions:
                unique_insertions.append(ins)
            else:
                non_unique_insertions.append(ins)
            return True


    if cont22 < cont11:
        ins = Insertion(al1.node, cont22, cont11, is_rc, al1.chrom, ref_breakpoint, abs(cont22 - cont21), abs(cont12 - cont11))
        if ins.size() > 500:
            print(ins.size())
            print("Insertion: " + str(al1) + " \t " + str(al2))
        if cont11 - cont22 > 5:
            if ins not in unique_insertions:
                unique_insertions.append(ins)
            else:
                non_unique_insertions.append(ins)
            return True
    return False


def AnalyzeCurrentSet(al1, al2):
    if BigAnchors(al1, al2):
        if IsInsertion(al1, al2):
            pass
            #print(str(al1) + " \t " +  str(al2))
        if IsDeletion(al1, al2):
            pass
            #print(str(al1) + " \t " +  str(al2))


class ContigFilter(object):

    def __init__(self):
        self.coverage_cutoff = 0.0
        self.length_cutoff = 0

    def ParseCoverage(self, contig_name):
        m = re.search("cov_([\d\.]+)", contig_name)
        return float(m.group(1))

    def ParseLength(self, contig_name):
        m = re.search("length_([\d\.]+)_", contig_name)
        return int(m.group(1))


    def PassFilter(self, contig_name):
        return self.ParseCoverage(contig_name) > self.coverage_cutoff \
               and self.ParseLength(contig_name) > self.length_cutoff




class Alignment(object):
    def __init__(self, ref1, ref2, cont1, cont2, chrom, node):
        self.ref1 = ref1
        self.ref2 = ref2
        self.cont1 = cont1
        self.cont2 = cont2
        self.chrom = chrom
        self.node = node

    def __str__(self):
        return str(self.ref1) + "\t" + str(self.ref2) + "\t" + str(self.cont1) + "\t" + str(self.cont2) + "\t" + self.chrom + "\t" + self.node
i = 0

filter = ContigFilter()

with open(sys.argv[1], "r") as nucmer_file:
    lines = nucmer_file.readlines()
    for i in range(1, len(lines) - 1):
        line = lines[i]
        if line.startswith("local misassembly") or line.startswith("transloc") or line.startswith("relocation") or line.startswith("indel") or line.startswith("fake"):
            line1 = lines[i-1]
            line2 = lines[i+1]
            if line2.startswith("CON"):
                continue

            chrom = line1.split("\t")[4].strip()

            ref11 = int(line1.split("\t")[0].strip())
            ref12 = int(line1.split("\t")[1].strip())

            ref21 = int(line2.split("\t")[0].strip())
            ref22 = int(line2.split("\t")[1].strip())

            cont11 = int(line1.split("\t")[2].strip())
            cont12 = int(line1.split("\t")[3].strip())

            cont21 = int(line2.split("\t")[2].strip())
            cont22 = int(line2.split("\t")[3].strip())

            contig_name = line1.split("\t")[5].strip()

            if not filter.PassFilter(contig_name):
                continue

            i += 1
            if i > 100000:
                pass
                #break

            al1 = Alignment(ref11, ref12, cont11, cont12, chrom, contig_name)
            al2 = Alignment(ref21, ref22, cont21, cont22, chrom, contig_name)

            AnalyzeCurrentSet(al1, al2)


records = list(SeqIO.parse(sys.argv[2], "fasta"))
record_dict = {}
for record in records:
    if record.id not in record_dict.keys():
        record_dict[record.id] = record
with open(sys.argv[3], "w") as vcf, open("insertions_with_anchors.fasta", "w") as fasta:
    for insertion in unique_insertions:
        if insertion.size() >= 300:
            fasta.write(">" + insertion.node + "\n")
            fasta.write(str(record_dict[insertion.node].seq))
            fasta.write("\n")

        ins_seq = str(record_dict[insertion.node].seq[insertion.pos1 : insertion.pos2] if not insertion.rc else record_dict[insertion.node].seq[insertion.pos1 : insertion.pos2].reverse_complement())
        vcf.write(str(insertion.ref_id) + "\t" + str(insertion.ref_pos) + "\t" + "." + "\t" + str(record_dict[insertion.node].seq[insertion.pos1]) + "\t" + ins_seq + "\t" + "." + "\t" + "PASS" + "\t" + "DP=100" + "\t" + insertion.node + "\t" + str(insertion.anchor1) + "\t" + str(insertion.anchor2) + "\t" + str(insertion.pos1) + "\t" + str(insertion.pos2) + "\n")
#    for deletion in unique_deletions:
#        del_seq = str(record_dict[deletion.node].seq[deletion.pos1 : deletion.pos2] if not insertion.rc else record_dict[deletion.node].seq[deletion.pos1 : deletion.pos2].reverse_complement())
#        vcf.write(str(deletion.ref_id) + "\t" + str(deletion.ref_pos1) + "\t" + "." + "\t" + str(record_dict[deletion.node].seq[deletion.pos1]) + "\t" + del_seq + "\t" + "." + "\t" + "PASS" + "\t" + "DP=100" + "\t" + deletion.node + "\t" + str(deletion.anchor1) + "\t" + str(deletion.anchor2) + "\t" + str(deletion.pos1) + "\t" + str(deletion.pos2) + "\n")