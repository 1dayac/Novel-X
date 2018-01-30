import sys
from Bio import SeqIO



current_set = []
current_node = ""
current_chrom = ""


unique = []

non_unique = []


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
        if self.node != other.node:
            return False
        return (abs(self.pos1 - other.pos1) < 1000 or abs(self.pos1 - other.pos2) < 1000 or abs(self.pos2 - other.pos1) < 1000 or abs(self.pos2 - other.pos2) < 1000)

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
    return abs(al1.cont1 - al1.cont2) > 500 or abs(al2.cont1 - al2.cont2) > 500

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
    ref_breakpoint = 0
    if is_rc:
        ref_breakpoint = al1.ref1
    else:
        ref_breakpoint = al1.ref2
    if cont12 < cont21:
        if cont21 - cont12 > 5:
            ins = Insertion(al1.node, cont12, cont21, is_rc, al1.chrom, ref_breakpoint, abs(cont22 - cont21), abs(cont12 - cont11))
            unique.append(ins)
            return True

    if cont22 < cont11:
        ins = Insertion(al1.node, cont22, cont11, is_rc, al1.chrom, ref_breakpoint, abs(cont22 - cont21), abs(cont12 - cont11))
        if ins in non_unique:
            return False
        if cont11 - cont22 > 5:
            unique.append(ins)
            return True
    return False


def AnalyzeCurrentSet(al1, al2):
    if BigAnchors(al1, al2):
        if IsInsertion(al1, al2):
            print(str(al1) + " \t " +  str(al2))



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
with open(sys.argv[1], "r") as nucmer_file:
    lines = nucmer_file.readlines()
    for i in range(1, len(lines) - 1):
        line = lines[i]
        if line.startswith("local misassembly") or line.startswith("relocation"):
            line1 = lines[i-1]
            line2 = lines[i+1]

            chrom = line1.split("\t")[4].strip()
            if chrom == "chrX" and line1.startswith("9913"):
                pass
            ref11 = int(line1.split("\t")[0].strip())
            ref12 = int(line1.split("\t")[1].strip())

            ref21 = int(line2.split("\t")[0].strip())
            ref22 = int(line2.split("\t")[1].strip())

            cont11 = int(line1.split("\t")[2].strip())
            cont12 = int(line1.split("\t")[3].strip())

            cont21 = int(line2.split("\t")[2].strip())
            cont22 = int(line2.split("\t")[3].strip())

            contig_name = line1.split("\t")[5].strip()

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
    if record not in record_dict:
        record_dict[record.id] = record
with open(sys.argv[3], "w") as vcf:
    for insertion in unique:
        ins_seq = str(record_dict[insertion.node].seq[insertion.pos1 : insertion.pos2] if not insertion.rc else record_dict[insertion.node].seq[insertion.pos1 : insertion.pos2].reverse_complement())
        vcf.write(str(insertion.ref_id) + "\t" + str(insertion.ref_pos) + "\t" + "." + "\t" + str(record_dict[insertion.node].seq[insertion.pos1]) + "\t" + ins_seq + "\t" + "." + "\t" + "PASS" + "\t" + "DP=100" + "\t" + insertion.node + "\t" + str(insertion.anchor1) + "\t" + str(insertion.anchor2) + "\n")