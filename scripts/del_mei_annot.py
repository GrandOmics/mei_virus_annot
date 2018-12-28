"""
Mobile Element Insertion Annotation for DELs.
50% overlap reciprocally
"""
import sys
from collections import Counter

import sv_vcf
from intervaltree import Interval, IntervalTree

def format_rmsk(rmsk_file, rmsk_out):
    out = open(rmsk_out, "w")
    with open(rmsk_file, "r") as io:
        header = io.readline()
        header = "#rmsk_id" + "\t" + header.strip() 
        print(header, file=out)
        n = 0
        for line in io:
            line = "rmsk_{}".format(n) + "\t" + line.strip()
            print(line, file=out)
            n += 1        
    out.close()

class rmsk(object):
    def __init__(self, record):
        """rmsk_id bin, swScore, milliDiv, milliDel, milliIns, genoName,
        genoStart, genoEnd, genoLeft, strand, repName, repClass, repFamily,
        repStart, repEnd, repLeft, id"""
        self.fields = record.strip().split("\t")
        self.chr = self.fields[6]
        self.start = int(self.fields[7])
        self.end = int(self.fields[8])
        self.id = self.fields[0]

    @property
    def mei_type(self):
        if self.fields[11][:2] == "L1":
            return "L1"
        elif self.fields[11][:3] == "Alu":
            return "Alu"
        elif self.fields[11][:3] == "SVA":
            return "SVA"
        else:
            return "Other"

class rmsk_db(object):
    def __init__(self):
        self.db = dict()

    def loadDB(self, rmsk_formated):
        with open(rmsk_formated, "r") as io:
            io.readline() # remove header
            for line in io:
                rmsk_record = rmsk(line)
                if rmsk_record.chr not in self.db:
                    self.db[rmsk_record.chr] = IntervalTree()
                    self.db[rmsk_record.chr].addi(rmsk_record.start,
                        rmsk_record.end, rmsk_record)
                else:
                    self.db[rmsk_record.chr].addi(rmsk_record.start,
                        rmsk_record.end, rmsk_record)

    def search(self, chrom, start, end):
        if chrom not in self.db:
            return 0
        overlaps = self.db[chrom][start:end]
        if overlaps != set(): # not a empty set
            for _iterval in overlaps:
                overlap = min([_iterval.end, end]) - max([_iterval.begin,
                    start])
                if overlap/(end - start) >= 0.5 and overlap/(_iterval.end -
                    _iterval.begin) >= 0.5:
                    return _iterval.data
                else:
                    return 0
        else:
            return 0

def main():
    if len(sys.argv) < 4:
        print("Usage: python3 {} <vcf> <rmsk> <outfile>")
        sys.exit(1)
    
    _rmsk_db = rmsk_db()
    _rmsk_db.loadDB(sys.argv[2])

    out = open(sys.argv[3],"w")

    summary = []

    with open(sys.argv[1],"r") as io:
        for line in io:
            if line[0] == "#":
                continue
            sv = sv_vcf.sv_vcf_record(line)
            if (sv.svtype == "DEL" and sv.svlen != "NA" and int(sv.svlen) > 100
                and int(sv.svlen) < 11000):
                if "chr" not in sv.chrom1:
                    sv.chrom1 = "chr" + sv.chrom1
                if int(sv.pos1) > int(sv.pos2):
                    raise RuntimeError("pos1 > pos2 for DEL {}".format(sv.id))
                target = _rmsk_db.search(sv.chrom1, int(sv.pos1), int(sv.pos2))
                if target != 0:
                    summary.append(target)
                    print("{}\t{}\t{}\t{}".format(sv.id,
                        sv.svtype, target.mei_type, target.id), file=out)
    out.close()

    out = open(sys.argv[3]+".summary","w")
    _c = Counter(summary)
    for i in _c:
        print("{}\t{}".format(i,_c[i]),file=out)
    out.close()

if __name__ == "__main__":
    main()


