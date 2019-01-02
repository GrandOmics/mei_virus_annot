"""
Mobile Element Insertion Annotation for DELs.
50% overlap reciprocally
"""
import sys

import sv_vcf
import bin_index

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
        """rmsk_id, bin, swScore, milliDiv, milliDel, milliIns, genoName,
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
                _interval = bin_index.interval(rmsk_record.start,
                    rmsk_record.end, rmsk_record.mei_type)
                if rmsk_record.chr not in self.db:
                    self.db[rmsk_record.chr] = bin_index.BinIndex()
                    self.db[rmsk_record.chr].add_interval(_interval)
                else:
                    self.db[rmsk_record.chr].add_interval(_interval)

    def search(self, chrom, query_interval):
        if chrom not in self.db:
            return 0
        overlaps = self.db[chrom].get_overlap(query_interval)
        if overlaps:
            for _interval in overlaps:
                overlap = min(_interval.end, query_interval.end) - max(
                    _interval.start, query_interval.start) + 1
                if (overlap/query_interval.size >= 0.5 and
                    overlap/_interval.size) >= 0.5:
                    if _interval.data == "Alu":
                        if (abs(query_interval.start - _interval.start) <= 20
                            and abs(query_interval.end - _interval.end) <= 20):
                            return _interval.data # first match
                    elif _interval.data == "L1" or _interval.data == "SVA":
                        if (abs(query_interval.start - _interval.start) <= 200
                            and abs(query_interval.end - _interval.end) <= 200):
                            return _interval.data # first match
            return 0 # no match
        else:
            return 0

class del_reader(object):
    def __init__(self, vcf):
        self.vcf = vcf
    
    def __iter__(self):
        with open(self.vcf,"r") as io:
            for line in io:
                if line[0] == "#":
                    continue
                sv = sv_vcf.sv_vcf_record(line)
                if (sv.svtype == "DEL" and sv.svlen != "NA" and
                    abs(sv.svlen) > 100 and abs(sv.svlen) < 11000):
                    if "chr" not in sv.chrom1:
                        sv.chrom1 = "chr" + sv.chrom1
                    if int(sv.pos1) > int(sv.pos2):
                        raise RuntimeError("pos1 > pos2 for DEL {}".format(
                            sv.id))
                    yield sv.chrom1, bin_index.interval(int(sv.pos1),
                        int(sv.pos2), sv.id)       

def main():
    if len(sys.argv) < 4:
        print("Usage: python3 {} <vcf> <rmsk> <outfile>")
        sys.exit(1)
    
    _rmsk_db = rmsk_db()
    _rmsk_db.loadDB(sys.argv[2])

    out = open(sys.argv[3],"w")
    for chrom, q_interval in del_reader(sys.argv[1]):
        result = _rmsk_db.search(chrom, q_interval)
        if result:
            print("{}\t{}".format(q_interval.data, result), file=out)
    out.close()

if __name__ == "__main__":
    main()