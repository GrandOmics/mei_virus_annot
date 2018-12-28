"""
author: Archieyoung <yangqi2@grandomics.com>
"""
import sys
import subprocess
import os
import argparse


def run_get_ins_seq(prog, vcf, bam, output):
    out_fp = open(output, "w")
    subprocess.run([prog, vcf, bam], stdout=out_fp)
    out_fp.close()


def run_blast(prog, query, db, nthread, output):
    out_fp = open(output, "w")
    outfmt = ("6 qaccver saccver qlen slen qstart qend sstart send sstrand "
        "pident length mismatch gapopen evalue bitscore")
    subprocess.run([prog, "-query", query, "-db", db, "-outfmt", outfmt,
        "-num_threads", nthread], stdout=out_fp)
    out_fp.close()


class hsp(object):
    def __init__(self, record_fields):
        (self.qaccver, self.saccver, self.qlen, self.slen, self.qstart,
            self.qend, self.sstart, self.send, self.sstrand, self.pident,
            self.length, self.mismatch, self.gapopen, self.evalue,
            self.bitscore) = record_fields

        int_vals = (self.qlen, self.slen, self.qstart, self.qend,
            self.sstart, self.send, self.length, self.mismatch,
            self.gapopen)
        (self.qlen, self.slen, self.qstart, self.qend,
            self.sstart, self.send, self.length, self.mismatch,
            self.gapopen) = [int(i) for i in int_vals]
        
        float_vals = (self.pident, self.evalue, self.bitscore)
        (self.pident, self.evalue, self.bitscore) = [float(i)
            for i in float_vals]

        if self.send < self.sstart:
            tmp = self.send
            self.send = self.sstart
            self.sstart = tmp


def overlap(a_start, a_end, b_start, b_end):
    if a_end < b_start or a_start > b_end:
        return 0
    else:
        return min(a_end, b_end) - max(a_start, b_start) + 1            


class hit(object):
    def __init__(self, hsps):
        self.hsps = sorted(hsps, key = lambda x: x.qstart)
        self.qaccver = hsps[0].qaccver
        self.saccver = hsps[0].saccver
        self.qlen = hsps[0].qend
        self.slen = hsps[0].slen
    
    @property
    def mean_pident(self):
        return sum([i.pident for i in self.hsps])/len(self.hsps)

    @property
    def query_cov(self):
        if len(self.hsps) == 1:
            return (self.hsps[0].qend - self.hsps[0].qstart +
                1)/self.hsps[0].qlen
        
        accumulat_len = 0
        pre = None
        for i in self.hsps:
            if pre == None:
                accumulat_len = i.qend - i.qstart + 1
            else:
                _overlap = overlap(i.qstart, i.qend, pre.qstart, pre.qend)
                if _overlap:
                    if i.qend > pre.qend:
                        accumulat_len += i.qend - pre.qend
                else:
                    accumulat_len += i.qend - i.qstart + 1
            pre = i
        return accumulat_len/self.hsps[0].qlen
    
    @property
    def subject_cov(self):
        if len(self.hsps) == 1:
            return (self.hsps[0].send - self.hsps[0].sstart +
                1)/self.hsps[0].slen
        
        hsps_tmp = sorted(self.hsps, key = lambda x: x.sstart)

        accumulat_len = 0
        pre = None
        for i in hsps_tmp:
            if pre == None:
                accumulat_len = i.send - i.sstart + 1
            else:
                _overlap = overlap(i.sstart, i.send, pre.sstart, pre.send)
                if _overlap:
                    if i.send > pre.send:
                        accumulat_len+= i.send - pre.send
                else:
                    accumulat_len += i.send - i.sstart + 1
            pre = i
        return accumulat_len/self.hsps[0].slen


class report_reader(object):
    def __init__(self, file):
        self.file = file
    
    def __iter__(self):
        qaccver_pre = None
        records = []
        with open(self.file, "r") as io:
            for line in io:
                fields = line.strip().split("\t")
                _hsp = hsp(fields)
                if fields[0] != qaccver_pre:
                    if qaccver_pre != None:
                        yield records
                    records = [_hsp]
                else:
                    records.append(_hsp)
                qaccver_pre = fields[0]
            yield records


class hit_iter(object):
    def __init__(self, records):
        self.records = records
    
    def __iter__(self):
        saccver_pre = None
        hsps = []
        for i in self.records:
            if i.send < i.sstart:
                print("Error!")
            if i.saccver != saccver_pre:
                if saccver_pre != None:
                    yield hit(hsps)
                hsps = [i]
            else:
                hsps.append(i)
            saccver_pre = i.saccver
        yield hit(hsps)


def get_best_hit(_hit_iter):
    """
    max coverage on query and coverage on subject
    """
    _max = 0
    _best_hit = None
    for hit in _hit_iter:
        cov_product = hit.query_cov * hit.subject_cov
        if cov_product > _max:
            _max = cov_product
            _best_hit = hit
    return _best_hit

def run_ins_annot(blast_report, outfile):
    reader = report_reader(blast_report)
    with open(outfile, "w") as io:
        print("#id\ttarget\ttarget_type\tquery_length\ttarget_lenght\t"
            "query_cov\ttarget_cov\tmean_identity", file=io)
        for report in reader:
            _best_hit = get_best_hit(hit_iter(report))
            if "." in _best_hit.saccver:
                target_type = _best_hit.saccver.split(".")[0]
            else:
                target_type = _best_hit.saccver
            fields = [_best_hit.qaccver, _best_hit.saccver, target_type,
                _best_hit.qlen, _best_hit.slen, _best_hit.query_cov,
                _best_hit.subject_cov, _best_hit.mean_pident]
            line = "\t".join([str(i) for i in fields])
            print(line, file = io)

def get_args():
    parser = argparse.ArgumentParser(description="Insertion sequence annotation"
        " for sniffles vcf", usage="%(prog)s [options]")
    parser.add_argument("--vcf", help="sniffles vcf file"
        " [default: %(default)s]", metavar="FILE")
    parser.add_argument("--bam", help="bam file which must match "
        "with sniffles vcf file [default: %(default)s]", metavar="FILE")
    parser.add_argument("--prefix", help="Output file prefix"
        " [default: %(default)s]", metavar="STR")

    if len(sys.argv) <= 1:
        parser.print_help()
        exit()
    return parser.parse_args()

            
def main():
    args = get_args()

    sv_ins_seq_prog = os.path.join(os.path.dirname(os.path.abspath(__file__)),
        "../build/bin/sv_ins_seq")
    
    ins_seq_fasta = args.prefix+".ins.fasta"
    run_get_ins_seq(sv_ins_seq_prog, args.vcf, args.bam,
        ins_seq_fasta)
    
    blast_db = os.path.join(os.path.dirname(os.path.abspath(__file__)),
        "../database/Homo_sapiens.mei_virus.db.fasta")
    blast_report = args.prefix+".ins.blast.txt"
    run_blast("blastn", ins_seq_fasta, blast_db, "1",
        blast_report)
    
    run_ins_annot(blast_report, args.prefix+".ins.annot.txt")

if __name__ == "__main__":
    main()