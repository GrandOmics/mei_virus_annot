"""
Mobile Element Insertion(MEI) Annotation for DELs. 50% overlap reciprocally
MEI Annotation and Virus Genome Sequence Annotation for INSs.
Archieyoung yangqi2@grandomics.com
"""

import sys
import os
import argparse
from collections import Counter

import sv_vcf
import del_mei_annot
import ins_seq_annot


# del mei annot worker
def run_del_mei_annot(vcf, rmsk_db_file):
    _rmsk_db = del_mei_annot.rmsk_db()
    _rmsk_db.loadDB(rmsk_db_file)

    _del_reader = del_mei_annot.del_reader(vcf)

    for chrom, query_interval in _del_reader:
        results = _rmsk_db.search(chrom, query_interval)
        if results:
            yield query_interval.data, results


# ins seq annot worker
def run_ins_annot(vcf, bam, blast_db, tmp, prog):
    
    ins_seq_fasta = tmp+".ins.fasta"
    ins_seq_annot.run_get_ins_seq_bam(prog, vcf, bam, ins_seq_fasta)
    
    blast_report = tmp+".ins.blast.txt"
    ins_seq_annot.run_blast("blastn", ins_seq_fasta, blast_db, "1",
        blast_report)
    
    for i in ins_seq_annot.ins_annot_iter(blast_report):
        yield i[0], i[2]


def vcf2bed_pe(sv_record, annot):
    return "\t".join(str(i) for i in [sv_record.chrom1, sv_record.pos1,
        sv_record.chrom2, sv_record.pos2, sv_record.svtype, sv_record.id,
        sv_record.svlen, annot])


def reporter(vcf, del_annot_dict, ins_annot_dict, outfile):
    del_num_total = 0
    ins_num_total = 0
    out_fp = open(outfile, "w")
    with open(vcf, "r") as io:
            for line in io:
                line = line.strip()
                if line[0] == "#":
                    continue
                 
                sv = sv_vcf.sv_vcf_record(line)
                if sv.svtype == "DEL":
                    del_num_total += 1
                    if sv.id in del_annot_dict:
                        print(vcf2bed_pe(sv, del_annot_dict[sv.id]),
                            file = out_fp)
                    else:
                        print(vcf2bed_pe(sv, "NA"), file = out_fp)
                if sv.svtype == "INS":
                    ins_num_total += 1
                    if sv.id in ins_annot_dict:
                        print(vcf2bed_pe(sv, ins_annot_dict[sv.id]),
                            file = out_fp)
                    else:
                        print(vcf2bed_pe(sv, "NA"), file = out_fp)
    out_fp.close()

    out_fp = open(outfile+".summary", "w")
    del_counter = Counter(del_annot_dict.values())
    ins_counter = Counter(ins_annot_dict.values())
    print("#DEL annot summay: \n#sv_type\tannot_type\tnumber", file = out_fp)
    for i in del_counter:
        print("DEL\t{}\t{}".format(i, del_counter[i]), file = out_fp)
    print("DEL\ttotal\t{}".format(del_num_total), file = out_fp)
    print("#INS annot summay: \n#sv_type\tannot_type\tnumber", file = out_fp)
    for i in ins_counter:
        print("INS\t{}\t{}".format(i, ins_counter[i]), file = out_fp)
    print("INS\ttotal\t{}".format(ins_num_total), file = out_fp)
    out_fp.close()


def get_args():
    parser = argparse.ArgumentParser(description="Insertion sequence annotation"
        " for sniffles vcf", usage="%(prog)s [options]")
    parser.add_argument("-v", "--vcf", help="sniffles vcf file"
        " [default: %(default)s]", metavar="FILE")
    parser.add_argument("-b", "--bam", help="bam file"
        " [default: %(default)s]", metavar="FILE")
    parser.add_argument("-o", "--outfile", help="Output file prefix"
        " [default: %(default)s]", metavar="STR")

    if len(sys.argv) <= 1:
        parser.print_help()
        exit()
    return parser.parse_args()


def main():
    args = get_args()

    rmsk_db_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
        "../database/rmsk.db")

    del_annot_dict = dict()
    for svid, annot in run_del_mei_annot(args.vcf, rmsk_db_file):
        if annot != "Other":
            del_annot_dict[svid] = annot

    blast_db = os.path.join(os.path.dirname(os.path.abspath(__file__)),
        "../database/Homo_sapiens.mei_virus.db.fasta")
    
    c_sv_ins_seq = os.path.join(os.path.dirname(os.path.abspath(__file__)),
        "../build/bin/sv_ins_seq")
    ins_annot_dict = dict()
    for svid, annot in run_ins_annot(args.vcf, args.bam, blast_db,
        args.outfile+".tmp", c_sv_ins_seq):
        ins_annot_dict[svid] = annot
    
    reporter(args.vcf, del_annot_dict, ins_annot_dict, args.outfile)

if __name__ == "__main__":
    main()
