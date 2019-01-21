#include "ins.h"

int main(int argc, char const *argv[]) {
    if (argc < 3) {
        std::cerr << 
            "Usage: sv_ins_seq <vcf> <bam> 1> ins.seq.fasta 2> error.log"
            << std::endl;
        std::exit(1);
    }

    vcfFile *fp_v = vcf_open(argv[1], "r");
    bcf_hdr_t *h_v = bcf_hdr_read(fp_v);
    bcf1_t *v = bcf_init();

    samFile *fp_b = sam_open(argv[2], "r");
    bam_hdr_t *h_b = sam_hdr_read(fp_b);
    hts_idx_t *idx_b = sam_index_load(fp_b, argv[2]);
    if (idx_b == nullptr) {
        std::cerr << "[sv_ins_seq::Error]: Can not load bam index!"
            << std::endl;
        std::exit(1);
    }

    while (vcf_read1(fp_v, h_v, v) >= 0) {
        // const char *chrom1 = bcf_hdr_id2name(h_v, v->rid);
        fprintf(stderr, "%d\n", v->d.id);
        int n1 = 0;
        char *svtype = NULL;
        int ret1 = bcf_get_info_string(h_v,v,"SVTYPE",&svtype,&n1);
        if (ret1 < 0) {
            fprintf(stderr, "[sv_ins_seq::Error]: Can not get SVTYPE in vcf "
                "info field for %s\n", v->d.id);
            std::exit(1);
        }
        std::string svtype_string = svtype;
        free(svtype);
        if (svtype_string == "INS") {
            ins _ins(h_v, v, fp_b, h_b, idx_b);
            if (_ins.consensus.begin() != _ins.consensus.end()) {
                std::cout << ">" << _ins.id << std::endl;
                std::cout << _ins.consensus << std::endl; 
            }
        }
    }

    sam_close(fp_b);
    bcf_destroy(v);
    bcf_hdr_destroy(h_v);
    vcf_close(fp_v);
}
