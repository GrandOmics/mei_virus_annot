#ifndef INS_H
#define INS_H

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <memory>

#include "htslib/vcf.h"
#include "htslib/sam.h"
#include "spoa/spoa.hpp"

#include "yutils.h"


struct region_t {
    std::string chrom;
    int32_t start;
    int32_t end;
};

class seq {
    public:
    seq(const std::string &n, const std::string &s): name(n), sequence(s){};
    std::string name;
    std::string sequence;
    // void print() {
    //     std::cout << ">" << name << std::endl;
    //     std::cout << sequence << std::endl;
    // }
};

class ins {
    public:
        ins() = default;
        ins(bcf_hdr_t *vcf_h, bcf1_t *vcf_v, samFile *sam_fp,
            bam_hdr_t *bam_h, hts_idx_t *bam_idx);

        std::string chrom;
        uint32_t pos;
        int32_t size;
        std::vector<std::string> qnames;
        std::string id;
        
        std::vector<std::shared_ptr<seq>> ins_sequences;
        std::string consensus;

        // void print() {
        //     std::cout << chrom << "\t" << pos+1 << "\t" << size << "\t" <<
        //         region.start << "\t" << region.end << "\t";
        //     for (const auto &s : qnames) {
        //         std::cout << s << "\t";
        //     }
        //     std::cout << std::endl;
        // }

    private:
        region_t region;
        std::shared_ptr<seq> get_ins_sequence(bam_hdr_t *h, bam1_t *b);
        std::vector<std::shared_ptr<seq>> get_ins_sequences(samFile *fp,
            bam_hdr_t *h, const hts_idx_t *bam_idx);
        std::string get_consensus();
};

#endif // INS_H
