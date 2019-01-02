#include "ins.h"

ins::ins(bcf_hdr_t *vcf_h, bcf1_t *vcf_v, samFile *sam_fp,
    bam_hdr_t *bam_h, hts_idx_t *bam_idx)
{
    chrom = bcf_hdr_id2name(vcf_h, vcf_v->rid);
    pos = vcf_v->pos;

    int n = 0;
    int32_t *tmp_i = NULL;
    bcf_get_info_int32(vcf_h, vcf_v, "SVLEN", &tmp_i, &n);
    size = *tmp_i;
    free(tmp_i);

    char *RNAMES = NULL;
    bcf_get_info_string(vcf_h, vcf_v, "RNAMES", &RNAMES, &n);
    std::string tmp_s = RNAMES;
    Tokenize(tmp_s, qnames, ',');
    free(RNAMES);

    region.chrom = chrom;
    if (pos > 1000) {
        region.start = pos - 1000;
    } else {
        region.start = 1;
    }

    bcf_idpair_t *ctg;
    ctg = vcf_h->id[BCF_DT_CTG];
    uint32_t ctg_len = ctg->val->info[0];
    if ( (pos + 1000) <= ctg_len) {
        region.end = pos + 1000;
    } else {
        region.end = ctg_len;
    }

    id = vcf_v->d.id;

    ins_sequences = get_ins_sequences(sam_fp, bam_h, bam_idx);
    consensus = get_consensus();
}

std::shared_ptr<seq> ins::get_ins_sequence(bam_hdr_t *h, bam1_t *b)
{
    int min_sv_size = 30;
    // std::cout << ":" << bam_get_qname(b) << std::endl;

    int32_t current_ref_pos = b->core.pos;
    uint32_t current_query_pos = 0;

    std::string query_seq;
    uint8_t *seq_prt = bam_get_seq(b);
    for (int32_t i = 0; i < b->core.l_qseq; ++i) {
        query_seq.push_back(seq_nt16_str[bam_seqi(seq_prt, i)]);
    }

    std::vector<uint32_t> ref_pos_vec;
    std::vector<uint32_t> query_pos_vec;
    std::vector<std::string> ins_seq_vec;

    uint32_t *ca_prt = bam_get_cigar(b);
    uint32_t ca_len = b->core.n_cigar;

    uint8_t c_op = 0;
    int32_t c_oplen = 0;
    int c_type = 0;
    for (uint32_t i = 0; i < ca_len; ++i) {
        c_op = bam_cigar_op(*(ca_prt+i));
        c_oplen = bam_cigar_oplen(*(ca_prt+i));
        c_type = bam_cigar_type(c_op);
        
        // bit 2 means the cigar operation consumes the reference
        if (c_type&2) {
            current_ref_pos += c_oplen;
        }

        // consumes query
        if(c_type&1) {
            current_query_pos += c_oplen;
        }

        if (c_op == BAM_CINS && c_oplen >= min_sv_size &&
                current_ref_pos >= region.start &&
                current_ref_pos <= region.end && c_oplen/float(size) >= 0.75 &&
                c_oplen/float(size) <= 1.33)
        {
            uint32_t ins_id = current_query_pos - c_oplen;
            std::string ins_seq = query_seq.substr(current_query_pos-c_oplen,
                c_oplen);
            ref_pos_vec.push_back(current_ref_pos);
            query_pos_vec.push_back(ins_id);
            ins_seq_vec.push_back(ins_seq);
        }
    }

    // select best match
    auto n = ref_pos_vec.size();
    if ( n == 0) {
        return nullptr;
    }
    std::size_t closest_idx = 0;
    if (n > 1) { 
        int min_dist = std::abs(ref_pos_vec[0] - pos);
        int dist = std::abs(ref_pos_vec[0] - pos);;
        for (std::size_t i = 1; i < n; ++i) {
            dist = std::abs(ref_pos_vec[i] - pos);
            if (dist < min_dist) {
                min_dist = dist;
                closest_idx = i;
            }
        }
    }
    std::string seq_name = id + "/" + std::string(bam_get_qname(b)) + "/"
        + std::string(h->target_name[b->core.tid]) + "/"
        + std::to_string(ref_pos_vec[closest_idx]) + "/"
        + std::to_string(query_pos_vec[closest_idx]);
    std::string sequence = ins_seq_vec[closest_idx];
    return std::shared_ptr<seq> (new seq(seq_name, sequence));
}

std::vector<std::shared_ptr<seq>> ins::get_ins_sequences(samFile *fp,
    bam_hdr_t *h, const hts_idx_t *bam_idx)
{
    bam1_t *b = bam_init1();
    int tid = bam_name2id(h, region.chrom.c_str());
    hts_itr_t *iter = sam_itr_queryi(bam_idx, tid, region.start, region.end);

    std::vector<std::shared_ptr<seq>> seqs;

    int ret;
    while ((ret = sam_itr_next(fp, iter, b)) >= 0) {
        std::shared_ptr<seq> ins_seq = get_ins_sequence(h, b);
        if (ins_seq != nullptr) {
            seqs.push_back(ins_seq);
        }
    }
    
    sam_itr_destroy(iter);
    bam_destroy1(b);
    if (seqs.size() < 1) {
        fprintf(stderr, 
            "[sv_ins_seq::Warning]: Can not find proper insert sequence for"
            " INS: %s\n", id.c_str());
    }
    return seqs;
}

std::string ins::get_consensus() {
    if (ins_sequences.size() > 1) {
        auto alignment_engine = spoa::createAlignmentEngine(
            static_cast<spoa::AlignmentType>(0), 5, -4, -8);

        auto graph = spoa::createGraph();

        for (const auto& it: ins_sequences) {
            auto alignment = alignment_engine->align_sequence_with_graph(
                it->sequence, graph);
            graph->add_alignment(alignment, it->sequence);
        }

        std::string consensus = graph->generate_consensus();

        // std::vector<std::string> msa;
        // graph->generate_multiple_sequence_alignment(msa);

        return consensus;
    } else if (ins_sequences.size() == 1) {
        return ins_sequences[0]->sequence;
    } else {
        return "";
    }
}