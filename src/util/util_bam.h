#pragma once

#include <string>
#include <sam.h>
#include <vector>
#include <htslib/faidx.h>
#include <set>
#include <string>
#include <sstream>

using namespace std;
std::string get_seq_str_bam(bam1_t *b);
void replace_cigar(bam1_t *b, int n, uint32_t *cigar);
uint64_t get_total_reads(std::string bam_file);
uint32_t combine_genome_chr_pos(bam_header_t *header, int chromID, int32_t position);
std::string get_right_neighbor_sequence_nib(std::string chrom, int32_t pos_1based, int length);
std::string get_left_neighbor_sequence_nib(std::string chrom, int32_t pos_1based, int length);

/// soft clip overlap part of an alignment of pair
bool soft_clip_overlap(bam1_t *b, bam1_t *b_copy);
std::string chromID2ChrName(int refID);