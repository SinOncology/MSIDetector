#pragma once

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <getopt.h>
#include <math.h>
#include <sstream>
#include <algorithm>
#include <htslib-1.3.1/htslib/sam.h>
#include <sam.h>

#include "BamAlignment.h"
#include "installdir.h"
#include "nibtools.h"
#include <set>
#include "InsertMole.h"
#include "util/util_string.h"
#include "util/util.h"

using namespace std;

string msi_help = "      SYNOPSIS \n \t MSIdetector <options>\n\n \
     DESCRIPTION\n \
     \t -h -? -help \t help\n \
     \t -i        \t tumor alignment file (sorted and indexed bam-file)\n \
     \t -b      \t whether to take baseline file as reference, if set means on,default is off [off]\n \
     \t -o      \t output name (prefix only, required)\n \
     \t -l      \t msi loci partition file absolute path( required)\n\
     \t -ref    \t baseline file (prefix only, required)\n\
     \t -plot   \t plot switch, if set means on,default is off [off]\n\
      \n";

typedef struct
{
  string chr;
  long range_start = -1;
  long range_end = -1;
  int32_t msi_exact_start; //左闭右闭的1Based区间
  int32_t msi_exact_end;  //左闭右闭的1Based区间
  string msi_marker_name;
  string rep_unit;
} msi_part_data;

typedef struct// single alignment support rep len info
{
  bool support = false;
  int rep_len = -1;
  bool exact = false;
  bool scope = false;
  char relation = 'n'; // enum{'n', 'e', 'g', 'l'} n: no relation, e:equal, g: greater than(>=), l: less than(<=)
  string overlapCigar = "";
  string overlapOC = "";
  bool overlap_bwa_soft = false;
  int32_t overlap_start = -1;
  int32_t overlap_end = -1;
  bool has_insert = false;
  string insert_base = "";
} msi_rep_len;

typedef struct
{
  string msi_marker = "";
  vector<int> repeatCountSet;
  vector<int> modelIndexSet;
  string modelIndexStr;
  double sd = -1.0;
  double mean = -1.0;
  vector<float> normalValueSet;
  string normalValueStr;
  double diffThreshold = 0.001;
  double value_thres = -1;
  double pvalue_thres = -1;
  string valid_direction = "";
  int baselineSampleSize = 0;
} loci_baseline;

typedef struct
{
  string msi_marker_name = "";
  vector<int> modelRepeatCounts;
  double modelMean = -1.0;
  double modelSd = -1.0;
  double sampleModelValue = -1.0;
  double sampleModelPvalue = -1.0;
  double sampleLogPvalue = -1.0;
  double modelPvalueThreshold = 0.0;
  double modelValueThreshold = 0.0;
  double logPvalueThreshold = 0.0;
  string loci_label = "";
} sample_loci_model_data;


void process_msi_partition_table(vector<msi_part_data> &msi_part, const string msi_file);

void
getRepeatLengthsForEachMsiLoc(vector<msi_rep_len> &msi_loci_rep_lengths, const msi_part_data &msi_part, samfile_t *fp,
                              bam_index_t *idx_bam, int readLen, int min_coverage, map<int, float> &quality_map);

void
computeRepeatLenPerAlignment(const msi_part_data &msi_part, InsertMole& insert_molecule, msi_rep_len& tmp_rep_len);

bool check_baseline_file(string base_file);
bool cmp_pos(BamAlignment &a1, BamAlignment &a2);

void writeDiffStats(string diff_stats_file, map<string, sample_loci_model_data> &loci_model_data_map,
                    vector<msi_part_data> &msi_parts);
void newWriteDiffStats(string diff_stats_file, map<string, sample_loci_model_data> &loci_model_data_map,
                    vector<msi_part_data> &msi_parts);

bool msiLociOverlap(InsertMole &mole, const msi_part_data &msi_part, msi_rep_len &tmp_rep_len);

int count_char(const string &test_string, char c);

void getExactRepLensFreqForEachMsiLoc(const vector<msi_rep_len> &tmp_repeat_lengths, const msi_part_data &msi_part,
                                      map<int, int> &tmp_exact_rep_lens_freq);

void process_baseline_file(string baseline_file, map<string, loci_baseline> &loci_baseline_map);

void computeDiffFromBaseline(map<string, loci_baseline> &loci_baseline_map,
                             map<string, map<int, float>> &msi_locs_freq_rate_map,
                             map<string, double> &loci_difference_map, string out_name,
                             map<string, string> &loci_label_map,
                             map<string, sample_loci_model_data> &loci_model_data_map);

void newProcessBaselineFile(string baseline_file, map<string, loci_baseline> &loci_baseline_map);

void newComputeDiffFromBaseline(map<string, loci_baseline> &loci_baseline_map,
                             map<string, map<int, float>> &msi_locs_freq_rate_map,
                             map<string, double> &loci_difference_map, string out_name,
                             map<string, string> &loci_label_map,
                             map<string, sample_loci_model_data> &loci_model_data_map);

void normalizeExactFrequencyForEachMsi(const msi_part_data &msi_part, map<int, int> &frequency_map,
                                       map<int, float> &frequency_rate_map);

void
writeLocsExactLensSupportReadsCountMatrix(map<string, vector<msi_rep_len>> &msi_locs_rep_lens_map,
                                          string total_reads_matrix_file);

void writeLocsExactLensFreqMatrix(map<string, map<int, int>>& msi_locs_exact_lens_freq, string out_freq_matrix);
void
writeSampleLocisFreqRateMatrix(map<string, map<int, float>> &msi_locs_exact_lens_freq_rate, string out_freq_rate_matrix,
                               string loci_plot_area_file, map<string, vector<int>> &loci_plot_area);

bool cmp_read_names(BamAlignment a1, BamAlignment a2);

int get_rep_count(CigarRoller &roller, const msi_part_data &msi_part, const string &molecule_seq,
                  int32_t align_start_pos);

void writeLocsSupportReadsMappingQualMatrix(map<string,map<int, float> >&loci_map_qual_map, string file_to_write);

