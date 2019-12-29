//
// Created by jinlf on 1/18/18.
//
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
#include <regex>


class InsertMole
{
public:
  InsertMole();
  
  void Set(BamAlignment &origin_left, BamAlignment &origin_right);
  void Set(BamAlignment &single_end);
  
  InsertMole(BamAlignment &single_end);
  
  ~InsertMole();
  
  void reset();
  string getLeftCigarString();
  string getRightCigarString() ;
  string getSingleCigarString();
  
  int getSumMappingQual();
  
  string getChrName();
  
  long getAlignmentLeftStart();
  
  long getAlignmentLeftEnd();
  
  long getAlignmentRightStart();
  
  long getAlignmentRightEnd();
  
  long getAlignmentSingleStart();
  
  long getAlignmentSingleEnd();
  int get_align_record_num();
  string getReadName();
  bool isGapped();
  bool isOverlapped();
  string getSumReadSeq();

  void getSumCigar(CigarRoller& result);

  string getLeftReadSeq();
  string getRightReadSeq();
  string getSingleReadSeq();
  void getLeftCigarRoller(CigarRoller& left_read_cigarRoller);
  void getRightCigarRoller(CigarRoller& left_read_cigarRoller);
  void getSingleCigarRoller(CigarRoller &left_read_cigarRoller);

private:
  CigarRoller left_CigarRoller, right_CigarRoller, single_CigarRoller, sum_CigarRoller;
  long left_start_pos, left_end_pos, right_end_pos, right_start_pos,
    single_start_pos, single_end_pos;
  string left_seq, right_seq, single_seq;
  int derive_from;
  bool blank;
  float left_mapping_qual, right_mapping_qual, single_mapping_qual;
  string ChrName;
  long gap_start, gap_end;
  long overlap_start, overlap_end;
  bool overlapped, gapped;
  string readName;
  string sum_read_seq;
  BamAlignment left_read,right_read,single_read;
  
};


