//
// Created by jinlf on 1/18/18.
//

#include "InsertMole.h"


InsertMole::InsertMole()
{
  left_read.reset();
  right_read.reset();
  single_read.reset();
  left_CigarRoller.Set("");
  right_CigarRoller.Set("");
  single_CigarRoller.Set("");
  sum_CigarRoller.Set("");
  
  left_start_pos = 0;
  right_start_pos = 0;
  single_start_pos = 0;
  left_end_pos = 0;
  right_end_pos = 0;
  single_end_pos = 0;
  
  left_seq = "";
  right_seq = "";
  single_seq = "";
  
  derive_from = 0;
  blank = true;
  
  left_mapping_qual = 0;
  right_mapping_qual = 0;
  single_mapping_qual = 0;
  
  ChrName = "";
  
  gap_start = 0;
  gap_end = 0;
  overlap_start = 0;
  overlap_end = 0;
  overlapped = false;
  gapped = false;
  readName = "";
  sum_read_seq = "";
  
}

void InsertMole::Set(BamAlignment &origin_left, BamAlignment &origin_right)
{
  readName = origin_right.getReadName();
  single_seq = "";
  derive_from = 2;
  single_CigarRoller.Set("");
  blank = false;
  ChrName = origin_right.getChrName();
  CigarRoller tmp;
  string sum_cigar = "";
  std::smatch sm;
  single_read.reset();
  
  if (origin_left.getAlignmentStart() > origin_right.getAlignmentStart())
  {///change read rank
    left_read = origin_right;
    right_read = origin_left;
    
    left_start_pos = origin_right.getAlignmentStart();
    right_start_pos = origin_left.getAlignmentStart();
    single_start_pos = 0;
    left_seq = origin_right.getReadSeq();
    right_seq = origin_left.getReadSeq();
    right_end_pos = origin_left.getAlignmentEnd();
    left_mapping_qual = origin_right.getMappingQual();
    right_mapping_qual = origin_left.getMappingQual();
    string oc = origin_right.originCigar();
    tmp.Set(origin_left.getCigarString().c_str());
    right_CigarRoller = tmp;
    
    if (oc.empty())
    {///oc is empty
      tmp.Set(origin_right.getCigarString().c_str());
      left_CigarRoller = tmp;
      left_end_pos = origin_right.getAlignmentEnd();
      ///judge the two ends is overlap or not
      if (left_end_pos >= right_start_pos)
      {
        overlapped = true;
        gapped = false;
        gap_start = 0;
        gap_end = 0;
        overlap_start = right_start_pos;
        overlap_end = left_end_pos;
        
        ///sum_cigar
        int left_cigar_str_len = left_CigarRoller.getExpandedCigarIndexFromRefPos(right_start_pos,
                                                                                  (int32_t) (left_start_pos));
        string left_part_expanded = left_CigarRoller.getExpandedString().substr(0, left_cigar_str_len);
        string left_part_condensed = condense_repeat_string(left_part_expanded);
        int right_start_cigar_index = right_CigarRoller.getExpandedCigarIndexFromRefPos(right_start_pos,
                                                                                        right_start_pos);
        int right_cigar_str_len = right_CigarRoller.getExpandedString().size() - right_start_cigar_index;
        string right_part_expanded = right_CigarRoller.getExpandedString().substr(right_start_cigar_index,
                                                                                  right_cigar_str_len);
        string right_part_condensed = condense_repeat_string(right_part_expanded);
        sum_cigar = link_cigar_string(left_part_condensed, right_part_condensed);
        tmp.Set(sum_cigar.c_str());
        sum_CigarRoller = tmp;
        ///sum seq
        int left_part_len = origin_right.getQueryIndex(right_start_pos - 1, left_start_pos);
        int right_part_start_index = origin_left.getQueryIndex(right_start_pos, right_start_pos);
        int right_part_len = right_seq.size() - right_part_start_index;
        sum_read_seq = left_seq.substr(0, left_part_len) +
                       right_seq.substr(right_part_start_index, right_part_len);
        
      }
      else
      {
        overlapped = false;
        gapped = true;
        gap_start = left_end_pos;
        gap_end = right_start_pos;
        overlap_start = 0;
        overlap_end = 0;
        sum_CigarRoller.Set("");
        sum_read_seq = "";
      }
    }
    else
    {/// oc is not blank
      overlapped = true;
      gapped = false;
      gap_start = 0;
      gap_end = 0;
      
      tmp.Set(oc.c_str());
      left_CigarRoller = tmp;
      left_end_pos = tmp.getAlignmentEnd((uint32_t) left_start_pos);
      tmp.Set(origin_left.getCigarString().c_str());
      right_end_pos = origin_left.getAlignmentEnd();
      right_CigarRoller = tmp;
      overlap_start = right_start_pos;
      overlap_end = left_end_pos;
      ///sum_cigar
      int left_cigar_str_len = left_CigarRoller.getExpandedCigarIndexFromRefPos(right_start_pos,
                                                                                (int32_t) (left_start_pos));
      string left_part_expanded = left_CigarRoller.getExpandedString().substr(0, left_cigar_str_len);
      string left_part_condensed = condense_repeat_string(left_part_expanded);
      int right_start_cigar_index = right_CigarRoller.getExpandedCigarIndexFromRefPos(right_start_pos, right_start_pos);
      int right_cigar_str_len = right_CigarRoller.getExpandedString().size() - right_start_cigar_index;
      string right_part_expanded = right_CigarRoller.getExpandedString().substr(right_start_cigar_index,
                                                                                right_cigar_str_len);
      string right_part_condensed = condense_repeat_string(right_part_expanded);
      sum_cigar = link_cigar_string(left_part_condensed, right_part_condensed);
      tmp.Set(sum_cigar.c_str());
      sum_CigarRoller = tmp;
      
      int left_part_len = origin_right.getQueryIndex(right_start_pos - 1, left_start_pos);
      int right_part_start_index = origin_left.getQueryIndex(right_start_pos, right_start_pos);
      int right_part_len = right_seq.size() - right_part_start_index;
      sum_read_seq = left_seq.substr(0, left_part_len) +
                     right_seq.substr(right_part_start_index, right_part_len);
      ///oc is not blank end
    }
  }
  else
  { ///don't change read rank
    left_read = origin_left;
    right_read = origin_right;
    
    left_start_pos = origin_left.getAlignmentStart();
    right_start_pos = origin_right.getAlignmentStart();
    single_start_pos = 0;
    left_seq = origin_left.getReadSeq();
    right_seq = origin_right.getReadSeq();
    
    right_end_pos = origin_right.getAlignmentEnd();
    
    left_mapping_qual = origin_left.getMappingQual();
    right_mapping_qual = origin_right.getMappingQual();
    string oc = origin_left.originCigar();
    tmp.Set(origin_right.getCigarString().c_str());
    right_CigarRoller = tmp;
    
    if (oc.empty())
    {///oc is empty
      tmp.Set(origin_left.getCigarString().c_str());
      left_CigarRoller = tmp;
      left_end_pos = origin_left.getAlignmentEnd();
      ///judge the two ends is overlap or not
      if (left_end_pos >= right_start_pos)
      {
        overlapped = true;
        gapped = false;
        gap_start = 0;
        gap_end = 0;
        overlap_start = right_start_pos;
        overlap_end = left_end_pos;
        
        ///sum_cigar
        int left_cigar_str_len = left_CigarRoller.getExpandedCigarIndexFromRefPos(right_start_pos,
                                                                                  (int32_t) (left_start_pos));
        string left_part_expanded = left_CigarRoller.getExpandedString().substr(0, left_cigar_str_len);
        string left_part_condensed = condense_repeat_string(left_part_expanded);
        int right_start_cigar_index = right_CigarRoller.getExpandedCigarIndexFromRefPos(right_start_pos,
                                                                                        right_start_pos);
        int right_cigar_str_len = right_CigarRoller.getExpandedString().size() - right_start_cigar_index;
        string right_part_expanded = right_CigarRoller.getExpandedString().substr(right_start_cigar_index,
                                                                                  right_cigar_str_len);
        string right_part_condensed = condense_repeat_string(right_part_expanded);
        sum_cigar = link_cigar_string(left_part_condensed, right_part_condensed);
        tmp.Set(sum_cigar.c_str());
        sum_CigarRoller = tmp;
        ///sum seq
        int left_part_len = origin_left.getQueryIndex(right_start_pos - 1, left_start_pos);
        int right_part_start_index = origin_right.getQueryIndex(right_start_pos, right_start_pos);
        int right_part_len = right_seq.size() - right_part_start_index;
        sum_read_seq = left_seq.substr(0, left_part_len) +
                       right_seq.substr(right_part_start_index, right_part_len);
        
      }
      else
      {
        overlapped = false;
        gapped = true;
        gap_start = left_end_pos;
        gap_end = right_start_pos;
        overlap_start = 0;
        overlap_end = 0;
        sum_CigarRoller.Set("");
        sum_read_seq = "";
      }
    }
    else
    {/// oc is not blank
      overlapped = true;
      gapped = false;
      gap_start = 0;
      gap_end = 0;
      tmp.Set(oc.c_str());
      left_CigarRoller = tmp;
      left_end_pos = tmp.getAlignmentEnd((uint32_t) left_start_pos);
      overlap_start = right_start_pos;
      overlap_end = left_end_pos;
      tmp.Set(origin_right.getCigarString().c_str());
      right_CigarRoller = tmp;
      ///sum_cigar
      int left_cigar_str_len = left_CigarRoller.getExpandedCigarIndexFromRefPos(right_start_pos,
                                                                                (int32_t) (left_start_pos));
      string left_part_expanded = left_CigarRoller.getExpandedString().substr(0, left_cigar_str_len);
      string left_part_condensed = condense_repeat_string(left_part_expanded);
      int right_start_cigar_index = right_CigarRoller.getExpandedCigarIndexFromRefPos(right_start_pos, right_start_pos);
      int right_cigar_str_len = right_CigarRoller.getExpandedString().size() - right_start_cigar_index;
      string right_part_expanded = right_CigarRoller.getExpandedString().substr(right_start_cigar_index,
                                                                                right_cigar_str_len);
      string right_part_condensed = condense_repeat_string(right_part_expanded);
      sum_cigar = link_cigar_string(left_part_condensed, right_part_condensed);
      tmp.Set(sum_cigar.c_str());
      sum_CigarRoller = tmp;
      
      int left_part_len = origin_left.getQueryIndex(right_start_pos - 1, left_start_pos);
      int right_part_start_index = origin_right.getQueryIndex(right_start_pos, right_start_pos);
      int right_part_len = right_seq.size() - right_part_start_index;
      sum_read_seq = left_seq.substr(0, left_part_len) +
                     right_seq.substr(right_part_start_index, right_part_len);
      ///oc is not blank end
    }
  }
}

void InsertMole::Set(BamAlignment &single_end)
{
  left_read.reset();
  right_read.reset();
  single_read = single_end;
  left_CigarRoller.Set("");
  right_CigarRoller.Set("");
  CigarRoller tmp;
  left_start_pos = 0;
  right_start_pos = 0;
  single_start_pos = single_end.getAlignmentStart();
  
  left_end_pos = 0;
  right_end_pos = 0;
  string oc = single_end.originCigar();
  if (!oc.empty())
  {
    tmp.Set(oc.c_str());
    single_end_pos = tmp.getAlignmentEnd(single_start_pos);
    single_CigarRoller = tmp;
  }
  else
  {
    tmp.Set(single_end.getCigarString().c_str());
    single_CigarRoller = tmp;
    single_end_pos = single_end.getAlignmentEnd();
  }
  sum_CigarRoller = single_CigarRoller;
  left_seq = "";
  right_seq = "";
  single_seq = single_end.getReadSeq();
  derive_from = 1;
  blank = false;
  left_mapping_qual = 0;
  right_mapping_qual = 0;
  single_mapping_qual = single_end.getMappingQual();
  ChrName = single_end.getChrName();
  gap_start = 0;
  gap_end = 0;
  overlap_start = 0;
  overlap_end = 0;
  overlapped = false;
  gapped = false;
  readName = single_end.getReadName();
  sum_read_seq = single_seq;
}

InsertMole::~InsertMole()
{
  left_read.reset();
  right_read.reset();
  single_read.reset();
  left_CigarRoller.Set("");
  right_CigarRoller.Set("");
  single_CigarRoller.Set("");
  sum_CigarRoller.Set("");
  
  left_start_pos = 0;
  right_start_pos = 0;
  single_start_pos = 0;
  left_end_pos = 0;
  right_end_pos = 0;
  single_end_pos = 0;
  
  left_seq = "";
  right_seq = "";
  single_seq = "";
  
  derive_from = 0;
  blank = true;
  
  left_mapping_qual = 0;
  right_mapping_qual = 0;
  single_mapping_qual = 0;
  
  ChrName = "";
  
  gap_start = 0;
  gap_end = 0;
  overlap_start = 0;
  overlap_end = 0;
  overlapped = false;
  gapped = false;
  sum_read_seq = "";
}

void InsertMole::reset()
{
  left_read.reset();
  right_read.reset();
  single_read.reset();
  left_CigarRoller.Set("");
  right_CigarRoller.Set("");
  single_CigarRoller.Set("");
  sum_CigarRoller.Set("");
  
  left_start_pos = 0;
  right_start_pos = 0;
  single_start_pos = 0;
  left_end_pos = 0;
  right_end_pos = 0;
  single_end_pos = 0;
  
  left_seq = "";
  right_seq = "";
  single_seq = "";
  
  derive_from = 0;
  blank = true;
  
  left_mapping_qual = 0;
  right_mapping_qual = 0;
  single_mapping_qual = 0;
  
  ChrName = "";
  
  gap_start = 0;
  gap_end = 0;
  overlap_start = 0;
  overlap_end = 0;
  overlapped = false;
  gapped = false;
  sum_read_seq = "";
}

string InsertMole::getLeftCigarString()
{
  if (derive_from == 2)
    return right_CigarRoller.getExpandedString();
  else
  {
    cerr << "invalid call, because the molecule does not contain 2 reads\n";
    exit(1);
  }
}

string InsertMole::getRightCigarString()
{
  if (derive_from == 2)
    return right_CigarRoller.getExpandedString();
  else
  {
    cerr << "invalid call, because the molecule does not contain 2 reads\n";
    exit(1);
  }
}

string InsertMole::getSingleCigarString()
{
  if (derive_from == 2)
    return single_CigarRoller.getExpandedString();
  else
  {
    cerr << "invalid call, because the molecule does not contain 2 reads\n";
    exit(1);
  }
}


int InsertMole::getSumMappingQual()
{
  if (!blank)
  {
    if (derive_from == 0)
      return 0;
    if (derive_from == 1)
      return single_mapping_qual;
    if (derive_from == 2)
      return min(left_mapping_qual, right_mapping_qual);
  }
  else
    return 0;
}

string InsertMole::getChrName()
{
  return ChrName;
}

long InsertMole::getAlignmentLeftStart()
{
  
  if (derive_from == 2)
    return left_start_pos;
  else
  {
    cerr << "invalid call, because the molecule does not contain 2 reads\n";
    exit(1);
  }
}


long InsertMole::getAlignmentLeftEnd()
{
  if (derive_from == 2)
    return left_end_pos;
  else
  {
    cerr << "invalid call, because the molecule does not contain 2 reads\n";
    exit(1);
  }
}

long InsertMole::getAlignmentRightStart()
{
  if (derive_from == 2)
    return right_start_pos;
  else
  {
    cerr << "invalid call, because the molecule does not contain 2 reads\n";
    exit(1);
  }
}

long InsertMole::getAlignmentRightEnd()
{
  if (derive_from == 2)
    return right_end_pos;
  else
  {
    cerr << "invalid call, because the molecule does not contain 2 reads\n";
    exit(1);
  }
}

long InsertMole::getAlignmentSingleStart()
{
  if (derive_from == 1)
    return single_start_pos;
  else
  {
    cerr << "invalid call, because the molecule does not contain exact 1 reads\n";
    exit(1);
  }
  
}

long InsertMole::getAlignmentSingleEnd()
{
  if (derive_from == 1)
    return single_end_pos;
  else
  {
    cerr << "invalid call, because the molecule does not contain exact 1 reads\n";
    exit(1);
  }
  
}

int InsertMole::get_align_record_num()
{
  return derive_from;
}

string InsertMole::getReadName()
{
  return readName;
}

bool InsertMole::isGapped()
{
  return gapped;
}

bool InsertMole::isOverlapped()
{
  return overlapped;
}

string InsertMole::getSumReadSeq()
{
  
  if (derive_from == 2)
    return sum_read_seq;
  else
  {
    cerr << "invalid call, because the molecule does not contain 2 reads\n";
    exit(1);
  }
}


void InsertMole::getSumCigar(CigarRoller &result)
{
  if (derive_from == 2)
    result = sum_CigarRoller;
  else
  {
    cerr << "invalid call, because the molecule does not contain 2 reads\n";
    exit(1);
  }
  
}

string InsertMole::getLeftReadSeq()
{
  if (derive_from == 2)
    return left_seq;
  else
  {
    cerr << "invalid call, because the molecule does not contain 2 reads\n";
    exit(1);
  }
}

string InsertMole::getRightReadSeq()
{
  if (derive_from == 2)
    return right_seq;
  else
  {
    cerr << "invalid call, because the molecule does not contain 2 reads\n";
    exit(1);
  }
}

string InsertMole::getSingleReadSeq()
{
  if (derive_from == 1)
    return single_seq;
  else
  {
    cerr << "invalid call, because the molecule does not contain exact 1 reads\n";
    exit(1);
  }
}

void InsertMole::getLeftCigarRoller(CigarRoller &read_cigarRoller)
{
  if (derive_from == 2)
    read_cigarRoller = left_CigarRoller;
  else
  {
    cerr << "invalid call, because the molecule does not contain 2 reads\n";
    exit(1);
  }
}

void InsertMole::getRightCigarRoller(CigarRoller &read_cigarRoller)
{
  if (derive_from == 2)
    read_cigarRoller = right_CigarRoller;
  else
  {
    cerr << "invalid call, because the molecule does not contain 2 reads\n";
    exit(1);
  }
  
}

void InsertMole::getSingleCigarRoller(CigarRoller &read_cigarRoller)
{
  if (derive_from == 1)
    read_cigarRoller = single_CigarRoller;
  else
  {
    cerr << "invalid call, because the molecule does not contain exact 1 reads\n";
    exit(1);
  }
  
}


