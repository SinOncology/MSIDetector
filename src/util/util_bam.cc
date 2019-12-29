#include <string>
#include <sstream>
#include <sam.h>
#include "util_bam.h"
#include "../installdir.h"
#include "../nibtools.h"


/**
  * @brief Compute a text-based sequence string from a BAM structure
  * @param b The BAM structure pointer
  *
  * @returns A string containing the sequence
  */
std::string get_seq_str_bam(bam1_t *b)
{
  bam1_core_t *c = &b->core;
  std::stringstream line;
  
  uint8_t *ss = bam1_seq(b); //pointer to the sequence
  //Get the sequence
  line.str("");
  line.clear();
  for (int j = 0; j < c->l_qseq; ++j)
  {
    line << bam_nt16_rev_table[bam1_seqi(ss, j)];
  }
  
  return line.str();
}


/**
 * @brief Change the CIGAR string for soft-clipped reads because partially overlapped with its mate (insert size < 2 * read length)
 *
 *        When soft-clipping occurs, update the CIGAR string of the first mate.
 *        Required to write the updated CIGAR in BAM format.
 *        Use memmove and memcpy for doing that
 *
 * @param b a pointer to an alignment structure
 * @param n the number of CIGAR operations
 * @param cigar a pointer to the cigar data
 */
void replace_cigar(bam1_t *b, int n, uint32_t *cigar)
{
  if (n != b->core.n_cigar)
  {
    int o = b->core.l_qname + b->core.n_cigar * 4;
    if (b->data_len + (n - b->core.n_cigar) * 4 > b->m_data)
    {
      b->m_data = b->data_len + (n - b->core.n_cigar) * 4;
      kroundup32(b->m_data);
      b->data = (uint8_t *) realloc(b->data, b->m_data);
    }
    memmove(b->data + b->core.l_qname + n * 4, b->data + o, b->data_len - o);
    memcpy(b->data + b->core.l_qname, cigar, n * 4);
    b->data_len += (n - b->core.n_cigar) * 4;
    b->core.n_cigar = n;
  }
  
  else
  {
    memcpy(b->data + b->core.l_qname, cigar, n * 4);
  }
}


/**
 * @brief get total reads (secondary alignments are included)
 * @param bam_file_name input bam file
 * @return total reads
 */
uint64_t get_total_reads(std::string input_bam)
{
  uint64_t total_reads = 0;
  samfile_t *fp = NULL;
  /// read the bam file
  fp = samopen(input_bam.c_str(), "rb", 0);
  if (fp == NULL)
  {
    std::cerr << "Error: can not open bam-file: " << input_bam << std::endl;
    exit(1);
  }
  bam_index_t *idx_bam;
  //hts_idx_t* idx_bam;
  /// load the bam index .bai file
  idx_bam = bam_index_load(input_bam.c_str());
  if (idx_bam == NULL)
  {
    std::cerr << "Error: please index bam file first:\t" << input_bam << std::endl;
    exit(1);
  }
  bam_header_t *header = fp->header;
  for (int i = 0; i < header->n_targets; ++i)
  {
    /// Now fetch info about it from the meta bin
    uint64_t u, v;
    hts_idx_get_stat(idx_bam, i, &u, &v);
    total_reads += u;
    total_reads += v;
  }
  /// information about unmapped reads
  total_reads += hts_idx_get_n_no_coor(idx_bam);
  hts_idx_destroy(idx_bam);
  samclose(fp);
  return total_reads;
}


/**
 * @brief given a 0-based position in a chromosome, return a 0-based genome-wide position
 * @param header bam file header
 * @param chromID 0-based chromsome ID
 * @param position 0-based position
 * @return 0-based genome-wide position
 */
uint32_t combine_genome_chr_pos(bam_header_t *header, int chromID, int32_t position)
{
  /// chr_pos is 0-based
  uint32_t chr_pos = 0;
  for (int i = 0; i < chromID; ++i)
  {
    chr_pos += header->target_len[i];
  }
  /// position is 0-based, chr_pos is also 0-based
  chr_pos += position;
  return chr_pos;
}


/**
 * @brief given a position, get the bases (length) in the right neighbor (do not include the base in the given position)
 * @param chrom chrom name
 * @param pos_1based 1-based reference position, for example, from vcf format
 * @param length, the length of the wanted bases
 * @return
 */
std::string get_right_neighbor_sequence_nib(std::string chrom, int32_t pos_1based, int length)
{
  std::string seq = "";
  //nib access
  nib nibObj;
  std::string build = "hg19";
  
  //open genome access (nib)
  std::string fn = (std::string) nib_folder + "/" + build + "_" + chrom + ".nib";
  int stat = nibObj.open(fn);
  char base;
  for (int32_t i = pos_1based; i < pos_1based + length; i++)
  {
    stat = nibObj.getBase(&base, i);
    seq += base;
  }
  nibObj.close();
  return seq;
}

/**
 * @brief given a positition, get the bases (length) in the left neighbor (do not include the base in the given position)
 * @param chrom chrom name
 * @param pos_1based 1-based reference position, for example, from vcf format
 * @param length, the length of the wanted bases
 * @return
 */
std::string get_left_neighbor_sequence_nib(std::string chrom, int32_t pos_1based, int length)
{
  std::string seq = "";
  //nib access
  nib nibObj;
  std::string build = "hg19";
  //open genome access (nib)
  std::string fn = (std::string) nib_folder + "/" + build + "_" + chrom + ".nib";
  int stat = nibObj.open(fn);
  char base;
  for (int32_t i = pos_1based - length; i < pos_1based; i++)
  {
    stat = nibObj.getBase(&base, i - 1);
    seq += base;
  }
  nibObj.close();
  return seq;
}


/**
 * @brief soft clip overlap part of an alignment of pair
 * @param b
 * @param b_copy
 */
bool soft_clip_overlap(bam1_t *b, bam1_t *b_copy)
{
  bool is_overlap = false;
  bam1_core_t *c = &b->core;
  /// cleaning the overlap
  /// the insert size lies within a certain interval and both mates map against the same chromosome
  /// exclude secondary alignment
  if (c->n_cigar > 0 && c->isize > 0 && c->isize <= 3 * c->l_qseq && !(c->flag & BAM_FUNMAP) &&
      !(c->flag & BAM_FMUNMAP) && c->mtid == c->tid && !(c->flag & BAM_FSECONDARY))
  {
    
    /// Combination of the two vectors 'cigar_pos' and 'cigar_type' encodes the CIGAR string of a read
    /// 'cigar_pos': Stores number of bases having the same CIGAR letter
    /// 'cigar_type': Letters of CIGAR string
    std::vector<int> cigar_pos;
    std::vector<char> cigar_type;
    int s_cigar, s_remain;
    int offset = 0;
    std::string orig_cigar;
    std::string new_cigar;
    int cigar_pos_tmp;
    
    
    uint32_t *cigar = bam1_cigar(b);
    
    std::stringstream line;
    
    cigar_pos.clear();
    cigar_type.clear();
    if (c->n_cigar == 0)
    {
      orig_cigar = "*";
    }
    else
    {
      line.str("");
      line.clear();
      /// Go through all elements of the CIGAR string
      for (int i = 0; i < c->n_cigar; ++i)
      {
        /// Collect length of constant CIGAR letter
        cigar_pos_tmp = bam1_cigar(b)[i] >> BAM_CIGAR_SHIFT;
        line << cigar_pos_tmp;
        cigar_pos.push_back(cigar_pos_tmp);
        /// Collect current cigar letter
        line << BAM_CIGAR_STR[bam1_cigar(b)[i] & BAM_CIGAR_MASK];
        cigar_type.push_back(BAM_CIGAR_STR[bam1_cigar(b)[i] & BAM_CIGAR_MASK]);
      }
      orig_cigar = line.str();
      //orig_cigar = get_cigar_str(cigar,c->n_cigar);
    }
    new_cigar = orig_cigar;
    
    /// a new cigar to store soft clip bases if overlap found
    uint32_t *cigar_new;
    cigar_new = bam1_cigar(b);
    /// to keep track the length of cigar operators
    int n_cigar_new = 0;
    std::string change_cigar;
    
    int i = 0;
    
    
    if (cigar_type[0] == 'S')
      offset = cigar_pos[0];
    else
      offset = 0;
    
    int j = 0;
    s_cigar = c->pos - offset + 1;
    for (int i = 0; i < c->n_cigar; ++i)
    {
      if (cigar_type[i] != 'D')
        s_cigar += cigar_pos[i];
      
      if (s_cigar < c->mpos + 2 || cigar_type[i] != 'M')
        j++;
      else
        break;
    }
    
    /// if R2 overlap with R1
    if (j < c->n_cigar && c->mpos - c->pos + offset < c->l_qseq)
    {
      b_copy = bam_copy1(b_copy, b);
      /// determining new cigar string
      line.str("");
      line.clear();
      /// keep the left cigar string
      for (int i = 0; i < j; ++i)
      {
        line << cigar_pos[i] << cigar_type[i];
        int op = bam_cigar_op(cigar[i]);
        int ol = bam_cigar_oplen(cigar[i]);
        cigar_new[i] = bam_cigar_gen(ol, op);
        n_cigar_new = n_cigar_new + 1;
      }
      
      if (line.str() != orig_cigar)
      {
        
        //cout<<orig_cigar<<endl;
        s_remain = 0;
        for (i = j; i < c->n_cigar; ++i)
        {
          if (cigar_type[i] != 'D')
            s_remain += cigar_pos[i];
        }
        
        int d = c->l_qseq - (c->mpos - c->pos + offset);
        if (d <= 0 || s_remain < d)
        {
          /// keep the cigar
          new_cigar = orig_cigar;
        }
        else
        {
          if (s_remain == d)
          {
            line << s_remain << "S";
            cigar_new[j] = bam_cigar_gen(s_remain, BAM_CSOFT_CLIP);
            n_cigar_new = n_cigar_new + 1;
          }
          else if (cigar_pos[j] >= s_remain - d)
          {
            line << s_remain - d << "M" << d << "S";
            /// modify the cigar array value
            cigar_new[j] = bam_cigar_gen(s_remain - d, BAM_CMATCH);
            cigar_new[j + 1] = bam_cigar_gen(d, BAM_CSOFT_CLIP);
            n_cigar_new = n_cigar_new + 2;
          }
          else
          {
            line << cigar_pos[j] << "M" << s_remain - cigar_pos[j] << "S";
            cigar_new[j] = bam_cigar_gen(cigar_pos[j], BAM_CMATCH);
            cigar_new[j + 1] = bam_cigar_gen(s_remain - cigar_pos[j], BAM_CSOFT_CLIP);
            n_cigar_new = n_cigar_new + 2;
          }
          new_cigar = line.str();
          if (orig_cigar != new_cigar)
          {
            
            /// replace old cigar with new cigar
            /// be aware of the length of cigar operaters also change to n_cigare_new
            replace_cigar(b_copy, n_cigar_new, cigar_new);
            bam_aux_append(b_copy, "OC", 'Z', strlen(orig_cigar.c_str()) + 1, (uint8_t *) orig_cigar.c_str());
            /// update the bin value
            b_copy->core.bin = bam_reg2bin(b_copy->core.pos, bam_endpos(b_copy));
            is_overlap = true;
          }
          
        }
      }
    }
  }
  return is_overlap;
  
  
}


// refID is 0-based in bam file
std::string chromID2ChrName(int refID)
{
  std::string chromName = "";
  std::ostringstream os;
  if (refID == 23)
    chromName = "chrY";
  else if (refID == 22)
    chromName = "chrX";
  else if (refID >= 0 && refID < 22)
  {
    os << "chr" << refID + 1;
    chromName = os.str();
  }
  return (chromName);
}


