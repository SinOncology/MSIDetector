#include <string>
#include <stdexcept>
#include <map>
#include "iostream"
#include "util.h"
#include <algorithm>
#include <regex>
#include <set>


using namespace std;

int chrom2number(std::string chr)
{
  int i;
  if (chr == "chrX")
    i = 23;
  else if (chr == "chrY")
    i = 24;
  else
  {
    try
    {
      i = std::stoi((chr.substr(chr.find("chr") + 3)).c_str());
    }
    catch (exception e)
    {
      //std::cerr << "no valid chrom numbre " << chr << std::endl;
      return -1;
    }
  }
  return (i - 1);
}

/**
 * @brief Convert from digit to string
 *
 * With the following table:
 * 0 => A,
 * 1 => C,
 * 2 => G,
 * 3 => T,
 * all other entries => N
 *
 * @param base integer
 *
 * @return the mapped string
 */
std::string num2base(int base)
{
  if (base == 0)
    return ("A");
  else if (base == 1)
    return ("C");
  else if (base == 2)
    return ("G");
  else if (base == 3)
    return ("T");
  else
    return ("N");
}

/**
 * @brief Convert from base character to digit
 * 
 * With the following table:
 * A => 0,
 * C => 1,
 * G => 2,
 * T => 3,
 * N => 4
 * other => 5
 * 
 * @param base base character
 * 
 * @return the mapped number
 */
long base2num(char base)
{
  if (toupper(base) == 'A')
    return (0);
  else if (toupper(base) == 'C')
    return (1);
  else if (toupper(base) == 'G')
    return (2);
  else if (toupper(base) == 'T')
    return (3);
  else
    return (4);
}

/**
 * @brief Convert from base character to digit
 * 
 * With the following table:
 * A => 0,
 * C => 1,
 * G => 2,
 * T => 3,
 * N => 4
 * other => 5
 * 
 * @param base base character
 * 
 * @return the mapped number
 */

long dsbase2num(char base)
{
  if (toupper(base) == 'A')
    return (0);
  else if (toupper(base) == 'C')
    return (1);
  else if (toupper(base) == 'G')
    return (2);
  else if (toupper(base) == 'T')
    return (3);
  else if (toupper(base) == 'N')
    return (4);
  else
    return (5);
}

//note that due to the number space of long int 
//index captures a resulution of 10bp
long mk_index_ds(std::string chr, long pos)
{
  return (chrom2number(chr) * 50000000 + pos / 10);
}

// FIXME both functions are identical...
long mk_index(std::string chr, long pos)
{
  return (chrom2number(chr) * 50000000 + pos / 10);
}


//in silico translation
std::string codon2protein(std::string codon)
{
  std::string amino;
  switch (codon[0])
  {
  case 'T':
    switch (codon[1])
    {
    case 'T':
      switch (codon[2])
      {
      case 'T':
        amino = "F";
        break;
      case 'C':
        amino = "F";
        break;
      case 'A':
        amino = "L";
        break;
      case 'G':
        amino = "L";
      }
      break;
    case 'C':
      amino = "S";
      break;
    case 'A':
      switch (codon[2])
      {
      case 'T':
        amino = "Y";
        break;
      case 'C':
        amino = "Y";
        break;
      case 'A':
        amino = "*";
        break; //stop codon
      case 'G':
        amino = "*"; //stop codon
      }
      break;
    case 'G':
      switch (codon[2])
      {
      case 'T':
        amino = "C";
        break;
      case 'C':
        amino = "C";
        break;
      case 'G':
        amino = "W";
        break;
      case 'A':
        amino = "*"; //stop codon
      }
    }
    break;
  case 'C':
    switch (codon[1])
    {
    case 'T':
      amino = "L";
      break;
    case 'C':
      amino = "P";
      break;
    case 'A':
      switch (codon[2])
      {
      case 'T':
        amino = "H";
        break;
      case 'C':
        amino = "H";
        break;
      case 'A':
        amino = "Q";
        break;
      case 'G':
        amino = "Q";
      }
      break;
    case 'G':
      amino = "R";
      break;
    }
    break;
  case 'A':
    switch (codon[1])
    {
    case 'T':
      switch (codon[2])
      {
      case 'G':
        amino = "M";
        break;
      default:
        amino = "I";
      }
      break;
    case 'C':
      amino = "T";
      break;
    case 'A':
      switch (codon[2])
      {
      case 'T':
        amino = "N";
        break;
      case 'C':
        amino = "N";
        break;
      case 'A':
        amino = "K";
        break;
      case 'G':
        amino = "K";
      }
      break;
    case 'G':
      switch (codon[2])
      {
      case 'T':
        amino = "S";
        break;
      case 'C':
        amino = "S";
        break;
      case 'A':
        amino = "R";
        break;
      case 'G':
        amino = "R";
      }
    }
    break;
  case 'G':
    switch (codon[1])
    {
    case 'T':
      amino = "V";
      break;
    case 'C':
      amino = "A";
      break;
    case 'A':
      switch (codon[2])
      {
      case 'T':
        amino = "D";
        break;
      case 'C':
        amino = "D";
        break;
      case 'A':
        amino = "E";
        break;
      case 'G':
        amino = "E";
      }
      break;
    case 'G':
      amino = "G";
      break;
    }
    break;
  default:
    amino = "";
  }
  return (amino);
}

/**
 * @brief combine chromsome ID and position into a 64 bit int
 * @param chromID
 * @param position
 * @return
 */

uint64_t combineChromPos(int32_t chromID, int32_t position)
{
  return (((uint64_t) chromID << 32) | (position & 0xFFFFFFFF));
}

void read_hotspot_fusion_file(const std::string hotspot_file, map<string, int> &hotspot_fusions)
{
  hotspot_fusions.clear();
  ifstream infile;
  infile.open(hotspot_file.c_str());
  if (!infile.is_open())
  {
    std::cerr << "Error: cannot open \t" << hotspot_file << std::endl;
    exit(1);
  }
  std::string line;
  std::string record1, record2, record3, record4;
  std::stringstream linestream;
  getline(infile, line, '\n');
  
  while (getline(infile, line, '\n'))
  {
    linestream.str("");
    linestream.clear();
    linestream << line;
    linestream >> record1;
    linestream >> record2;
    linestream >> record3;
    linestream >> record4;
    hotspot_fusions[record1 + ":" + record2 + "-" + record3 + ":" + record4] = 1;
  }
  infile.close();
}


void read_cosmic_fusion_file(const string cosmic_file,
                             map<string, string> &cosmic_fusions)
{
  cosmic_fusions.clear();
  ifstream infile;
  infile.open(cosmic_file.c_str());
  if (!infile.is_open())
  {
    std::cerr << "Error: cannot open \t" << cosmic_file << std::endl;
    exit(1);
  }
  std::string line;
  std::string record1, record2, record3, record4, tmp, cosmic_id;
  std::stringstream linestream;
  getline(infile, line, '\n');
  
  while (getline(infile, line, '\n'))
  {
    linestream.str("");
    linestream.clear();
    linestream << line;
    linestream >> cosmic_id;
    linestream >> record1;
    linestream >> record2;
    linestream >> record3;
    linestream >> record4;
    cosmic_fusions[record1 + ":" + record2 + "-" + record3 + ":" + record4] = cosmic_id;
  }
  infile.close();
  
}

void read_fusion_anno_file(const string &anno_file, vector<map<string, string>> &anno_fusion_info)
{
  map<string, string> fusion_5_3_map, fusion_5_map, fusion_3_map, fusion_no_map;
  map<string, string>::iterator fusion_map_it;
  fusion_5_3_map.clear();
  fusion_5_map.clear();
  fusion_3_map.clear();
  fusion_no_map.clear();
  ifstream infile;
  infile.open(anno_file.c_str());
  if (!infile.is_open())
  {
    std::cerr << "Error: cannot open \t" << anno_file << std::endl;
    exit(1);
  }
  std::string line;
  std::string gene1, exon1, gene2, exon2, id;
  std::stringstream linestream;
  getline(infile, line, '\n');
  
  while (getline(infile, line, '\n'))
  {
    linestream.str("");
    linestream.clear();
    linestream << line;
    linestream >> id;
    linestream >> gene1;
    linestream >> exon1;
    linestream >> gene2;
    linestream >> exon2;
    if (exon1 == ".")
    {
      if (exon2 == ".")
      {
        fusion_no_map[id] = gene1 + "-" + gene2;
      }
      else
      {
        fusion_3_map[id] = gene1 + "-" + gene2 + ":" + exon2;
      }
    }
    else
    {
      if (exon2 == ".")
      {
        fusion_5_map[id] = gene1 + ":" + exon1 + "-" + gene2;
      }
      else
      {
        fusion_5_3_map[id] = gene1 + ":" + exon1 + "-" + gene2 + ":" + exon2;
      }
    }
  }
  infile.close();
  anno_fusion_info.clear();
  anno_fusion_info.push_back(fusion_5_3_map);
  anno_fusion_info.push_back(fusion_5_map);
  anno_fusion_info.push_back(fusion_3_map);
  anno_fusion_info.push_back(fusion_no_map);
  
}

string number2chrom(int num)
{
  // 0 is chr1, 21 is chr22, and so on
  if (num == 22) return "chrX";
  if (num == 23) return "chrY";
  return "chr" + std::to_string(num + 1);
}

