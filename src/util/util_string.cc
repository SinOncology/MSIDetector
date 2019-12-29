#include <string>
#include <vector>
#include <algorithm>
#include "util_string.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <set>


/**
 * @brief Performs an in-place sequence reversal
 * 
 * @param[in,out] seq the sequence to reverse
 */
void reverseSequence(std::string &seq)
{
  reverse(seq.begin(), seq.end());
}

/**
 *@brief Performs an in-place reverse complement conversion
 * 
 * @param[in,out] seq the sequence to reverse complement
 */
void reverseComplement(std::string &seq)
{
  
  // reverse the sequence
  reverseSequence(seq);
  
  // swap the bases
  for (unsigned int i = 0; i < seq.length(); i++)
  {
    switch (seq[i])
    {
    case 'A':
      seq[i] = 'T';
      break;
    case 'C':
      seq[i] = 'G';
      break;
    case 'G':
      seq[i] = 'C';
      break;
    case 'T':
      seq[i] = 'A';
      break;
    case 'a':
      seq[i] = 't';
      break;
    case 'c':
      seq[i] = 'g';
      break;
    case 'g':
      seq[i] = 'c';
      break;
    case 't':
      seq[i] = 'a';
      break;
    default:
      break;
    }
  }
}


std::vector<std::string> split_string(const std::string &s,
                                      const std::string &delim)
{
  bool keep_empty = false;
  std::vector<std::string> result;
  result.empty();
  if (delim.empty())
  {
    result.push_back(s);
    return result;
  }
  
  std::string::const_iterator substart = s.begin(), subend;
  while (true)
  {
    subend = search(substart, s.end(), delim.begin(), delim.end());
    std::string temp(substart, subend);
    if (keep_empty || !temp.empty())
    {
      result.push_back(temp);
    }
    if (subend == s.end())
    {
      break;
    }
    substart = subend + delim.size();
  }
  return result;
}

int find_longest_repeat_substring(const std::string &str)
{
  vector<repeat_str> repeat_vec;
  uint16_t len = str.length();
  uint16_t start_index = 0;
  repeat_str repeat_tmp;
  while (start_index < len)
  {
    uint16_t sub_len = 0;
    char begin_char = str.substr(start_index, 1).c_str()[0];
    bool same = true;
    while (same)
    {
      sub_len++;
      char current;
      current = str.substr(start_index + sub_len, 1).c_str()[0];
      same = (current == begin_char);
    }
    repeat_tmp.length = sub_len;
//    repeat_tmp.max_index = min_index+sub_len;
    repeat_tmp.start_index = start_index;
    repeat_tmp.rawstring = str;
    repeat_tmp.sub_str = str.substr(start_index, sub_len);
    repeat_vec.push_back(repeat_tmp);
    start_index = start_index + sub_len;
  }
  uint16_t max_len = 0;
  int max_repeat_index = 0;
  for (int i = 0; i < repeat_vec.size(); ++i)
  {
    if (repeat_vec[i].length > max_len)
    {
      max_len = repeat_vec[i].length;
      max_repeat_index = i;
    }
  }
  return repeat_vec[max_repeat_index].length;
}

//int find_substring() judge the string is composed by repeat
pair<int, string> strRepeatKmp(const string &str)
{
  bool repeat = false;
  vector<string> suffixArray;
  string s1, s2, s3, s4, s5, s6;
  
  int len = (int) str.size();
  for (int i = 0; i < len; ++i)
  {
    suffixArray.push_back(str.substr(i));
  }
  int count = 1;
  int maxCount = 1;
  string suffix;
  for (unsigned long i = 0; i < len; ++i)
  {
    for (unsigned long j = i + 1; j < len; ++j)
    {
      count = 1;
      s1 = suffixArray[i].substr(0, j - i);
      s2 = suffixArray[j].substr(0, j - i);
      s3 = suffixArray[j];
      if (suffixArray[i].substr(0, j - i) == suffixArray[j].substr(0, j - i))
      {
        ++count;
        for (unsigned long k = j + j - i; k < len; k += j - 1)
        {
          s4 = suffixArray[k];
          s5 = suffixArray[i].substr(0, j - i);
          s6 = suffixArray[k].substr(0, j - i);
          if (suffixArray[i].substr(0, j - i) == suffixArray[k].substr(0, j - i))
          {
            ++count;
          }
          else
          {
            break;
          }
        }
        if (count > maxCount)
        {
          maxCount = count;
          suffix = suffixArray[i].substr(0, j - i);
        }
      }
    }
  }
  return make_pair(maxCount, suffix);
  
}


/*
 * usage: judge the input string is fully composed by repeat unit
 */
pair<int, string> strRepeatEasy(const string &str)
{
  //violent window
  int repeatCount;
  bool repeat = false;
  std::string sub = "", current_sub;
  unsigned long len = str.size();
  pair<int, string> result_pair;
  for (int i = 1; i < str.size() / 2; ++i)
  {
    //i : window size
    if (len % i != 0) continue;//if window size is not exact division of string size : illegal
    sub = str.substr(0, i);
    repeatCount = 1;
    for (int j = 1; j <= str.size() / i; ++j)
    {
      //j: repeat unit count
      current_sub = str.substr(j * i, i);
      if (current_sub != sub)
      {
        break;
      }
      else
      {
        repeatCount++;
      }
    }
    if (repeatCount == str.size() / i)
    {
      break;
    }
    else
    {
      repeatCount = 0;
      sub = "";
    }
  }
  return make_pair(repeatCount, sub);
}

int countCharInStr(string str, char i)
{
  return 0;
}


/**判断str1是否以str2开头
 * 如果是返回1
 * 不是返回0
 * 出错返回-1
 * */
int is_begin_with(const char *str1, char *str2)
{
  if (str1 == NULL || str2 == NULL)
    return -1;
  size_t len1 = strlen(str1);
  size_t len2 = strlen(str2);
  if ((len1 < len2) || (len1 == 0 || len2 == 0))
    return -1;
  char *p = str2;
  int i = 0;
  while (*p != '\0')
  {
    if (*p != str1[i])
      return 0;
    p++;
    i++;
  }
  return 1;
}


/**
 * Brute Force search algorithms: to find tandem repeats
 * @param shorter_str
 * @param longer_str
 * @return
 */
bool has_common_rep_unit(const string shorter_str, const string longer_str)
{
  //cerr << "Now we compare str1: " << shorter_str << ", str2: " << longer_str << std::endl;
  bool result = false;
  int match = 0;
  std::set<string> shorter_str_rep_units;
  std::set<string> longer_str_rep_units;
  std::set<string>::iterator it_str;
  string rep_unit;
  shorter_str_rep_units.clear();
  longer_str_rep_units.clear();
  if (shorter_str.size() >= 2)
  {
    pair<int, string> short_str_rep_pair = strRepeatEasy(shorter_str);
    if (short_str_rep_pair.first != 0)
      rep_unit = short_str_rep_pair.second;
    else
      rep_unit = shorter_str;
  }
  else
    rep_unit = shorter_str;
  
  string regex_str;
  regex regex_single("(" + rep_unit + "){5,}[A-Za-z]*");
  regex regex_multi("(" + rep_unit + "){3,}[A-Za-z]*");
  
  
  if ((rep_unit.size() == 1 && regex_match(longer_str, regex_single)) ||
      (rep_unit.size() > 1 && regex_match(longer_str, regex_multi)))
  {
    result = true;
    //cout << "yes, they have the common unit:" << min_rep_unit << "\n";
  }
  return result;
}

void get_rep_units(std::set<string> &rep_unit, const string query_str)
{
  int match = 0;
  regex regex2("[A-Za-z]{1,}");
  bool find = false;
  for (int i = 1; i <= query_str.size() / 2; ++i)
  {
    for (int j = 0; j < query_str.size() - i; ++j)
    {
//      cout << "i: " << i << std::endl;
//      cout << "j: " << j << std::endl;
      if (query_str[j] == query_str[j + i])
      {
//        cout << "enter" << std::endl;
//        cout << query_str[j - 1] << std::endl;
//        cout << query_str[j + i - 1] << std::endl;
        match++;
        if (match >= i && regex_match(query_str.substr(j - i + 1, i), regex2))
        {
          rep_unit.insert(query_str.substr(j - i + 1, i));
          find = true;
//          cout << "repeat: " << query_str.substr(j - i + 1, i) << std::endl;
        }
      }
      else
        match = 0;
    }
    if (find)
      break;
  }
}

string condense_repeat_string(string expanded_string)
{
  string condense = "";
  int count = 1;
  char c_c;
  char c_n;
  char c_p = expanded_string[0];
  for (int i = 1; i < expanded_string.size();)
  {
    c_c = expanded_string[i];
    if (c_c == c_p)
    {
      c_p = c_c;
      count++;
      i++;
    }
    else
    {
      condense += to_string(count) + c_p;
      count = 1;
      c_p = c_c;
      i++;
    }
    
  }
  condense += to_string(count) + c_p;
  return condense;
}


int countRepeatUnit(const string &src_str, const string &rep_unit, const string part)
{
  int count = 0;
  if (part == "right")
  {
//    string expr = "\\w*(" + rep_unit + "{1,})$";
//    regex regex2(expr);
//    std::smatch sm;
//    std::regex_match(src_str, sm, regex2);
    for (int i = src_str.size() - rep_unit.size(); i >= 0;)
    {
      string sub = src_str.substr(i, rep_unit.size());
      if (sub == rep_unit)
      {
        count++;
        i -= rep_unit.size();
        continue;
      }
      else
      {
        break;
      }
    }
    return count;
  }
  
  if (part == "left")
  {
    count = 0;
//    string expr = "^(" + rep_unit + "{1,})\\w*$";
//    regex regex2(expr);
//    std::smatch sm;
//    std::regex_match(src_str, sm, regex2);
    for (int i = 0; i <= src_str.size() - rep_unit.size();)
    {
      string sub = src_str.substr(i, rep_unit.size());
      if (sub == rep_unit)
      {
        count++;
        i += rep_unit.size();
        continue;
      }
      else
      {
        break;
      }
    }
    return count;
  }
}