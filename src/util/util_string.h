#pragma once

#include <string>
#include <iostream>
#include <stdexcept>
#include <map>
#include "iostream"
#include <algorithm>
#include <regex>
#include <set>

#include <regex>
using namespace std;

void reverseSequence(std::string &seq);

void reverseComplement(std::string &seq);

std::vector<std::string> split_string(const std::string &s,
                                      const std::string &delim);

int find_longest_repeat_substring(const std::string &s);

typedef struct
{
  std::string rawstring;
  uint16_t start_index;
  std::string sub_str;
  uint16_t length;
} repeat_str;

pair<int, string> strRepeatKmp(const string &str);

pair<int, string> strRepeatEasy(const string &str);

int countCharInStr(string str, char i);

int is_begin_with(const char *str1, char *str2);

bool has_common_rep_unit(const string shorter_str, const string longer_str);

void get_rep_units(std::set<string> &rep_units, const string query_str);
string condense_repeat_string(string expanded_string);

int countRepeatUnit(const string &src_str, const string &rep_unit, const string part);
