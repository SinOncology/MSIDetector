#pragma once

#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>

using namespace std;
int chrom2number(std::string chr);

std::string num2base(int base);

long base2num(char base);

long dsbase2num(char base);

long mk_index_ds(std::string chr, long pos);

long mk_index(std::string chr, long pos);

std::string codon2protein(std::string codon);


uint64_t combineChromPos(int32_t chromID, int32_t position);


void read_fusion_anno_file(const string &anno_file, vector<map<string, string>> &anno_fusion_info);

string number2chrom(int num);


