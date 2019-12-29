
#pragma once

#include "msi.h"


//main_program
int main(int argc, char **argv)
{
  int longindex, opt;
  // option definition
  
  static struct option longopts[] = {
    {"help", 0, 0, 'h'},
    {"help", 0, 0, '?'},
    {"i",    1, 0, 1},
    {"o",    1, 0, 2},
    {"b",    0, 0, 3},
    {"plot", 0, 0, 4},
    {"l",    1, 0, 5},
    {"ref",  1, 0, 6},
    {0,      0, 0, 0}
  };
  
  string bam_file = "";
  string out_name = "";
  string build = "hg19";
  ofstream out;
  ifstream in;
  long i, j;
  string tmp;
  bool baseline = false;
  string base_file = "";
  bool plot = false;
  string msi_file = "";
  optind = 0;
  //parse command line arguments
  while ((opt = getopt_long_only(argc, argv, "h?", longopts, &longindex)) != -1)
  {
    switch (opt)
    {
    case 'h':
      cerr << msi_help;
      break;
    case '?':
      cerr << msi_help;
      break;
    case 1:
      bam_file = (string) optarg;
      break;
    case 2:
      out_name = (string) optarg;
      break;
    case 3:
      baseline = true;
      break;
    case 4:
      plot = true;
      break;
    case 5:
      msi_file = (string) optarg;
      break;
    case 6:
      base_file = (string) optarg;
      break;
    default:
      cerr << "Error: cannot parse arguments.\n";
      exit(1);
    }
  }
  
  if ((bam_file == "" && msi_file == "") || out_name == "")
  {
    cerr << "Please specifiy input/output parameters\n";
    cerr << msi_help;
    cerr << endl;
    exit(1);
  }
 
  
  if (baseline && base_file == "")
  {
    cerr << "Please specifiy baseline file path  when you choose turn baseline on\n";
    
    cerr << msi_help;
    cerr << endl;
    exit(1);
  }
  if (base_file != "")
  {
    if (!check_baseline_file(base_file))
    {
      cerr << "your baseline file is illegal\n";
      
      cerr << msi_help;
      cerr << endl;
      exit(1);
    }
  }
  vector<msi_part_data> msi_parts;
  //读入msi bed 格式的文件  形成 msi_part_data的数组
  process_msi_partition_table(msi_parts, msi_file);
  
  if (msi_parts.size() == 0)
  {
    cerr << "Error: msi-part data empty\n";
    exit(1);
  }
  
  map<string, vector<msi_rep_len>> msi_locs_rep_lens_map; //
  map<string, vector<msi_rep_len>>::iterator msi_locs_rep_lens_iterator;
  vector<msi_rep_len> tmp_repeat_lengths;
  map<int, int> tmp_exact_rep_lens_freq;
  map<int, int>::iterator tmp_exact_rep_lens_freq_it;
  map<int, float> tmp_exact_rep_lens_freq_rate;
  map<int, float>::iterator tmp_exact_rep_lens_freq_rate_it;
  map<string, map<int, int >> msi_locs_exact_lens_freq;
  map<string, map<int, int >>::iterator msi_locs_exact_lens_freq_iterator;
  map<string, map<int, float >> msi_locs_exact_lens_freq_rate;
  map<string, map<int, float >>::iterator msi_locs_exact_lens_freq_rate_iterator;
  
  ///打开bam文件
  samfile_t *fp = samopen(bam_file.c_str(), "rb", 0);
  bam_index_t *idx_bam;
  if (fp == NULL)
  {
    std::cerr << "Error: can not open input bam-file:\t" << bam_file << std::endl;
    exit(1);
  }
  idx_bam = bam_index_load(bam_file.c_str());
  if (idx_bam == NULL)
  {
    std::cerr << "Error: please index bam-file first:\t" << bam_file << std::endl;
    exit(1);
  }
  int readLen = 150;
  int min_coverage = 50;
  map<string, map<int, float>> loci_map_qual_map;
  map<int, float> single_loci_map_qual_map;
  for (int k = 0; k < msi_parts.size(); ++k)
  {
    //cout << "Now it's " << msi_parts[k].msi_marker_name << "'s  turn.\n ";

//    if (msi_parts[k].msi_marker_name == "MONO−27")
//    {
    //遍历每一个msi loci main compute module 计算此loci的所有支持reads数 装入msi_rep_len的数组里
    tmp_repeat_lengths.clear();
    tmp_exact_rep_lens_freq.clear();
    single_loci_map_qual_map.clear();
    //TODO take multi-threading into consideration
    getRepeatLengthsForEachMsiLoc(tmp_repeat_lengths, msi_parts[k], fp, idx_bam, readLen, min_coverage,
                                  single_loci_map_qual_map);
    //
    //cout << "get the repeat length is completed.\n";
    
    if (tmp_repeat_lengths.size() != 0)
    {
      loci_map_qual_map[msi_parts[k].msi_marker_name] = single_loci_map_qual_map;
      msi_locs_rep_lens_iterator = msi_locs_rep_lens_map.find(msi_parts[k].msi_marker_name);
      if (msi_locs_rep_lens_iterator == msi_locs_rep_lens_map.end())
      {
        msi_locs_rep_lens_map[msi_parts[k].msi_marker_name] = tmp_repeat_lengths; // store the loci rep length data
      }
      
      getExactRepLensFreqForEachMsiLoc(tmp_repeat_lengths, msi_parts[k], tmp_exact_rep_lens_freq);
      
      ////map: 0~50
      for (int len = 0; len < 51; ++len)
      {
        tmp_exact_rep_lens_freq_it = tmp_exact_rep_lens_freq.find(len);
        if (tmp_exact_rep_lens_freq_it == tmp_exact_rep_lens_freq.end())
        {
          tmp_exact_rep_lens_freq[len] = 0;
        }
      }
      
      msi_locs_exact_lens_freq_iterator = msi_locs_exact_lens_freq.find(msi_parts[k].msi_marker_name);
      if (msi_locs_exact_lens_freq_iterator == msi_locs_exact_lens_freq.end())
      {
        msi_locs_exact_lens_freq[msi_parts[k].msi_marker_name] = tmp_exact_rep_lens_freq; // store the loci rep length data
      }
      
      //cout << "get the exact repeat length frequency  is completed.\n";
      normalizeExactFrequencyForEachMsi(msi_parts[k], tmp_exact_rep_lens_freq, tmp_exact_rep_lens_freq_rate);
      
      ////map: 0~50
      for (int len = 0; len < 51; ++len)
      {
        tmp_exact_rep_lens_freq_rate_it = tmp_exact_rep_lens_freq_rate.find(len);
        if (tmp_exact_rep_lens_freq_rate_it == tmp_exact_rep_lens_freq_rate.end())
        {
          tmp_exact_rep_lens_freq_rate[len] = 0.0;
        }
      }
      
      msi_locs_exact_lens_freq_rate_iterator = msi_locs_exact_lens_freq_rate.find(msi_parts[k].msi_marker_name);
      if (msi_locs_exact_lens_freq_rate_iterator == msi_locs_exact_lens_freq_rate.end())
      {
        msi_locs_exact_lens_freq_rate[msi_parts[k].msi_marker_name] = tmp_exact_rep_lens_freq_rate; // store the loci rep length data
      }
      
    }
    
    tmp_repeat_lengths.clear();

//    }
  }
  samclose(fp);
  //write sample distribution statistics file
  vector<string> msi_rep_lens_files;
  vector<string> msi_rep_exact_lens_freq_files;
  vector<string> msi_rep_exact_lens_freq_rate_files;
  vector<string> msi_rep_lens_plots;
  vector<string> msi_rep_lens_rate_plots;
  string msi_rep_lens_plots_file = out_name + "_len_line_charts_integrated.pdf";
  string msi_rep_lens_rate_plots_file = out_name + "_rate_line_charts_integrated.pdf";
  map<string, vector<int>> loci_plot_area;
  string total_reads_matrix_file = out_name + "_msi_support_reads_matrix.txt";
  string mapping_qual_matrix_file = out_name + "_msi_support_reads_mapping_quality_matrix.txt";
  
  
  writeLocsExactLensSupportReadsCountMatrix(msi_locs_rep_lens_map, total_reads_matrix_file);
  writeLocsSupportReadsMappingQualMatrix(loci_map_qual_map, mapping_qual_matrix_file);
  
  string loci_plot_area_file = out_name + "_msi_loci_plot_area.txt";
  string out_freq_rate_matrix = out_name + "_msi_locis_matrix.txt";
  string out_freq_matrix = out_name + "_msi_locis_freq_matrix.txt";
  string locs_distribution_file = out_name + "_msi_distribution_plots.pdf";
  string loci_distribution_file = out_name + "_msi_distribution_plots";
  writeSampleLocisFreqRateMatrix(msi_locs_exact_lens_freq_rate, out_freq_rate_matrix, loci_plot_area_file,
                                 loci_plot_area);
  writeLocsExactLensFreqMatrix(msi_locs_exact_lens_freq, out_freq_matrix);
  
  //string plot_r = (string) INSTALLDIR + "/scripts/plot_msi_distribution.R";
  //string cmd = (string) Rscript_BIN + "  " + plot_r + " " + out_freq_rate_matrix + "  " + loci_plot_area_file + "  " +
   //            locs_distribution_file;
  if (plot)
  {
    //system(cmd.c_str());
    //system(("rm -f " + loci_plot_area_file).c_str());
  }
  string diff_stats_file = out_name + "_msi_status.txt";
  map<string, loci_baseline> loci_baseline_map;
  map<string, double> loci_difference_map;
  map<string, string> loci_label_map;
  map<string, sample_loci_model_data> loci_model_data_map;
  if (baseline)
  {
    newProcessBaselineFile(base_file, loci_baseline_map);
    newComputeDiffFromBaseline(loci_baseline_map, msi_locs_exact_lens_freq_rate, loci_difference_map, out_name,
                            loci_label_map,
                            loci_model_data_map);
    
    newWriteDiffStats(diff_stats_file, loci_model_data_map, msi_parts);
  }
}

void writeLocsSupportReadsMappingQualMatrix(map<string, map<int, float> > &loci_map_qual_map, string file_to_write)
{
  map<string, map<int, float> >::iterator iterator1;
  map<int, float>::iterator iterator2;
  set<int> qual_set;
  set<int>::iterator qual_set_it;
  vector<int> qual_vec;
  for (iterator1 = loci_map_qual_map.begin(); iterator1 != loci_map_qual_map.end(); iterator1++)
  {
    for (iterator2 = iterator1->second.begin(); iterator2 != iterator1->second.end(); iterator2++)
    {
      qual_set.insert(iterator2->first);
    }
  }
  string head = "loci\t";
  
  for (qual_set_it = qual_set.begin(); qual_set_it != qual_set.end(); ++qual_set_it)
  {
    qual_vec.push_back(*qual_set_it);
  }
  sort(qual_vec.begin(), qual_vec.end());
  for (int i = 0; i < qual_vec.size(); ++i)
  {
    head += to_string(qual_vec[i]) + "\t";
  }
  head += "\n";
  
  ofstream out;
  string tmp_str;
  out.open(file_to_write.c_str());
  out << head;
  for (iterator1 = loci_map_qual_map.begin(); iterator1 != loci_map_qual_map.end(); iterator1++)
  {
    tmp_str = iterator1->first + "\t";
    for (int i = 0; i < qual_vec.size(); ++i)
    {
      string value = "0.00";
      iterator2 = iterator1->second.find(qual_vec[i]);
      if (iterator2 != iterator1->second.end())
        value = to_string(iterator2->second);
      tmp_str += value + "\t";
    }
    
    tmp_str += "\n";
    out << tmp_str;
  }
  out.close();
}

bool compare_index(int a, int b)
{
  return a < b;
}

void writeLocsExactLensFreqMatrix(map<string, map<int, int>> &freq_map, string out_freq_matrix)
{
  map<string, map<int, int >>::iterator it;
  map<int, int> tmp_freqs;
  map<int, int>::iterator tmp_freqs_it;
  ofstream out;
  string tmp_str;
  int tmp_start, tmp_end;
  vector<int> non_zero_index;
  vector<int> area;
  out.open((out_freq_matrix).c_str());
  out << "loci\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\t17\t18\t19\t20\t"
      << "21\t22\t23\t24\t25\t26\t27\t28\t29\t30\t31\t32\t33\t34\t35\t36\t37\t38\t39\t40\t"
      << "41\t42\t43\t44\t45\t46\t47\t48\t49\t50\n";
  for (it = freq_map.begin(); it != freq_map.end(); ++it)
  {
    
    tmp_str = "";
    non_zero_index.clear();
    for (int i = 0; i < 51; ++i)
    {
      tmp_freqs_it = it->second.find(i);
      if (tmp_freqs_it != it->second.end())
      {
        tmp_str = tmp_str + "\t" + to_string(tmp_freqs_it->second);
      }
      else
      {
        tmp_str = tmp_str + "\t" + to_string(-1.0); //exception
        cerr << "the frequency rate map is abnormal, lacking the index: " << i << std::endl;
      }
    }
    out << it->first << tmp_str << std::endl;
    
  }
  out.close();
  
}

void
writeSampleLocisFreqRateMatrix(map<string, map<int, float>> &freq_rate_map, string out_freq_rate_matrix,
                               string loci_plot_area_file, map<string, vector<int>> &loci_plot_area)
{
  map<string, map<int, float >>::iterator it;
  map<int, float> tmp_freqs;
  map<int, float>::iterator tmp_freqs_it;
  ofstream out;
  string tmp_str;
  int tmp_start, tmp_end;
  vector<int> non_zero_index;
  vector<int> area;
  out.open((out_freq_rate_matrix).c_str());
  out << "loci\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\t17\t18\t19\t20\t"
      << "21\t22\t23\t24\t25\t26\t27\t28\t29\t30\t31\t32\t33\t34\t35\t36\t37\t38\t39\t40\t"
      << "41\t42\t43\t44\t45\t46\t47\t48\t49\t50\n";
  for (it = freq_rate_map.begin(); it != freq_rate_map.end(); ++it)
  {
    
    tmp_str = "";
    non_zero_index.clear();
    for (int i = 0; i < 51; ++i)
    {
      tmp_freqs_it = it->second.find(i);
      if (tmp_freqs_it != it->second.end())
      {
        tmp_str = tmp_str + "\t" + to_string(tmp_freqs_it->second);
        if (tmp_freqs_it->second != 0)
        {
          non_zero_index.push_back(tmp_freqs_it->first);
        }
      }
      else
      {
        tmp_str = tmp_str + "\t" + to_string(-1.0); //exception
        cerr << "the frequency rate map is abnormal, lacking the index: " << i << std::endl;
      }
    }
    out << it->first << tmp_str << std::endl;
    sort(non_zero_index.begin(), non_zero_index.end(), compare_index);
    if (non_zero_index.size() == 0)
    {
      cerr << "abnormal distribution." << std::endl;
      tmp_start = 0;
      tmp_end = 50;
    }
    else
    {
      tmp_start = non_zero_index[0] - 5;
      if (tmp_start < 0)
        tmp_start = 0;
      tmp_end = non_zero_index[non_zero_index.size() - 1] + 5;
      if (tmp_end > 50)
        tmp_end = 50;
    }
    
    area.clear();
    area.push_back(tmp_start);
    area.push_back(tmp_end);
    loci_plot_area[it->first] = area;
  }
  out.close();
  
  out.open((loci_plot_area_file).c_str());
  out << "loci\tareaStart\tareaEnd" << std::endl;
  map<string, vector<int>>::iterator area_it;
  for (area_it = loci_plot_area.begin(); area_it != loci_plot_area.end(); ++area_it)
  {
    out << area_it->first << "\t" << area_it->second[0] << "\t" << area_it->second[1] << std::endl;
  }
  out.close();
  
}

void computeDiffFromBaseline(map<string, loci_baseline> &loci_baseline_map,
                             map<string, map<int, float>> &msi_locs_freq_rate_map,
                             map<string, double> &loci_difference_map, string out_name,
                             map<string, string> &loci_label_map,
                             map<string, sample_loci_model_data> &loci_model_data_map)
{
  map<string, double>::iterator loci_difference_map_it;
  map<string, map<int, float>>::iterator msi_locs_freq_rate_map_it;
  
  map<string, loci_baseline>::iterator loci_baseline_map_it;
  map<string, sample_loci_model_data>::iterator loci_model_data_map_it;
  sample_loci_model_data tmp_sample_loci_model;
  double tmp_base_mean, tmp_base_sd, tmp_base_diff_threshold;
  double tmp_model_value, tmp_p_value, tmp_log_pvalue;
  int countIndex;
  stringstream line;
  string tmp_out_file = out_name + "_tmp_pnorm_result.txt";
  ifstream in;
  string tmp;
  vector<string> tmp_substr;
  string tmp_label;
  map<string, string>::iterator loci_label_map_it;
  
  
  for (msi_locs_freq_rate_map_it = msi_locs_freq_rate_map.begin();
    msi_locs_freq_rate_map_it != msi_locs_freq_rate_map.end(); ++msi_locs_freq_rate_map_it)
  {
    tmp_p_value = -1.0;
    tmp_model_value = 0.0;
    tmp_log_pvalue = -1.0;
    tmp_label = "";
    
    
    loci_baseline_map_it = loci_baseline_map.find(msi_locs_freq_rate_map_it->first);
    if (loci_baseline_map_it != loci_baseline_map.end())
    {
      tmp_base_mean = loci_baseline_map_it->second.mean;
      tmp_base_sd = loci_baseline_map_it->second.sd;
      tmp_base_diff_threshold = loci_baseline_map_it->second.diffThreshold;
      
      tmp_sample_loci_model.msi_marker_name = msi_locs_freq_rate_map_it->first;
      tmp_sample_loci_model.modelMean = tmp_base_mean;
      tmp_sample_loci_model.modelSd = tmp_base_sd;
      tmp_sample_loci_model.modelPvalueThreshold = tmp_base_diff_threshold;
      tmp_sample_loci_model.modelRepeatCounts = loci_baseline_map_it->second.repeatCountSet;
      tmp_sample_loci_model.logPvalueThreshold = -log10(tmp_base_diff_threshold);
      
      for (int i = 0; i < loci_baseline_map_it->second.repeatCountSet.size(); ++i)
      {
        countIndex = loci_baseline_map_it->second.repeatCountSet[i];
        tmp_model_value += double((msi_locs_freq_rate_map_it->second)[countIndex]);
      }
      
      tmp_sample_loci_model.sampleModelValue = tmp_model_value;
      
      
      line.str("");
      line.clear();
      line << Rscript_BIN << "  " << (string) INSTALLDIR + "/R/pnorm_rscript.R" << "  " << tmp_model_value
           << "   " << tmp_base_mean << "   " << tmp_base_sd << "  >  " << tmp_out_file;
      
      system(line.str().c_str());
      in.open(tmp_out_file.c_str());
      if (!in.is_open())
      {
        cerr << "Error: cannot open msi-partition table:\t" << tmp_out_file << endl;
        exit(1);
      }
      tmp_substr.clear();
      getline(in, tmp, '\n');
      tmp_substr = split_string(tmp, " ");
      if (tmp_substr.size() == 2)
      {
        tmp_p_value = stod(tmp_substr[1]);
      }
      else
      {
        tmp_p_value = -1.0;
      }
      tmp_sample_loci_model.sampleModelPvalue = tmp_p_value;
      
      tmp_substr.clear();
      getline(in, tmp, '\n');
      tmp_substr = split_string(tmp, " ");
      if (tmp_substr.size() == 2)
      {
        if (tmp_substr[1] != "Inf")
        {
          tmp_log_pvalue = stod(tmp_substr[1]);
        }
        else
        {
          tmp_log_pvalue = -1.0;
        }
      }
      
      tmp_sample_loci_model.sampleLogPvalue = tmp_log_pvalue;
      
      in.clear();
      in.close();
      loci_difference_map_it = loci_difference_map.find(msi_locs_freq_rate_map_it->first);
      if (loci_difference_map_it == loci_difference_map.end())
      {
        loci_difference_map[msi_locs_freq_rate_map_it->first] = tmp_p_value;
      }
      
      tmp_label = (tmp_p_value >= 0 && tmp_p_value < tmp_base_diff_threshold) ? "unstable" : "stable";
      
      tmp_sample_loci_model.loci_label = tmp_label;
      cout << "loci: " << msi_locs_freq_rate_map_it->first << " label: " << tmp_label << std::endl;
      loci_label_map_it = loci_label_map.find(msi_locs_freq_rate_map_it->first);
      if (loci_label_map_it == loci_label_map.end())
      {
        loci_label_map[msi_locs_freq_rate_map_it->first] = tmp_label;
      }
      
      loci_model_data_map_it = loci_model_data_map.find(msi_locs_freq_rate_map_it->first);
      if (loci_model_data_map_it == loci_model_data_map.end())
      {
        loci_model_data_map[msi_locs_freq_rate_map_it->first] = tmp_sample_loci_model;
      }
      
    }
    else
    {
      cerr << "the baseline has no info about loci: " << msi_locs_freq_rate_map_it->first << std::endl;
    }
    
  }
  system(("rm -rf " + tmp_out_file).c_str());
}

/**
 *
 * @param loci_baseline_map BASELINE FILE INFO:
 * @param msi_locs_freq_rate_map
 * @param loci_difference_map
 * @param out_name
 * @param loci_label_map
 * @param loci_model_data_map
 */

void newComputeDiffFromBaseline(map<string, loci_baseline> &loci_baseline_map,
                                map<string, map<int, float>> &msi_locs_freq_rate_map,
                                map<string, double> &loci_difference_map, string out_name,
                                map<string, string> &loci_label_map,
                                map<string, sample_loci_model_data> &loci_model_data_map)
{
  map<string, double>::iterator loci_difference_map_it;
  map<string, map<int, float>>::iterator msi_locs_freq_rate_map_it;
  
  map<string, loci_baseline>::iterator loci_baseline_map_it;
  map<string, sample_loci_model_data>::iterator loci_model_data_map_it;
  sample_loci_model_data tmp_sample_loci_model;
  
  int countIndex;
  stringstream line;
  
  ifstream in;
  string tmp;
  vector<string> tmp_substr;
  string tmp_label;
  map<string, string>::iterator loci_label_map_it;
  
  
  for (msi_locs_freq_rate_map_it = msi_locs_freq_rate_map.begin();
    msi_locs_freq_rate_map_it != msi_locs_freq_rate_map.end(); ++msi_locs_freq_rate_map_it)
  {
    
    loci_baseline_map_it = loci_baseline_map.find(msi_locs_freq_rate_map_it->first);
    if (loci_baseline_map_it != loci_baseline_map.end())
    {
      tmp_sample_loci_model.msi_marker_name = msi_locs_freq_rate_map_it->first;
      tmp_sample_loci_model.modelMean = loci_baseline_map_it->second.mean;
      tmp_sample_loci_model.modelSd = loci_baseline_map_it->second.sd;
      tmp_sample_loci_model.modelPvalueThreshold = loci_baseline_map_it->second.pvalue_thres;
      tmp_sample_loci_model.modelRepeatCounts = loci_baseline_map_it->second.modelIndexSet;
   
      
      tmp_sample_loci_model.sampleModelValue = 0.0;
      for (int i = 0; i < loci_baseline_map_it->second.modelIndexSet.size(); ++i)
      {
        countIndex = loci_baseline_map_it->second.modelIndexSet[i];
        tmp_sample_loci_model.sampleModelValue += double((msi_locs_freq_rate_map_it->second)[countIndex]);
      }
      
      tmp_sample_loci_model.sampleModelPvalue = getProbOfNDvalue(tmp_sample_loci_model.sampleModelValue,
                                                                 tmp_sample_loci_model.modelMean,
                                                                 tmp_sample_loci_model.modelSd);
      if (loci_baseline_map_it->second.valid_direction == "g")
      {
        if (tmp_sample_loci_model.sampleModelValue > loci_baseline_map_it->second.value_thres)
        {
          tmp_sample_loci_model.loci_label = "stable";
        }
        else
        {
          tmp_sample_loci_model.loci_label = "unstable";
        }
      }
      else if (loci_baseline_map_it->second.valid_direction == "l")
      {
        if (tmp_sample_loci_model.sampleModelValue < loci_baseline_map_it->second.value_thres)
        {
          tmp_sample_loci_model.loci_label = "stable";
        }
        else
        {
          tmp_sample_loci_model.loci_label = "unstable";
        }
      }
      else
      {
        cerr << "illegal model valid direction(must be \'g\' or \'l\')\n";
        exit(-1);
      }
  
      cout << "loci: " << msi_locs_freq_rate_map_it->first << " label: " << tmp_sample_loci_model.loci_label << std::endl;
      loci_label_map_it = loci_label_map.find(msi_locs_freq_rate_map_it->first);
      if (loci_label_map_it == loci_label_map.end())
      {
        loci_label_map[msi_locs_freq_rate_map_it->first] = tmp_sample_loci_model.loci_label;
      }
  
      loci_model_data_map_it = loci_model_data_map.find(msi_locs_freq_rate_map_it->first);
      if (loci_model_data_map_it == loci_model_data_map.end())
      {
        loci_model_data_map[msi_locs_freq_rate_map_it->first] = tmp_sample_loci_model;
      }
    }
    else
    {
      cerr << "the baseline has no info about loci: " << msi_locs_freq_rate_map_it->first << std::endl;
    }
  }
}

void process_baseline_file(string baseline_file, map<string, loci_baseline> &loci_baseline_map)
{
  loci_baseline tmp_baseline;
  map<string, loci_baseline>::iterator loci_baseline_map_it;
  ifstream in;
  in.open(baseline_file.c_str());
  string tmp;
  stringstream line;
  vector<string> tmp_value_set;
  vector<int> tmp_repeat_count_set;
  vector<float> tmp_base_rate_set;
  double tmp_mean, tmp_sum, tmp_square_sum;
  if (!in.is_open())
  {
    cerr << "Error: cannot open msi-baseline table:\t" << baseline_file << endl;
    exit(1);
  }
  getline(in, tmp, '\n');
  while (getline(in, tmp, '\n'))
  {
    line.str();
    line.clear();
    line << tmp;
    getline(line, tmp, '\t');
    tmp_baseline.msi_marker = tmp;
    getline(line, tmp, '\t');
    tmp_baseline.diffThreshold = stof(tmp);
    
    getline(line, tmp, '\t');
    getline(line, tmp, '\t');
    tmp_value_set = split_string(tmp, ",");
    tmp_repeat_count_set.clear();
    for (int i = 0; i < tmp_value_set.size(); ++i)
    {
      tmp_repeat_count_set.push_back(stoi(tmp_value_set[i]));
    }
    tmp_baseline.repeatCountSet = tmp_repeat_count_set;
    
    getline(line, tmp, '\t');
    tmp_value_set = split_string(tmp, ",");
    tmp_base_rate_set.clear();
    for (int i = 0; i < tmp_value_set.size(); ++i)
    {
      tmp_base_rate_set.push_back(stof(tmp_value_set[i].c_str()));
    }
    tmp_baseline.normalValueSet = tmp_base_rate_set;
    
    tmp_sum = 0;
    for (int j = 0; j < tmp_base_rate_set.size(); ++j)
    {
      tmp_sum += tmp_base_rate_set[j];
    }
    tmp_mean = tmp_sum / tmp_base_rate_set.size();
    tmp_square_sum = 0;
    for (int j = 0; j < tmp_base_rate_set.size(); ++j)
    {
      tmp_square_sum += (tmp_base_rate_set[j] - tmp_mean) * (tmp_base_rate_set[j] - tmp_mean);
    }
    tmp_baseline.mean = tmp_mean;
    tmp_baseline.sd = sqrt(tmp_square_sum / tmp_base_rate_set.size());
    loci_baseline_map_it = loci_baseline_map.find(tmp_baseline.msi_marker);
    if (loci_baseline_map_it == loci_baseline_map.end())
    {
      loci_baseline_map[tmp_baseline.msi_marker] = tmp_baseline;
    }
  }
}

void newProcessBaselineFile(string baseline_file, map<string, loci_baseline> &loci_baseline_map)
{
  loci_baseline tmp_baseline;
  map<string, loci_baseline>::iterator loci_baseline_map_it;
  ifstream in;
  in.open(baseline_file.c_str());
  string tmp;
  stringstream line;
  vector<string> tmp_value_set, tmp_model_index_set;
  vector<int> tmp_repeat_count_set;
  vector<float> tmp_base_rate_set;
  double tmp_mean, tmp_sum, tmp_square_sum;
  if (!in.is_open())
  {
    cerr << "Error: cannot open msi-baseline table:\t" << baseline_file << endl;
    exit(1);
  }
  getline(in, tmp, '\n');
  while (getline(in, tmp, '\n'))
  {
    line.str();
    line.clear();
    line << tmp;
    getline(line, tmp, '\t');
    tmp_baseline.msi_marker = tmp;
    getline(line, tmp, '\t');
    tmp_baseline.value_thres = stof(tmp);
    getline(line, tmp, '\t');
    tmp_baseline.pvalue_thres = stof(tmp);
    getline(line, tmp, '\t');
    tmp_baseline.valid_direction = tmp;
    getline(line, tmp, '\t');
    tmp_baseline.mean = stof(tmp);
    getline(line, tmp, '\t');
    tmp_baseline.sd = stof(tmp);
    
    getline(line, tmp, '\t');
    tmp_baseline.modelIndexStr = tmp;
    tmp_model_index_set = split_string(tmp, ",");
    tmp_baseline.modelIndexSet.clear();
    for (int i = 0; i < tmp_model_index_set.size(); ++i)
    {
      tmp_baseline.modelIndexSet.push_back(stoi(tmp_model_index_set[i]));
    }
    
    getline(line, tmp, '\t');
    tmp_baseline.normalValueStr = tmp;
    tmp_value_set = split_string(tmp, ",");
    tmp_baseline.normalValueSet.clear();
    for (int i = 0; i < tmp_value_set.size(); ++i)
    {
      tmp_baseline.normalValueSet.push_back(stof(tmp_value_set[i]));
    }
    getline(line, tmp, '\t');
    tmp_baseline.baselineSampleSize = stoi(tmp);
    
    loci_baseline_map_it = loci_baseline_map.find(tmp_baseline.msi_marker);
    if (loci_baseline_map_it == loci_baseline_map.end())
    {
      loci_baseline_map[tmp_baseline.msi_marker] = tmp_baseline;
    }
  }
}

void normalizeExactFrequencyForEachMsi(const msi_part_data &msi_part, map<int, int> &frequency_map,
                                       map<int, float> &frequency_rate_map)
{
  frequency_rate_map.clear();
  
  if (frequency_map.size() == 0)
  {
    return;
  }
  map<int, float>::iterator freq_rate_it;
  map<int, int>::iterator freq_it;
  int sum = 0;
  for (freq_it = frequency_map.begin(); freq_it != frequency_map.end(); ++freq_it)
  {
    sum += freq_it->second;
  }
  for (freq_it = frequency_map.begin(); freq_it != frequency_map.end(); ++freq_it)
  {
    freq_rate_it = frequency_rate_map.find(freq_it->first);
    if (freq_rate_it == frequency_rate_map.end())
    {
      frequency_rate_map[freq_it->first] = (float) freq_it->second / (float) sum;
    }
  }
  return;
}

void getExactRepLensFreqForEachMsiLoc(const vector<msi_rep_len> &tmp_repeat_lengths, const msi_part_data &msi_part,
                                      map<int, int> &tmp_exact_rep_lens_freq)
{
  //此函数已完成
  tmp_exact_rep_lens_freq.clear();
  map<int, int>::iterator iterator1;
  if (tmp_repeat_lengths.size() == 0)
  {
    tmp_exact_rep_lens_freq.clear();
    return;
  }
  
  for (int i = 0; i < tmp_repeat_lengths.size(); ++i)
  {
    if (tmp_repeat_lengths[i].exact && tmp_repeat_lengths[i].rep_len != -1)
    {
      iterator1 = tmp_exact_rep_lens_freq.find(tmp_repeat_lengths[i].rep_len);
      if (iterator1 == tmp_exact_rep_lens_freq.end())
      {
        tmp_exact_rep_lens_freq[tmp_repeat_lengths[i].rep_len] = 0;
      }
    }
  }
  
  for (int i = 0; i < tmp_repeat_lengths.size(); ++i)
  {
    if (tmp_repeat_lengths[i].exact && tmp_repeat_lengths[i].rep_len != -1)
    {
      iterator1 = tmp_exact_rep_lens_freq.find(tmp_repeat_lengths[i].rep_len);
      if (iterator1 == tmp_exact_rep_lens_freq.end())
      {
        cerr << "abnormal repeat counts: " << tmp_repeat_lengths[i].rep_len << std::endl;
        exit(-1);
      }
      tmp_exact_rep_lens_freq[tmp_repeat_lengths[i].rep_len]++;
    }
  }
  
}

void writeDiffStats(string diff_stats_file, map<string, sample_loci_model_data> &loci_model_data_map,
                    vector<msi_part_data> &msi_parts)
{
  ofstream out;
  string tmp_str;
  stringstream ss;
  int sum = loci_model_data_map.size();
  int unstable_count = 0;
  double unstable_rate = 0.0;
  map<string, sample_loci_model_data>::iterator it;
  out.open(diff_stats_file.c_str());
  out << "loci\tsampleModelValue\tsampleModelPvalue\tsampleMinusLog10Pvalue\t"
      << "modelRepeatCountIndexes\tmodelMean\tmodelSd\tmodelPvalueDiffThreshold\tmodelMinusLog10PvalueDiffThreshold\tsampleLociLabel\n";
  for (it = loci_model_data_map.begin(); it != loci_model_data_map.end(); ++it)
  {
    
    ////start to write the file
    
    out << it->first << "\t";
    out << it->second.sampleModelValue << "\t";
    out << it->second.sampleModelPvalue << "\t";
    out << it->second.sampleLogPvalue << "\t";
    
    ss.str("");
    ss.clear();
    for (int j = 0; j < it->second.modelRepeatCounts.size() - 1; ++j)
    {
      ss << it->second.modelRepeatCounts[j] << ",";
    }
    ss << it->second.modelRepeatCounts[(it->second.modelRepeatCounts.size() - 1)];
    out << ss.str() << "\t";
    
    out << it->second.modelMean << "\t";
    out << it->second.modelSd << "\t";
    out << it->second.modelPvalueThreshold << "\t";
    out << -log10(it->second.modelPvalueThreshold) << "\t";
    out << it->second.loci_label;
    out << std::endl;
    
    if (it->second.loci_label == "unstable")
      unstable_count++;
    
  }
  if (sum == 0)
  {
    out << "conclusion: the sample has no msi loci which contain valued data.\n";
    out.close();
    return;
  }
  
  tmp_str = "not been decided yet";
  unstable_rate = double(unstable_count) / double(sum);
  if (unstable_count >= 0 && unstable_rate < 0.2)
  {
    tmp_str = "MSI-Stable";
  }
  if (unstable_rate >= 0.2 && unstable_rate < 0.4)
  {
    tmp_str = "MSI-Low";
  }
  if (unstable_rate <= 1 && unstable_rate >= 0.4)
  {
    tmp_str = "MSI-High";
  }
  
  out
    << "#Judge Criterion: \n#  Unstable Rate is in [0,0.2): MSI-Stable\n#  Unstable Rate is in [0.2,0.4): MSI-Low\n#  Unstable Rate is in [0.4,1]: MSI-High\n";
  out << "#Note: The -1 in \"sampleMinusLog10Pvalue\" column means infinite value.\n";
  out << "#Conclusion:\n#  Unstable loci count is: " << unstable_count << ", Total Loci Count is: " << sum
      << "\n#  The unstable rate is: " << unstable_rate << "\n#  The sample MSI feature is: " << tmp_str
      << std::endl;
  
  out.close();
  
}

bool check_baseline_file(string base_file)
{
  return true;
}

bool cmp_pos(BamAlignment &a1, BamAlignment &a2)
{
  return (a1.getAlignmentStart() < a1.getAlignmentStart());
}

bool cmp_read_names(BamAlignment a1, BamAlignment a2)
{
  return (a1.getReadName() < a2.getReadName());
}

void
getRepeatLengthsForEachMsiLoc(vector<msi_rep_len> &msi_loci_rep_lengths, const msi_part_data &msi_part, samfile_t *fp,
                              bam_index_t *idx_bam, int readLen, int min_coverage, map<int, float> &quality_map)
{
  bam1_t *b = bam_init1();
  BamAlignment src_align;
  msi_rep_len tmp_rep_len;
  BamAlignment r1_align, r2_align;
  InsertMole molecule, current_insert;
  vector<InsertMole> actual_evident_reads;
  std::string chromosome = msi_part.chr;
  int tid_t;
  int beg_t, end_t;
  bam_iter_t iter_bam;
  int32_t left_pos = msi_part.msi_exact_start;
  int32_t right_pos = msi_part.msi_exact_end;
  
  bam_parse_region(fp->header, chromosome.c_str(), &tid_t, &beg_t, &end_t);
  vector<BamAlignment> src_aligns;
  int readCount = 0;
  iter_bam = bam_iter_query(idx_bam, tid_t, left_pos - 2 * readLen, right_pos + 2 * readLen);
  bool covered;
  CigarRoller tmp_cigar;
  /// load all the reads overlap with the region
  while (bam_iter_read(fp->x.bam, iter_bam, b) >= 0)
  {
    if (b->core.flag & (BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FUNMAP | BAM_FSECONDARY))
    {
      continue;
    }
    if (!b->core.flag & BAM_FPROPER_PAIR)
    {
      continue;
    }
    src_aligns.push_back(BamAlignment());
    
  }
  
  iter_bam = bam_iter_query(idx_bam, tid_t, left_pos - 2 * readLen, right_pos + 2 * readLen);
  int i = 0;
  /// load all the reads overlap with the region
  while (bam_iter_read(fp->x.bam, iter_bam, b) >= 0)
  {
    if (b->core.flag & (BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FUNMAP | BAM_FSECONDARY))
    {
      continue;
    }
    if (!b->core.flag & BAM_FPROPER_PAIR)
    {
      continue;
    }
    src_aligns[i].duplicateBamRecord(b);
    i++;
  }
  
  actual_evident_reads.clear();
  string r1_read_name, r2_read_name;
  map<string, vector<BamAlignment>> src_aligns_map;
  map<string, vector<BamAlignment>>::iterator src_aligns_map_it;
  vector<BamAlignment> tmp;
  BamAlignment tmp_align;
  string read_name;
  for (int j = 0; j < src_aligns.size(); ++j)
  {
    if (src_aligns[j].getExpandedCigarString().rfind('M') !=string::npos)
    {
      read_name = src_aligns[j].getReadName();
      src_aligns_map_it = src_aligns_map.find(read_name);
      if (src_aligns_map_it == src_aligns_map.end())
      {
        tmp_align.readBamRecord(src_aligns[j].getBamRecord());
        src_aligns_map[read_name] = vector<BamAlignment>(1, tmp_align);
        tmp_align.reset();
      }
      else
      {
        tmp_align.readBamRecord(src_aligns[j].getBamRecord());
        src_aligns_map[read_name].push_back(tmp_align);
        tmp_align.reset();
      }
    }
 
  }
  
  for (src_aligns_map_it = src_aligns_map.begin(); src_aligns_map_it != src_aligns_map.end(); ++src_aligns_map_it)
  {
    int align_num = src_aligns_map_it->second.size();
    if (align_num > 2)
    {
      cerr << "read: " << src_aligns_map_it->first << " has " << align_num << " reads, (more than 2 reads)\n";
//      exit(-1);
    }
    else
    {
      if (align_num == 2)
      {
        molecule.Set(src_aligns_map_it->second[0], src_aligns_map_it->second[1]);
        actual_evident_reads.push_back(molecule);
        continue;
      }
      if (align_num == 1)
      {
        molecule.Set(src_aligns_map_it->second[0]);
        actual_evident_reads.push_back(molecule);
        continue;
      }
    }
  }
  
  cout << "total " << actual_evident_reads.size() << " molecule spanning the loci: " << msi_part.msi_marker_name
       << std::endl;
  msi_loci_rep_lengths.clear();
  vector<int> map_qual_vec;
  for (int i = 0; i < actual_evident_reads.size(); ++i)
  {
    current_insert = actual_evident_reads[i];
    tmp_rep_len.support = false;
    tmp_rep_len.rep_len = -1;
    tmp_rep_len.exact = false;
    tmp_rep_len.scope = false;
    tmp_rep_len.relation = 'n'; // enum{'n', 'e', 'g', 'l'} n: no relation, e:equal, g: greater than(>=), l: less than(<=)
    tmp_rep_len.overlapCigar = "";
    tmp_rep_len.overlapOC = "";
    tmp_rep_len.overlap_bwa_soft = false;
    tmp_rep_len.overlap_start = -1;
    tmp_rep_len.overlap_end = -1;
    tmp_rep_len.has_insert = false;
    tmp_rep_len.insert_base = "";
    //    if (current_align.getReadName() == "ST-E00161:435:H7K37CCXY:5:2107:29792:44450")
//    {
    if (current_insert.getSumMappingQual() >= 10 && msiLociOverlap(current_insert, msi_part, tmp_rep_len))
    {
      readCount++;
//      if (current_insert.getReadName() == "ST-E00161:435:H7K37CCXY:5:2222:19015:64562")
//      {
      computeRepeatLenPerAlignment(msi_part, current_insert, tmp_rep_len);
      map_qual_vec.push_back(current_insert.getSumMappingQual());
      
      msi_loci_rep_lengths.push_back(tmp_rep_len);
      if (tmp_rep_len.rep_len > 35)
      {
        cout << "abnormal length: " << tmp_rep_len.rep_len
             << ".\nread name: " << current_insert.getReadName()
             << ".\nthe msi loci is: " << msi_part.msi_marker_name << std::endl;
      }

//      }
      
    }
    else
    {
      continue;
    }
//    }
  }
  map<int, int> map_qual_count_map;
  map<int, int>::iterator map_qual_count_map_it;
  
  for (int k = 0; k < map_qual_vec.size(); ++k)
  {
    map_qual_count_map[map_qual_vec[k]] = 0;
  }
  for (int k = 0; k < map_qual_vec.size(); ++k)
  {
    map_qual_count_map[map_qual_vec[k]]++;
  }
  for (map_qual_count_map_it = map_qual_count_map.begin();
    map_qual_count_map_it != map_qual_count_map.end(); map_qual_count_map_it++)
  {
    quality_map[map_qual_count_map_it->first] = (float) map_qual_count_map_it->second / (float) map_qual_vec.size();
  }


//  float average_coverage = float(readCount) / float(right_pos - left_pos);
  cout << "the loci: " << msi_part.msi_marker_name << " has  exact span reads coverage: " << readCount << std::endl;
  if (readCount < min_coverage)
  {
    cout << "the loci: " << msi_part.msi_marker_name << " has no enough exact span reads coverage ." << std::endl;
    msi_loci_rep_lengths.clear();
    return;
  }
}


/**
 *
 * @param mole
 * @param msi_part
 * @param tmp_rep_len
 * @return
 */

bool msiLociOverlap(InsertMole &mole, const msi_part_data &msi_part, msi_rep_len &tmp_rep_len)
{
  bool covered, isolate, special_first;
  string align_chr;
  align_chr = mole.getChrName();
  if (align_chr != msi_part.chr)
    return false;
  int record_num = mole.get_align_record_num();
  if (record_num == 2)
  {
    if (mole.isGapped())
    {
      covered = ((mole.getAlignmentLeftStart() <= msi_part.msi_exact_start) &&
                 (mole.getAlignmentLeftEnd() >= msi_part.msi_exact_end)) ||
                ((mole.getAlignmentRightStart() <= msi_part.msi_exact_start) &&
                 (mole.getAlignmentRightEnd() >= msi_part.msi_exact_end));
      
    }
    if (mole.isOverlapped())
    {
      covered = (mole.getAlignmentLeftStart() <= msi_part.msi_exact_start) &&
                (mole.getAlignmentRightEnd() >= msi_part.msi_exact_end);
    }
  }
  else if (record_num == 1)
  {
    covered = (mole.getAlignmentSingleStart() <= msi_part.msi_exact_start) &&
              (mole.getAlignmentSingleEnd() >= msi_part.msi_exact_end);
  }
  else
  {
    cerr << "the insert molecule:" << mole.getReadName() << " has " << record_num << " alignments" << std::endl;
    exit(-1);
  }
  
  if (covered)
  {
    tmp_rep_len.exact = true;
    tmp_rep_len.support = true;
    return true;
  }


//  isolate = (mole.getAlignmentRightStart() < msi_part.msi_exact_start) ||
//            (mole.getAlignmentStart() > msi_part.msi_exact_end);
//  if (!(isolate || covered))
//  {
//    tmp_rep_len.scope = true;
//    tmp_rep_len.support = true;
//    return true;
//  }
  return false;
}

void computeRepeatLenPerAlignment(const msi_part_data &msi_part,
                                  InsertMole &insert_molecule,
                                  msi_rep_len &tmp_rep_len)
{
  
  string expand_cigar = "";
  string insert_bases = "";
//  string readSeq = insert_molecule.getSumReadSeq();
  CigarRoller cigar_roller;
  string molecule_seq;
  int32_t align_start_pos;
  if (tmp_rep_len.exact)
  {
    tmp_rep_len.relation = 'e';
    if (insert_molecule.get_align_record_num() == 2)
    {
      if (insert_molecule.isOverlapped())
      {
        tmp_rep_len.rep_len = 0;
        insert_molecule.getSumCigar(cigar_roller);
        molecule_seq = insert_molecule.getSumReadSeq();
        align_start_pos = insert_molecule.getAlignmentLeftStart();
        tmp_rep_len.rep_len += get_rep_count(cigar_roller, msi_part, molecule_seq, align_start_pos);
        return;
      }
      if (insert_molecule.isGapped())
      {
        if (insert_molecule.getAlignmentLeftStart() <= msi_part.msi_exact_start &&
            insert_molecule.getAlignmentLeftEnd() >= msi_part.msi_exact_end)
        {
          tmp_rep_len.rep_len = 0;
          insert_molecule.getLeftCigarRoller(cigar_roller);
          molecule_seq = insert_molecule.getLeftReadSeq();
          align_start_pos = insert_molecule.getAlignmentLeftStart();
          tmp_rep_len.rep_len += get_rep_count(cigar_roller, msi_part, molecule_seq, align_start_pos);
          return;
        }
        if (insert_molecule.getAlignmentRightStart() <= msi_part.msi_exact_start &&
            insert_molecule.getAlignmentRightEnd() >= msi_part.msi_exact_end)
        {
          tmp_rep_len.rep_len = 0;
          insert_molecule.getRightCigarRoller(cigar_roller);
          molecule_seq = insert_molecule.getRightReadSeq();
          align_start_pos = insert_molecule.getAlignmentRightStart();
          tmp_rep_len.rep_len += get_rep_count(cigar_roller, msi_part, molecule_seq, align_start_pos);
          return;
        }
      }
    }
    else if (insert_molecule.get_align_record_num() == 1)
    {
      tmp_rep_len.rep_len = 0;
      insert_molecule.getSingleCigarRoller(cigar_roller);
      molecule_seq = insert_molecule.getSingleReadSeq();
      align_start_pos = insert_molecule.getAlignmentSingleStart();
      tmp_rep_len.rep_len += get_rep_count(cigar_roller, msi_part, molecule_seq, align_start_pos);
      return;
    }
    else
    {
      cerr << "the insert molecule:" << insert_molecule.getReadName() << " has "
           << insert_molecule.get_align_record_num()
           << " alignments" << std::endl;
      exit(-1);
    }
    
  }
  
}

int get_rep_count(CigarRoller &roller, const msi_part_data &msi_part, const string &molecule_seq,
                  int32_t align_start_pos)
{
  int count = 0;
  int overlap_seq_start_index, overlap_seq_end_index, target_seq_start_index, target_seq_end_index;
  int overlap_cigar_start_index, overlap_cigar_end_index, target_cigar_start_index, target_cigar_end_index;
  int32_t align_end_pos = roller.getAlignmentEnd(align_start_pos);
  string cigar_str, target_cigar_str;
  string before_start_seq = "null";
  string after_end_seq = "null";
  if (align_start_pos == msi_part.msi_exact_start)
  {
    before_start_seq = "";
  }
  
  if (align_end_pos == msi_part.msi_exact_start)
  {
    after_end_seq = "";
  }
  
  if (roller.getCigarCharOpFromRefPos(msi_part.msi_exact_start - 1, align_start_pos) == 'D' ||
      roller.getCigarCharOpFromRefPos(msi_part.msi_exact_start - 1, align_start_pos) == 'S')
  {
    
    before_start_seq = "";
  }
  if (roller.getCigarCharOpFromRefPos(msi_part.msi_exact_end + 1, align_start_pos) == 'D' ||
      roller.getCigarCharOpFromRefPos(msi_part.msi_exact_end + 1, align_start_pos) == 'S')
  {
    
    after_end_seq = "";
  }
  overlap_seq_start_index = roller.getQueryIndex(msi_part.msi_exact_start - 1, align_start_pos) + 1;
  overlap_seq_end_index = roller.getQueryIndex(msi_part.msi_exact_end + 1, align_start_pos) - 1;
  target_seq_start_index = roller.getQueryIndex(msi_part.msi_exact_start, align_start_pos);
  target_seq_end_index = roller.getQueryIndex(msi_part.msi_exact_end, align_start_pos);
  overlap_cigar_start_index = roller.getExpandedCigarIndexFromRefPos(msi_part.msi_exact_start - 1, align_start_pos) + 1;
  overlap_cigar_end_index = roller.getExpandedCigarIndexFromRefPos(msi_part.msi_exact_end + 1, align_start_pos) - 1;
  target_cigar_start_index = roller.getExpandedCigarIndexFromRefPos(msi_part.msi_exact_start, align_start_pos);
  target_cigar_end_index = roller.getExpandedCigarIndexFromRefPos(msi_part.msi_exact_end, align_start_pos);
  
  cigar_str = roller.getExpandedString();
  target_cigar_str = cigar_str.substr(target_cigar_start_index, target_cigar_end_index - target_cigar_start_index + 1);
  if (before_start_seq == "null")
  {
    if (overlap_cigar_start_index < target_cigar_start_index)
    {
      if (overlap_seq_start_index != -1 && target_seq_start_index != -1)
        before_start_seq = molecule_seq.substr(overlap_seq_start_index,
                                               target_seq_start_index - overlap_seq_start_index);
      else
        cerr << "invalid seq index\n";
    }
    else
    {
      before_start_seq = "";
    }
  }
  
  if (after_end_seq == "null")
  {
    if (overlap_cigar_end_index > target_cigar_end_index)
    {
      if (overlap_cigar_end_index != -1 && target_cigar_end_index != -1)
        after_end_seq = molecule_seq.substr(target_seq_end_index + 1, overlap_seq_end_index - target_seq_end_index);
      else
        cerr << "invalid seq index\n";
    }
    else
    {
      after_end_seq = "";
    }
  }
  
  
  if (before_start_seq != "")
    count += countRepeatUnit(before_start_seq, msi_part.rep_unit, "right");
  if (after_end_seq != "")
    count += countRepeatUnit(after_end_seq, msi_part.rep_unit, "left");
  
  count += count_char(target_cigar_str, 'M');
  if (target_cigar_str.find('I') != std::string::npos)
  {
    vector<int> insert_cigar_indexes;
    for (int i = 0; i < target_cigar_str.size(); ++i)
    {
      if (target_cigar_str[i] == 'I')
        insert_cigar_indexes.push_back(i + target_cigar_start_index);
    }
    
    for (int j = 0; j < insert_cigar_indexes.size(); ++j)
    {
      for (int i = target_seq_start_index; i <= target_seq_end_index; ++i)
      {
        if (roller.getExpandedCigarIndexFromQueryIndex(i) == insert_cigar_indexes[j] &&
            molecule_seq[i] == msi_part.rep_unit[0])
          count++;
      }
    }
  }
  return count;
}


int count_char(const string &test_string, char c)
{
  
  int num = 0;
  for (int i = 0; i < test_string.size(); ++i)
  {
    if (test_string[i] == c)
    {
      num++;
    }
  }
  return num;
}

void process_msi_partition_table(vector<msi_part_data> &msi_part, const string msi_file)
{
  msi_part_data tmp_msi_part;
  ifstream in;
  string tmp;
  stringstream line;
  in.open(msi_file.c_str());
  
  if (!in.is_open())
  {
    cerr << "Error: cannot open msi-partition table:\t" << msi_file << endl;
    exit(1);
  }
  getline(in, tmp, '\n'); //throw away header line
  while (getline(in, tmp, '\n'))
  {
    line.str("");
    line.clear();
    line << tmp;
    line >> tmp_msi_part.chr;
    line >> tmp_msi_part.msi_exact_start;
    line >> tmp_msi_part.msi_exact_end;
    line >> tmp_msi_part.msi_marker_name;
    line >> tmp_msi_part.rep_unit;
    tmp_msi_part.msi_marker_name =
      tmp_msi_part.chr + "-" + to_string(tmp_msi_part.msi_exact_start) + "-" +
      to_string(tmp_msi_part.msi_exact_end) + "-" + tmp_msi_part.msi_marker_name;
    tmp_msi_part.msi_exact_start++;//将bed格式的左闭右开+0Based  变为  左闭右闭的1Based区间
    msi_part.push_back(tmp_msi_part);
  }
  in.close();
}

void writeLocsExactLensSupportReadsCountMatrix(map<string, vector<msi_rep_len>> &msi_locs_rep_lens_map,
                                               string total_reads_matrix_file)
{
  map<string, vector<msi_rep_len>>::iterator it;
  string new_file;
  ofstream out;
  out.open((total_reads_matrix_file).c_str());
  out << "loci\ttotalSupportReads\n";
  for (it = msi_locs_rep_lens_map.begin(); it != msi_locs_rep_lens_map.end(); ++it)
  {
    out << it->first << "\t" << it->second.size() << std::endl;
  }
  out.close();
}

void newWriteDiffStats(string diff_stats_file, map<string, sample_loci_model_data> &loci_model_data_map,
                       vector<msi_part_data> &msi_parts)
{
  ofstream out;
  string tmp_str;
  stringstream ss;
  int sum = loci_model_data_map.size();
  int unstable_count = 0;
  double unstable_rate = 0.0;
  map<string, sample_loci_model_data>::iterator it;
  out.open(diff_stats_file.c_str());
  out << "loci\tsampleModelValue\tsampleModelPvalue\t"
      << "modelRepeatCountIndexes\tmodelMean\tmodelSd\tmodelValueThreshold\tmodelPvalueThreshold\t"<<
   "sampleLociLabel\n";
  for (it = loci_model_data_map.begin(); it != loci_model_data_map.end(); ++it)
  {
    ////start to write the file
    
    out << it->first << "\t";
    out << it->second.sampleModelValue << "\t";
    out << it->second.sampleModelPvalue << "\t";
    
    ss.str("");
    ss.clear();
    for (int j = 0; j < it->second.modelRepeatCounts.size() - 1; ++j)
    {
      ss << it->second.modelRepeatCounts[j] << ",";
    }
    ss << it->second.modelRepeatCounts[(it->second.modelRepeatCounts.size() - 1)];
    out << ss.str() << "\t";
    
    out << it->second.modelMean << "\t";
    out << it->second.modelSd << "\t";
    out << it->second.modelValueThreshold << "\t";
    out << it->second.modelPvalueThreshold << "\t";
    
    out << it->second.loci_label;
    out << std::endl;
    
    if (it->second.loci_label == "unstable")
      unstable_count++;
    
  }
  if (sum == 0)
  {
    out << "conclusion: the sample has no msi loci which contain valued data.\n";
    out.close();
    return;
  }
  
  tmp_str = "not been decided yet";
  unstable_rate = double(unstable_count) / double(sum);
  if (unstable_count >= 0 && unstable_rate < 0.2)
  {
    tmp_str = "MSI-Stable";
  }
  if (unstable_rate >= 0.2 && unstable_rate < 0.4)
  {
    tmp_str = "MSI-Low";
  }
  if (unstable_rate <= 1 && unstable_rate >= 0.4)
  {
    tmp_str = "MSI-High";
  }
  
  out
    << "#Judge Criterion: \n#  Unstable Rate is in [0,0.2): MSI-Stable\n#  Unstable Rate is in [0.2,0.4): MSI-Low\n#  Unstable Rate is in [0.4,1]: MSI-High\n";
  out << "#Note: The -1 in \"sampleMinusLog10Pvalue\" column means infinite value.\n";
  out << "#Conclusion:\n#  Unstable loci count is: " << unstable_count << ", Total Loci Count is: " << sum
      << "\n#  The unstable rate is: " << unstable_rate << "\n#  The sample MSI feature is: " << tmp_str
      << std::endl;
  
  out.close();
}



