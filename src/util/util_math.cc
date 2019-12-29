#include <vector>
#include <cmath>
#include <algorithm>
#include "util_math.h"

/*
 * @brief: Function to compute the median of a set of input values
 *
 * The function is called by 'median_smooth'
 *
 * @param data (std::vector<double>&): input data
 * @return: median of the set of input values
 */
double median(std::vector<double> &data)
{
  long i, j, n = data.size();
  std::vector<double> tmp;
  double median;
  
  if (data.size() > 0)
  {
    tmp.clear();
    for (i = 0; i < n; i++)
      tmp.push_back(data[i]);
    
    sort(tmp.begin(), tmp.end());
    j = tmp.size();
    if (j - 2 * (long) (j / 2) == 0)
    {
      median = 0.5 * (tmp[j / 2 - 1] + tmp[j / 2]);
    }
    else
    {
      median = tmp[j / 2];
    }
  }
  else
    median = nan("");
  
  return (median);
}


//note that the actual bandwidth of the smoother is 2*w+1
void median_smooth(std::vector<double> &data, long w)
{
  long i, j, k, l;
  long n = data.size();
  std::vector<double> sd(n, 0);
  std::vector<double> med_v(2 * w + 1, 0);
  
  for (i = 0; i < n; i++)
  {
    l = 0;
    for (j = -w; j <= w; j++)
    {
      k = i + j;
      //boundary treatment: copy end points over
      if (k < 0)
        k = (-1) * k;
      if (k >= n)
        k = n - 2 - (k - n);
      //copy to med_v for median computation
      med_v[l] = data[k];
      l++;
    }
    sd[i] = median(med_v);
  }
  
  //replace input data by median smoothed values
  data = sd;
}
