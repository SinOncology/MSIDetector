#pragma once

#include <vector>

#include "math/gamma_q.h" // gamma_q()

double median(std::vector<double> &data);

void median_smooth(std::vector<double> &data, long w);
