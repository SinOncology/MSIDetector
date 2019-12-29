/*
 * @file
 * 
 * This file provides the normalized upper incomplete gamma function gamma_q.
 * This code is extracted and adapted from the Boost library special functions.
 * See: http://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/math_toolkit/sf_gamma/igamma.html
 * The Boost licence can be found here: http://www.boost.org/LICENSE_1_0.txt
 * 
 *  Boost Software License - Version 1.0 - August 17th, 2003
 *
 *  Permission is hereby granted, free of charge, to any person or organization
 *  obtaining a copy of the software and accompanying documentation covered by
 *  this license (the "Software") to use, reproduce, display, distribute,
 *  execute, and transmit the Software, and to prepare derivative works of the
 *  Software, and to permit third-parties to whom the Software is furnished to
 *  do so, all subject to the following:
 * 
 *  The copyright notices in the Software and this entire statement, including
 *  the above license grant, this restriction and the following disclaimer,
 *  must be included in all copies of the Software, in whole or in part, and
 *  all derivative works of the Software, unless such copies or derivative
 *  works are solely in the form of machine-executable object code generated by
 *  a source language processor.
 * 
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 *  SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 *  FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 *  ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 *  DEALINGS IN THE SOFTWARE.
 * 
 */

#include <cmath>
#include <cstdint> // std::uintmax_t
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <numeric>

const double g_pi = 3.141592653589793238462643383279502884e+00;
//const float g_pi_f = 3.141592653589793238462643383279502884;
const double g_log_max_value = 709.0;
const double g_log_min_value = -708.0;
const double g_epsilon = std::numeric_limits<double>::epsilon();
const int g_digits = std::numeric_limits<double>::digits;
const std::uintmax_t g_max_iter = 1000000;
/*==============================================================================
    CLASS DECLARATION
==============================================================================*/
/*------------------------------------------------------------------------------
    lower_incomplete_gamma_series
------------------------------------------------------------------------------*/
class lower_incomplete_gamma_series
{
public :
  lower_incomplete_gamma_series(double a, double z) : a_(a), z_(z), result_(1)
  {
  }
  
  double operator()()
  {
    double r = result_;
    a_ += 1.0;
    result_ *= (z_ / a_);
    return r;
  }

private :
  double a_;
  double z_;
  double result_;
};

/*------------------------------------------------------------------------------
    small_gamma2_series
------------------------------------------------------------------------------*/
class small_gamma2_series
{
public :
  small_gamma2_series(double a, double x) : result_(-x), x_(-x), apn_(a + 1.0), n_(1)
  {
  }
  
  double operator()()
  {
    double r = (result_ / (apn_));
    result_ *= x_;
    result_ /= ++n_;
    apn_ += 1;
    return r;
  }

private :
  double result_;
  double x_;
  double apn_;
  int n_;
};

/*------------------------------------------------------------------------------
    upper_incomplete_gamma_frac
------------------------------------------------------------------------------*/
class upper_incomplete_gamma_frac
{
public :
  upper_incomplete_gamma_frac(double a, double z) : a_(a), z_(z - a + 1), k_(0)
  {
  }
  
  std::pair<double, double> operator()()
  {
    ++k_;
    z_ += 2.0;
    return std::pair<double, double>(k_ * (a_ - k_), z_);
  }

private :
  double a_;
  double z_;
  int k_;
};

/*==============================================================================
    FUNCTION DECLARATION
==============================================================================*/
double gamma_q(double a, double x);

double continued_fraction_a(class upper_incomplete_gamma_frac &g, const double &factor);

double evaluate_polynomial(const double *poly, const double &z, std::size_t count);

double finite_gamma_q(double a, double x);

double finite_half_gamma_q(double a, double x);

double lower_gamma_series(double a, double z, double init_value = 0);

double igamma_temme_large(double a, double x);

double regularised_gamma_prefix(double a, double z);

double root_epsilon();

template<class Functor>
double sum_series(Functor &func, const double &factor, std::uintmax_t &max_terms, const double &init_value);

double tgamma_small_upper_part(double a, double x, double *pgam = 0, bool invert = false);

double upper_gamma_fraction(double a, double z, double epsilon);

/*------------------------------------------------------------------------------
    gamma_q()
------------------------------------------------------------------------------*/
double gamma_q(double a, double x)
{
  bool invert = true;
  double result = 0;
  
  if (a <= 0)
  {
    std::cerr << "Error : a <= 0\n";
    exit(1);
  }
  if (x < 0)
  {
    std::cerr << "Error : x < 0\n";
    exit(1);
  }
  
  bool a_is_int = false;
  bool a_is_half_int = false;
  bool a_is_small = ((a < 30.0) && (a <= (x + 1.0)) && (x < g_log_max_value));
  
  if (a_is_small)
  {
    double floor_a = std::floor(a);
    if (floor_a == a)
    {
      a_is_int = true;
    }
    else
    {
      a_is_half_int = std::fabs(floor_a - a) == 0.5;
    }
  }
  else
  {
    a_is_int = a_is_half_int = false;
  }
  
  int method = 0;
  
  if (a_is_int && (x > 0.6))
  {
    invert = false;
    method = 0;
  }
  else if (a_is_half_int && (x > 0.2))
  {
    invert = false;
    method = 1;
  }
  else if ((x < root_epsilon()) && (a > 1.0))
  {
    method = 6;
  }
  else if (x < 0.5)
  {
    if ((-0.4 / std::log(x)) < a)
    {
      method = 2;
    }
    else
    {
      method = 3;
    }
  }
  else if (x < 1.1)
  {
    if ((x * 0.75) < a)
    {
      method = 2;
    }
    else
    {
      method = 3;
    }
  }
  else
  {
    bool use_temme = false;
    if (a > 20.0)
    {
      double sigma = std::fabs((x - a) / a);
      if ((a > 200.0) && (g_digits <= 113))
      {
        if ((20.0 / a) > (sigma * sigma))
        {
          use_temme = true;
        }
      }
      else if (g_digits <= 64)
      {
        if (sigma < 0.4)
        {
          use_temme = true;
        }
      }
    }
    if (use_temme == true)
    {
      method = 5;
    }
    else
    {
      if ((x - (1.0 / (3.0 * x))) < a)
      {
        method = 2;
      }
      else
      {
        invert = false;
        method = 4;
      }
    }
  }
  
  switch (method)
  {
  case 0:
  {
    result = finite_gamma_q(a, x);
    break;
  }
  case 1:
  {
    result = finite_half_gamma_q(a, x);
  }
  case 2:
  {
    result = regularised_gamma_prefix(a, x);
    if (result != 0)
    {
      double init_value = 0;
      bool optimised_invert = false;
      if (invert)
      {
        init_value = 1.0;
        init_value /= result;
        init_value *= -a;
        optimised_invert = true;
      }
      result *= (lower_gamma_series(a, x, init_value) / a);
      if (optimised_invert)
      {
        invert = false;
        result = -result;
      }
    }
    break;
  }
  case 3:
  {
    invert = false;
    double g = 0;
    result = tgamma_small_upper_part(a, x, &g, invert);
    result /= g;
    break;
  }
  case 4:
  {
    result = regularised_gamma_prefix(a, x);
    if (result != 0)
    {
      result *= upper_gamma_fraction(a, x, g_epsilon);
    }
    break;
  }
  case 5:
  {
    result = igamma_temme_large(a, x);
    if (x >= a)
    {
      invert = false;
    }
    break;
  }
  case 6:
  {
    result = (std::pow(x, a) / std::tgamma(a + 1.0));
    result *= (1.0 - ((a * x) / (a + 1.0)));
  }
  }
  if (result > 1.0)
  {
    result = 1.0;
  }
  if (invert)
  {
    result = (1.0 - result);
  }
  
  return result;
}

/*------------------------------------------------------------------------------
    continued_fraction_a()
------------------------------------------------------------------------------*/
double continued_fraction_a(class upper_incomplete_gamma_frac &g, const double &factor)
{
  double tiny = std::numeric_limits<double>::min();
  
  auto v = g();
  
  double f = 0;
  double C = 0;
  double D = 0;
  double delta = 0;
  double a0 = 0;
  
  f = v.second;
  a0 = v.first;
  if (f == 0)
  {
    f = tiny;
  }
  C = f;
  D = 0;
  
  auto counter = std::numeric_limits<std::uintmax_t>::max();
  
  do
  {
    v = g();
    D = (v.second + (v.first * D));
    if (D == 0)
    {
      D = tiny;
    }
    C = (v.second + (v.first / C));
    if (C == 0)
    {
      C = tiny;
    }
    D = (1 / D);
    delta = (C * D);
    f = (f * delta);
  } while ((std::fabs(delta - 1.0) > factor) && --counter);
  
  return (a0 / f);
}

/*------------------------------------------------------------------------------
    evaluate_polynomial()
------------------------------------------------------------------------------*/
double evaluate_polynomial(const double *poly, const double &z, std::size_t count)
{
  double sum = poly[count - 1];
  for (int i = (static_cast<int>(count) - 2); i >= 0; --i)
  {
    sum *= z;
    sum += poly[i];
  }
  
  return sum;
}

/*------------------------------------------------------------------------------
    finite_gamma_q()
------------------------------------------------------------------------------*/
double finite_gamma_q(double a, double x)
{
  double e = std::exp(-x);
  double sum = e;
  if (sum != 0)
  {
    double term = sum;
    for (unsigned n = 1; n < a; ++n)
    {
      term /= n;
      term *= x;
      sum += term;
    }
  }
  
  return sum;
}

/*------------------------------------------------------------------------------
    finite_half_gamma_q()
------------------------------------------------------------------------------*/
double finite_half_gamma_q(double a, double x)
{
  double e = std::erfc(std::sqrt(x));
  if ((e != 0) && (a > 1.0))
  {
    double term = (std::exp(-x) / std::sqrt(g_pi * x));
    term *= x;
    term /= 0.5;
    double sum = term;
    for (unsigned n = 2; n < a; ++n)
    {
      term /= (n - 0.5);
      term *= x;
      sum += term;
    }
    e += sum;
  }
  
  return e;
}

/*------------------------------------------------------------------------------
    lower_gamma_series()
------------------------------------------------------------------------------*/
double lower_gamma_series(double a, double z, double init_value)
{
  class lower_incomplete_gamma_series s(a, z);
  std::uintmax_t max_iter = g_max_iter;
  double factor = g_epsilon;
  double result = sum_series(s, factor, max_iter, init_value);
  if (max_iter >= g_max_iter)
  {
    std::cerr << "max_iter >= g_max_iter\n";
    exit(1);
  }
  
  return result;
}

/*------------------------------------------------------------------------------
    igamma_temme_large()
------------------------------------------------------------------------------*/
double igamma_temme_large(double a, double x)
{
  double sigma = ((x - a) / a);
  double phi = -(std::log(sigma + 1.0) - sigma);
  double y = (a * phi);
  double z = std::sqrt(2.0 * phi);
  if (x < a)
  {
    z = -z;
  }
  
  double workspace[10];
  
  static const double C0[15] =
    {
      static_cast<double>(-0.33333333333333333L   ),
      static_cast<double>( 0.083333333333333333L  ),
      static_cast<double>(-0.014814814814814815L  ),
      static_cast<double>( 0.0011574074074074074L ),
      static_cast<double>( 0.0003527336860670194L ),
      static_cast<double>(-0.00017875514403292181L),
      static_cast<double>( 0.39192631785224378e-4L),
      static_cast<double>(-0.21854485106799922e-5L),
      static_cast<double>(-0.185406221071516e-5L  ),
      static_cast<double>( 0.8296711340953086e-6L ),
      static_cast<double>(-0.17665952736826079e-6L),
      static_cast<double>( 0.67078535434014986e-8L),
      static_cast<double>( 0.10261809784240308e-7L),
      static_cast<double>(-0.43820360184533532e-8L),
      static_cast<double>( 0.91476995822367902e-9L)
    };
  workspace[0] = evaluate_polynomial(C0, z, 15);
  
  
  static const double C1[13] =
    {
      static_cast<double>(-0.0018518518518518519L ),
      static_cast<double>(-0.0034722222222222222L ),
      static_cast<double>( 0.0026455026455026455L ),
      static_cast<double>(-0.00099022633744855967L),
      static_cast<double>( 0.00020576131687242798L),
      static_cast<double>(-0.40187757201646091e-6L),
      static_cast<double>(-0.18098550334489978e-4L),
      static_cast<double>( 0.76491609160811101e-5L),
      static_cast<double>(-0.16120900894563446e-5L),
      static_cast<double>( 0.46471278028074343e-8L),
      static_cast<double>( 0.1378633446915721e-6L ),
      static_cast<double>(-0.5752545603517705e-7L ),
      static_cast<double>( 0.11951628599778147e-7L),
    };
  workspace[1] = evaluate_polynomial(C1, z, 13);
  
  static const double C2[11] =
    {
      static_cast<double>( 0.0041335978835978836L ),
      static_cast<double>(-0.0026813271604938272L ),
      static_cast<double>( 0.00077160493827160494L),
      static_cast<double>( 0.20093878600823045e-5L),
      static_cast<double>(-0.00010736653226365161L),
      static_cast<double>( 0.52923448829120125e-4L),
      static_cast<double>(-0.12760635188618728e-4L),
      static_cast<double>( 0.34235787340961381e-7L),
      static_cast<double>( 0.13721957309062933e-5L),
      static_cast<double>(-0.6298992138380055e-6L ),
      static_cast<double>( 0.14280614206064242e-6L)
    };
  workspace[2] = evaluate_polynomial(C2, z, 11);
  
  static const double C3[9] =
    {
      static_cast<double>( 0.00064943415637860082L),
      static_cast<double>( 0.00022947209362139918L),
      static_cast<double>(-0.00046918949439525571L),
      static_cast<double>( 0.00026772063206283885L),
      static_cast<double>(-0.75618016718839764e-4L),
      static_cast<double>(-0.23965051138672967e-6L),
      static_cast<double>( 0.11082654115347302e-4L),
      static_cast<double>(-0.56749528269915966e-5L),
      static_cast<double>( 0.14230900732435884e-5L)
    };
  workspace[3] = evaluate_polynomial(C3, z, 9);
  
  static const double C4[7] =
    {
      static_cast<double>(-0.0008618882909167117L ),
      static_cast<double>( 0.00078403922172006663L),
      static_cast<double>(-0.00029907248030319018L),
      static_cast<double>(-0.14638452578843418e-5L),
      static_cast<double>( 0.66414982154651222e-4L),
      static_cast<double>(-0.39683650471794347e-4L),
      static_cast<double>( 0.11375726970678419e-4L)
    };
  workspace[4] = evaluate_polynomial(C4, z, 7);
  
  static const double C5[9] =
    {
      static_cast<double>(-0.00033679855336635815L),
      static_cast<double>(-0.69728137583658578e-4L),
      static_cast<double>( 0.00027727532449593921L),
      static_cast<double>(-0.00019932570516188848L),
      static_cast<double>( 0.67977804779372078e-4L),
      static_cast<double>( 0.1419062920643967e-6L ),
      static_cast<double>(-0.13594048189768693e-4L),
      static_cast<double>( 0.80184702563342015e-5L),
      static_cast<double>(-0.22914811765080952e-5L)
    };
  workspace[5] = evaluate_polynomial(C5, z, 9);
  
  static const double C6[7] =
    {
      static_cast<double>( 0.00053130793646399222L),
      static_cast<double>(-0.00059216643735369388L),
      static_cast<double>( 0.00027087820967180448L),
      static_cast<double>( 0.79023532326603279e-6L),
      static_cast<double>(-0.81539693675619688e-4L),
      static_cast<double>( 0.56116827531062497e-4L),
      static_cast<double>(-0.18329116582843376e-4L)
    };
  workspace[6] = evaluate_polynomial(C6, z, 7);
  
  static const double C7[5] =
    {
      static_cast<double>( 0.00034436760689237767L),
      static_cast<double>( 0.51717909082605922e-4L),
      static_cast<double>(-0.00033493161081142236L),
      static_cast<double>( 0.0002812695154763237L ),
      static_cast<double>(-0.00010976582244684731L)
    };
  workspace[7] = evaluate_polynomial(C7, z, 5);
  
  static const double C8[3] =
    {
      static_cast<double>(-0.00065262391859530942L),
      static_cast<double>( 0.00083949872067208728L),
      static_cast<double>(-0.00043829709854172101L),
    };
  workspace[8] = evaluate_polynomial(C8, z, 3);
  workspace[9] = static_cast<double>(-0.00059676129019274625L);
  
  double result = evaluate_polynomial(workspace, (1.0 / a), 10.0);
  result *= (std::exp(-y) / std::sqrt(2.0 * g_pi * a));
  if (x < a)
  {
    result = -result;
  }
  result += (std::erfc(std::sqrt(y)) / 2.0);
  
  return result;
}

/*------------------------------------------------------------------------------
    regularised_gamma_prefix()
------------------------------------------------------------------------------*/
double regularised_gamma_prefix(double a, double z)
{
  double result = 0;
  
  double limit = (std::max)(10.0, a);
  double sum = (lower_gamma_series(a, limit) / a);
  sum += upper_gamma_fraction(a, limit, g_epsilon);
  
  if (a < 10.0)
  {
    double prefix = std::pow((z / 10.0), a);
    prefix *= std::exp(10.0 - z);
    if (prefix == 0)
    {
      prefix = std::pow(((z * exp((10.0 - z) / a)) / 10.0), a);
    }
    prefix /= sum;
    
    result = prefix;
  }
  else
  {
    double zoa = (z / a);
    double amz = (a - z);
    double alzoa = (a * std::log(zoa));
    double prefix = 0;
    if (((std::min)(alzoa, amz) <= g_log_min_value) || ((std::max)(alzoa, amz) >= g_log_max_value))
    {
      double amza = (amz / a);
      if ((amza <= g_log_min_value) || (amza >= g_log_max_value))
      {
        prefix = std::exp(alzoa + amz);
      }
      else
      {
        prefix = std::pow(zoa * std::exp(amza), a);
      }
    }
    else
    {
      prefix = (std::pow(zoa, a) * std::exp(amz));
    }
    prefix /= sum;
    
    result = prefix;
  }
  
  return result;
}

/*------------------------------------------------------------------------------
    root_epsilon()
------------------------------------------------------------------------------*/
double root_epsilon()
{
  double result = 0;
  
  if (g_digits == 24)
  {
    result = static_cast<double>(0.00034526698300124390839884978618400831996329879769945L);
  }
  else if (g_digits == 53)
  {
    result = static_cast<double>(0.1490116119384765625e-7L);
  }
  else if (g_digits == 64)
  {
    result = static_cast<double>(0.32927225399135962333569506281281311031656150598474e-9L);
  }
  else if (g_digits == 113)
  {
    result = static_cast<double>(0.1387778780781445675529539585113525390625e-16L);
  }
  
  return result;
}

/*------------------------------------------------------------------------------
    sum_series()
------------------------------------------------------------------------------*/
template<class Functor>
double sum_series(Functor &func, const double &factor, std::uintmax_t &max_terms, const double &init_value)
{
  auto counter = max_terms;
  
  double result = init_value;
  double next_term = 0;
  
  do
  {
    next_term = func();
    result += next_term;
  } while ((std::fabs(factor * result) < std::fabs(next_term)) && --counter);
  
  max_terms = (max_terms - counter);
  
  return result;
}

/*------------------------------------------------------------------------------
    tgamma_small_upper_part()
------------------------------------------------------------------------------*/
double tgamma_small_upper_part(double a, double x, double *pgam, bool invert)
{
  double result = (std::tgamma(a + 1.0) - 1.0);
  if (pgam)
  {
    *pgam = ((result + 1.0) / a);
  }
  double p = (std::pow(x, a) - 1.0);
  result -= p;
  result /= a;
  class small_gamma2_series s(a, x);
  std::uintmax_t max_iter = (g_max_iter - 10);
  p += 1.0;
  double init_value = (invert ? *pgam : 0);
  result = (-p * sum_series(s, g_epsilon, max_iter, ((init_value - result) / p)));
  if (max_iter >= g_max_iter)
  {
    std::cerr << "max_iter >= g_max_iter\n";
    exit(1);
  }
  if (invert)
  {
    result = -result;
  }
  
  return result;
}

/*------------------------------------------------------------------------------
    upper_gamma_fraction()
------------------------------------------------------------------------------*/
double upper_gamma_fraction(double a, double z, double epsilon)
{
  class upper_incomplete_gamma_frac f(a, z);
  double result = (1.0 / (z - a + 1.0 + continued_fraction_a(f, epsilon)));
  
  return result;
}
