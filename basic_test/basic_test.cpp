///  \file basic_test.cpp
///  Basic unit tests for the standard stable distribution
/// \author Joseph Dunn
/// \copyright 2017 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <iomanip>
using std::setw;
using std::setprecision;
using std::scientific;
using std::fixed;
using std::right;
#include <fstream>
using std::ofstream;
#include <sstream>
using std::stringstream;
#include <random>
using std::mt19937;
using std::uniform_real_distribution;
#include <algorithm>
using std::sort;

#include <boost/timer/timer.hpp>
using boost::timer::auto_cpu_timer;

#include <boost/filesystem.hpp>
#include <boost/math/special_functions/airy.hpp>

#include "stable_config.h"

#define MPREAL
#define CPP_BIN_FLOAT
#define MPFR_FLOAT
#include "stable_distribution.h"
#include "kolmogorov.h"
#define LIBRARY
#include "zolotarev.h"
#undef LIBRARY
#include "bracket_and_solve.h"

using namespace stable_distribution;

template<typename myFloat>
myFloat epsdiff(myFloat r1, myFloat r2) {
  if (r1 == r2)
    return 0;
  else if (!boost::math::isfinite(r1) || !boost::math::isfinite(r2))
    return std::numeric_limits<myFloat>::infinity();
  else {
    myFloat one{1};
    return fabs(r1-r2)/(max(one, max(fabs(r1), fabs(r2)))*std::numeric_limits<myFloat>::epsilon());
  }
}

template<typename myFloat>
string fmt_eps(myFloat eps) {
  stringstream ss;
  if (eps < 10000)
    ss << setw(10) << setprecision(1) << fixed << eps;
  else
    ss << setw(10) << right << "*******" ;
  return ss.str();
}

template<typename myFloat>
myFloat stable_levy_cdf(myFloat x, bool lower_tail=true, bool log_p=false) {
  // pm = S1 representation
  myFloat c= 1;
  x = (x)/c;
  if (x <=0)
    return  lower_tail ? (log_p ? StandardStableDistribution<myFloat>::NegInf : 0)
    : (log_p ? 0 : 1);
  if (x == StandardStableDistribution<myFloat>::PosInf)
    return lower_tail ? (log_p ? 0 : 1)
    : (log_p ? StandardStableDistribution<myFloat>::NegInf : 0);
  myFloat u = 1/sqrt(2*x);
  if(log_p) {
    if(lower_tail)
      return log(my_erfc(u));
    else
      return log(erf(u));
  } else {
    if(lower_tail)
      return my_erfc(u);
    else
      return erf(u);
  }
}

template<typename myFloat>
int test_stable_cdf(ostream& out, Controllers<myFloat> ctls) {
  auto_cpu_timer timer(out);
  int verbose = 0;
  
  bool log_p=true, lower_tail=false;
  myFloat alpha = .5;
  myFloat beta = 1;
  
  vector<myFloat> xs; // this is x - zeta
  xs.push_back(0.);
  
  for (int i=1; i<100; i++)
    xs.push_back(i/1000.);
  for (int i=1; i<=100; i++) {
    xs.push_back(i/10.);
  }
  for (int i=2; i<100; i++) {
    xs.push_back(pow(10.,i));
  }
  ;
  xs.push_back(pow(10.,120));
  xs.push_back(pow(10.,150));
  xs.push_back(pow(10.,200));
  xs.push_back(std::numeric_limits<myFloat>::infinity());
  
  Fmt<myFloat> fmt;
  out << endl;
  out << "Comparison of cdf to Levy formula for alpha = .5, beta = 1, digits10 = "
  << fmt.digits10 << endl << endl;
  StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctls, verbose);
  Zolotarev<myFloat> zol(alpha, beta, &ctls.controller, verbose);
  out << setw(8) << right << "alpha"
  << setw(8) << right << "beta"
  << setw(13) << right << "x"
  << setw(8) << right << "tail"
  << setw(fmt.width) << right << "log(pLevy)"
  << setw(fmt.width) << right << "log(cdf)"
  << setw(10) << right << "epsdiff"
  << setw(15) << right << "abserr"
  << setw(4) << right << "tc"
  << setw(7) << right << "neval"
  << setw(fmt.width) << right << "log(zol_cdf)"
  << setw(10) << right << "epsdiff"
  << setw(2) << right << "T"
  << endl << endl;
  bool pass = true;
  for (auto x : xs) {
    lower_tail = true;
    myFloat r_Levy = stable_levy_cdf<myFloat>(x, lower_tail, log_p);
    if (r_Levy > (log_p ? log(.5) : .5)) {
      lower_tail = false;
      r_Levy = stable_levy_cdf<myFloat>(x, lower_tail, log_p);
    }
    myFloat r = std_stable_dist.cdf(x, lower_tail, log_p, S1);
    myFloat eps = epsdiff(r, r_Levy);
    myFloat r_zol = log(zol.cdf(x, lower_tail, S1));
    myFloat eps_zol = epsdiff(r_zol, r_Levy);
    bool pass1 = !boost::math::isnan(eps) && boost::math::isfinite(eps)
                 && eps < 100;

    out << setw(8) << setprecision(3) << fixed << myFloat(alpha)
    << setw(8) << setprecision(3) << fixed << myFloat(beta)
    << setw(13) << setprecision(3) << scientific << myFloat(x)
    << setw(8) << right << ((lower_tail) ? "lower" : "upper")
    << fmt << r_Levy
    << fmt << r
    << setw(10) << setprecision(1) << fixed << eps
    << setw(15) << setprecision(5) << scientific << std_stable_dist.abserr
    << setw(4) << std_stable_dist.termination_code
    << setw(7) << std_stable_dist.neval
    << fmt << r_zol
    << setw(10) << fmt_eps(eps_zol)
    << setw(2) << (zol.result_type==asymptotic ? "A" : "C")
    << (pass1 ? "" : " FAIL") << endl;
    pass = pass && pass1;
  }
  out << endl << (pass?"Test Passed":"Test Failed") << endl;
  return !pass;
}

template <typename myFloat>
myFloat Ai(myFloat x) {
  if (x < .1) {
    const myFloat two{2}, three{3};
    const myFloat ai_0 = 1/(pow(three,two/three)*tgamma(two/three));
    const myFloat ai_prime_0 = -1/(pow(three, 1/three) * tgamma(1/three));
    myFloat term1 = ai_0;
    myFloat term2 = ai_prime_0 * x;
    myFloat sum1 = term1;
    myFloat sum2 = term2;
    myFloat tol = std::numeric_limits<myFloat>::epsilon();
    myFloat x_cubed = x*x*x;
    for (int k = 1; fabs(term1)> tol*fabs(sum1) || fabs(term2)>tol*fabs(sum2); ++k) {
      term1 *= x_cubed/((3*k)*(3*k-1));
      sum1 += term1;
      term2 *= x_cubed/((3*k+1)*(3*k));
      sum2 += term2;
    }
    return sum1+sum2;
  } else {
    return boost::math::airy_ai(x);
  }
}

template <>
mpreal Ai<mpreal>(mpreal x) {
  // Boost special function determine precision at compile time
  // MpfrFloat precision is defined at compile time. mpreal is not.
  MpfrFloat tmp{x.mpfr_ptr()};
  return mpreal(Ai(tmp).backend().data());
}

template <typename myFloat>
myFloat Ai_prime(myFloat x) {
  if (fabs(x) < .1) {
    const myFloat two{2}, three{3};
    const myFloat ai_0 = 1/(pow(three,two/three)*tgamma(two/three));
    const myFloat ai_prime_0 = -1/(pow(three, 1/three) * tgamma(1/three));
    myFloat term1 = ai_prime_0;
    myFloat term2 = ai_0 * x * x/2;
    myFloat sum1 = term1;
    myFloat sum2 = term2;
    myFloat tol = std::numeric_limits<myFloat>::epsilon();
    myFloat x_cubed = x*x*x;
    for (int k = 1; fabs(term1)> tol*fabs(sum1) || fabs(term2)>tol*fabs(sum2); ++k) {
      term1 *= x_cubed/((3*k)*(3*k-2));
      sum1 += term1;
      term2 *= x_cubed/((3*k+2)*(3*k));
      sum2 += term2;
    }
    return sum1+sum2;

  } else {
    return boost::math::airy_ai_prime(x);
  }
}

template <>
mpreal Ai_prime<mpreal>(mpreal x) {
  MpfrFloat tmp{x.mpfr_ptr()};
  return mpreal(Ai_prime(tmp).backend().data());
}

template<typename myFloat>
myFloat stable_taleb_pdf(myFloat x, bool log_flag=false){
  // stable for alpha=3/2, beta=1, pm=1
  // From Taleb Silent Risk page 50 with slight modification (x=-x)
  const myFloat one{1}, two{2}, three{3}, four{4};
  const myFloat denom1 = pow(three,four/three)*pow(two, two/three);
  const myFloat factor1 = -pow(two, one/three) * pow(three , -four/three);
  const myFloat factor2 = -pow(two, two/three) * pow(three, -two/three);
  myFloat arg = pow(x,2)/denom1;
  if (x < 100) {
    myFloat res = exp(pow(x,3)/27)*(factor1 * x * Ai(arg)+ factor2 * Ai_prime(arg));
    return log_flag ? log(res) : res;
  } else {
    // For x large both Ai and Ai_prime contain of factor exp(-(x^3)/27)
    // Which underflows at the same time that the exp((x^3)/27) overflows
    // We'll go directly to the asymptotic expansion of Ai and Ai_prime
    myFloat airy_zeta = pow(x,3)/27;
    myFloat pi = const_pi<myFloat>();
    myFloat u = 1;
    myFloat v = 1;
    myFloat sum_Ai = 0;
    myFloat sum_Ai_prime = 0;
    myFloat zeta_k = 1;
    myFloat old_term=1;
    for (int k = 1; k<100; ++k) {
      u = u * (6*k-5)*(6*k-3)*(6*k-1)/((2*k-1)*216*k);
      v = u * (6*k+1)/(1-6*k);
      zeta_k *= -1/airy_zeta;
      myFloat term = u * zeta_k;
      if (fabs(term) > fabs(old_term) || fabs(term) < std::numeric_limits<myFloat>::min()) break;
      sum_Ai += term;
      sum_Ai_prime += v * zeta_k;
      old_term = term;
    }
    myFloat res = factor1 * x / (2 * sqrt(pi) * pow(arg,one/four))*sum_Ai
                 -factor2 * pow(arg,one/four)/(2 * sqrt(pi))*sum_Ai_prime;
    return log_flag ? log(res) : res;
  }
}

template<typename myFloat>
myFloat stable_levy_pdf(myFloat x, bool log_flag=false) {
  // pm=1 representation, zeta = -beta*tan(pi/4) = -1
  myFloat c = 1;
  x = x;
  myFloat pi = const_pi<myFloat>();
  // ensure f(0) = 0 {not NaN}:
  if (x <= 0 || x==std::numeric_limits<myFloat>::infinity())
    return log_flag ? -std::numeric_limits<myFloat>::infinity() : 0;
  else return log_flag ? (log(c/(2*pi)) + -c/x - 3*log(x))/2
    : sqrt(c/(2*pi)) * exp(-c/(2*x)) / pow(x,(1.5));
}

template<typename myFloat>
int test_stable_pdf(ostream& out, Controllers<myFloat> ctls) {
  auto_cpu_timer timer(out);
  myFloat alpha = .5;
  myFloat beta = 1;
  vector<myFloat> xs;
  xs.push_back(0.);
  for (int i=1; i<100; i++) {
    xs.push_back(i/1000.);
  }
  for (int i=1; i<=100; i++) {
    xs.push_back(i/10.);
  }
  for (int i=2; i<100; i++) {
    xs.push_back(pow(10.,i));
  }
  
  xs.push_back(pow(10.,120));
  xs.push_back(pow(10.,150));
  xs.push_back(pow(10.,200));
  xs.push_back(std::numeric_limits<myFloat>::infinity());
  
  int verbose = 0;
  Fmt<myFloat> fmt;
  out << endl;
  out << "Comparison of pdf to Levy formula for alpha = .5, beta = 1, digits10 = "
  << fmt.digits10 << endl << endl;
  StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctls, verbose);
  Zolotarev<myFloat> zol(alpha, beta, &ctls.controller, verbose);
  int log_flag=1;
  out << setw(13) << right << "alpha"
  << setw(13) << right << "beta"
  << setw(13) << right << "x"
  << setw(fmt.width) << right << "log(dLevy)"
  << setw(fmt.width) << right << "log(pdf_myFloat)"
  << setw(10) << right << "epsdiff"
  << setw(15) << right << "abserr"
  << setw(4) << right << "TC"
  << setw(7) << right << "neval"
  << setw(fmt.width) << right << "log(zol_pdf)"
  << setw(10) << right << "epsdiff"
  << setw(2) << right << "T"
  << endl << endl;
  bool pass = true;
  for (auto x : xs) {
    myFloat r_Levy = stable_levy_pdf<myFloat>(x, log_flag);
    myFloat r = std_stable_dist.pdf(x, log_flag, S1);
    myFloat eps = epsdiff(r, r_Levy);
    myFloat r_zol = log(zol.pdf(x, S1));
    myFloat eps_zol = epsdiff(r_zol, r_Levy);
    bool pass1 = !boost::math::isnan(eps) && boost::math::isfinite(eps)
    && eps < 100;

    out << setw(13) << setprecision(5) << fixed << myFloat(alpha)
    << setw(13) << setprecision(5) << fixed << myFloat(beta)
    << setw(13) << setprecision(3) << scientific << myFloat(x)
    << fmt << r_Levy
    << fmt << r
    << setw(10) << setprecision(1) << fixed << eps
    << setw(15) << setprecision(5) << scientific << std_stable_dist.abserr
    << setw(4) << std_stable_dist.termination_code
    << setw(7) << std_stable_dist.neval
    << fmt << r_zol
    << setw(10) << fmt_eps(eps_zol)
    << setw(2) << (zol.result_type==asymptotic ? "A" : "C")
    << (pass1 ? "" : " FAIL") << endl;
    pass = pass && pass1;

  }
  out << endl;
  out << "Comparison of pdf to Taleb formula for alpha = 1.5, beta = 1, digits10 = "
  << fmt.digits10 << endl << endl;
  alpha = 1.5;
  StandardStableDistribution<myFloat> stable_1_pt_5(alpha, beta, ctls, verbose);
  Zolotarev<myFloat> zol_1_pt_5(alpha, beta, &ctls.controller, verbose);
  out << setw(13) << right << "alpha"
  << setw(13) << right << "beta"
  << setw(13) << right << "x"
  << setw(fmt.width) << right << "log(pdf_Taleb)"
  << setw(fmt.width) << right << "log(pdf_myFloat)"
  << setw(10) << right << "epsdiff"
  << setw(15) << right << "abserr"
  << setw(4) << right << "TC"
  << setw(7) << right << "neval"
  << setw(fmt.width) << right << "log(zol_pdf)"
  << setw(10) << right << "epsdiff"
  << setw(2) << right << "T"
  << endl << endl;
  vector<myFloat> xs2;
  for (auto x : xs) {
    if (x <= 1e100) {
      xs2.push_back(-x);
      xs2.push_back(x);
    }
  }
  sort(xs2.begin(), xs2.end());
  
  for (auto x : xs2) {
    myFloat r_taleb = stable_taleb_pdf<myFloat>(x, log_flag);
    myFloat r = stable_1_pt_5.pdf(x, log_flag, S1);
    myFloat eps = epsdiff(r, r_taleb);
    myFloat r_zol = log(zol_1_pt_5.pdf(x, S1));
    myFloat eps_zol = epsdiff(r_zol, r_taleb);
    bool pass1 = !boost::math::isnan(eps) && boost::math::isfinite(eps)
    && eps < 300;
    
    out << setw(13) << setprecision(5) << fixed << myFloat(alpha)
    << setw(13) << setprecision(5) << fixed << myFloat(beta)
    << setw(13) << setprecision(3) << scientific << myFloat(x)
    << fmt << r_taleb
    << fmt << r
    << setw(10) << setprecision(1) << fixed << eps
    << setw(15) << setprecision(5) << scientific << stable_1_pt_5.abserr
    << setw(4) << stable_1_pt_5.termination_code
    << setw(7) << stable_1_pt_5.neval
    << fmt << r_zol
    << setw(10) << fmt_eps(eps_zol)
    << setw(2) << (zol_1_pt_5.result_type==asymptotic ? "A" : "C")
    << (pass1 ? "" : " FAIL") << endl;
    pass = pass && pass1;
    
  }
  out << endl << (pass?"Test Passed":"Test Failed") << endl;
  return !pass;
} //test_stable_pdf

template<typename myFloat>
myFloat stable_levy_ddx_pdf(myFloat x) {
  // pm=1 representation
  myFloat c = 1;
  x = x;
  myFloat pi = const_pi<myFloat>();
  // ensure f(0) = 0 {not NaN}:
  if (x <= 0 || x==std::numeric_limits<myFloat>::infinity())
    return 0;
  else
    return sqrt(c/(2*pi)) * exp(-c/(2*x)) * (c/(2*x)-myFloat(3)/2)/pow(x,2.5);
}

template<typename myFloat>
int test_stable_ddx_pdf(ostream& out, Controllers<myFloat> ctls) {
  auto_cpu_timer timer(out);
  myFloat alpha = .5;
  myFloat beta = 1;
  vector<myFloat> xs;
  xs.push_back(0.);
  for (int i=1; i<100; i++) {
    xs.push_back(i/1000.);
  }
  for (int i=1; i<=100; i++) {
    xs.push_back(i/10.);
  }
  for (int i=2; i<100; i++) {
    xs.push_back(pow(10.,i));
  }
  
  xs.push_back(pow(10.,120));
  xs.push_back(pow(10.,150));
  xs.push_back(pow(10.,200));
  xs.push_back(std::numeric_limits<myFloat>::infinity());
  
  int verbose = 0;
  Fmt<myFloat> fmt;
  out << endl;
  out << "Comparison of ddx_pdf to Levy formula for alpha = .5, beta = 1, digits10 = "
  << fmt.digits10 << endl << endl;
  StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctls, verbose);
  Zolotarev<myFloat> zol(alpha, beta, &ctls.controller, verbose);
  out << setw(13) << right << "alpha"
  << setw(13) << right << "beta"
  << setw(13) << right << "x"
  << setw(fmt.width) << right << "Levy_ddx_pdf"
  << setw(fmt.width) << right << "ddx_pdf_myFloat"
  << setw(10) << right << "epsdiff"
  << setw(15) << right << "abserr"
  << setw(4) << right << "termination_code"
  << setw(7) << right << "neval"
  << setw(fmt.width) << right << "zol_ddx_pdf)"
  << setw(10) << right << "epsdiff"
  << setw(2) << right << "T"
  << endl << endl;
  bool pass = true;
  for (auto x : xs) {
    myFloat r_Levy = stable_levy_ddx_pdf<myFloat>(x);
    myFloat r = std_stable_dist.ddx_pdf(x, S1);
    myFloat eps = epsdiff(log(fabs(r)), log(fabs(r_Levy)));
    myFloat r_zol = zol.ddx_pdf(x,S1);
    myFloat eps_zol = epsdiff(log(fabs(r_zol)), log(fabs(r_Levy)));
    bool pass1 = !boost::math::isnan(eps) && boost::math::isfinite(eps)
    && eps < 100;

    out << setw(13) << setprecision(5) << fixed << myFloat(alpha)
    << setw(13) << setprecision(5) << fixed << myFloat(beta)
    << setw(13) << setprecision(3) << scientific << myFloat(x)
    << fmt << r_Levy
    << fmt << r
    << setw(10) << setprecision(1) << fixed << eps
    << setw(15) << setprecision(5) << scientific << std_stable_dist.abserr
    << setw(4) << std_stable_dist.termination_code
    << setw(7) << std_stable_dist.neval
    << fmt << r_zol
    << setw(10) << fmt_eps(eps_zol)
    << setw(2) << (zol.result_type==asymptotic ? "A" : "C")
    << (pass1 ? "" : " FAIL") << endl;
    pass = pass && pass1;
  }
  out << endl << (pass?"Test Passed":"Test Failed") << endl;
  return !pass;
} //test_stable_ddx_pdf

template<typename myFloat>
myFloat stable_levy_quantile(myFloat pp, bool lower_tail, bool log_p) {
  // pm = S1 representation
  myFloat c= 1;
  myFloat p = log_p ? exp(pp) : pp;
  if (p == 0)
    return  lower_tail ? StandardStableDistribution<myFloat>::NegInf
    : StandardStableDistribution<myFloat>::PosInf;
  if (p == 1)
    return lower_tail ? StandardStableDistribution<myFloat>::PosInf
    : StandardStableDistribution<myFloat>::NegInf;
  
  myFloat u = lower_tail ? erfc_inv(p)
                        : erf_inv(p);
  return c/(2*u*u);
}

template<typename myFloat>
int test_stable_quantile(ostream& out, Controllers<myFloat> ctls) {
  auto_cpu_timer timer(out);
  int verbose = 0;
  
  bool log_p=false;
  myFloat q_tol = 64 * std::numeric_limits<myFloat>::epsilon();
  myFloat alpha = .5;
  myFloat beta = 1;
  
  vector<myFloat> ps;
  vector<int> lower_tails;
  ps.push_back(0);
  lower_tails.push_back(1);
  for (int i = -6; i<-1; ++i) {
    ps.push_back(pow(static_cast<myFloat>(10),i));
    lower_tails.push_back(1);
  }
  for (int i = 1; i<=5; ++i) {
    ps.push_back(static_cast<myFloat>(i)/10);
    lower_tails.push_back(1);
  }
  for (int i = 4; i>=1; --i) {
    ps.push_back(static_cast<myFloat>(i)/10);
    lower_tails.push_back(0);
  }
  for (int i = -2; i>=-6; --i) {
    ps.push_back(pow(static_cast<myFloat>(10),i));
    lower_tails.push_back(0);
  }
  ps.push_back(0);
  lower_tails.push_back(0);
  
  out << endl;
  Fmt<myFloat> fmt;
  out << "Comparison of quantile to Levy formula for alpha = .5, beta = 1, digits10 = "
  << fmt.digits10 << endl << endl;
  StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctls, verbose);
  out << setw(8) << right << "alpha"
  << setw(8) << right << "beta"
  << setw(13) << right << "p"
  << setw(8) << right << "tail"
  << setw(fmt.width) << right << "Levy_quantile"
  << setw(fmt.width) << right << "quantile"
  << setw(10) << right << "epsdiff"
  << setw(7) << right << "iters"
  << setw(7) << right << "neval"
  << endl << endl;
  bool pass = true;
  for (int i=0; i<ps.size(); ++i) {
    myFloat p = ps.at(i);
    int lower_tail = lower_tails.at(i);
    myFloat r_Levy = stable_levy_quantile<myFloat>(p, lower_tail, log_p);
    myFloat r = std_stable_dist.quantile(p, lower_tail, log_p,
                                         q_tol, S1);
    myFloat eps = epsdiff(r, r_Levy);
    bool pass1 = !boost::math::isnan(eps) && boost::math::isfinite(eps)
    && eps < 100;

    out << setw(8) << setprecision(3) << fixed << myFloat(alpha)
    << setw(8) << setprecision(3) << fixed << myFloat(beta)
    << setw(13) << setprecision(3) << scientific << myFloat(p)
    << setw(8) << right << ((lower_tail) ? "lower" : "upper")
    << fmt << r_Levy
    << fmt << r
    << setw(10) << setprecision(1) << fixed << eps
    << setw(7) << std_stable_dist.num_iter
    << setw(7) << std_stable_dist.neval
    << (pass1 ? "" : " FAIL") << endl;
    pass = pass && pass1;
  }
  out << endl << (pass?"Test Passed":"Test Failed") << endl;
  return !pass;
}

template<typename myFloat>
void print_stable_mode_heading(ostream& os, const vector<myFloat>& betas) {
  os << setw(65) << right << "beta" << endl
  << setw(20) << right << "alpha";
  for (auto beta : betas)
    os << setw(30) << setprecision(11) << fixed <<  beta;
  os << endl;
  
}

template<typename myFloat>
int test_stable_mode(ostream& out, Controllers<myFloat> ctls) {
  auto_cpu_timer timer(out);
  myFloat pi2 = const_pi<myFloat>()/2;
  out << endl;
  out << "Test of stable_mode" << endl << endl;
  
  vector<myFloat> alphas;
  for (int i=8; i>=2; i--)
    alphas.push_back(pow(static_cast<myFloat>(10), -i));
  for (int i=1; i<10; ++i)
    alphas.push_back(i/static_cast<myFloat>(10));
  for (int i=2; i<=4; i+=1)
    alphas.push_back(1-pow(static_cast<myFloat>(10),-i));
  alphas.push_back(1.);
  for (int i=4; i>=2; --i)
    alphas.push_back(1+pow(static_cast<myFloat>(10),-i));
  for (int i=1; i<10; i++)
    alphas.push_back(1+(i/static_cast<myFloat>(10)));
  for (int i=2; i<=4; i++)
    alphas.push_back(2 - pow(static_cast<myFloat>(10), -i));
  alphas.push_back(2);
  
  vector<myFloat> betas;
  for (int i=0; i<=2; i++) betas.push_back(i/static_cast<myFloat>(2));
  myFloat ddx_tol = 64*std::numeric_limits<myFloat>::epsilon();
  //vector<myFloat> out(alphas.size()*betas.size());
  
  out << setw(75) << right
  << "Position of Mode in the S1 Parameterization" << endl << endl;
  print_stable_mode_heading(out, betas);
  
  int verbose_mode = 0;
  int verbose = 0;
  auto pm = S1;
  vector<myFloat> mode_S1;
  vector<myFloat> pdf_at_mode;
  myFloat first_good_alpha = 2;
  bool pass = true;
  for (auto alpha : alphas) {
    out << setw(20) << setprecision(14) << fixed << alpha;
    if (lgamma(1+1/alpha) > log(std::numeric_limits<myFloat>::max())) {
      out << setw(30) << "alpha is too small." << endl;
      continue;
    } else {
      first_good_alpha = min(first_good_alpha, alpha);
    }
    bool pass1 = true;
    for (auto beta : betas) {
      StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctls, verbose);
      std::pair<myFloat, myFloat> mode = std_stable_dist.mode(ddx_tol, verbose_mode, pm);
      out << setw(30) << setprecision(6) << scientific << mode.first;
      mode_S1.push_back(mode.first);
      pdf_at_mode.push_back(mode.second);
      pass1 = pass1 && !boost::math::isnan(mode.first) && !boost::math::isnan(mode.second);
    }
    out << (pass1 ? "" : " FAIL") << endl;
    pass = pass && pass1;
  }
  out << endl << setw(75) << right
  << "Position of Mode in the S0 Parameterization" << endl << endl;
  print_stable_mode_heading(out, betas);
  
  int j = 0;
  for (myFloat alpha : alphas) {
    out << setw(20) << setprecision(14) << fixed << alpha;
    if (alpha < first_good_alpha) {
      out << setw(30) << "alpha is too small." << endl;
      continue;
    }
    for (myFloat beta : betas) {
      myFloat zeta = (alpha!=1) ? -beta * tan(pi2 * alpha) : 0;
      out << setw(30) << setprecision(6) << scientific << mode_S1.at(j++)+zeta;
    }
    out << endl;
  }
  
  out << endl << setw(65) << right
  << "Probability Density Function at Mode"<< endl << endl;
  print_stable_mode_heading(out, betas);
  j = 0;
  for (myFloat alpha : alphas) {
    out << setw(20) << setprecision(14) << fixed << alpha;
    if (alpha < first_good_alpha) {
      out << setw(30) << "alpha is too small." << endl;
      continue;
    }
    for (int i = 0; i< betas.size(); ++i)
      out << setw(30) << setprecision(6) << scientific << pdf_at_mode.at(j++);
    out << endl;
  }
  out << endl << (pass?"Test Passed":"Test Failed") << endl;
  return !pass;
}

template<typename myFloat>
/// Calculated the D for Kolomogorov Smirnov Test
myFloat D(vector<myFloat>& data, StandardStableDistribution<myFloat>& dist) {
  sort(data.begin(), data.end());
  int n = int(data.size());
  myFloat ret{0};
  for (int i=0; i<n; ++i) {
    int lower_tail{1};
    int log_flag{0};
    myFloat F = dist.cdf(data[i], lower_tail, log_flag);
    ret= std::max(ret, std::max(F -myFloat(i)/n, (myFloat(i)+1)/n -F));
  }
  return ret;
}

template<typename myFloat>
int test_stable_random (ostream& out, int n, myFloat alpha, myFloat beta,
                    Controllers<myFloat> ctls) {
  auto_cpu_timer timer(out);
  mt19937 gen(200);
  uniform_real_distribution<> dis;
  vector<myFloat> r(n);
  for (int i=0; i<n; i++) {
    myFloat u1 = dis(gen);
    myFloat u2 = dis(gen);
    r.at(i) = random_stable(alpha, beta, u1, u2);
  }
  sort(r.begin(),r.end());
  
  int verbose = 0;
  
  StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctls, verbose);
  
  int log_p=0;
  
  myFloat sum_neginf = 0, sum_posinf = 0, sum_zero = 0, sum_nan = 0;
  for (typename vector<myFloat>::iterator pr=r.begin(); pr < r.end(); pr++) {
    sum_neginf+=(*pr < -1e300);
    sum_posinf+=(*pr > 1e300);
    sum_zero+=(fabs(*pr)<.01);
    sum_nan+=boost::math::isnan(*pr);
  }
  
  out << endl;
  int digits10 = static_cast<int>(-log(std::numeric_limits<myFloat>::epsilon())/log(myFloat(10)));
  out << "Call to random with n = " << n << ", alpha = " << alpha
  << ", beta = " << beta << ", digits10 = " << digits10
  << endl << endl;
  out << "p_neginf = " << setw(15) << sum_neginf/n
      << ", cdf = "<<  setw(15) << std_stable_dist.cdf(-1e300,true,log_p) << endl
      << "p_posinf = " << setw(15) << sum_posinf/n
      << ", cdf = " << setw(15) << std_stable_dist.cdf(1e300,false,log_p) << endl
      << "p_zero   = "  << setw(15) << sum_zero/n
      << ", cdf = "<< setw(15) << std_stable_dist.cdf(.01,true, log_p)-std_stable_dist.cdf(-.01,true, log_p) << endl
      << "p_nan    = " << setw(15) << sum_nan/n << endl
      << "p_other  = " << setw(15) << (n-sum_neginf-sum_posinf-sum_zero-sum_nan)/n << endl;
  
  out << "Kolmogorov Smirnov Test" << endl << endl;
  
  double d = static_cast<double>(D(r, std_stable_dist));
  
  double one_minus_k = 1-kolmogorov_cdf(n, d);
  bool pass = one_minus_k > .01;
  out << "D = " << d << ", KS Probability = " << kolmogorov_asymptotic_cdf(d * sqrt(n))
      << ", 1 - K(n,d) = " << one_minus_k << (pass ? "" : " FAIL") << endl;
  out << endl << (pass?"Test Passed":"Test Failed") << endl;
  return !pass;
}

static void show_usage (string name){
  boost::filesystem::path p(name);
  cerr << "Usage: " << p.filename().string() << "test_name" << endl;
}

int main(int argc, const char * argv[]) {
  // Check the number of parameters
  if (argc != 2) {
    // Tell the user how to run the program
    show_usage(string(argv[0]));
    return 1;
  }
  
  string test_name = string(argv[1]);
  // First create some high precision coeficients for
  mpreal::set_default_prec(128);
  Kronrod<mpreal> k_big(10);
  mpreal::set_default_prec(96);
  int noext = 1;
  double epsabs_double = 0;
  double epsrel_double = 64 * std::numeric_limits<double>::epsilon();
  int limit = 1000;
  int verbose_integration = 0;
  IntegrationController<double> ctl_double(noext, k_big,
                                           epsabs_double, epsrel_double,
                                           limit, verbose_integration);
  if (verbose_integration)
    cout << "ctl_double:" << endl << ctl_double << endl;
  Controllers<double> ctls_double(ctl_double, ctl_double);
  
#ifdef MPREAL
  mpreal epsabs_mpreal = 0;
  mpreal epsrel_mpreal = 64 * std::numeric_limits<mpreal>::epsilon();
  IntegrationController<mpreal> ctl_mpreal(noext, k_big,
                                           epsabs_mpreal, epsrel_mpreal,
                                           limit, verbose_integration);
  if (verbose_integration)
    cout << "ctl_mpreal:" << endl << ctl_mpreal << endl;
  Controllers<mpreal> ctls_mpreal(ctl_mpreal, ctl_double);
#endif
  
#ifdef CPP_BIN_FLOAT
  Kronrod<BigCppBinFloat> k_big_cpp_bin_float(10);
  CppBinFloat epsabs_cpp_bin_float = 0;
  CppBinFloat epsrel_cpp_bin_float = 64 * std::numeric_limits<CppBinFloat>::epsilon();
  IntegrationController<CppBinFloat> ctl_cpp_bin_float(noext, k_big_cpp_bin_float,
                                           epsabs_cpp_bin_float, epsrel_cpp_bin_float,
                                           limit, verbose_integration);
  if (verbose_integration)
    cout << "ctl_cpp_bin_float:" << endl << ctl_cpp_bin_float << endl;
  Controllers<CppBinFloat> ctls_cpp_bin_float(ctl_cpp_bin_float, ctl_double);
#endif
  
#ifdef MPFR_FLOAT
  Kronrod<BigMpfrFloat> k_big_mpfr_float(10);
  MpfrFloat epsabs_mpfr_float = 0;
  MpfrFloat epsrel_mpfr_float = 64 * std::numeric_limits<MpfrFloat>::epsilon();
  IntegrationController<MpfrFloat> ctl_mpfr_float(noext, k_big_mpfr_float,
                                                       epsabs_mpfr_float, epsrel_mpfr_float,
                                                       limit, verbose_integration);
  if (verbose_integration)
    cout << "ctl_mpfr_float:" << endl << ctl_mpfr_float << endl;
  Controllers<MpfrFloat> ctls_mpfr_float(ctl_mpfr_float, ctl_double);
#endif

  string out_dir = string(OUT_DIR);
  if (!boost::filesystem::is_directory(out_dir))
    boost::filesystem::create_directory(out_dir);
  bool fail = false;
  if (test_name == "cdf_double") {
    string outfile = out_dir + "/test_stable_cdf_double.out"; 
    cout << "Writing output to " << outfile << endl;
    ofstream cdf_double(outfile);
    fail = test_stable_cdf<double>(cdf_double, ctls_double);
  } else if (test_name == "pdf_double") {
    string outfile = out_dir + "/test_stable_pdf_double.out"; 
    cout << "Writing output to " << outfile << endl;
    ofstream pdf_double(outfile);
    fail = test_stable_pdf<double>(pdf_double, ctls_double);
  } else if (test_name == "ddx_pdf_double") {
    string outfile = out_dir + "/test_stable_ddx_pdf_double.out"; 
    cout << "Writing output to " << outfile << endl;
    ofstream ddx_pdf_double(outfile);
    fail = test_stable_ddx_pdf<double>(ddx_pdf_double, ctls_double);
  } else if (test_name == "quantile_double") {
    string outfile = out_dir + "/test_stable_quantile_double.out"; 
    cout << "Writing output to " << outfile << endl;
    ofstream quantile_double(outfile);
    fail = test_stable_quantile<double>(quantile_double, ctls_double);
  } else if (test_name == "mode_double") {
    string outfile = out_dir + "/test_stable_mode_double.out"; 
    cout << "Writing output to " << outfile << endl;
    ofstream mode_double(outfile);
    fail = test_stable_mode<double>(mode_double, ctls_double);
  } else if (test_name == "random_double") {
    string outfile = out_dir + "/test_stable_random_double.out"; 
    cout << "Writing output to " << outfile << endl;
    ofstream random_double(outfile);
    fail = test_stable_random<double>(random_double, 10000, .5, .5, ctls_double);
    fail = test_stable_random<double>(random_double, 10000, 1.-1./128., .5, ctls_double) || fail;
    fail = test_stable_random<double>(random_double, 10000, 1, .5, ctls_double) || fail;
    fail = test_stable_random<double>(random_double, 10000, 1.+1./128., .5, ctls_double) || fail;
    fail = test_stable_random<double>(random_double, 10000, 1.5, .5, ctls_double) || fail;
  }
#ifdef MPREAL
  else if (test_name == "cdf_mpreal") {
    string outfile = out_dir + "/test_stable_cdf_mpreal.out"; 
    cout << "Writing output to " << outfile << endl;
    ofstream cdf_mpreal(outfile);
    fail = test_stable_cdf<mpreal>(cdf_mpreal,  ctls_mpreal);
  } else if (test_name == "pdf_mpreal") {
    string outfile = out_dir + "/test_stable_pdf_mpreal.out"; 
    cout << "Writing output to " << outfile << endl;
    ofstream pdf_mpreal(outfile);
    fail = test_stable_pdf<mpreal>(pdf_mpreal,  ctls_mpreal);
  } else if (test_name == "ddx_pdf_mpreal") {
    string outfile = out_dir + "/test_stable_ddx_pdf_mpreal.out"; 
    cout << "Writing output to " << outfile << endl;
    ofstream ddx_pdf_mpreal(outfile);
    fail = test_stable_ddx_pdf<mpreal>(ddx_pdf_mpreal,  ctls_mpreal);
  } else if (test_name == "quantile_mpreal") {
    string outfile = out_dir + "/test_stable_quantile_mpreal.out"; 
    cout << "Writing output to " << outfile << endl;
    ofstream quantile_mpreal(outfile);
    fail = test_stable_quantile<mpreal>(quantile_mpreal, ctls_mpreal);
  } else if (test_name == "mode_mpreal") {
    string outfile = out_dir + "/test_stable_mode_mpreal.out"; 
    cout << "Writing output to " << outfile << endl;
    ofstream mode_mpreal(outfile);
    fail = test_stable_mode<mpreal>(mode_mpreal,  ctls_mpreal);
  } else if (test_name == "random_mpreal") {
    string outfile = out_dir + "/test_stable_random_mpreal.out"; 
    cout << "Writing output to " << outfile << endl;
    ofstream random_mpreal(outfile);
    fail = test_stable_random<mpreal>(random_mpreal, 10000, 1.5, .5,  ctls_mpreal);
  }
#endif //MPREAL
#ifdef CPP_BIN_FLOAT
  else if (test_name == "cdf_cpp_bin_float") {
    string outfile = out_dir + "/test_stable_cdf_cpp_bin_float.out"; 
    cout << "Writing output to " << outfile << endl;
    ofstream cdf_cpp_bin_float(outfile);
    fail = test_stable_cdf<CppBinFloat>(cdf_cpp_bin_float,  ctls_cpp_bin_float);
  } else if (test_name == "pdf_cpp_bin_float") {
    string outfile = out_dir + "/test_stable_pdf_cpp_bin_float.out"; 
    cout << "Writing output to " << outfile << endl;
    ofstream pdf_cpp_bin_float(outfile);
    fail = test_stable_pdf<CppBinFloat>(pdf_cpp_bin_float,  ctls_cpp_bin_float);
  } else if (test_name == "ddx_pdf_cpp_bin_float") {
    string outfile = out_dir + "/test_stable_ddx_pdf_cpp_bin_float.out"; 
    cout << "Writing output to " << outfile << endl;
    ofstream ddx_pdf_cpp_bin_float(outfile);
    fail = test_stable_ddx_pdf<CppBinFloat>(ddx_pdf_cpp_bin_float,  ctls_cpp_bin_float);
  } else if (test_name == "quantile_cpp_bin_float") {
    string outfile = out_dir + "/test_stable_quantile_cpp_bin_float.out"; 
    cout << "Writing output to " << outfile << endl;
    ofstream quantile_cpp_bin_float(outfile);
    fail = test_stable_quantile<CppBinFloat>(quantile_cpp_bin_float, ctls_cpp_bin_float);
  } else if (test_name == "mode_cpp_bin_float") {
    string outfile = out_dir + "/test_stable_mode_cpp_bin_float.out"; 
    cout << "Writing output to " << outfile << endl;
    ofstream mode_cpp_bin_float(outfile);
    fail = test_stable_mode<CppBinFloat>(mode_cpp_bin_float,  ctls_cpp_bin_float);
  } else if (test_name == "random_cpp_bin_float") {
    string outfile = out_dir + "/test_stable_random_cpp_bin_float.out"; 
    cout << "Writing output to " << outfile << endl;
    ofstream random_cpp_bin_float(outfile);
    fail = test_stable_random<CppBinFloat>(random_cpp_bin_float, 10000, 1.5, .5,  ctls_cpp_bin_float);
  }
#endif //CPP_BIN_FLOAT
#ifdef MPFR_FLOAT
  else if (test_name == "cdf_mpfr_float") {
    string outfile = out_dir + "/test_stable_cdf_mpfr_float.out"; 
    cout << "Writing output to " << outfile << endl;
    ofstream cdf_mpfr_float(outfile);
    fail = test_stable_cdf<MpfrFloat>(cdf_mpfr_float,  ctls_mpfr_float);
  } else if (test_name == "pdf_mpfr_float") {
    string outfile = out_dir + "/test_stable_pdf_mpfr_float.out"; 
    cout << "Writing output to " << outfile << endl;
    ofstream pdf_mpfr_float(outfile);
    fail = test_stable_pdf<MpfrFloat>(pdf_mpfr_float,  ctls_mpfr_float);
  } else if (test_name == "ddx_pdf_mpfr_float") {
    string outfile = out_dir + "/test_stable_ddx_pdf_mpfr_float.out"; 
    cout << "Writing output to " << outfile << endl;
    ofstream ddx_pdf_mpfr_float(outfile);
    fail = test_stable_ddx_pdf<MpfrFloat>(ddx_pdf_mpfr_float,  ctls_mpfr_float);
  } else if (test_name == "quantile_mpfr_float") {
    string outfile = out_dir + "/test_stable_quantile_mpfr_float.out"; 
    cout << "Writing output to " << outfile << endl;
    ofstream quantile_mpfr_float(outfile);
    fail = test_stable_quantile<MpfrFloat>(quantile_mpfr_float, ctls_mpfr_float);
  } else if (test_name == "mode_mpfr_float") {
    string outfile = out_dir + "/test_stable_mode_mpfr_float.out"; 
    cout << "Writing output to " << outfile << endl;
    ofstream mode_mpfr_float(outfile);
    fail = test_stable_mode<MpfrFloat>(mode_mpfr_float,  ctls_mpfr_float);
  } else if (test_name == "random_mpfr_float") {
    string outfile = out_dir + "/test_stable_random_mpfr_float.out"; 
    cout << "Writing output to " << outfile << endl;
    ofstream random_mpfr_float(outfile);
    fail = test_stable_random<MpfrFloat>(random_mpfr_float, 10000, 1.5, .5,  ctls_mpfr_float);
  }
#endif //MPFR_FLOAT
  else {
    cerr << "Improper test name: " << test_name << endl;
    return 1;
  }
  
  return fail;
}
