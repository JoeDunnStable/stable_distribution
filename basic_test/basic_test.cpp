//
//  basic_test.cpp
//  stable_distribution_v3
//
//  Created by Joseph Dunn on 9/25/17.
//  Copyright © 2017 Joseph Dunn. All rights reserved.
//

#include <iostream>
using std::cout;
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
#include <boost/timer/timer.hpp>
using boost::timer::auto_cpu_timer;
#include <boost/filesystem.hpp>


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
  else
    return fabs(r1-r2)/(max(fabs(r1), fabs(r2))*std::numeric_limits<myFloat>::epsilon());
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
      return log(erfc(u));
    else
      return log(erf(u));
  } else {
    if(lower_tail)
      return erfc(u);
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
  xs.push_back(pow(10.,300));
  xs.push_back(std::numeric_limits<myFloat>::infinity());
  
  out << endl;
  out << "Comparison of cdf to Levy formula for alpha = .5, beta = 1, digits10 = "
  << std::numeric_limits<myFloat>::digits10 << endl << endl;
  StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctls, verbose);
  Zolotarev<myFloat> zol(alpha, beta, &ctls.controller, verbose);
  out << setw(8) << right << "alpha"
  << setw(8) << right << "beta"
  << setw(13) << right << "x"
  << setw(8) << right << "tail"
  << setw(35) << right << "log(pLevy)"
  << setw(35) << right << "log(cdf)"
  << setw(10) << right << "epsdiff"
  << setw(15) << right << "abserr"
  << setw(4) << right << "tc"
  << setw(7) << right << "neval"
  << setw(35) << right << "log(zol_cdf)"
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
    << setw(35) << setprecision(26) << scientific << r_Levy
    << setw(35) << setprecision(26) << scientific << r
    << setw(10) << setprecision(1) << fixed << eps
    << setw(15) << setprecision(5) << scientific << std_stable_dist.abserr
    << setw(4) << std_stable_dist.termination_code
    << setw(7) << std_stable_dist.neval
    << setw(35) << setprecision(26) << scientific << r_zol
    << setw(10) << fmt_eps(eps_zol)
    << setw(2) << (zol.result_type==asymptotic ? "A" : "C")
    << (pass1 ? "" : " FAIL") << endl;
    pass = pass && pass1;
  }
  out << endl;
  return !pass;
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
  xs.push_back(pow(10.,300));
  xs.push_back(std::numeric_limits<myFloat>::infinity());
  
  int verbose = 0;
  out << endl;
  out << "Comparison of pdf to Levy formula for alpha = .5, beta = 1, digits10 = "
  << std::numeric_limits<myFloat>::digits10 << endl << endl;
  StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctls, verbose);
  Zolotarev<myFloat> zol(alpha, beta, &ctls.controller, verbose);
  int log_flag=1;
  out << setw(13) << right << "alpha"
  << setw(13) << right << "beta"
  << setw(13) << right << "x"
  << setw(35) << right << "log(dLevy)"
  << setw(35) << right << "log(pdf_myFloat)"
  << setw(10) << right << "epsdiff"
  << setw(15) << right << "abserr"
  << setw(4) << right << "TC"
  << setw(7) << right << "neval"
  << setw(35) << right << "log(zol_pdf)"
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
    << setw(35) << setprecision(26) << scientific << r_Levy
    << setw(35) << setprecision(26) << scientific << r
    << setw(10) << setprecision(1) << fixed << eps
    << setw(15) << setprecision(5) << scientific << std_stable_dist.abserr
    << setw(4) << std_stable_dist.termination_code
    << setw(7) << std_stable_dist.neval
    << setw(35) << setprecision(26) << scientific << r_zol
    << setw(10) << fmt_eps(eps_zol)
    << setw(2) << (zol.result_type==asymptotic ? "A" : "C")
    << (pass1 ? "" : " FAIL") << endl;
    pass = pass && pass1;

  }
  out << endl;
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
  xs.push_back(pow(10.,300));
  xs.push_back(std::numeric_limits<myFloat>::infinity());
  
  int verbose = 0;
  out << endl;
  out << "Comparison of ddx_pdf to Levy formula for alpha = .5, beta = 1, digits10 = "
  << std::numeric_limits<myFloat>::digits10 << endl << endl;
  StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctls, verbose);
  Zolotarev<myFloat> zol(alpha, beta, &ctls.controller, verbose);
  out << setw(13) << right << "alpha"
  << setw(13) << right << "beta"
  << setw(13) << right << "x"
  << setw(35) << right << "Levy_ddx_pdf"
  << setw(35) << right << "ddx_pdf_myFloat"
  << setw(10) << right << "epsdiff"
  << setw(15) << right << "abserr"
  << setw(4) << right << "termination_code"
  << setw(7) << right << "neval"
  << setw(35) << right << "zol_ddx_pdf)"
  << setw(10) << right << "epsdiff"
  << setw(2) << right << "T"
  << endl << endl;
  bool pass = true;
  for (auto x : xs) {
    myFloat r_Levy = stable_levy_ddx_pdf<myFloat>(x);
    myFloat r = std_stable_dist.ddx_pdf(x, S1);
    myFloat eps = epsdiff(r, r_Levy);
    myFloat r_zol = zol.ddx_pdf(x,S1);
    myFloat eps_zol = epsdiff(r_zol, r_Levy);
    bool pass1 = !boost::math::isnan(eps) && boost::math::isfinite(eps)
    && eps < 500;

    out << setw(13) << setprecision(5) << fixed << myFloat(alpha)
    << setw(13) << setprecision(5) << fixed << myFloat(beta)
    << setw(13) << setprecision(3) << scientific << myFloat(x)
    << setw(35) << setprecision(26) << scientific << r_Levy
    << setw(35) << setprecision(26) << scientific << r
    << setw(10) << setprecision(1) << fixed << eps
    << setw(15) << setprecision(5) << scientific << std_stable_dist.abserr
    << setw(4) << std_stable_dist.termination_code
    << setw(7) << std_stable_dist.neval
    << setw(35) << setprecision(26) << scientific << r_zol
    << setw(10) << fmt_eps(eps_zol)
    << setw(2) << (zol.result_type==asymptotic ? "A" : "C")
    << (pass1 ? "" : " FAIL") << endl;
    pass = pass && pass1;
  }
  out << endl;
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
    ps.push_back(pow<myFloat>(10,i));
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
    ps.push_back(pow<myFloat>(10,i));
    lower_tails.push_back(0);
  }
  ps.push_back(0);
  lower_tails.push_back(0);
  
  out << endl;
  out << "Comparison of quantile to Levy formula for alpha = .5, beta = 1, digits10 = "
  << std::numeric_limits<myFloat>::digits10 << endl << endl;
  StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctls, verbose);
  out << setw(8) << right << "alpha"
  << setw(8) << right << "beta"
  << setw(13) << right << "p"
  << setw(8) << right << "tail"
  << setw(35) << right << "Levy_quantile"
  << setw(35) << right << "quantile"
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
    << setw(35) << setprecision(26) << scientific << r_Levy
    << setw(35) << setprecision(26) << scientific << r
    << setw(10) << setprecision(1) << fixed << eps
    << setw(7) << std_stable_dist.num_iter
    << setw(7) << std_stable_dist.neval
    << (pass1 ? "" : " FAIL") << endl;
    pass = pass && pass1;
  }
  out << endl;
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
    if (tgamma(1+1/alpha) > std::numeric_limits<myFloat>::max()) {
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
    for (myFloat beta : betas)
      out << setw(30) << setprecision(6) << scientific << pdf_at_mode.at(j++);
    out << endl;
  }
  out << endl;
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
  out << "Call to random with n = " << n << ", alpha = " << alpha
  << ", beta = " << beta << ", digits10 = " << std::numeric_limits<myFloat>::digits10
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
  out << endl;
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
  mpreal eps_mpreal = std::numeric_limits<mpreal>::epsilon();
  double epsabs_double = 0;
  double epsrel_double = 64 * std::numeric_limits<double>::epsilon();
  int limit = 1000;
  int verbose_integration = 0;
  IntegrationController<double> ctl_double(noext, k_big,
                                           epsabs_double, epsrel_double,
                                           limit, verbose_integration);
  if (verbose_integration)
    cout << "ctl_double:" << endl << ctl_double << endl;
  
  mpreal epsabs_mpreal = 0;
  mpreal epsrel_mpreal = 64 * std::numeric_limits<mpreal>::epsilon();
  IntegrationController<mpreal> ctl_mpreal(noext, k_big,
                                           epsabs_mpreal, epsrel_mpreal,
                                           limit, verbose_integration);
  if (verbose_integration)
    cout << "ctl_mpreal:" << endl << ctl_mpreal << endl;

  bool fail = false;
  if (test_name == "cdf_double") {
    cout << "Writing output to ../output/test_stable_cdf_double.out" << endl;
    ofstream cdf_double("../output/test_stable_cdf_double.out");
    fail = test_stable_cdf<double>(cdf_double, Controllers<double>(ctl_double, ctl_double));
  } else if (test_name == "pdf_double") {
    cout << "Writing output to ../output/test_stable_pdf_double.out" << endl;
    ofstream pdf_double("../output/test_stable_pdf_double.out");
    fail = test_stable_pdf<double>(pdf_double, Controllers<double>(ctl_double, ctl_double));
  } else if (test_name == "ddx_pdf_double") {
    cout << "Writing output to ../output/test_stable_ddx_pdf_double.out" << endl;
    ofstream ddx_pdf_double("../output/test_stable_ddx_pdf_double.out");
    fail = test_stable_ddx_pdf<double>(ddx_pdf_double, Controllers<double>(ctl_double, ctl_double));
  } else if (test_name == "quantile_double") {
    cout << "Writing output to ../output/test_stable_quantile_double.out" << endl;
    ofstream quantile_double("../output/test_stable_quantile_double.out");
    fail = test_stable_quantile<double>(quantile_double, Controllers<double>(ctl_double, ctl_double));
  } else if (test_name == "mode_double") {
    cout << "Writing output to ../output/test_stable_mode_double.out" << endl;
    ofstream mode_double("../output/test_stable_mode_double.out");
    fail = test_stable_mode<double>(mode_double, Controllers<double>(ctl_double, ctl_double));
  } else if (test_name == "random_double") {
    cout << "Writing output to ../output/test_stable_random_double.out" << endl;
    ofstream random_double("../output/test_stable_random_double.out");
    fail = test_stable_random<double>(random_double, 10000, .5, .5, Controllers<double>(ctl_double, ctl_double));
    fail = test_stable_random<double>(random_double, 10000, 1.-1./128., .5, Controllers<double>(ctl_double, ctl_double)) || fail;
    fail = test_stable_random<double>(random_double, 10000, 1, .5, Controllers<double>(ctl_double, ctl_double)) || fail;
    fail = test_stable_random<double>(random_double, 10000, 1.+1./128., .5, Controllers<double>(ctl_double, ctl_double)) || fail;
    fail = test_stable_random<double>(random_double, 10000, 1.5, .5, Controllers<double>(ctl_double, ctl_double)) || fail;
  } else if (test_name == "cdf_mpreal") {
    cout << "Writing output to ../output/test_stable_cdf_mpreal.out" << endl;
    ofstream cdf_mpreal("../output/test_stable_cdf_mpreal.out");
    fail = test_stable_cdf<mpreal>(cdf_mpreal,  Controllers<mpreal>(ctl_mpreal, ctl_double));
  } else if (test_name == "pdf_mpreal") {
    cout << "Writing output to ../output/test_stable_pdf_mpreal.out" << endl;
    ofstream pdf_mpreal("../output/test_stable_pdf_mpreal.out");
    fail = test_stable_pdf<mpreal>(pdf_mpreal,  Controllers<mpreal>(ctl_mpreal, ctl_double));
  } else if (test_name == "ddx_pdf_mpreal") {
    cout << "Writing output to ../output/test_stable_ddx_pdf_mpreal.out" << endl;
    ofstream ddx_pdf_mpreal("../output/test_stable_ddx_pdf_mpreal.out");
    fail = test_stable_ddx_pdf<mpreal>(ddx_pdf_mpreal,  Controllers<mpreal>(ctl_mpreal, ctl_double));
  } else if (test_name == "quantile_mpreal") {
    cout << "Writing output to ../output/test_stable_quantile_mpreal.out" << endl;
    ofstream quantile_mpreal("../output/test_stable_quantile_mpreal.out");
    fail = test_stable_quantile<mpreal>(quantile_mpreal, Controllers<mpreal>(ctl_mpreal, ctl_double));
  } else if (test_name == "mode_mpreal") {
    cout << "Writing output to ../output/test_stable_mode_mpreal.out" << endl;
    ofstream mode_mpreal("../output/test_stable_mode_mpreal.out");
    fail = test_stable_mode<mpreal>(mode_mpreal,  Controllers<mpreal>(ctl_mpreal, ctl_double));
  } else if (test_name == "random_mpreal") {
    cout << "Writing output to ../output/test_stable_random_mpreal.out" << endl;
    ofstream random_mpreal("../output/test_stable_random_mpreal.out");
    fail = test_stable_random<mpreal>(random_mpreal, 10000, 1.5, .5,  Controllers<mpreal>(ctl_mpreal, ctl_double));
  } else {
    cerr << "Improper test name: " << test_name << endl;
    return 1;
  }
  
  return fail;
}
