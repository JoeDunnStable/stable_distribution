//
//  main.cpp
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
#include <random>
using std::mt19937;
using std::uniform_real_distribution;

#include "stable_distribution.h"
#include "kolmogorov.h"

using namespace stable_distribution;

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
void test_stable_cdf(ostream& out, IntegrationController<myFloat>& ctl) {
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
  
  out << "Comparison of cdf to Levy formula for alpha = .5, beta = 1, digits10 = "
  << std::numeric_limits<myFloat>::digits10 << endl << endl;
  StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctl, verbose);
  out << setw(8) << right << "alpha"
  << setw(8) << right << "beta"
  << setw(13) << right << "x"
  << setw(8) << right << "tail"
  << setw(35) << right << "log(pLevy)"
  << setw(35) << right << "log(cdf)"
  << setw(15) << right << "absdiff"
  << setw(15) << right << "abserr"
  << setw(5) << right << "code"
  << setw(7) << right << "neval"
  << endl << endl;
  for (auto x : xs) {
    lower_tail = true;
    myFloat r_Levy = stable_levy_cdf<myFloat>(x, lower_tail, log_p);
    if (r_Levy > (log_p ? log(.5) : .5)) {
      lower_tail = false;
      r_Levy = stable_levy_cdf<myFloat>(x, lower_tail, log_p);
    }
    myFloat r = std_stable_dist.cdf(x, lower_tail, log_p, StandardStableDistribution<myFloat>::Parameterization::S1);
    
    out << setw(8) << setprecision(3) << fixed << myFloat(alpha)
    << setw(8) << setprecision(3) << fixed << myFloat(beta)
    << setw(13) << setprecision(3) << scientific << myFloat(x)
    << setw(8) << right << ((lower_tail) ? "lower" : "upper")
    << setw(35) << setprecision(26) << scientific << r_Levy
    << setw(35) << setprecision(26) << scientific << r
    << setw(15) << setprecision(5) << scientific << fabs((r-r_Levy))
    << setw(15) << setprecision(5) << scientific << std_stable_dist.abserr
    << setw(4) << std_stable_dist.termination_code
    << setw(7) << std_stable_dist.neval << endl;
  }
  out << endl << endl;
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
void test_stable_pdf(ostream& out, IntegrationController<myFloat>& ctl) {
  
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
  out << "Comparison of pdf to Levy formula for alpha = .5, beta = 1, digits10 = "
  << std::numeric_limits<myFloat>::digits10 << endl << endl;
  StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctl, verbose);
  int log_flag=1;
  out << setw(13) << right << "alpha"
  << setw(13) << right << "beta"
  << setw(13) << right << "x"
  << setw(35) << right << "log(dLevy)"
  << setw(35) << right << "log(pdf_myFloat)"
  << setw(15) << right << "absdiff"
  << setw(15) << right << "abserr"
  << setw(4) << right << "termination_code"
  << setw(7) << right << "neval"
  << endl << endl;
  for (auto x : xs) {
    myFloat r = std_stable_dist.pdf(x, log_flag, StandardStableDistribution<myFloat>::Parameterization::S1);
    myFloat r_Levy = stable_levy_pdf<myFloat>(x, log_flag);
    
    out << setw(13) << setprecision(5) << fixed << myFloat(alpha)
    << setw(13) << setprecision(5) << fixed << myFloat(beta)
    << setw(13) << setprecision(3) << scientific << myFloat(x)
    << setw(35) << setprecision(26) << scientific << r_Levy
    << setw(35) << setprecision(26) << scientific << r
    << setw(15) << setprecision(5) << scientific << fabs((r-r_Levy))
    << setw(15) << setprecision(5) << scientific << std_stable_dist.abserr
    << setw(4) << std_stable_dist.termination_code
    << setw(7) << std_stable_dist.neval << endl;
  }
  out << endl << endl;
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
void test_stable_ddx_pdf(ostream& out, IntegrationController<myFloat>& ctl) {
  
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
  out << "Comparison of ddx_pdf to Levy formula for alpha = .5, beta = 1, digits10 = "
  << std::numeric_limits<myFloat>::digits10 << endl << endl;
  StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctl, verbose);
  out << setw(13) << right << "alpha"
  << setw(13) << right << "beta"
  << setw(13) << right << "x"
  << setw(35) << right << "Levy_ddx_pdf"
  << setw(35) << right << "ddx_pdf_myFloat"
  << setw(15) << right << "absdiff"
  << setw(15) << right << "abserr"
  << setw(4) << right << "termination_code"
  << setw(7) << right << "neval"
  << endl << endl;
  for (auto x : xs) {
    myFloat r = std_stable_dist.ddx_pdf(x, StandardStableDistribution<myFloat>::Parameterization::S1);
    myFloat r_Levy = stable_levy_ddx_pdf<myFloat>(x);
    
    out << setw(13) << setprecision(5) << fixed << myFloat(alpha)
    << setw(13) << setprecision(5) << fixed << myFloat(beta)
    << setw(13) << setprecision(3) << scientific << myFloat(x)
    << setw(35) << setprecision(26) << scientific << r_Levy
    << setw(35) << setprecision(26) << scientific << r
    << setw(15) << setprecision(5) << scientific << fabs((r-r_Levy))
    << setw(15) << setprecision(5) << scientific << std_stable_dist.abserr
    << setw(4) << std_stable_dist.termination_code
    << setw(7) << std_stable_dist.neval << endl;
  }
  out << endl << endl;
} //test_stable_ddx_pdf

template<typename myFloat>
class erf_solve {
  myFloat target;
  bool lower_tail;
public:
  erf_solve(myFloat target, bool lower_tail) : target(target), lower_tail(lower_tail){}
  myFloat operator() (myFloat x) {
    return (lower_tail ? erfc(x) : erf(x)) - target;
  }
};

template<typename myFloat>
myFloat stable_levy_quantile(myFloat pp, bool lower_tail=true, bool log_p=false) {
  // pm = S1 representation
  myFloat c= 1;
  myFloat p = log_p ? exp(pp) : pp;
  if (p == 0)
    return  lower_tail ? StandardStableDistribution<myFloat>::NegInf
                       : StandardStableDistribution<myFloat>::PosInf;
  if (p == 1)
    return lower_tail ? StandardStableDistribution<myFloat>::PosInf
                      : StandardStableDistribution<myFloat>::NegInf;
  
  double u = lower_tail ? erfc_inv(static_cast<double>(p))
                        : erf_inv(static_cast<double>(p));
  myFloat guess = c/(2*u*u);
  erf_solve<myFloat> erf_s(p, lower_tail);
  pair<myFloat,myFloat> r;
  RelativeComparisonTolerance<myFloat> tol(16*std::numeric_limits<myFloat>::epsilon());
  myFloat factor = max<myFloat>(static_cast<myFloat>(1),static_cast<myFloat>(.1)*fabs(guess));
  bool rising = !lower_tail;
  boost::uintmax_t maxiter = 1000;
  r=boost::math::tools::bracket_and_solve_root2(erf_s,guess,factor,rising,tol,maxiter);
  myFloat u2 =(r.first + r.second)/2;
  return c/(2*u2*u2);
}

template<>
double stable_levy_quantile<double>(double pp, bool lower_tail, bool log_p) {
  // pm = S1 representation
  double c= 1;
  double p = log_p ? exp(pp) : pp;
  if (p == 0)
    return  lower_tail ? StandardStableDistribution<double>::NegInf
    : StandardStableDistribution<double>::PosInf;
  if (p == 1)
    return lower_tail ? StandardStableDistribution<double>::PosInf
    : StandardStableDistribution<double>::NegInf;
  
  double u = lower_tail ? erfc_inv(p)
                        : erf_inv(p);
  return c/(2*u*u);
}

template<typename myFloat>
void test_stable_quantile(ostream& out, IntegrationController<myFloat>& ctl) {
  int verbose = 0;
  
  bool log_p=false;
  myFloat q_tol = 64 * std::numeric_limits<myFloat>::epsilon();
  myFloat alpha = .5;
  myFloat beta = 1;
  
  vector<myFloat> ps;
  vector<int> lower_tails;
  ps.push_back(0);
  lower_tails.push_back(1);
  for (int i = -4; i<0; ++i) {
    ps.push_back(pow<myFloat>(10,i));
    lower_tails.push_back(1);
  }
  ps.push_back(.5);
  lower_tails.push_back(1);
  for (int i = -1; i>-5; --i) {
    ps.push_back(pow<myFloat>(10,i));
    lower_tails.push_back(0);
  }
  ps.push_back(0);
  lower_tails.push_back(0);
  
  out << "Comparison of quantile to Levy formula for alpha = .5, beta = 1, digits10 = "
  << std::numeric_limits<myFloat>::digits10 << endl << endl;
  StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctl, verbose);
  out << setw(8) << right << "alpha"
  << setw(8) << right << "beta"
  << setw(13) << right << "p"
  << setw(8) << right << "tail"
  << setw(35) << right << "Levy_quantile"
  << setw(35) << right << "quantile"
  << setw(15) << right << "reldiff"
  << setw(7) << right << "iters"
  << setw(7) << right << "neval"
  << endl << endl;
  for (int i=0; i<ps.size(); ++i) {
    myFloat p = ps.at(i);
    int lower_tail = lower_tails.at(i);
    myFloat r_Levy = stable_levy_quantile<myFloat>(p, lower_tail, log_p);
    myFloat r = std_stable_dist.quantile(p, lower_tail, log_p,
                                         q_tol, StandardStableDistribution<myFloat>::Parameterization::S1);
    
    out << setw(8) << setprecision(3) << fixed << myFloat(alpha)
    << setw(8) << setprecision(3) << fixed << myFloat(beta)
    << setw(13) << setprecision(3) << scientific << myFloat(p)
    << setw(8) << right << ((lower_tail) ? "lower" : "upper")
    << setw(35) << setprecision(26) << scientific << r_Levy
    << setw(35) << setprecision(26) << scientific << r
    << setw(15) << setprecision(5) << scientific << fabs((r-r_Levy)/r_Levy)
    << setw(7) << std_stable_dist.num_iter
    << setw(7) << std_stable_dist.neval << endl;
  }
  out << endl << endl;
}

template<typename myFloat>
void print_stable_mode_heading(ostream& os, const vector<myFloat>& betas) {
  os << setw(65) << right << "beta" << endl
  << setw(20) << right << "alpha";
  for (auto beta : betas)
    os << setw(15) << setprecision(11) << fixed <<  beta;
  os << endl;
  
}

template<typename myFloat>
void test_stable_mode(ostream& out, IntegrationController<myFloat>& ctl) {
  myFloat pi2 = const_pi<myFloat>()/2;
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
  for (int i=-2; i<=2; i++) betas.push_back(i/static_cast<myFloat>(2));
  myFloat ddx_tol = 64*std::numeric_limits<myFloat>::epsilon();
  //vector<myFloat> out(alphas.size()*betas.size());
  
  out << setw(75) << right
  << "Position of Mode in the S1 Parameterization" << endl << endl;
  print_stable_mode_heading(out, betas);
  
  int verbose_mode = 0;
  int verbose = 0;
  auto pm = StandardStableDistribution<myFloat>::Parameterization::S1;
  vector<myFloat> mode_S1;
  vector<myFloat> pdf_at_mode;
  myFloat first_good_alpha = 2;
  for (auto alpha : alphas) {
    out << setw(20) << setprecision(14) << fixed << alpha;
    if (tgamma(1+1/alpha) > std::numeric_limits<myFloat>::max()) {
      out << setw(30) << "alpha is too small." << endl;
      continue;
    } else {
      first_good_alpha = min(first_good_alpha, alpha);
    }
    for (auto beta : betas) {
      StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctl, verbose);
      std::pair<myFloat, myFloat> mode = std_stable_dist.mode(ddx_tol, verbose_mode, pm);
      out << setw(15) << setprecision(6) << scientific << mode.first;
      mode_S1.push_back(mode.first);
      pdf_at_mode.push_back(mode.second);
    }
    out << endl;
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
      out << setw(15) << setprecision(6) << scientific << mode_S1.at(j++)+zeta;
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
      out << setw(15) << setprecision(6) << scientific << pdf_at_mode.at(j++);
    out << endl;
  }
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
void test_stable_random (ostream& out, int n, myFloat alpha, myFloat beta,
                    IntegrationController<myFloat>& ctl) {
  out << "n = " << n
      << ", alpha = " << alpha
      << ", beta = " << beta << endl << endl;
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
  
  StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctl, verbose);
  
  int log_p=0;
  
  myFloat sum_neginf = 0, sum_posinf = 0, sum_zero = 0, sum_nan = 0;
  for (typename vector<myFloat>::iterator pr=r.begin(); pr < r.end(); pr++) {
    sum_neginf+=(*pr < -1e300);
    sum_posinf+=(*pr > 1e300);
    sum_zero+=(fabs(*pr)<.01);
    sum_nan+=boost::math::isnan(*pr);
  }
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
  
  out << "D = " << d << ", KS Probability = " << kolmogorov_asymptotic_cdf(d * sqrt(n))
      << ", 1 - K(n,d) = " << 1-kolmogorov_cdf(n, d) << endl;
}

int main(int argc, const char * argv[]) {
  // First create some high precision coeficients for
  mpreal::set_default_prec(mpfr::digits2bits(60));
  Kronrod<mpreal> k_big(10);
  mpreal::set_default_prec(mpfr::digits2bits(50));
  int noext = 1;
  mpreal eps_mpreal = std::numeric_limits<mpreal>::epsilon();
  double epsabs_double = 0;
  double epsrel_double = 64 * std::numeric_limits<double>::epsilon();
  int limit = 3000;
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
/*
  test_stable_cdf<double>(cout, ctl_double);
  test_stable_pdf<double>(cout, ctl_double);
  test_stable_ddx_pdf<double>(cout, ctl_double);
  test_stable_quantile<double>(cout, ctl_double);
  test_stable_mode<double>(cout, ctl_double);
 */
  test_stable_random<double>(cout, 10000, 1.5, .5, ctl_double);
/*
  test_stable_cdf<mpreal>(cout, ctl_mpreal);
  test_stable_pdf<mpreal>(cout, ctl_mpreal);
  test_stable_ddx_pdf<mpreal>(cout, ctl_mpreal);
  test_stable_quantile<mpreal>(cout,ctl_mpreal);
  test_stable_mode<mpreal>(cout, ctl_mpreal);
*/
  test_stable_random<mpreal>(cout, 10000, 1.5, .5, ctl_mpreal);

  return 0;
}
