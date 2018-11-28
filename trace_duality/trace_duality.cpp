/// \file trace_duality.cpp
/// Checks that Zolotarev duality theorem is satisfied
/// \author Joseph Dunn
/// \copyright 2017, 2018 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include "stable_config.h"

#define MPREAL
#define MPFR_FLOAT
#define CPP_BIN_FLOAT
#include "stable_distribution.h"
using namespace stable_distribution;

#include <iomanip>
using std::setw;
using std::setprecision;
using std::right;
using std::scientific;
using std::fixed;

#include <string>
using std::string;
using std::getline;

#include <sstream>
using std::istringstream;
using std::ostringstream;

#include <boost/timer/timer.hpp>
using boost::timer::auto_cpu_timer;

#include <boost/filesystem.hpp>
using boost::filesystem::path;

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
  ostringstream ss;
  if (eps < 1e4)
    ss << setw(10) << setprecision(1) << fixed << eps;
  else
    ss << setw(10) << right << "*******" ;
  return ss.str();
}

template<typename myFloat>
void reset_prec(int digits) {}

template<>
void reset_prec<mpreal>(int digits) {
  mpreal::set_default_prec(digits);
}

template<typename myFloat>
void duality_check(double alpha_m_1, double beta, double x_prime,
                   Controllers<myFloat>& ctls, int verbose) {
  Fmt<myFloat> fmt;
  myFloat pi2 = const_pi<myFloat>()/2;
  myFloat inf = std::numeric_limits<myFloat>::infinity();
  int log_flag = 1;
  // Zolotarev Theorem 2.3.2
  if (!(alpha_m_1 > 0 && -1 <= beta && beta <=1)) {
    throw std::range_error("duality check parameter error");
  }
  myFloat alpha = myFloat(alpha_m_1)+1;
  /*
  bool near = abs(alpha_m_1) < StandardStableDistribution<myFloat>::threshhold_1*abs(beta);
   */
  myFloat alpha_prime, alpha_prime_m_1, zeta, Q, beta_star, D;
  if (true /* near */) {
    alpha_prime_m_1 = expm1(-log1p(myFloat(alpha_m_1)));
    alpha_prime = 1 + alpha_prime_m_1;
    zeta = beta/tan(pi2*alpha_m_1);
    Q = (zeta<0) ? atan(1/abs(zeta))/pi2
    : 2-atan(1/abs(zeta))/pi2;
    beta_star = (alpha==2 || beta==-1)
    ? 1
    :-tan(pi2*alpha_prime_m_1)/tan(pi2*Q*alpha_prime);
    D = pow(1+zeta*zeta,1/(2*myFloat(alpha)))
        /sin(pi2*Q*alpha_prime);
  } else {
    alpha_prime = myFloat(1)/alpha;
    alpha_prime_m_1 = alpha_prime - 1;
    zeta = -beta*tan(pi2*alpha);
    Q = 1-atan(-zeta)/pi2;
    beta_star = (alpha==2 || beta==-1)
    ? 1
    :1/(tan(pi2*alpha_prime)*tan(pi2*Q*alpha_prime));
    D = pow(1+pow(beta*tan(pi2*alpha),2),1/(2*myFloat(alpha)))
    /sin(pi2*Q*alpha_prime);
  }
  myFloat min_ln = log(std::numeric_limits<myFloat>::min());
  StandardStableDistribution<myFloat> dist1(AlphaMinusOne<myFloat>(alpha_prime_m_1), beta_star, ctls, verbose);
  StandardStableDistribution<myFloat> dist2(AlphaMinusOne<myFloat>(alpha_m_1), beta, ctls, verbose);
  myFloat x = D*pow(myFloat(x_prime),-alpha_prime);
  myFloat ln_pdf1 = dist1.pdf(x_prime, log_flag, S1);
  myFloat ln_tmp = dist2.pdf(x, log_flag, S1);
  myFloat ln_pdf2 = log(D) + (-1-alpha_prime)*log(myFloat(x_prime)) + ln_tmp;
  ;
  if (ln_pdf1 < min_ln) ln_pdf1 = -inf;
  if (ln_pdf2 < min_ln) ln_pdf2 = -inf;
  myFloat eps_diff = epsdiff(ln_pdf1, ln_pdf2);
  cout << "ln_pdf1  =         " << fmt << ln_pdf1 << endl
       << "D        =         " << fmt << D << endl
       << "log(D)   =         " << fmt << log(D) << endl
       << "log(x'^-1-alpha' = " << fmt << (-1-alpha_prime)*log(myFloat(x_prime)) << endl
       << "log(dist2.pdf) =   " << fmt << ln_tmp << endl
       << "ln_pdf2  =         " << fmt << ln_pdf2 << endl
       << "eps_diff = " << fixed << setprecision(1) << eps_diff << endl;
   return;
}

static void show_usage (path p){
  cerr << "Usage: " << p.filename().string() << " float_type alpha_m_1 beta x_prime"<< endl << endl
  << "         where float_type is double"
#ifdef MPREAL
  << ", mpreal"
#endif
#ifdef CPP_BIN_FLOAT
  << ", cpp_bin_float"
#endif
#ifdef MPFR_FLOAT
  << ", mpfr_float"
#endif
  << endl;
}

int main(int argc, char *argv[]) {
  
  // Check the number of parameters
  if (argc != 5) {
    // Tell the user how to run the program
    show_usage(string(argv[0]));
    return 1;
  }
  string name(argv[0]);
  path p((string(argv[0])));
  string float_type(argv[1]);
  istringstream ss2((argv[2])), ss3((argv[3])), ss4((argv[4]));
  double alpha_m_1;
  if (!(ss2>>alpha_m_1)) {
    cerr << "Unreadable alpha_m_1" << endl;
    show_usage(p);
    return 1;
  }
  double beta;
  if (!(ss3>>beta)) {
    cerr << "unreadable beta" << endl;
    show_usage(p);
    return 1;
  }
  double x_prime;
  if (!(ss4 >> x_prime)) {
    cerr << "Ureadable x_prime" << endl;
    show_usage(p);
    return 1;
  }
  int verbose = 4;
  
  string title ;
  Fmt<double> fmt_dbl;
  ostringstream oss(title);
  oss << p.filename().string() << ": alpha_m_1 = " << fmt_dbl << alpha_m_1
      << ", beta = " << fmt_dbl << beta
      << ", x_prime = " << fmt_dbl << x_prime << endl;
  

  
#ifdef MPREAL
  // First create some high precision coeficients for
  mpreal::set_default_prec(128);
  Kronrod<mpreal> k_big(10);
  mpreal::set_default_prec(96);
#define BIGFLOAT mpreal
#elif MPFR_FLOAT
  Kronrod<BigMpfrFloat> k_big(10);
#define BIGFLOAT BigMpfrFloat
#elif CPP_BIN_FLOAT
  Kronrod<BigCppBinFloat> k_big(10);
#define BIGFLOAT BigCppBinFloat
#else
  Kronrod<double> k_big(10);
#define BIGFLOAT double
#endif
  bool noext = true;
  double epsabs_dbl = 0;
  double epsrel_dbl = 64 * std::numeric_limits<double>::epsilon();
  int subdivisions = 1000;
  int verbose_integration = 0;
  IntegrationController<double> cntl_double(noext, k_big, epsabs_dbl, epsrel_dbl,
                                            subdivisions, verbose_integration);

  if (float_type == "double") {
    Controllers<double> ctls(cntl_double, cntl_double);
    cout << title << endl
       << "Machine epsilon = " << std::numeric_limits<double>::epsilon() << endl << endl;
    duality_check<double>(alpha_m_1, beta, x_prime, ctls, verbose);
  }
  
#ifdef MPREAL
  else if(float_type == "mpreal") {
    mpreal epsabs=0;
    mpreal epsrel=64*std::numeric_limits<mpreal>::epsilon();
    mpreal::set_default_prec(128);
    Kronrod<mpreal> g_k_big(10);
    mpreal::set_default_prec(96);
    StandardStableDistribution<mpreal>::initialize();

    IntegrationController<mpreal> cntl_mpreal(noext, g_k_big, epsabs,
                                              epsrel, subdivisions, verbose);
    Controllers<mpreal> ctls(cntl_mpreal, cntl_double);
    cout << title << endl
    << "Machine epsilon = " << std::numeric_limits<mpreal>::epsilon() << endl << endl;
    duality_check<mpreal>(alpha_m_1, beta, x_prime, ctls, verbose);
  }
#endif
  
#ifdef CPP_BIN_FLOAT
  else if (float_type == "cpp_bin_float") {
    CppBinFloat epsabs=0;
    CppBinFloat epsrel=64*std::numeric_limits<CppBinFloat>::epsilon();
    Kronrod<BigCppBinFloat> g_k_big(10);
    IntegrationController<CppBinFloat> cntl_cpp_bin_float(noext, g_k_big, epsabs,
                                              epsrel, subdivisions, verbose);
    Controllers<CppBinFloat> ctls(cntl_cpp_bin_float, cntl_double);
    cout << title << endl
    << "Machine epsilon = " << std::numeric_limits<CppBinFloat>::epsilon() << endl << endl;
    duality_check<CppBinFloat>(alpha_m_1, beta, x_prime, ctls, verbose);
  }
#endif
  
#ifdef MPFR_FLOAT
  else if (float_type == "mpfr_float") {
    MpfrFloat epsabs=0;
    MpfrFloat epsrel=64*std::numeric_limits<MpfrFloat>::epsilon();
    Kronrod<BigMpfrFloat> g_k_big(10);
    IntegrationController<MpfrFloat> cntl_mpfr_float(noext, g_k_big, epsabs,
                                                          epsrel, subdivisions, verbose);
    Controllers<MpfrFloat> ctls(cntl_mpfr_float, cntl_double);
    cout << title << endl
    << "Machine epsilon = " << std::numeric_limits<MpfrFloat>::epsilon() << endl << endl;
    duality_check<MpfrFloat>(alpha_m_1, beta, x_prime, ctls, verbose);
  }
#endif
  else {
    show_usage(p);
    return 1;
    
  } //
  
  return 0;
}




