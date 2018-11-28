/// \file trace_zolotarev.cpp
/// Traces for the routines in zolotarev.h
/// \author Joseph Dunn
/// \copyright 2017, 2018 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <iomanip>
#define MPREAL
#include "stable_distribution.h"
#define LIBRARY
#include "zolotarev.h"
#undef LIBRARY
#include <boost/math/constants/info.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/trim.hpp>
using boost::algorithm::trim;
#include <vector>

using std::setw;
using std::fixed;
using std::scientific;
using std::right;
using std::istringstream;
using std::setprecision;
using std::string;
using std::vector;

using namespace stable_distribution;

template<typename myFloat>
void calc_result(ostream& os, string type, double alpha, double beta, double x,
                 Parameterization pm,
                 Controllers<myFloat> ctls, int verbose, int verbose_integration) {
  os << "Calculating " << type << endl << endl;
  int log_flag=0, lower_tail=x<0;
  myFloat error_std, error_zol, error_series_small_x, error_series_large_x;
  myFloat res_zol, res_std, res_series_small_x, res_series_large_x;
  int n_series_large_x, n_series_small_x;
  Zolotarev<myFloat> zol(myFloat(alpha), myFloat(beta), &ctls.controller,
                        verbose, verbose_integration) ;
  StandardStableDistribution<myFloat> std_dist(myFloat(alpha), myFloat(beta), ctls, 2);
  if (type == "cdf_double" || type == "cdf_mpreal") {
    res_zol = zol.cdf(myFloat(x), lower_tail, pm) ;
    error_zol = zol.error;
    res_std = std_dist.cdf(myFloat(x), lower_tail, log_flag, pm);
    error_std= std_dist.abserr;
    res_series_small_x = std_dist.series_small_x_cdf(myFloat(x),lower_tail, pm);
    error_series_small_x = std_dist.error_series;
    n_series_small_x = std_dist.n_series;
    res_series_large_x = std_dist.series_large_x_cdf(myFloat(x), lower_tail, pm);
    error_series_large_x = std_dist.error_series;
    n_series_large_x = std_dist.n_series;
  } else if (type == "pdf_double" || type == "pdf_mpreal"){
    res_zol = zol.pdf(myFloat(x), pm) ;
    error_zol = zol.error;
    res_std = std_dist.pdf(myFloat(x), log_flag, pm);
    error_std= std_dist.abserr;
    res_series_small_x = std_dist.series_small_x_pdf(myFloat(x), pm);
    error_series_small_x = std_dist.error_series;
    n_series_small_x = std_dist.n_series;
    res_series_large_x = std_dist.series_large_x_pdf(myFloat(x), pm);
    error_series_large_x = std_dist.error_series;
    n_series_large_x = std_dist.n_series;
  } else {
    res_zol = zol.ddx_pdf(myFloat(x), pm) ;
    error_zol = zol.error;
    res_std = std_dist.ddx_pdf(myFloat(x), pm);
    error_std= std_dist.abserr;
    res_series_small_x = std_dist.series_small_x_ddx_pdf(myFloat(x), pm);
    error_series_small_x = std_dist.error_series;
    n_series_small_x = std_dist.n_series;
    res_series_large_x = std_dist.series_large_x_ddx_pdf(myFloat(x), pm);
    error_series_large_x = std_dist.error_series;
    n_series_large_x = std_dist.n_series;

  }
  Fmt<myFloat> fmt;
  os.unsetf(std::ios_base::floatfield);
  os << setw(25) << "alpha = " << setw(6)<< setprecision(3) << alpha << endl
     << setw(25) << "beta = " << setw(6)<<setprecision(3) << beta << endl
     << setw(25) << "x = " << setw(6) << setprecision(2) << x << endl
     << setw(25) << "Parameterization = " << ((pm == S0)?"S0":"S1") << endl
     << setw(25) << "good theta2 = "<< setw(2) << std_dist.good_theta2 << endl
     << setw(25) << "result stable = " << fmt  << res_std << endl
     << setw(25) << "error stable = " << fmt << error_std << endl
     << setw(25) << "result zol = " << fmt << res_zol << endl
     << setw(25) << "error zol = " << fmt << error_zol << endl
     << setw(25) << "n zol = " << setw(4) << zol.n << endl
     << setw(25) << "rel error zol= " << fmt << res_zol/res_std-1 << endl
     << setw(25) << "result series_small = " << fmt << res_series_small_x << endl
     << setw(25) << "error series small = " << fmt << error_series_small_x << endl
     << setw(25) << "n series_small = " << setw(4) << n_series_small_x << endl
     << setw(25) << "rel error series small x = " << fmt << res_series_small_x/res_std-1 << endl
     << setw(25) << "result series_large = " << fmt << res_series_large_x << endl
     << setw(25) << "error series large = " << fmt << error_series_large_x << endl
     << setw(25) << "n series_large = " << setw(4) << n_series_large_x << endl
  << setw(25) << "rel error series large = " << fmt << res_series_large_x/res_std-1 << endl;
}


static void show_usage (string name){
  boost::filesystem::path p(name);
  cerr << endl;
  cerr << "Usage: " << p.filename().string() << " test_name ..." << endl;
  cerr << "cdf_double alpha beta x [pm=0] [verbose=4] {verbose_integration=0]" << endl;
  cerr << "pdf_double alpha beta x [pm=0] [verbose=4] [verbose_integration=0]" << endl;
  cerr << "ddx_pdf_double alpha beta x [pm=0] [verbose=4] [verbose_integration=0]" << endl;
  cerr << "Similarly for ";
  cerr << "cdf_mpreal, etc." << endl;
}

int main(int argc, const char * argv[]) {
  
  // Check the number of parameters
  if (argc < 5 || argc > 8) {
    // Tell the user how to run the program
    show_usage(string(argv[0]));
    return 1;
  }
  
  string test_name = argv[1];
  if (test_name != "cdf_double" &&
      test_name != "pdf_double" &&
      test_name != "ddx_pdf_double" &&
      test_name != "cdf_mpreal" &&
      test_name != "pdf_mpreal" &&
      test_name != "ddx_pdf_mpreal") {show_usage(string(argv[0])); return 1;}
  
  cout << test_name << ":" << endl;
  istringstream ss2((string(argv[2]))), ss3((string(argv[3])));
  double alpha; ss2 >> alpha;
  cout << "alpha = " << alpha;

  double beta; ss3 >> beta;
  cout << ", beta = " << beta;
  
  string x_str{argv[4]};
  trim(x_str);
  double x;
  if (x_str == "-inf")
    x = -std::numeric_limits<double>::infinity();
  else if (x_str == "inf")
    x = std::numeric_limits<double>::infinity();
  else {
    istringstream ss4(x_str);
    ss4 >> x;
  }
  cout << ", x = " << x;

  Parameterization pm = S0;
  if (argc > 5) {
    istringstream ss5((string(argv[5])));
    int tmp;
    ss5 >> tmp;
    pm = tmp ? S1 : S0;
  }
  cout << ", pm = " << pm << endl;

  int verbose = 4;
  if (argc > 6) {
    istringstream ss6((string(argv[6])));
    ss6 >> verbose;
  }
  
  int verbose_integration = 0;
  if (argc > 7) {
    istringstream ss7((string(argv[7])));
    ss7 >> verbose_integration;
  }

  mpfr::mpreal::set_default_prec(128);
  int n_gauss = 10;
  Kronrod<mpreal> g_k_big(n_gauss);
  mpfr::mpreal::set_default_prec(96);
    
  bool noext = true;   //disable extrapolation
  mpreal epsabs = 0.;
  mpreal epsrel = 64*std::numeric_limits<mpreal>::epsilon();
  int subdivisions = 10000;
  int verbose_cntl = 0;
  
  IntegrationController<mpreal> controller(noext, g_k_big, epsabs, epsrel, subdivisions,
                                           verbose_cntl);
  IntegrationController<double> ctl_double(noext, g_k_big,
                                           static_cast<double>(epsabs),
                                           static_cast<double>(epsrel),
                                           subdivisions,
                                           verbose_cntl);
  Controllers<mpreal> ctls(controller, ctl_double);
  Controllers<double> ctls_double(ctl_double, ctl_double);
  
  if (test_name == "cdf_double" || test_name == "pdf_double" || test_name == "ddx_pdf_double")
    calc_result<double>(cout, test_name, alpha, beta, x, pm, ctls_double, verbose, verbose_integration);
  else
    calc_result<mpreal>(cout, test_name, alpha, beta, x, pm, ctls, verbose, verbose_integration);
  return 0;
  
}
