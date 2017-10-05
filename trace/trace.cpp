//
/// \file trace_stable_mode.cpp
/// \author Joseph Dunn
/// \copyright 2016 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include "stable_distribution.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <boost/filesystem.hpp>

using std::cout;
using std::endl;
using std::setw;
using std::setprecision;
using std::right;
using std::fixed;
using std::scientific;
using std::istringstream;

using std::vector;

using mpfr::mpreal;
using mpfr::bits2digits;

using namespace stable_distribution;

static void show_usage (string name){
  boost::filesystem::path p(name);
  cerr << endl;
  cerr << "Usage: " << p.filename().string() << " test_name ..." << endl;
  cerr << "cdf_double alpha beta x [pm=0] [lower_tail=1] [verbose = 4]" << endl;
  cerr << "pdf_double alpha beta x [pm=0] [verbose = 4]" << endl;
  cerr << "ddx_pdf_double alpha beta x [pm=0] [verbose=4]" << endl;
  cerr << "quantile_double alpha beta p [pm=0] [lower_tail=1] [verbose=4]" << endl;
  cerr << "mode_double alpha beta [pm=0] [verbose_mode = 4] [verbose_pdf = 0]" << endl;
  cerr << "Similarly for cdf_mpreal, etc." << endl;
}

int main(int argc, const char * argv[]) {
  // First create some high precision coeficients for
  mpreal::set_default_prec(128);
  Kronrod<mpreal> k_big(10);
  mpreal::set_default_prec(96);
  int digits = bits2digits(mpfr::mpreal::get_default_prec());
  int noext = 1;
  mpreal eps_mpreal = std::numeric_limits<mpreal>::epsilon();
  double epsabs_double = 0;
  double epsrel_double = 64 * std::numeric_limits<double>::epsilon();
  int limit = 1000;
  int verbose_integration = 4;
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
  
  Controllers<double> ctls_double(ctl_double, ctl_double);
  Controllers<mpreal> ctls_mpreal(ctl_mpreal, ctl_double);
  
  // Check the number of parameters
  if (argc < 4 || argc > 8) {
    // Tell the user how to run the program
    show_usage(string(argv[0]));
    return 1;
  }
  
  string test_name = argv[1];
  cout << test_name << ":" << endl;
  istringstream ss2((string(argv[2]))), ss3((string(argv[3])));
  double alpha; ss2 >> alpha;
  cout << "alpha = " << alpha;
  double beta; ss3 >> beta;
  cout << ", beta = " << beta;
  int log_p = 0;
  
  if (test_name == "cdf_double" || test_name == "cdf_mpreal" ||
      test_name == "quantile_double" || test_name == "quantile_mpreal") {
    if (argc < 5) { show_usage(string(argv[0])); return 1;}
    istringstream ss4((string(argv[4])));
    double x; ss4 >> x;
    cout << ", x = " << x;
    Parameterization pm = S0;
    if (argc > 5) {
      istringstream ss5((string(argv[5])));
      int tmp;
      ss5 >> tmp;
      pm = tmp ? S1 : S0;
    }
    cout << ", pm = " << pm;
    int lower_tail = 1;
    if (argc > 6) {
      istringstream ss6((string(argv[6])));
      ss6 >> lower_tail;
    }
    cout << ", lower_tail = " << lower_tail << endl;
    int verbose = 4;
    if (argc > 7) {
      istringstream ss7((string(argv[7])));
      ss7 >> verbose;
    }
    if (test_name == "cdf_double" || test_name == "quantile_double") {
      StandardStableDistribution<double> dist(alpha, beta, ctls_double, verbose);
      double result;
      cout.precision(std::numeric_limits<double>::digits10);
      if (test_name == "cdf_double") {
        result = dist.cdf(x, lower_tail, log_p, pm);
      } else {
        result = dist.quantile(x, lower_tail, log_p, pm);
      }
      cout << setw(15) << "Result = " << result << endl;
    } else {
      StandardStableDistribution<mpreal> dist(alpha, beta, ctls_mpreal, verbose);
      mpreal result;
      cout.precision(bits2digits(mpfr::mpreal::get_default_prec()));
      if (test_name == "cdf_mpreal") {
        result = dist.cdf(x, lower_tail, log_p, pm);
      } else {
        result = dist.quantile(x, lower_tail, log_p, pm);
      }
      cout << setw(15) << "Result = " << result << endl;
    }
  } else if (test_name == "pdf_double" || test_name == "pdf_mpreal" ||
             test_name == "ddx_pdf_double" || test_name == "ddx_pdf_mpreal") {
    if (argc < 5) { show_usage(string(argv[0])); return 1;}
    istringstream ss4((string(argv[4])));
    double x; ss4 >> x;
    cout << ", x = " << x;
    Parameterization pm = S0;
    if (argc > 5) {
      istringstream ss5((string(argv[5])));
      int tmp;
      ss5 >> tmp;
      pm = tmp ? S1 : S0;
    }
    cout << ", pm = " << pm;
    int verbose = 4;
    if (argc > 6) {
      istringstream ss6((string(argv[6])));
      ss6 >> verbose;
    }
    if (test_name == "pdf_double" || test_name == "ddx_pdf_double") {
      StandardStableDistribution<double> dist(alpha, beta, ctls_double, verbose);
      double result;
      cout.precision(std::numeric_limits<double>::digits10);
      if (test_name == "pdf_double") {
        result = dist.pdf(x, log_p, pm);
      } else {
        result = dist.ddx_pdf(x, pm);
      }
      cout << setw(15) << "Result = " << result << endl;
    } else {
      StandardStableDistribution<mpreal> dist(alpha, beta, ctls_mpreal, verbose);
      mpreal result;
      cout.precision(bits2digits(mpfr::mpreal::get_default_prec()));
      if (test_name == "pdf_mpreal") {
        result = dist.pdf(x, log_p, pm);
      } else {
        result = dist.ddx_pdf(x, pm);
      }
      cout << setw(15) << "Result = " << result << endl;
    }
  } else if (test_name == "mode_double" || test_name == "mode_mpreal") {
    if (argc < 4) { show_usage(string(argv[0])); return 1;}
    Parameterization pm = S0;
    if (argc > 4) {
      istringstream ss4((string(argv[4])));
      int tmp;
      ss4 >> tmp;
      pm = tmp ? S1 : S0;
    }
    cout << ", pm = " << pm;
    int verbose_mode = 4;
    if (argc > 5) {
      istringstream ss5((string(argv[5])));
      ss5 >> verbose_mode;
    }
    int verbose_pdf = 0;
    if (argc > 6) {
      istringstream ss6((string(argv[6])));
      ss6 >> verbose_pdf;
    }
    if (test_name == "mode_double") {
      StandardStableDistribution<double> dist(alpha, beta, ctls_double, verbose_pdf);
      std::pair<double, double> result;
      cout.precision(std::numeric_limits<double>::digits10);
      double mode_tol = 64*std::numeric_limits<double>::epsilon();
      result = dist.mode(mode_tol, verbose_mode, pm);
      cout << setw(15) << "Result = (" << result.first << ", " << result.second << ")"<< endl;
    } else {
      StandardStableDistribution<mpreal> dist(alpha, beta, ctls_mpreal, verbose_pdf);
      std::pair<mpreal, mpreal> result;
      cout.precision(bits2digits(mpfr::mpreal::get_default_prec()));
      mpreal mode_tol = 64 * std::numeric_limits<mpreal>::epsilon();
      result = dist.mode(mode_tol, verbose_mode, pm);
      cout << setw(15) << "Result = (" << result.first << ", " << result.second << ")"<< endl;

    }
  } else {
    cerr << "Incorrect test name: " << test_name << endl;
    show_usage(string(argv[0]));
    return 1;
  }
  

  cout << "tan(StandardStableDistibution<double>::pi/4) - 1= "<< setprecision(16) << tan(StandardStableDistribution<double>::pi/4)-1 << endl;
  cout << "tan(StandardStableDistibution<mpreal>::pi/4) - 1 = "<< setprecision(16) << tan(StandardStableDistribution<mpreal>::pi/4) - 1 << endl;
  return 0;
}
