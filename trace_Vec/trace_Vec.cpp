//
/// \file trace_Vec.cpp
/// \author Joseph Dunn
/// \copyright 2017 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#define MPREAL
#include "stable_distribution_Vec.h"
#include <iomanip>
#include <vector>
#include <boost/filesystem.hpp>

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
  cerr << "cdf_double alpha beta gamma delta x [pm=0] [lower_tail=1] [verbose = 4]" << endl;
  cerr << "pdf_double alpha beta gamma delta x [pm=0] [verbose = 4]" << endl;
  cerr << "ddx_pdf_double alpha beta gamma delta x [pm=0] [verbose=4]" << endl;
  cerr << "quantile_double alpha beta gamma delta p [pm=0] [lower_tail=1] [verbose=4]" << endl;
  cerr << "Similarly for cdf_mpreal, etc." << endl;
}

int main(int argc, const char * argv[]) {
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
  
  Controllers<double> ctls_double(ctl_double, ctl_double);
  Controllers<mpreal> ctls_mpreal(ctl_mpreal, ctl_double);
  
  // Check the number of parameters
  if (argc < 6 || argc > 10) {
    // Tell the user how to run the program
    show_usage(string(argv[0]));
    return 1;
  }
  
  string test_name = argv[1];
  cout << test_name << ":" << endl;
  istringstream ss2((string(argv[2]))), ss3((string(argv[3])));
  double alpha; ss2 >> alpha;
  cout << "alpha = " << alpha << endl;
  double beta; ss3 >> beta;
  cout << "beta = " << beta << endl;
  istringstream ss4((string(argv[4]))), ss5((string(argv[5])));
  double gamma; ss4 >> gamma;
  cout << "gamma = " << gamma << endl;
  Matrix<double, Dynamic, 1> gamma_double(1); gamma_double(0) = gamma;
  Matrix<mpreal, Dynamic, 1> gamma_mpreal(1); gamma_mpreal(0) = gamma;
  double delta; ss5 >> delta;
  cout << "delta = " << delta << endl;
  Matrix<double, Dynamic, 1> delta_double(1); delta_double(0) = delta;
  Matrix<mpreal, Dynamic, 1> delta_mpreal(1); delta_mpreal(0) = delta;

  int log_p = 0;
  
  if (test_name == "cdf_double" || test_name == "cdf_mpreal" ||
      test_name == "quantile_double" || test_name == "quantile_mpreal") {
    if (argc < 7) { show_usage(string(argv[0])); return 1;}
    istringstream ss6((string(argv[6])));
    double x; ss6 >> x;
    cout << "x = " << x << endl;
    Matrix<double, Dynamic, 1> x_double(1); x_double(0) = x;
    Matrix<mpreal, Dynamic, 1> x_mpreal(1); x_mpreal(0) = x;
    int pm = 0;
    if (argc > 7) {
      istringstream ss7((string(argv[7])));
      ss7 >> pm;
    }
    cout << "pm = " << pm << endl;
    int lower_tail = 1;
    if (argc > 8) {
      istringstream ss8((string(argv[8])));
      ss8 >> lower_tail;
    }
    cout << "lower_tail = " << lower_tail << endl;
    int verbose = 4;
    if (argc > 9) {
      istringstream ss9((string(argv[9])));
      ss9 >> verbose;
    }
    cout << "verbose = " << verbose << endl;
    if (test_name == "cdf_double" || test_name == "quantile_double") {
      Matrix<double, Dynamic, 1> result;
      cout.precision(std::numeric_limits<double>::digits10);
      if (test_name == "cdf_double") {
        result = cdf(x_double, alpha, beta, gamma_double, delta_double,
                     pm, lower_tail, log_p, ctls_double, verbose);
      } else {
        double q_tol = 64*std::numeric_limits<double>::epsilon();
        result = quantile(x_double, alpha, beta, gamma_double, delta_double,
                          pm, lower_tail, log_p, q_tol, ctls_double, verbose);
      }
      cout << endl << "Result = " << result << endl;
    } else {
      Matrix<mpreal, Dynamic, 1> result;
      cout.precision(bits2digits(mpfr::mpreal::get_default_prec()));
      if (test_name == "cdf_mpreal") {
        result = cdf(x_mpreal, mpreal(alpha), mpreal(beta), gamma_mpreal, delta_mpreal,
                     pm, lower_tail, log_p, ctls_mpreal, verbose);
      } else {
        mpreal q_tol = 64*std::numeric_limits<mpreal>::epsilon();
        result = quantile(x_mpreal, mpreal(alpha), mpreal(beta), gamma_mpreal, delta_mpreal,
                          pm, lower_tail, log_p, q_tol, ctls_mpreal, verbose);;
      }
      cout << endl << "Result = " << result << endl;
    }
  } else if (test_name == "pdf_double" || test_name == "pdf_mpreal" ||
             test_name == "ddx_pdf_double" || test_name == "ddx_pdf_mpreal") {
    if (argc < 7) { show_usage(string(argv[0])); return 1;}
    istringstream ss6((string(argv[6])));
    double x; ss6 >> x;
    cout << "x = " << x << endl;
    Matrix<double, Dynamic, 1> x_double(1); x_double(0) = x;
    Matrix<mpreal, Dynamic, 1> x_mpreal(1); x_mpreal(0) = x;
    int pm = 0;
    if (argc > 7) {
      istringstream ss7((string(argv[7])));
      ss7 >> pm;
    }
    cout << "pm = " << pm << endl;
    int verbose = 4;
    if (argc > 8) {
      istringstream ss8((string(argv[8])));
      ss8 >> verbose;
    }
    cout << "verbose = " << verbose << endl;
    if (test_name == "pdf_double" || test_name == "ddx_pdf_double") {
      Matrix<double, Dynamic, 1> result;
      cout.precision(std::numeric_limits<double>::digits10);
      if (test_name == "pdf_double") {
        result = pdf(x_double, alpha, beta, gamma_double, delta_double,
                     pm, log_p, ctls_double, verbose);
      } else {
        result = ddx_pdf(x_double, alpha, beta, gamma_double, delta_double,
                     pm, ctls_double, verbose);
      }
      cout << endl << "Result = " << result << endl;
    } else {
      Matrix<mpreal, Dynamic, 1> result;
      cout.precision(bits2digits(mpfr::mpreal::get_default_prec()));
      if (test_name == "pdf_mpreal") {
        result = pdf(x_mpreal, mpreal(alpha), mpreal(beta), gamma_mpreal, delta_mpreal,
                     pm, log_p, ctls_mpreal, verbose);
      } else {
        result = ddx_pdf(x_mpreal, mpreal(alpha), mpreal(beta), gamma_mpreal, delta_mpreal,
                         pm, ctls_mpreal, verbose);
      }
      cout << endl << "Result = " << result << endl;
    }
  } else {
    cerr << "Incorrect test name: " << test_name << endl;
    show_usage(string(argv[0]));
    return 1;
  }
  
  return 0;
}

