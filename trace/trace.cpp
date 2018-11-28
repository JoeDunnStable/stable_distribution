/// \file trace.cpp
/// Traces for the standard stable distribution
/// \author Joseph Dunn
/// \copyright 2017, 2018 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#define MPREAL
#define CPP_BIN_FLOAT
#define MPFR_FLOAT
#include "stable_distribution.h"
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
  cerr << "cdf_double alpha beta x [pm=0] [lower_tail=1] [verbose = 4]" << endl;
  cerr << "pdf_double alpha beta x [pm=0] [verbose = 4]" << endl;
  cerr << "ddx_pdf_double alpha beta x [pm=0] [verbose=4]" << endl;
  cerr << "quantile_double alpha beta p [pm=0] [lower_tail=1] [verbose=4]" << endl;
  cerr << "mode_double alpha beta [pm=0] [verbose_mode = 4] [verbose_pdf = 0]" << endl;
  cerr << "Similarly for ";
#ifdef MPREAL
  cerr << "cdf_mpreal, ";
#endif
#ifdef CPP_BIN_FLOAT
  cerr << "cdf_cpp_bin_float, ";
#endif
#ifdef MPFR_FLOAT
  cerr << "cdf_mpfr_float, ";
#endif
  cerr << " etc." << endl;
}

int main(int argc, const char * argv[]) {
  // First create some high precision coeficients for
  mpreal::set_default_prec(128);
  Kronrod<mpreal> k_big(10);
  mpreal::set_default_prec(96);
  int noext = 1;
  double epsabs_double = 0;
  double epsrel_double = 64 * std::numeric_limits<double>::epsilon();
  int limit = 1000;
  int verbose_integration = 4;
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
  
  if (test_name.substr(0,3) == "cdf" || test_name.substr(0,8) == "quantile") {
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
        double q_tol = 64 * std::numeric_limits<double>::epsilon();
        result = dist.quantile(x, lower_tail, log_p, q_tol, pm);
      }
      cout << setw(15) << "Result = " << Fmt<double>() << result << endl;
    }
#ifdef MPREAL
    else if (test_name == "cdf_mpreal" || test_name == "quantile_mpreal"){
      StandardStableDistribution<mpreal> dist(alpha, beta, ctls_mpreal, verbose);
      mpreal result;
      cout.precision(bits2digits(mpfr::mpreal::get_default_prec()));
      if (test_name == "cdf_mpreal") {
        result = dist.cdf(x, lower_tail, log_p, pm);
      } else {
        mpreal q_tol = 64 * std::numeric_limits<mpreal>::epsilon();
        result = dist.quantile(x, lower_tail, log_p, q_tol, pm);
      }
      cout << setw(15) << "Result = " << Fmt<mpreal>() << result << endl;
    }
#endif //MPREAL
#ifdef CPP_BIN_FLOAT
    else if (test_name == "cdf_cpp_bin_float" || test_name == "quantile_cpp_bin_float"){
      StandardStableDistribution<CppBinFloat> dist(alpha, beta, ctls_cpp_bin_float, verbose);
      CppBinFloat result;
      cout.precision(std::numeric_limits<CppBinFloat>::digits10);
      if (test_name == "cdf_cpp_bin_float") {
        result = dist.cdf(x, lower_tail, log_p, pm);
      } else {
        CppBinFloat q_tol = 64 * std::numeric_limits<CppBinFloat>::epsilon();
        result = dist.quantile(x, lower_tail, log_p, q_tol, pm);
      }
      cout << setw(15) << "Result = " << Fmt<CppBinFloat>() << result << endl;
    }
#endif //CPP_BIN_FLOAT
#ifdef MPFR_FLOAT
    else if (test_name == "cdf_mpfr_float" || test_name == "quantile_mpfr_float"){
      StandardStableDistribution<MpfrFloat> dist(alpha, beta, ctls_mpfr_float, verbose);
      MpfrFloat result;
      cout.precision(std::numeric_limits<MpfrFloat>::digits10);
      if (test_name == "cdf_mpfr_float") {
        result = dist.cdf(x, lower_tail, log_p, pm);
      } else {
        MpfrFloat q_tol = 64 * std::numeric_limits<MpfrFloat>::epsilon();
        result = dist.quantile(x, lower_tail, log_p, q_tol, pm);
      }
      cout << setw(15) << "Result = " << Fmt<MpfrFloat>() << result << endl;
    }
#endif //MPFR_FLOAT
  } else if (test_name.substr(0,3) == "pdf" || test_name.substr(0,7) == "ddx_pdf") {
    if (argc < 5) { show_usage(string(argv[0])); return 1;}
	size_t len = test_name.length();
    istringstream ss4((string(argv[4])));
    double x; ss4 >> x;
	if (len >= 6 && test_name.substr(len-6,6)=="double")
      cout << ", x = " << x;
/*
#ifdef MPREAL
	istringstream ss4_mpreal((string(argv[4])));
	mpreal x_mpreal; ss4_mpreal >> x_mpreal;
	if (len >=6 && test_name.substr(len - 6, 6) == "mpreal")
	  cout << ", x_mpreal = " << x_mpreal;
#endif
#ifdef MPFR_FLOAT
	istringstream ss4_mpfr_float((string(argv[4])));
	MpfrFloat x_mpfr_float; ss4_mpfr_float >> x_mpfr_float;
	if (len >= 10 && test_name.substr(len - 10, 10) == "mpfr_float")
		cout << ", x_mpfr_float = " << x_mpfr_float;
#endif
#ifdef CPP_BIN_FLOAT
	istringstream ss4_cpp_bin_float((string(argv[4])));
	CppBinFloat x_cpp_bin_float;  ss4_cpp_bin_float >> x_cpp_bin_float;
	if (len >= 13 && test_name.substr(len - 13, 13) == "cpp_bin_float")
		cout << ", x_cpp_bin_float = " << x_cpp_bin_float;
#endif
*/
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
    if (test_name == "pdf_double" || test_name == "ddx_pdf_double") {
      StandardStableDistribution<double> dist(alpha, beta, ctls_double, verbose);
      double result;
      cout.precision(std::numeric_limits<double>::digits10);
      if (test_name == "pdf_double") {
        result = dist.pdf(x, log_p, pm);
      } else {
        result = dist.ddx_pdf(x, pm);
      }
      cout << setw(15) << "Result = " << Fmt<double>() << result << endl;
    }
#ifdef MPREAL
    else if( test_name == "pdf_mpreal" || test_name == "ddx_pdf_mpreal") {
      StandardStableDistribution<mpreal> dist(alpha, beta, ctls_mpreal, verbose);
      mpreal result;
      cout.precision(bits2digits(mpfr::mpreal::get_default_prec()));
      if (test_name == "pdf_mpreal") {
        result = dist.pdf(x, log_p, pm);
      } else {
        result = dist.ddx_pdf(x, pm);
      }
      cout << setw(15) << "Result = " << Fmt<mpreal>() << result << endl;
    }
#endif // MPREAL
#ifdef CPP_BIN_FLOAT
    else if( test_name == "pdf_cpp_bin_float" || test_name == "ddx_pdf_cpp_bin_float") {
      StandardStableDistribution<CppBinFloat> dist(alpha, beta, ctls_cpp_bin_float, verbose);
      CppBinFloat result;
      cout.precision(std::numeric_limits<CppBinFloat>::digits10);
      if (test_name == "pdf_cpp_bin_float") {
        result = dist.pdf(x, log_p, pm);
      } else {
        result = dist.ddx_pdf(x, pm);
      }
      cout << setw(15) << "Result = " << Fmt<CppBinFloat>() << result << endl;
    }
#endif // CPP_BIN_FLOAT
#ifdef MPFR_FLOAT
    else if( test_name == "pdf_mpfr_float" || test_name == "ddx_pdf_mpfr_float") {
      StandardStableDistribution<MpfrFloat> dist(alpha, beta, ctls_mpfr_float, verbose);
      MpfrFloat result;
      cout.precision(std::numeric_limits<MpfrFloat>::digits10);
      if (test_name == "pdf_mpfr_float") {
        result = dist.pdf(x, log_p, pm);
      } else {
        result = dist.ddx_pdf(x, pm);
      }
      cout << setw(15) << "Result = " << Fmt<MpfrFloat>() << result << endl;
    }
#endif // MPFRFLOAT
    else {
      cerr << "Invalid test name: " << test_name << endl;
      show_usage(string(argv[0])); return 1;
    }
  } else if (test_name.substr(0,4) == "mode") {
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
      cout << setw(15) << "Result = (" << Fmt<double>() << result.first << ", " << result.second << ")"<< endl;
      cout << setw(15) << "ddx at result = " << Fmt<double>() << dist.ddx_pdf(result.first, pm) << endl;
    }
#ifdef MPREAL
    else if (test_name == "mode_mpreal") {
      StandardStableDistribution<mpreal> dist(alpha, beta, ctls_mpreal, verbose_pdf);
      std::pair<mpreal, mpreal> result;
      cout.precision(bits2digits(mpfr::mpreal::get_default_prec()));
      mpreal mode_tol = 64 * std::numeric_limits<mpreal>::epsilon();
      result = dist.mode(mode_tol, verbose_mode, pm);
      cout << setw(15) << "Result = (" << Fmt<mpreal>() << result.first << ", " << result.second << ")"<< endl;
      cout << setw(15) << "ddx at result = " << Fmt<mpreal>() << dist.ddx_pdf(result.first, pm) << endl;
    }
#endif // MPREAL
#ifdef CPP_BIN_FLOAT
    else if (test_name == "mode_cpp_bin_float") {
      StandardStableDistribution<CppBinFloat> dist(alpha, beta, ctls_cpp_bin_float, verbose_pdf);
      std::pair<CppBinFloat, CppBinFloat> result;
      cout.precision(std::numeric_limits<CppBinFloat>::digits10);
      CppBinFloat mode_tol = 64 * std::numeric_limits<CppBinFloat>::epsilon();
      result = dist.mode(mode_tol, verbose_mode, pm);
      cout << setw(15) << "Result = (" << Fmt<CppBinFloat>() << result.first << ", " << result.second << ")"<< endl;
      cout << setw(15) << "ddx at result = " << Fmt<CppBinFloat>() << dist.ddx_pdf(result.first, pm) << endl;
    }
#endif // CPP_BIN_FLOAT
#ifdef MPFR_FLOAT
    else if (test_name == "mode_mpreal") {
      StandardStableDistribution<MpfrFloat> dist(alpha, beta, ctls_mpfr_float, verbose_pdf);
      std::pair<MpfrFloat, MpfrFloat> result;
      cout.precision(std::numeric_limits<MpfrFloat>::digits10);
      MpfrFloat mode_tol = 64 * std::numeric_limits<MpfrFloat>::epsilon();
      result = dist.mode(mode_tol, verbose_mode, pm);
      cout << setw(15) << "Result = (" << Fmt<MpfrFloat>() << result.first << ", " << result.second << ")"<< endl;
      cout << setw(15) << "ddx at result = " << Fmt<MpfrFloat>() << dist.ddx_pdf(result.first, pm) << endl;
    }
#endif // MPFR_FLOAT
    else {
      cerr << "Invalid test name: " << test_name << endl;
      show_usage(string(argv[0]));
      return 1;
    }
  } else {
    cerr << "Invalid test name: " << test_name << endl;
    show_usage(string(argv[0]));
    return 1;
  }

  return 0;
}
