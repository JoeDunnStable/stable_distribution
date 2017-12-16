/// \file trace_zolotarev.cpp
/// Traces for the routines in zolotarev.h
/// \author Joseph Dunn
/// \copyright 2017 Joseph Dunn
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
                 Controllers<myFloat> ctls, int verbose, int verbose_integration) {
  os << "Calculating " << type << endl << endl;
  if (type == "cdf_double" || type == "pdf_double" || type == "ddx_pdf_double")
    os.precision(15);
  else
    os.precision(25);
  int log_flag=0, lower_tail=x<0;
  myFloat stable_error, zol_error;
  myFloat ret, res;
  Zolotarev<myFloat> zol(myFloat(alpha), myFloat(beta), &ctls.controller,
                        verbose, verbose_integration) ;
  StandardStableDistribution<myFloat> std_dist(myFloat(alpha), myFloat(beta), ctls, 1);
  if (type == "cdf_double" || type == "cdf_mpreal") {
    ret = zol.cdf(myFloat(x), lower_tail, S1) ;
    res = std_dist.cdf(myFloat(x), lower_tail, log_flag, S1);
    zol_error = zol.error;
    stable_error= std_dist.abserr;
  } else if (type == "pdf_double" || type == "pdf_mpreal"){
    ret = zol.pdf(myFloat(x), S1) ;
    res = std_dist.pdf(myFloat(x), log_flag, S1);
    zol_error = zol.error;
    stable_error= std_dist.abserr;
  } else {
    ret = zol.ddx_pdf(myFloat(x), S1) ;
    res = std_dist.ddx_pdf(myFloat(x), S1);
    zol_error = zol.error;
    stable_error= std_dist.abserr;

  }
  os.unsetf(std::ios_base::floatfield);
  os << setw(15) << "alpha = " << setw(6)<< setprecision(3) << alpha << endl
     << setw(15) << "beta = " << setw(6)<<setprecision(3) << beta << endl
     << setw(15) << "x = " << setw(6) << setprecision(2) << x << endl
     << setw(15) << "good theta2 = "<< setw(2) << std_dist.good_theta2 << endl
     << setw(15) << "result stable = " << setw(25) << setprecision(16) << scientific  << res << endl
     << setw(15) << "error stable = " << setw(12) << setprecision(5) << scientific  << stable_error << endl
     << setw(15) << "n zol = " << setw(4) << zol.n << endl
     << setw(15) << "result zol = " << setw(25) << scientific << setprecision(16) << ret << endl
     << setw(15) << "error zol = " << setw(12) << scientific << setprecision(5) << zol_error << endl
     << setw(15) << "rel error = " << setw(12) << scientific << setprecision(5) << ret/res-1 << endl;
}




static void show_usage (string name){
  boost::filesystem::path p(name);
  cerr << "Usage: " << p.filename().string() << "type[cdf_double, pdf_double, ddx_pdf_double, cdf_mpreal, pdf_mpreal, or ddx_pdf_mpreal]  alpha beta x [verbose = 1] [verbose_integration=0]" << endl;
}

int main(int argc, const char * argv[]) {
  
  // Check the number of parameters
  if (argc < 5 || argc > 7) {
    // Tell the user how to run the program
    show_usage(string(argv[0]));
    return 1;
  }
  
  istringstream ss1{string(argv[1])}, ss2{string(argv[2])}, ss3{string(argv[3])},
  ss4{string(argv[4])};
  string type; ss1 >> type;
  if (type != "cdf_double" && type != "pdf_double" && type != "ddx_pdf_double"
      && type != "cdf_mpreal" && type!= "pdf_mpreal" && type != "ddx_pdf_mpreal") {
    // Tell the user how to run the program
    show_usage(string(argv[0]));
    return 1;
  }
  double alpha; ss2 >> alpha;
  double beta; ss3 >> beta;
  double x; ss4 >> x;
  int verbose = 1;
  if (argc > 5) {
    istringstream ss5((string(argv[5])));
    ss5 >> verbose;
  }
  int verbose_integration = 0;
  if (argc > 6) {
    istringstream ss6((string(argv[6])));
    ss6 >> verbose_integration;
  }
  
  using namespace stable_distribution;
  
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
  
  if (type == "cdf_double" || type == "pdf_double" || type == "ddx_pdf_double")
    calc_result<double>(cout, type, alpha, beta, x, ctls_double, verbose, verbose_integration);
  else
    calc_result<mpreal>(cout, type, alpha, beta, x, ctls, verbose, verbose_integration);
  
  return 0;
  
}
