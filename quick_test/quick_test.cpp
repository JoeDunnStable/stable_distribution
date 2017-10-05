//
/// \file quick_test.cpp
/// \author Joseph Dunn
/// \copyright 2016, 2017 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <chrono>
#include "stable_distribution_fit.h"
using Vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

//#define BOOST_MATH_In_gaussSTRUMEn_gaussT
#include <boost/math/tools/toms748_solve.hpp>

using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::right;
using std::setprecision;
using std::fixed;
using std::scientific;
using std::ofstream;
using std::string;
using std::stringstream;
using std::sort;
using boost::math::tools::toms748_solve;
using std::chrono::high_resolution_clock;
using std::chrono::duration;

using namespace stable_distribution;

static void show_usage (string name){
  cerr << "Usage: " << name << "" << endl;
}

class quick_result {
public:
  double alpha;
  double beta;
  double p;
  double q;
  double d_exact;
  unsigned int reg;
  double pt;
  double d_quick;
  double delta;
  quick_result(double alpha, double beta, double p, double q, double d_exact,
               unsigned int reg, double pt, double d_quick) :
  alpha(alpha), beta(beta), p(p), q(q), d_exact(d_exact),
  reg(reg), pt(pt), d_quick(d_quick), delta(fabs(d_exact-d_quick)) {}
};

bool cmp_result(const quick_result& lhs, const quick_result& rhs) {
  return lhs.delta > rhs.delta;
}

int main(int argc, char *argv[]) {
  
  StandardStableDistribution<double>::initialize();
  
  // Check the number of parameters
  if (argc > 1) {
    // Tell the user how to run the program
    show_usage(string(argv[0]));
    return 1;
  }
  
  string out_file = "../output/quick_test.out";
  cout << "Writing output to " + out_file << endl;
  ofstream out(out_file);
  
  Vec alphas(6);    alphas << .1, .5, 1, 1.5, 1.99, 2;
  Vec betas(5);     betas << -1., -.5, 0, .5, 1;
  
  int n_gauss = 10;   // Use 10point Gauss, 21 point Kronrod rule
  mpreal::set_default_prec(64);
  Kronrod<mpreal> g_k_big(n_gauss);
  bool noext = true;   //disable extrapolation
  double epsabs = 0.;
  double epsrel = 64*std::numeric_limits<double>::epsilon();
  int subdivisions = 1000;
  IntegrationController<double> ctrl(noext, g_k_big , epsabs, epsrel, subdivisions, 0);
  int verbose = 0;
  Controllers<double> ctls(ctrl, ctrl);
  unsigned int n = 2000;
  Vec probs = (Vec::LinSpaced(n,0,n-1).array()+.5)/n;
  Vec gamma; gamma.setOnes(n);
  Vec delta; delta.setZero(n);
  Parameterization pm = S0;
  vector<quick_result> results;
  out << setw(7) << right << "alpha"
  << setw(7) << right << "beta"
  << setw(20) << right << "pdf_exact"
  << setw(20) << right << "pdf_quick"
  << setw(15) << right << "delta"
  << setw(12) << right << "t_exact ms"
  << setw(12) << right << "t_quick ms"
  << endl << endl;
  bool pass = true;
  for (int ia=0; ia<alphas.size(); ia++) {
    for (int ib=0; ib<betas.size(); ib++) {
      Vec q =quantile(probs , alphas(ia), betas(ib), gamma, delta, pm, true, false,
                      1e-9, ctls, verbose);
      StandardStableDistribution<double> std_stable_dist(alphas(ia), betas(ib), ctls, verbose);
      DstableQuick<double> std_stable_dist_quick(&std_stable_dist);
      high_resolution_clock::time_point t0 = high_resolution_clock::now();
      Vec d_exact = pdf(q, alphas(ia), betas(ib), gamma, delta, pm, true, ctls, verbose);
      high_resolution_clock::time_point t1 = high_resolution_clock::now();
      Vec d_quick = std_stable_dist_quick(q);
      high_resolution_clock::time_point t2 = high_resolution_clock::now();
      duration<double, std::milli> time_exact = (t1-t0);
      duration<double, std::milli> time_quick = (t2-t1);
      uVec reg = std_stable_dist_quick.region(q);
      Vec pt = std_stable_dist_quick.pt(q);
      double delta = fabs(d_exact.sum()-d_quick.sum())/n;
      bool pass1 = delta < .0001;
      pass = pass && pass1;
      for (int i=0; i<n; i++)
        results.push_back(quick_result(alphas(ia), betas(ib), probs(i), q(i), d_exact(i), reg(i), pt(i), d_quick(i)));
      out << setw(7) << setprecision(2) <<  fixed << alphas(ia)
           << setw(7) << setprecision(1) << betas(ib)
           << setw(20) << setprecision(10) << fixed << d_exact.sum()/n
           << setw(20) << setprecision(10) << d_quick.sum()/n
           << setw(15) << setprecision(10) << fabs(d_exact.sum()-d_quick.sum())/n
           << setw(12) << setprecision(2) << time_exact.count()
           << setw(12) << setprecision(2) << time_quick.count()
           << (pass1 ? "" : " FAIL") << endl;
    }
  }
  
  sort(results.begin(), results.end(), cmp_result);
  out << endl << "The hundred results with the biggest discrepancy" << endl << endl;
  out << setw(7) << right << "alpha"
  << setw(7) << right << "beta"
  << setw(20) << right << "p"
  << setw(20) << right << "q"
  << setw(20) << right << "pdf_exact"
  << setw(7) << right << "region"
  << setw(20) << right << "pt"
  << setw(20) << right << "pdf_quick"
  << setw(10) << right << "delta" << endl << endl;
  for (int i=0; i<100; i++)
    out << setw(7) << setprecision(2) <<  fixed << results.at(i).alpha
    << setw(7) << setprecision(1) << results.at(i).beta
    << setw(20) << setprecision(8) << results.at(i).p
    << setw(20) << setprecision(10) << scientific << results.at(i).q
    << setw(20) << setprecision(10) << fixed << results.at(i).d_exact
    << setw(7) << results.at(i).reg
    << setw(20) << setprecision(10) << results.at(i).pt
    << setw(20) << setprecision(10) << results.at(i).d_quick
    << setw(10) << setprecision(5) << results.at(i).delta << endl;
  
  return !pass;
}
