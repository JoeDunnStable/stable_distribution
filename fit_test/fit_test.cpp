//
/// \file fit_test.cpp
/// \author Joseph Dunn
/// \copyright 2016, 2017 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <sstream>
using std::stringstream;

#include <fstream>
using std::ofstream;

#include <string>
using std::string;

#include <random>
using std::mt19937;
using std::uniform_real_distribution;
#include <boost/filesystem.hpp>

#include <mpreal.h>
using mpfr::mpreal;
#include "stable_distribution_fit.h"
using Vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

using namespace stable_distribution;

static void show_usage (string name){
  boost::filesystem::path p(name);
  cerr << "Usage: " << p.filename().string() << " type(q, q_mle,  or mle)" << endl;
}

class NullBuffer : public std::streambuf
{
public:
    int overflow(int c) { return c; }
};

int main(int argc, char *argv[]) {
  
  StandardStableDistribution<double>::initialize();
  
  // Check the number of parameters
  if (argc != 2) {
    // Tell the user how to run the program
    show_usage(string(argv[0]));
    return 1;
  }
  stringstream ss1(argv[1]);
  string type;
  ss1 >> type;
  if (type != "q" && type != "q_mle" && type != "mle") {
    show_usage(string(argv[0]));
  }
  string out_file = "../output/test_stable_fit.out";
  cout << "Writing output to " + out_file << endl;
  ofstream out(out_file);

  NullBuffer null_buffer;
  ostream null_stream(&null_buffer);
  mt19937 gen(200);
  uniform_real_distribution<> dis;
  Vec alphas(6);    alphas << .5, .9, 1., 1.1, 1.5, 2.;
  Vec betas(5);     betas << -1., -.5, 0, .5, 1;
  
  mpreal::set_default_prec(96);
  int n_gauss = 10;   // Use 10point Gauss, 21 point Kronrod rule
  Kronrod<mpreal> g_k_big(n_gauss);  //high precision nodes and weights
  bool noext = true;   //disable extrapolation
  double epsabs = 0.;
  double epsrel = 64*std::numeric_limits<double>::epsilon();
  int subdivisions = 1000;
  IntegrationController<double> ctrl(noext, g_k_big, epsabs, epsrel, subdivisions, 0);
  Controllers<double> ctls(ctrl, ctrl);
  int verbose = 0;
  bool quick = true;
  int n = 10000;
  Vec r(n);
  Vec probs(5); probs << .05, .25, .5, .75, .95;
  result_heading(out);
  bool pass = true;
  for (int ia=0; ia<alphas.size(); ia++) {
    for (int ib=0; ib<betas.size(); ib++) {
      gen.seed(200);
      for (int i=0; i<n; i++) {
        double u1 = dis(gen);
        double u2 = dis(gen);
        r(i) = random_stable(alphas(ia), betas(ib), u1, u2);
      }
      double actual_value = capped_pdf(r, alphas(ia), betas(ib), 1. , 0.,
                                       quick, ctls, verbose);
      Vec q = quantile(r, probs);
      FitResult<double> fr_actual("Actual", alphas(ia), betas(ib), 1., 0., -actual_value, n, q, "NA", 0, 0.);
      out << fr_actual;
      vector<FitResult<double> > results;
      results = stable_fit(r, ctls, 1e-12, type, quick, 0 );
      for (int ir=0; ir<results.size(); ++ir) {
        string method = results.at(ir).method;
        double tol = (method == "Dunn" || method == "McCulloch") ? .05 : .02;
        bool pass1 = fabs(results.at(ir).two_ll_n + actual_value) < tol * fabs(actual_value);
        out << results.at(ir);
        pass = pass & pass1;
        if (!pass1) out << "FAIL" << endl;
      }
    }
    
  }
  return !pass;
  
}
