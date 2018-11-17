/// \file FMStable_test.cpp
/// Unit test for standard stable distribution compare to FMStable
/// \author Joseph Dunn
/// \copyright 2016, 2017 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include <mpreal.h>
using mpfr::mpreal;
#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
#include "stable_config.h"
#include "stable_distribution.h"
#include <iomanip>
#include <sstream>
#include <fstream>
#include <limits>
#include <boost/timer/timer.hpp>
using boost::timer::auto_cpu_timer;

#include <boost/filesystem.hpp>

// This is the only function that we'll use from FMStable
  void tailsMSS(int n,double x[],double d[],double logd[],double F[],
                double logF[],double cF[],double logcF[],
                double alpha,double oneminusalpha, double twominusalpha,
                double location,double logscale);

using std::string;
using std::ofstream;
using std::setw;
using std::setprecision;
using std::right;
using std::fixed;
using std::scientific;
using std::max;

using namespace stable_distribution;


class error_comp {
public:
  double x;
  double stable_error;
  error_comp() : x(NAN), stable_error(NAN){}
  error_comp(double x, double r_double, double r_FMStable, double r_Pareto) :x(x) {
    stable_error = fabs(r_double-r_FMStable);
  }
};

class result {
public:
  double alpha;
  double beta;
  double x;
  double r_double;
  double r_FMStable;
  double absdiff;
  double reldiff;
  double g_theta2_error;
  double abserr_double;
  bool good_theta2_double;
  result(double alpha, double beta, double x, double r_double, double r_FMStable,
         double absdiff, double reldiff, double g_theta2_error, double abserr_double,
         bool good_theta2_double) :
  alpha(alpha), beta(beta), x(x), r_double(r_double),
  r_FMStable(r_FMStable), absdiff(absdiff), reldiff(reldiff), g_theta2_error(g_theta2_error), abserr_double(abserr_double),
   good_theta2_double(good_theta2_double){}
};

bool comp_rel(const result lhs, const result rhs) {
  return lhs.reldiff > rhs.reldiff;
}
bool comp_abs(const result lhs, const result rhs) {
  return lhs.absdiff > rhs.absdiff;
}

class results {
  static vector<double> probs;
  string type;
  vector<result> abs_data;
  vector<result> rel_data;
  vector<double> rel_all;
  double abs_best = -std::numeric_limits<double>::infinity();
  double abs_worst = std::numeric_limits<double>::infinity();
  double rel_best = -std::numeric_limits<double>::infinity();
  double rel_worst = std::numeric_limits<double>::infinity();
  double abs_diff_sum = 0;
  double rel_diff_sum = 0;
  int count = 0;
  int max_size;
  vector<double> quantile() {
    double eps = 100 * std::numeric_limits<double>::epsilon();
    long np = probs.size();
    for (double prob : probs) {
      if (prob < -eps || prob > 1 + eps)
        throw std::range_error("quantile: 'probs' outside [0,1]");
      prob =std::max(static_cast<double>(0),std::min(static_cast<double>(1),prob));
    }
    vector<double> qs(np);
    vector<double>& x{rel_all};
    long n = x.size();
    if (n > 0 && np > 0) {
      std::sort(x.begin(), x.end());
      for (int j=0; j<np; ++j) {
        double index = (n - 1) * probs.at(j);
        int lo = static_cast<int>(floor(index));
        int hi = static_cast<int>(ceil(index));
        double h = index - lo;
        qs.at(j) = (1-h) * x.at(lo) + h * x.at(hi);
      }
      return qs;
    } else {
      throw std::range_error("quantile: Both x and prob must be of length greater than 0");
    }
  }  //quantile
public:
  results(string type, int max_size) : type(type), max_size(max_size) {
  }
  void add_result(ostream& out, double alpha, double beta, double x, double r_double, double r_FMStable,
                  double g_theta2_error, double abserr_double, bool good_theta2_double) {
    double dblmin = std::numeric_limits<double>::min()/std::numeric_limits<double>::epsilon();
    
    if (std::isnan(r_FMStable) || std::isinf(r_FMStable) || std::isnan(r_double) || std::isinf(r_double))
      out << type << ": " << alpha << " " << beta << " " << x << " " << r_double << " " << r_FMStable << endl;
    else {
      count++;
      double abs_diff = fabs(r_double-r_FMStable);
      abs_diff_sum += abs_diff;
      double rel_diff = abs_diff/(dblmin+max(fabs(r_double), fabs(r_FMStable)));
      rel_diff_sum += rel_diff;
      rel_all.push_back(rel_diff);
      if (abs_diff > abs_best)
      {
        abs_data.push_back(result(alpha, beta, x, r_double, r_FMStable,
                                  abs_diff, rel_diff, g_theta2_error, abserr_double, good_theta2_double));
        sort(abs_data.begin(), abs_data.end(), comp_abs);
        if (abs_data.size() > max_size) abs_data.pop_back();
        abs_best = abs_data.back().absdiff;
        abs_worst = abs_data.front().absdiff;
      }
      if (rel_diff > rel_best) {
        rel_data.push_back(result(alpha, beta, x, r_double, r_FMStable,
                                  abs_diff, rel_diff, g_theta2_error, abserr_double, good_theta2_double));
        sort(rel_data.begin(),rel_data.end(),comp_rel);
        if (rel_data.size() > max_size) {
          rel_data.pop_back();
        }
        rel_best = rel_data.back().reldiff;
        rel_worst = rel_data.front().reldiff;
      }
    }
  };
  friend ostream& operator<<(ostream& os, results r);
};

vector<double> results::probs = {.01, .05, .25, .5, .75, .95, .99, .999, .9999, 1};

ostream& operator<<(ostream& os, results r) {
  vector<double> qs=r.quantile();
  os << "Table of worst " << r.abs_data.size() << " out of " << r.count << " absolute differences for " << r.type << endl << endl;
  os  << setw(13) << right << "alpha"
  << setw(13) << right << "beta"
  << setw(13) << right << "x"
  << setw(30) << right << r.type+"_double"
  << setw(30) << right << r.type+"_FMStable"
  << setw(15) << right << "absdiff"
  << setw(15) << right << "reldiff"
  << setw(15) << right << "g_theta2_error"
  << setw(15) << right << "abserr_double"
  << setw(5) << right << "th2d"
  << endl << endl;
  for (vector<result>::iterator presult=r.abs_data.begin(); presult < r.abs_data.end(); presult++) {
    os << setw(13) << setprecision(5) << fixed << presult->alpha
    << setw(13) << setprecision(5) << presult->beta
    << setw(13) << setprecision(5) << scientific << presult->x
    << setw(30) << setprecision(20) << presult->r_double
    << setw(30) << setprecision(20) << presult->r_FMStable
    << setw(15) << setprecision(5) << presult->absdiff
    << setw(15) << setprecision(5) << presult->reldiff
    << setw(15) << setprecision(5) << presult->g_theta2_error
    << setw(15) << setprecision(5) << presult->abserr_double
    << setw(5) << right << presult->good_theta2_double
    << endl;
  }
  
  os << endl << setw(99) << "Average"
  << setw(15) << setprecision(5) << r.abs_diff_sum/r.count
  << setw(15) << setprecision(5) << r.rel_diff_sum/r.count << endl << endl;
  os << setw(99) << "Quantile" << endl << endl;
  for (int i =0; i<results::probs.size(); ++i)
    os << setw(98) << fixed << setprecision(0) << results::probs.at(i)*100 << "%"
    << setw(15) << " " << setw(15) << setprecision(5) << scientific<< qs.at(i) << endl;
  
  os << "Table of worst " << r.rel_data.size() << " out of " << r.count << " relative differences for " << r.type << endl << endl;
  os  << setw(13) << right << "alpha"
  << setw(13) << right << "beta"
  << setw(13) << right << "x"
  << setw(30) << right << r.type+"_double"
  << setw(30) << right << r.type+"_FMStable"
  << setw(15) << right << "absdiff"
  << setw(15) << right << "reldiff"
  << setw(15) << right << "g_theta2_error"
  << setw(15) << right << "abserr_double"
  << setw(5) << right << "th2d"
  << endl << endl;
  for (vector<result>::iterator presult=r.rel_data.begin(); presult < r.rel_data.end(); presult++) {
    os << setw(13) << setprecision(5) << fixed << presult->alpha
    << setw(13) << setprecision(5) << presult->beta
    << setw(13) << setprecision(5) << scientific << presult->x
    << setw(30) << setprecision(20) << presult->r_double
    << setw(30) << setprecision(20) << presult->r_FMStable
    << setw(15) << setprecision(5) << presult->absdiff
    << setw(15) << setprecision(5) << presult->reldiff
    << setw(15) << setprecision(5) << presult->g_theta2_error
    << setw(15) << setprecision(5) << presult->abserr_double
    << setw(5) << right << presult->good_theta2_double
    << endl;
  }
  
  os << endl << setw(99) << "Average"
  << setw(15) << setprecision(5) << r.abs_diff_sum/r.count
  << setw(15) << setprecision(5) << r.rel_diff_sum/r.count << endl << endl;
  os << setw(99) << "Quantile" << endl << endl;
  for (int i =0; i<results::probs.size(); ++i)
    os << setw(98) << fixed << setprecision(0) << results::probs.at(i)*100 << "%"
    << setw(15) << " " << setw(15) << setprecision(5) << scientific<< qs.at(i) << endl;
  
  
  return os;
}

static void show_usage (string name){
 boost::filesystem::path p(name);
 cerr << "Usage: " << p.filename().string() << endl;
}

int main(int argc, char *argv[]) {
  
  // Check the number of parameters
  if (argc != 1) {
    // Tell the user how to run the program
    show_usage(string(argv[0]));
    return 1;
  }
  string out_dir = string(OUT_DIR);
  if (!boost::filesystem::is_directory(out_dir))
    boost::filesystem::create_directory(out_dir);

  string out_file = out_dir + "/FMStable_test.out";
  cout << "Writing output to " + out_file << endl;
  ofstream out(out_file);
  auto_cpu_timer timer(out);
  
  StandardStableDistribution<double>::initialize();
  
  int n_gauss =10;
  mpreal::set_default_prec(96);
  Kronrod<mpreal> g_k_big(n_gauss);  //high precision nodes and weights
  // Parameters for calculation
  int log_flag=0;
  
  bool noext = true, lower_tail=true;
  int subdivisions = 1000;
  int verbose = 0;
  
  double epsabs = 0.;
  double epsrel = 64*std::numeric_limits<double>::epsilon();

  IntegrationController<double> ctl_dbl(noext, g_k_big, epsabs, epsrel, subdivisions, 0);
  Controllers<double> ctls(ctl_dbl, ctl_dbl);
  results r_cdf("cdf",100);
  results r_pdf("pdf",100);
  
  vector<double> alphas = {.01, .1, .2, .3, .4, .5, .6, .7, .8, .9, .99, 1.,
    1.01, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 1.99, 2.0};
  vector<double> xs;
  xs.push_back(0.);
  for (int i=1; i<=100; i++) {
    xs.push_back(static_cast<double>(i)/10);
    xs.push_back(static_cast<double>(-i)/10);
  }
  for (int i=2; i<100; i++) {
    xs.push_back(pow(10.,i));
    xs.push_back(-pow(10.,i));
  }
  ;
  xs.push_back(std::pow(10.,120));
  xs.push_back(-std::pow(10.,120));
  xs.push_back(std::pow(10.,150));
  xs.push_back(-std::pow(10.,150));
/* Eliminated because FMStable hangs
  xs.push_back(std::pow(10.,200));
  xs.push_back(-std::pow(10.,200));
  xs.push_back(std::pow(10.,300));
  xs.push_back(-std::pow(10.,300));
*/

  sort(xs.begin(),xs.end());

  bool pass = true;
  for (double alpha : alphas) {
    out << "alpha = " << fixed << alpha << endl;

    // FMStable uses the C parameterization for alpha < .5
    double angle = StandardStableDistribution<double>::pi2*alpha;
    double location = (alpha<.5) ? -tan(angle) : 0;
    double logscale = (alpha<.5) ? -log(cos(angle))/alpha : 0;

    for (double x : xs) {
      double cdf_FMStable, log_cdf_FMStable, compl_cdf_FMStable, log_compl_cdf_FMStable;
      double pdf_FMStable, log_pdf_FMStable;
      tailsMSS(1, &x, &pdf_FMStable, &log_pdf_FMStable, &cdf_FMStable,
                    &log_cdf_FMStable, &compl_cdf_FMStable, &log_compl_cdf_FMStable,
                    alpha, 1-alpha, 2-alpha,
                    location, logscale);
      StandardStableDistribution<double> std_stable_dist_double(alpha, 1, ctls, verbose);
      double cdf_double = std_stable_dist_double.cdf(x, lower_tail, log_flag);
      bool pass1 = fabs(cdf_double - cdf_FMStable) < 1e-13;
      bool good_theta2_double = std_stable_dist_double.good_theta2;
      double g_theta2_error = std_stable_dist_double.g_theta2_error;
      double cdf_abserr = std_stable_dist_double.abserr;
      r_cdf.add_result(out, alpha, 1, x, cdf_double, cdf_FMStable, g_theta2_error,
                       cdf_abserr, good_theta2_double);
      
      double pdf_double = std_stable_dist_double.pdf(x, log_flag);
      pass1 = pass1 && fabs(pdf_double - pdf_FMStable) < 1e-12;
      double pdf_abserr = std_stable_dist_double.abserr;
      r_pdf.add_result(out, alpha, 1, x, pdf_double, pdf_FMStable, g_theta2_error,
                       pdf_abserr, good_theta2_double);
      if (!pass1)
        out << "FAIL: alpha = " << alpha << ", x = " << scientific << x << endl;
      pass = pass && pass1;
    }
  }
  
  out << endl << endl;
  out << r_cdf << endl;
  out << r_pdf << endl;
  out << endl << (pass?"Test Passed":"Test Failed") << endl;
  
  return !pass;
}




