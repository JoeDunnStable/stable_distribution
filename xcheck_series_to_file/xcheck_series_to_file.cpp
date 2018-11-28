/// \file xcheck_series_to_file.cpp
/// Compare cdf, pdf and ddx_pdf of standard stable distribution to file values
/// \author Joseph Dunn
/// \copyright 2016, 2017, 2018 Joseph Dunn
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
#include <iomanip>
#include <sstream>
using std::stringstream;
#include <limits>
#include <fstream>
#include <boost/filesystem.hpp>

#include <boost/timer/timer.hpp>
using boost::timer::auto_cpu_timer;

using mpfr::mpreal;

#include <string>
using std::string;
using std::istringstream;
using std::setw;
using std::setprecision;
using std::right;
using std::fixed;
using std::scientific;
using std::ifstream;
using std::ofstream;
using std::stod;
using std::getline;
using std::stringstream;
using std::max;

using namespace stable_distribution;

class result {
public:
  double alpha;
  double beta;
  double x;
  double r_mpreal;
  double r_series;
  double absdiff;
  double reldiff;
  result(double alpha, double beta, double x, double r_mpreal, double r_series,
         double absdiff, double reldiff) :
  alpha(alpha), beta(beta), x(x),
  r_mpreal(r_mpreal), r_series(r_series), absdiff(absdiff), reldiff(reldiff) {}
};

bool comp_rel(const result lhs, const result rhs) {
  return lhs.reldiff > rhs.reldiff;
}
bool comp_abs(const result lhs, const result rhs) {
  return lhs.absdiff > rhs.absdiff;
}

class results {
public:
  static vector<mpreal> probs;
  string type;
  vector<result> abs_data;
  vector<result> rel_data;
  vector<mpreal> rel_all;
  double abs_best = -std::numeric_limits<double>::infinity();
  double abs_worst = std::numeric_limits<double>::infinity();
  double rel_best = -std::numeric_limits<double>::infinity();
  double rel_worst = std::numeric_limits<double>::infinity();
  double abs_diff_sum = 0;
  double rel_diff_sum = 0;
  int count = 0;
  int max_size;
  vector<mpreal> quantile() {
    mpreal eps = 100 * std::numeric_limits<mpreal>::epsilon();
    long np = probs.size();
    for (mpreal prob : probs) {
      if (prob < -eps || prob > 1 + eps)
        throw std::range_error("quantile: 'probs' outside [0,1]");
      prob =std::max(static_cast<mpreal>(0),std::min(static_cast<mpreal>(1),prob));
    }
    vector<mpreal> qs(np);
    vector<mpreal>& x{rel_all};
    long n = x.size();
    if (n > 0 && np > 0) {
      std::sort(x.begin(), x.end());
      for (int j=0; j<np; ++j) {
        mpreal index = (n - 1) * probs.at(j);
        int lo = static_cast<int>(floor(index));
        int hi = static_cast<int>(ceil(index));
        mpreal h = index - lo;
        qs.at(j) = (1-h) * x.at(lo) + h * x.at(hi);
      }
      return qs;
    } else {
      throw std::range_error("quantile: Both x and prob must be of length greater than 0");
    }
  }  //quantile
  results(string type, int max_size) : type(type), max_size(max_size) {
  }
  void add_result(double alpha, double beta, double x, double r_mpreal, double r_series_large_x, double r_series_small_x) {
    double dblmin = std::numeric_limits<double>::min()/std::numeric_limits<double>::epsilon();
    
    bool large_is_bad = (boost::math::isnan(r_series_large_x) || boost::math::isinf(r_series_large_x));
    bool small_is_bad = (boost::math::isnan(r_series_small_x) || boost::math::isinf(r_series_small_x));
    if (boost::math::isnan(r_mpreal) || boost::math::isinf(r_mpreal)
        || (large_is_bad && small_is_bad ))
      cout << type << ": " << alpha << " " << beta << " " << x << " " << r_mpreal
           << " " << r_series_large_x << " " << r_series_small_x << endl;
    else {
      count++;
      double abs_diff_large = fabs(r_series_large_x-r_mpreal);
      double abs_diff_small = fabs(r_series_small_x-r_mpreal);
      double r_series;
      double abs_diff;
      if (small_is_bad || abs_diff_large < abs_diff_small) {
        abs_diff = abs_diff_large;
        r_series = r_series_large_x;
      } else {
        abs_diff = abs_diff_small;
        r_series = r_series_small_x;
      }
      abs_diff_sum += abs_diff;
      double rel_diff = abs_diff/(dblmin+max(fabs(r_series), fabs(r_mpreal)));
      rel_diff_sum += rel_diff;
      rel_all.push_back(rel_diff);
      if (abs_diff > abs_best)
      {
        abs_data.push_back(result(alpha, beta, x, r_mpreal, r_series, abs_diff, rel_diff));
        sort(abs_data.begin(), abs_data.end(), comp_abs);
        if (abs_data.size() > max_size) abs_data.pop_back();
        abs_best = abs_data.back().absdiff;
        abs_worst = abs_data.front().absdiff;
      }
      if (rel_diff > rel_best) {
        rel_data.push_back(result(alpha, beta, x, r_mpreal, r_series, abs_diff, rel_diff));
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

vector<mpreal> results::probs = {.01, .05, .25, .5, .75, .95, .99, .999, .9999, 1};

ostream& operator<<(ostream& os, results r) {
  vector<mpreal> qs = r.quantile();
  os << "Table of worst " << r.abs_data.size() << " out of " << r.count << " absolute differences for " << r.type << endl << endl;
  os  << setw(13) << right << "alpha"
      << setw(13) << right << "beta"
      << setw(13) << right << "x"
      << setw(30) << right << r.type+"_series"
      << setw(30) << right << r.type+"_file"
      << setw(15) << right << "absdiff"
      << setw(15) << right << "reldiff"
      << endl << endl;
  for (vector<result>::iterator presult=r.abs_data.begin(); presult < r.abs_data.end(); presult++) {
    os << setw(13) << setprecision(5) << fixed << presult->alpha
       << setw(13) << setprecision(5) << presult->beta
       << setw(13) << setprecision(5) << scientific << presult->x
       << setw(30) << setprecision(20) << presult->r_series
       << setw(30) << setprecision(20) << presult->r_mpreal
       << setw(15) << setprecision(5) << presult->absdiff
       << setw(15) << setprecision(5) << presult->reldiff
       << endl;
  }
  
  os << endl << setw(99) << "Average"
  << setw(15) << setprecision(5) << r.abs_diff_sum/r.count
  << setw(15) << setprecision(5) << r.rel_diff_sum/r.count << endl << endl;
  os << setw(99) << "Quantile" << endl << endl;
  for (int i =0; i<results::probs.size(); ++i)
    os << setw(98) << fixed << setprecision(2) << results::probs.at(i)*100 << "%"
    << setw(15) << " " << setw(15) << setprecision(5) << scientific<< qs.at(i) << endl;
  
  os << "Table of worst " << r.rel_data.size() << " out of " << r.count << " relative differences for " << r.type << endl << endl;
  os  << setw(13) << right << "alpha"
  << setw(13) << right << "beta"
  << setw(13) << right << "x"
  << setw(30) << right << r.type+"_series"
  << setw(30) << right << r.type+"_file"
  << setw(15) << right << "absdiff"
  << setw(15) << right << "reldiff"
  << endl << endl;
  for (vector<result>::iterator presult=r.rel_data.begin(); presult < r.rel_data.end(); presult++) {
    os << setw(13) << setprecision(5) << fixed << presult->alpha
       << setw(13) << setprecision(5) << presult->beta
       << setw(13) << setprecision(5) << scientific << presult->x
       << setw(30) << setprecision(20) << presult->r_series
       << setw(30) << setprecision(20) << presult->r_mpreal
       << setw(15) << setprecision(5) << presult->absdiff
       << setw(15) << setprecision(5) << presult->reldiff
       << endl;
  }
  
  os << endl << setw(99) << "Average"
  << setw(15) << setprecision(5) << r.abs_diff_sum/r.count
  << setw(15) << setprecision(5) << r.rel_diff_sum/r.count << endl << endl;
  os << setw(99)<< "Quantile" << endl << endl;
  for (int i =0; i<results::probs.size(); ++i)
    os << setw(98) << fixed << setprecision(2) << results::probs.at(i)*100 << "%"
    << setw(15) << " " << setw(15) << setprecision(5) << scientific<< qs.at(i) << endl;
  
  return os;
}

static void show_usage (string name){
  boost::filesystem::path p(name);
  cerr << "Usage: " << p.filename().string() << "float_type" << endl;
  cerr << "Where float_type can be double, "
#ifdef MPREAL
  << "mpreal, "
#endif
#ifdef MPFR_FLOAT
  << "mpfr_float, "
#endif
#ifdef CPP_BIN_FLOAT
  << "cpp_bin "
#endif
  << endl;
}

int main(int argc, char *argv[]) {
  
  // Check the number of parameters
  if (argc != 2) {
    // Tell the user how to run the program
    show_usage(string(argv[0]));
    return 1;
  }
  
  stringstream ss1((string(argv[1])));
  string float_type_str; ss1 >> float_type_str;
  if (!(float_type_str == "double"
        || float_type_str == "mpreal"
        || float_type_str == "mpfr_float"
        || float_type_str == "cpp_bin_float")) {
    cerr << "Unknown float type" << endl;
    show_usage(string(argv[0]));
    return 1;
  }

  mpreal::set_default_prec(96);
  
  string out_dir = string(OUT_DIR);
  if (!boost::filesystem::is_directory(out_dir))
    boost::filesystem::create_directory(out_dir);

  string in_file = out_dir + "/stable_" + float_type_str +".out";
  cout << "Reading from " << in_file << endl;
  ifstream in(in_file);
  if (!in) {
    cerr << "Unable to open input file." << endl;
    return 1;
  }
  
  string out_file = out_dir + "/xcheck_series_to_" + float_type_str + ".out";
  cout << "Writing output to " + out_file << endl;
  ofstream out(out_file);
  auto_cpu_timer timer(out);
  out << stable_config << endl;
  
  StandardStableDistribution<double>::initialize();
  int N =10;
  
  // Parameters for calculation
  int log_flag=0;
  
  bool noext = true, lower_tail=true;
  int subdivisions = 10000;
  int verbose = 0;
  
  double dblepsabs = 0.;
  double dblepsrel = 64*std::numeric_limits<double>::epsilon();

  IntegrationController<double> ctl_dbl(noext, N, dblepsabs, dblepsrel, subdivisions, 0);
  Controllers<double> ctls(ctl_dbl, ctl_dbl);
  
  results r_cdf("cdf",100);
  results r_pdf("pdf",100);
  results r_ddx_pdf("ddx_pdf",100);
  
  mpreal p_mpreal, d_mpreal, ddx_d_mpreal;
  mpreal p_abserr_mpreal, d_abserr_mpreal, ddx_d_abserr_mpreal;
  bool good_theta2_mpreal;
  mpreal g_dd_theta2_mpreal;
  string x_str;
  double alpha, beta, x;
  double alpha_old=-1;
  double beta_old=-10;
  
  while (in >> alpha >> beta >> x_str >> p_mpreal >> d_mpreal >> ddx_d_mpreal
         >> p_abserr_mpreal >> d_abserr_mpreal >> ddx_d_abserr_mpreal
         >> good_theta2_mpreal >> g_dd_theta2_mpreal) {
    if (x_str=="inf")
      x = std::numeric_limits<double>::infinity();
    else if (x_str=="-inf")
      x = -std::numeric_limits<double>::infinity();
    else {
      stringstream ss(x_str);
      ss >> x;
    }
    if (alpha != alpha_old || beta != beta_old) {
      cout << "alpha = " << alpha << ", beta = " << beta << endl;
      alpha_old = alpha;
      beta_old = beta;
    }
    double cdf_mpreal = double(p_mpreal);
    double pdf_mpreal = double(d_mpreal);
    double ddx_pdf_mpreal = double(ddx_d_mpreal);
    StandardStableDistribution<double> std_stable_dist_double(alpha, beta, ctls, verbose);
    std_stable_dist_double.cdf(x, lower_tail, log_flag);
    if (!std_stable_dist_double.use_series()) {
      // special case or middle range where integral is best
      continue;
    }
    double series_large_x_cdf = std_stable_dist_double.series_large_x_cdf(x, lower_tail);
    double series_small_x_cdf = std_stable_dist_double.series_small_x_cdf(x, lower_tail);
    r_cdf.add_result(alpha, beta, x, cdf_mpreal, series_large_x_cdf, series_small_x_cdf);

    std_stable_dist_double.pdf(x, log_flag);
    double series_large_x_pdf = std_stable_dist_double.series_large_x_pdf(x);
    double series_small_x_pdf = std_stable_dist_double.series_small_x_pdf(x);
    r_pdf.add_result(alpha, beta, x, pdf_mpreal, series_large_x_pdf, series_small_x_pdf);

    std_stable_dist_double.ddx_pdf(x);
    double series_large_x_ddx_pdf = std_stable_dist_double.series_large_x_ddx_pdf(x);
    double series_small_x_ddx_pdf = std_stable_dist_double.series_small_x_ddx_pdf(x);
    r_ddx_pdf.add_result(alpha, beta, x, ddx_pdf_mpreal, series_large_x_ddx_pdf, series_small_x_ddx_pdf);
  }
  
  out << endl << endl;
  bool pass = true;

  pass = pass && r_cdf.abs_worst < 1e-13 && r_cdf.rel_worst < 1e-10;
  pass = pass && (r_cdf.abs_diff_sum/r_cdf.count < 1e-16) && (r_cdf.rel_diff_sum/r_cdf.count) < 1e-14;
  out << r_cdf << endl;

  pass = pass && r_pdf.rel_worst < 1e-10;
  pass = pass && (r_pdf.rel_diff_sum/r_pdf.count) < 1e-13;
  out << r_pdf << endl;

  pass = pass && r_ddx_pdf.abs_worst < 1e-10 && r_ddx_pdf.rel_worst < 1e-2;
  pass = pass && (r_ddx_pdf.abs_diff_sum/r_ddx_pdf.count < 1e-15) && (r_ddx_pdf.rel_diff_sum/r_ddx_pdf.count) < 1e-7;
  out << r_ddx_pdf << endl;
  out << endl << (pass?"Test Passed":"Test Failed") << endl;
  return !pass;
}




