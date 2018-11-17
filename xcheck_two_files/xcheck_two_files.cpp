/// \file xcheck_two_files.cpp
/// Compare values for cdf, pdf and ddx_pdf of standard stable distribution
/// from two files
/// \author Joseph Dunn
/// \copyright 2018 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::istream;
using std::ostream;
#include "stable_config.h"
#define MPREAL
#define MPFR_FLOAT
#define CPP_BIN_FLOAT
#include "stable_distribution.h"
#include <iomanip>
using std::setw;
using std::setprecision;
using std::right;
using std::fixed;
using std::scientific;
#include <sstream>
using std::stringstream;
using std::istringstream;
#include <limits>
#include <fstream>
using std::ofstream;
#include <boost/filesystem.hpp>
using boost::filesystem::path;
using boost::filesystem::is_regular_file;
using boost::filesystem::is_directory;
#include <boost/filesystem/fstream.hpp>
using boost::filesystem::ifstream;
#include <boost/timer/timer.hpp>
using boost::timer::auto_cpu_timer;

using mpfr::mpreal;

#include <string>
using std::string;
using std::stod;
using std::getline;
using std::max;

using namespace stable_distribution;

struct record {
  double alpha;
  double beta;
  double x;
  mpreal cdf;
  mpreal pdf;
  mpreal ddx_pdf;
  mpreal cdf_abserr;
  mpreal pdf_abserr;
  mpreal ddx_pdf_abserr;
  bool good_theta2;
  mpreal g_dd_theta2;
  friend bool operator> (const record& lhs, const record& rhs) {
    if (lhs.alpha > rhs.alpha) return true;
    else if (lhs.alpha < rhs.alpha) return false;
    else if (lhs.beta > rhs.beta) return true;
    else if (lhs.beta < rhs.beta) return false;
    else if (lhs.x > rhs.x) return true;
    else return false;
  }
  friend istream& operator>> (istream& is, record& r) {
    string buf;
    getline(is, buf);
    if (!is) return is;
    string x_str;
    istringstream ss(buf);
    is >> r.alpha >> r.beta >> x_str >> r.cdf >> r.pdf >> r.ddx_pdf
       >> r.cdf_abserr >> r.pdf_abserr >> r.ddx_pdf_abserr
       >> r.good_theta2 >> r.g_dd_theta2;
    if (x_str == "-inf")
      r.x = -std::numeric_limits<double>::infinity();
    else if (x_str == "inf")
      r.x = std::numeric_limits<double>::infinity();
    else {
      istringstream x_ss(x_str);
      x_ss >> r.x;
      if (!x_ss) {
        is.setstate(std::ios_base::failbit);
      }
    }
    return is;
  }
};

class result {
public:
  double alpha;
  double beta;
  double x;
  mpreal r_file1;
  mpreal r_file2;
  mpreal absdiff;
  mpreal reldiff;
  mpreal abserr_file1;
  mpreal abserr_file2;
  bool good_theta2_file1;
  bool good_theta2_file2;
  result(double alpha, double beta, double x, mpreal r_file1, mpreal r_file2,
         mpreal absdiff, mpreal reldiff, mpreal abserr_file1, mpreal abserr_file2,
         bool good_theta2_file1, bool good_theta2_file2) :
  alpha(alpha), beta(beta), x(x), r_file1(r_file1),
  r_file2(r_file2), absdiff(absdiff), reldiff(reldiff), abserr_file1(abserr_file1),
  abserr_file2(abserr_file2), good_theta2_file1(good_theta2_file1), good_theta2_file2(good_theta2_file2){}
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
  mpreal abs_best = -std::numeric_limits<mpreal>::infinity();
  mpreal abs_worst = std::numeric_limits<mpreal>::infinity();
  mpreal rel_best = -std::numeric_limits<mpreal>::infinity();
  mpreal rel_worst = std::numeric_limits<mpreal>::infinity();
  mpreal abs_diff_sum = 0;
  mpreal rel_diff_sum = 0;
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
  void add_result(double alpha, double beta, double x, mpreal r_file1, mpreal r_file2,
                  mpreal abserr_file1, mpreal abserr_file2, bool good_theta2_file1, bool good_theta2_file2) {
    
    if (boost::math::isnan(r_file2) || boost::math::isinf(r_file2)
        || boost::math::isnan(r_file1) || boost::math::isinf(r_file1))
      cout << type << ": " << alpha << " " << beta << " " << x << " " << r_file1 << " " << r_file2 << endl;
    else {
      count++;
      mpreal abs_diff = fabs(r_file1-r_file2);
      abs_diff_sum += abs_diff;
      mpreal rel_diff = (r_file1 !=0 || r_file2 != 0)
                         ? abs_diff/(max(fabs(r_file1), fabs(r_file2)))
                         :0;
      rel_diff_sum += rel_diff;
      rel_all.push_back(rel_diff);
      if (abs_diff > abs_best)
      {
        abs_data.push_back(result(alpha, beta, x, r_file1, r_file2,
                                  abs_diff, rel_diff, abserr_file1, abserr_file2, good_theta2_file1, good_theta2_file2));
        sort(abs_data.begin(), abs_data.end(), comp_abs);
        if (abs_data.size() > max_size) abs_data.pop_back();
        abs_best = abs_data.back().absdiff;
        abs_worst = abs_data.front().absdiff;
      }
      if (rel_diff > rel_best) {
        rel_data.push_back(result(alpha, beta, x, r_file1, r_file2,
                                  abs_diff, rel_diff, abserr_file1, abserr_file2, good_theta2_file1, good_theta2_file2));
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
  << setw(30) << right << r.type+"_file1"
  << setw(30) << right << r.type+"_file2"
  << setw(15) << right << "absdiff"
  << setw(15) << right << "reldiff"
  << setw(15) << right << "abserr_file1"
  << setw(15) << right << "abserr_file2"
  << setw(5) << right << "th2_1"
  << setw(5) << right << "th2_2"
  << endl << endl;
  for (vector<result>::iterator presult=r.abs_data.begin(); presult < r.abs_data.end(); presult++) {
    os << setw(13) << setprecision(5) << fixed << presult->alpha
    << setw(13) << setprecision(5) << presult->beta
    << setw(13) << setprecision(5) << scientific << presult->x
    << setw(30) << setprecision(20) << presult->r_file1
    << setw(30) << setprecision(20) << presult->r_file2
    << setw(15) << setprecision(5) << presult->absdiff
    << setw(15) << setprecision(5) << presult->reldiff
    << setw(15) << setprecision(5) << presult->abserr_file1
    << setw(15) << setprecision(5) << presult->abserr_file2
    << setw(5) << right << presult->good_theta2_file1
    << setw(5) << right << presult->good_theta2_file2 << endl;
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
  << setw(30) << right << r.type+"_file1"
  << setw(30) << right << r.type+"_file2"
  << setw(15) << right << "absdiff"
  << setw(15) << right << "reldiff"
  << setw(15) << right << "abserr_file1"
  << setw(15) << right << "abserr_file2"
  << setw(5) << right << "th2_1"
  << setw(5) << right << "th2_2"
  << endl << endl;
  for (vector<result>::iterator presult=r.rel_data.begin(); presult < r.rel_data.end(); presult++) {
    os << setw(13) << setprecision(5) << fixed << presult->alpha
    << setw(13) << setprecision(5) << presult->beta
    << setw(13) << setprecision(5) << scientific << presult->x
    << setw(30) << setprecision(20) << presult->r_file1
    << setw(30) << setprecision(20) << presult->r_file2
    << setw(15) << setprecision(5) << presult->absdiff
    << setw(15) << setprecision(5) << presult->reldiff
    << setw(15) << setprecision(5) << presult->abserr_file1
    << setw(15) << setprecision(5) << presult->abserr_file2
    << setw(5) << right << presult->good_theta2_file1
    << setw(5) << right << presult->good_theta2_file2 << endl;
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
  path p(name);
  cerr << "Usage: " << p.filename().string() << "label file1_path file2_path" << endl;
}

int main(int argc, char *argv[]) {
  
  // Check the number of parameters
  if (argc != 4) {
    // Tell the user how to run the program
    show_usage(string(argv[0]));
    return 1;
  }
  
  path path_file1((string(argv[2])));
  if (!is_regular_file(path_file1)) {
    cerr << "Bad file path for file1: " << argv[2] << endl;
    show_usage(string(argv[0]));
    return 1;
  }
  cout << "Opening " << path_file1 << " as file1" << endl;
  ifstream is1(path_file1);
  
  path path_file2((string(argv[3])));
  if (!is_regular_file(path_file2)) {
    cerr << "Bad file path for file2: " << argv[3] << endl;
    show_usage(string(argv[0]));
    return 1;
  }
  cout << "Opening " << path_file2 << " as file2" << endl;
  ifstream is2(path_file2);
  
  mpreal::set_default_prec(96);
  
  string out_dir = string(OUT_DIR);
  
  if (!boost::filesystem::is_directory(out_dir))
    boost::filesystem::create_directory(out_dir);
  
  string out_label{argv[1]};

  string out_file = out_dir + "/xcheck_" + out_label + ".out";
  cout << "Writing output to " + out_file << endl;
  ofstream out(out_file);
  out << "Opening " << path_file1 << " as file1" << endl;
  out << "Opening " << path_file2 << " as file2" << endl << endl;
  auto_cpu_timer timer(out);
  
  results r_cdf("cdf",100);
  results r_pdf("pdf",100);
  results r_ddx_pdf("ddx_pdf",100);
  
  double alpha_old{-999}, beta_old{-999};
  record rec1, rec2;
  bool need1{true}, need2{true};
  
  while (is1 && is2) {
    if (need1) {
      if (!(is1 >> rec1)) break;
      need1 = false;
    }
    if (need2) {
      if (!(is2 >> rec2)) break;
      need2 = false;
    }
    if (rec2 > rec1) {
      need1 = true;
      continue;
    } else if (rec1 > rec2) {
      need2 = true;
      continue;
    }
    
    // rec1 and rec2 now have the same alpha, beta and x;
    if (rec1.alpha != alpha_old || rec1.beta != beta_old) {
      cout << "alpha = " << fixed << setprecision(2) << rec1.alpha
           << ", beta = " << fixed << setprecision(2) << rec1.beta << endl;
      alpha_old = rec1.alpha;
      beta_old = rec1.beta;
    }

    r_cdf.add_result(rec1.alpha, rec1.beta, rec1.x,
                     rec1.cdf, rec2.cdf,
                     rec1.cdf_abserr, rec2.cdf_abserr,
                     rec1.good_theta2, rec2.good_theta2);

    r_pdf.add_result(rec1.alpha, rec1.beta, rec1.x,
                     rec1.pdf, rec2.pdf,
                     rec1.pdf_abserr, rec2.pdf_abserr,
                     rec1.good_theta2, rec2.good_theta2);
 
    r_ddx_pdf.add_result(rec1.alpha, rec1.beta, rec1.x,
                         rec1.ddx_pdf, rec2.ddx_pdf,
                         rec1.ddx_pdf_abserr, rec2.ddx_pdf_abserr,
                         rec1.good_theta2, rec2.good_theta2);
    need1=true;
    need2=true;
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




