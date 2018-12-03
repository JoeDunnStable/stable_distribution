/// \file xcheck_to_Zolotarev.cpp
/// Unit tests for routines in zolotarev.h
/// \author Joseph Dunn
/// \copyright 2017, 2018 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <iomanip>
#include <sstream>
#include <fstream>
#define MPREAL
#include "stable_distribution.h"
#include "stable_config.h"
#include <boost/math/constants/info.hpp>
#include <boost/filesystem.hpp>
#include <vector>
#include <algorithm>

#define LIBRARY
#include "zolotarev.h"
#undef LIBRARY

using std::setw;
using std::fixed;
using std::scientific;
using std::right;
using std::istringstream;
using std::ostream;
using std::ofstream;
using std::setprecision;
using std::string;
using std::vector;

using namespace stable_distribution;

template<typename myFloat>
class error_comp {
public:
  myFloat x;
  myFloat stable_error;
  error_comp() : x(NAN), stable_error(NAN){}
  error_comp(myFloat x, myFloat r_myFloat, myFloat, myFloat r_Zolotarev) :x(x) {
    stable_error = fabs(r_myFloat-r_Zolotarev);
  }
};

template<typename myFloat>
class Result {
public:
  myFloat alpha;
  myFloat beta;
  myFloat x;
  myFloat r_myFloat;
  myFloat r_Zolotarev;
  myFloat absdiff;
  myFloat reldiff;
  myFloat g_theta2_error;
  myFloat abserr_myFloat;
  myFloat error_Zolotarev;
  bool good_theta2_myFloat;
  ResultType result_type;
  Result(myFloat alpha, myFloat beta, myFloat x, myFloat r_myFloat, myFloat r_Zolotarev,
         myFloat absdiff, myFloat reldiff, myFloat g_theta2_error, myFloat abserr_myFloat,
         myFloat error_Zolotarev, bool good_theta2_myFloat,
         ResultType result_type) :
  alpha(alpha), beta(beta), x(x), r_myFloat(r_myFloat),
  r_Zolotarev(r_Zolotarev), absdiff(absdiff), reldiff(reldiff), g_theta2_error(g_theta2_error), abserr_myFloat(abserr_myFloat),
  error_Zolotarev(error_Zolotarev), good_theta2_myFloat(good_theta2_myFloat),
  result_type(result_type){}
};

template<typename myFloat>
bool comp_rel(const Result<myFloat> lhs, const Result<myFloat> rhs) {
  return lhs.reldiff > rhs.reldiff;
}
template<typename myFloat>
bool comp_abs(const Result<myFloat> lhs, const Result<myFloat> rhs) {
  return lhs.absdiff > rhs.absdiff;
}
template<typename myFloat>
bool comp_ab(const Result<myFloat> lhs, const Result<myFloat> rhs) {
  return lhs.alpha > rhs.alpha || (lhs.alpha == rhs.alpha && lhs.beta > rhs.beta);
}

template<typename myFloat>
class Results {
  static vector<myFloat> probs;
  string type;
  vector<Result<myFloat> > abs_data;
  vector<Result<myFloat> > rel_data;
  vector<Result<myFloat> > raw_data;
  vector<myFloat> rel_all;
  myFloat abs_best = -std::numeric_limits<myFloat>::infinity();
  myFloat abs_worst = std::numeric_limits<myFloat>::infinity();
  myFloat rel_best = -std::numeric_limits<myFloat>::infinity();
  myFloat rel_worst = std::numeric_limits<myFloat>::infinity();
  myFloat abs_diff_sum = 0;
  myFloat rel_diff_sum = 0;
  int count = 0;
  int max_size;
  
  vector<myFloat> quantile() {
    myFloat eps = 100 * std::numeric_limits<myFloat>::epsilon();
    long np = probs.size();
    for (myFloat prob : probs) {
      if (prob < -eps || prob > 1 + eps)
        throw std::range_error("quantile: 'probs' outside [0,1]");
      prob =std::max(static_cast<myFloat>(0),std::min(static_cast<myFloat>(1),prob));
    }
    vector<myFloat> qs(np);
    vector<myFloat>& x{rel_all};
    long n = x.size();
    if (n > 0 && np > 0) {
      std::sort(x.begin(), x.end());
      for (int j=0; j<np; ++j) {
        myFloat index = (n - 1) * probs.at(j);
        int lo = static_cast<int>(floor(index));
        int hi = static_cast<int>(ceil(index));
        myFloat h = index - lo;
        qs.at(j) = (1-h) * x.at(lo) + h * x.at(hi);
      }
      return qs;
    } else {
      throw std::range_error("quantile: Both x and prob must be of length greater than 0");
    }
  }  //quantile
public:
  Results(string type, int max_size) : type(type), max_size(max_size) {
  }
  void add_result(ostream& out, myFloat alpha, myFloat beta, myFloat x, myFloat r_myFloat,
                  myFloat r_Zolotarev, myFloat g_theta2_error, myFloat abserr_myFloat,
                  myFloat error_Zolotarev, bool good_theta2_myFloat,
                  ResultType result_type) {
    myFloat dblmin = std::numeric_limits<myFloat>::min()/   std::numeric_limits<myFloat>::epsilon();
    
    if (boost::math::isnan(r_Zolotarev) || boost::math::isinf(r_Zolotarev) || boost::math::isnan(r_myFloat) || boost::math::isinf(r_myFloat))
      out << type << ": " << alpha << " " << beta << " " << x << " " << r_myFloat << " " << r_Zolotarev << endl;
    else {
      count++;
      myFloat abs_diff = fabs(r_myFloat-r_Zolotarev);
      abs_diff_sum += abs_diff;
      myFloat rel_diff = abs_diff/(dblmin+std::max(fabs(r_myFloat), fabs(r_Zolotarev)));
      rel_diff_sum += rel_diff;
      rel_all.push_back(rel_diff);
      raw_data.push_back(Result<myFloat>(alpha, beta, x, r_myFloat, r_Zolotarev,
                                         abs_diff, rel_diff, g_theta2_error, abserr_myFloat,
                                         error_Zolotarev, good_theta2_myFloat,
                                         result_type));
      if (abs_data.size() < max_size || abs_diff >= abs_best)
      {
        abs_data.push_back(Result<myFloat>(alpha, beta, x, r_myFloat, r_Zolotarev,
                                  abs_diff, rel_diff, g_theta2_error, abserr_myFloat,
                                  error_Zolotarev, good_theta2_myFloat,
                                  result_type));
        sort(abs_data.begin(), abs_data.end(), comp_abs<myFloat>);
        if (abs_data.size() > max_size) abs_data.pop_back();
        abs_best = abs_data.back().absdiff;
        abs_worst = abs_data.front().absdiff;
      }
      if (rel_data.size() < max_size || rel_diff >= rel_best) {
        rel_data.push_back(Result<myFloat>(alpha, beta, x, r_myFloat, r_Zolotarev,
                                  abs_diff, rel_diff, g_theta2_error, abserr_myFloat,
                                  error_Zolotarev, good_theta2_myFloat,
                                  result_type));
        sort(rel_data.begin(),rel_data.end(),comp_rel<myFloat>);
        if (rel_data.size() > max_size) {
          rel_data.pop_back();
        }
        rel_best = rel_data.back().reldiff;
        rel_worst = rel_data.front().reldiff;
      }
    }
  };
  void print(ostream& os);
};

template<typename myFloat>
vector<myFloat> Results<myFloat>::probs = {.01, .05, .25, .5, .75, .95, .99, .999, .9999, 1};

template<typename myFloat>
void Results<myFloat>::print(ostream& os){
  vector<myFloat> qs = quantile();
  os << "Table of all " << count << " " << type << " raw results" << endl << endl;
  os  << setw(13) << right << "alpha"
  << setw(13) << right << "beta"
  << setw(13) << right << "x"
  << setw(30) << right << type+"_myFloat"
  << setw(30) << right << type+"_Zolotarev"
  << setw(15) << right << "absdiff"
  << setw(15) << right << "reldiff"
  << setw(15) << right << "g_theta2_error"
  << setw(15) << right << "abserr"
  << setw(15) << right << "err_Zolotarev"
  << setw(5) << right << "th2d"
  << setw(5) << right << "type"
  << endl << endl;
  for (auto result : raw_data) {
    os << setw(13) << setprecision(5) << fixed << result.alpha
    << setw(13) << setprecision(5) << result.beta
    << setw(13) << setprecision(5) << scientific << result.x
    << setw(30) << setprecision(20) << result.r_myFloat
    << setw(30) << setprecision(20) << result.r_Zolotarev
    << setw(15) << setprecision(5) << result.absdiff
    << setw(15) << setprecision(5) << result.reldiff
    << setw(15) << setprecision(5) << result.g_theta2_error
    << setw(15) << setprecision(5) << result.abserr_myFloat
    << setw(15) << setprecision(5) << result.error_Zolotarev
    << setw(5) << right << result.good_theta2_myFloat
    << setw(5) << right << (result.result_type==asymptotic?"A":"C")
    << endl;
  }
  
  os << endl << setw(99) << "Average"
  << setw(15) << setprecision(5) << abs_diff_sum/count
  << setw(15) << setprecision(5) << rel_diff_sum/count << endl << endl;
  os << setw(99) << "Quantile" << endl << endl;
  for (int i =0; i<Results<myFloat>::probs.size(); ++i)
    os << setw(98) << fixed << setprecision(2) << Results<myFloat>::probs.at(i)*100 << "%"
    << setw(15) << " " << setw(15) << setprecision(5) << scientific<< qs.at(i) << endl;
  
  os << "Table of worst " << abs_data.size() << " out of " << count << " absolute differences for " << type << endl << endl;
  os  << setw(13) << right << "alpha"
  << setw(13) << right << "beta"
  << setw(13) << right << "x"
  << setw(30) << right << type+"_myFloat"
  << setw(30) << right << type+"_Zolotarev"
  << setw(15) << right << "absdiff"
  << setw(15) << right << "reldiff"
  << setw(15) << right << "g_theta2_error"
  << setw(15) << right << "abserr"
  << setw(15) << right << "err_Zolotarev"
  << setw(5) << right << "th2d"
  << setw(5) << right << "type"
  << endl << endl;
  for (auto result : abs_data) {
    os << setw(13) << setprecision(5) << fixed << result.alpha
    << setw(13) << setprecision(5) << result.beta
    << setw(13) << setprecision(5) << scientific << result.x
    << setw(30) << setprecision(20) << result.r_myFloat
    << setw(30) << setprecision(20) << result.r_Zolotarev
    << setw(15) << setprecision(5) << result.absdiff
    << setw(15) << setprecision(5) << result.reldiff
    << setw(15) << setprecision(5) << result.g_theta2_error
    << setw(15) << setprecision(5) << result.abserr_myFloat
    << setw(15) << setprecision(5) << result.error_Zolotarev
    << setw(5) << right << result.good_theta2_myFloat
    << setw(5) << right << (result.result_type==asymptotic?"A":"C")
    << endl;
  }
  
  os << endl << setw(99) << "Average"
  << setw(15) << setprecision(5) << abs_diff_sum/count
  << setw(15) << setprecision(5) << rel_diff_sum/count << endl << endl;
  os << setw(99) << "Quantile" << endl << endl;
    for (int i =0; i<Results<myFloat>::probs.size(); ++i)
    os << setw(98) << fixed << setprecision(2) << Results<myFloat>::probs.at(i)*100 << "%"
    << setw(15) << " " << setw(15) << setprecision(5) << scientific<< qs.at(i) << endl;
  
  os << "Table of worst " << rel_data.size() << " out of " << count << " relative differences for " << type << endl << endl;
  os  << setw(13) << right << "alpha"
  << setw(13) << right << "beta"
  << setw(13) << right << "x"
  << setw(30) << right << type+"_myFloat"
  << setw(30) << right << type+"_Zolotarev"
  << setw(15) << right << "absdiff"
  << setw(15) << right << "reldiff"
  << setw(15) << right << "g_theta2_error"
  << setw(15) << right << "abserr"
  << setw(15) << right << "err_Zolotarev"
  << setw(5) << right << "th2d"
  << setw(5) << right << "type"
  << endl << endl;
  for (auto result : rel_data) {
    os << setw(13) << setprecision(5) << fixed << result.alpha
    << setw(13) << setprecision(5) << result.beta
    << setw(13) << setprecision(5) << scientific << result.x
    << setw(30) << setprecision(20) << result.r_myFloat
    << setw(30) << setprecision(20) << result.r_Zolotarev
    << setw(15) << setprecision(5) << result.absdiff
    << setw(15) << setprecision(5) << result.reldiff
    << setw(15) << setprecision(5) << result.g_theta2_error
    << setw(15) << setprecision(5) << result.abserr_myFloat
    << setw(15) << setprecision(5) << result.error_Zolotarev
    << setw(5) << right << result.good_theta2_myFloat
    << setw(5) << right << (result.result_type==asymptotic?"A":"C")
    << endl;
  }
  
  os << endl << setw(99) << "Average"
  << setw(15) << setprecision(5) << abs_diff_sum/count
  << setw(15) << setprecision(5) << rel_diff_sum/count << endl << endl;
  os << setw(99) << "Quantile" << endl << endl;
  for (int i =0; i<Results<myFloat>::probs.size(); ++i)
    os << setw(98) << fixed << setprecision(2) << Results<myFloat>::probs.at(i)*100 << "%"
    << setw(15) << " " << setw(15) << setprecision(5) << scientific << qs.at(i) << endl;
  
}

static void show_usage (string name){
  boost::filesystem::path p(name);
  cerr << "Usage: " << p.filename().string() << "[verbose = 0]" << endl;
}

template<typename myFloat>
void zolotarev_comparison(string float_type, Controllers<myFloat> ctls, int verbose) {
  
  string out_dir = string(OUT_DIR);
  if (!boost::filesystem::is_directory(out_dir))
    boost::filesystem::create_directory(out_dir);
  string outfile = out_dir + "/xcheck_to_Zolotarev_" + float_type + ".out";
  cout << "Writing output to " << outfile << endl;
  ofstream out(outfile);
  out << stable_config << endl;  

 vector<double> alphas {.01, .1, .5, .9, .99, 1,
    1.01, 1.1, 1.5, 1.9, 1.99};
  
  vector<double> betas {-1, -.5, 0, .5, 1.0};
  
  vector<double> xs {-1e300, -1e150, -1e100, -1e50, -1e20, -1e12, -1e9, -1e6, -1000, -100, -10, -5, -1, -.5, -.1, 0, .1, .5, 1, 5, 10, 100, 1000,
    1e6, 1e9, 1e12, 1e20, 1e50, 1e100, 1e150, 1e300};
  
  Results<myFloat> r_cdf("cdf",100);
  Results<myFloat> r_pdf("pdf",100);
  Results<myFloat> r_ddx_pdf("ddx_pdf",100);
  
  for (double alpha : alphas) {
    for (double beta : betas) {
      for (double x : xs) {
        Zolotarev<myFloat> zol(myFloat(alpha), myFloat(beta), &ctls.controller, verbose) ;
        StandardStableDistribution<myFloat> std_dist(myFloat(alpha), myFloat(beta), ctls, 0);
        int lower_tail= 1;
        myFloat r_Zolotarev, error_Zolotarev;
        r_Zolotarev = zol.cdf(myFloat(x), lower_tail, S1);
        error_Zolotarev = zol.error;
        ResultType result_type = zol.result_type;
        int log_flag=0;
        myFloat r_myFloat, r_abserr;
        r_myFloat = std_dist.cdf(myFloat(x), lower_tail, log_flag, S1);
        r_abserr = std_dist.abserr;
        bool good_theta2 = std_dist.good_theta2;
        myFloat g_theta2_error = std_dist.g_theta2_error;
        r_cdf.add_result(out, myFloat(alpha), myFloat(beta), myFloat(x),
                         r_myFloat, r_Zolotarev, g_theta2_error, r_abserr,
                         error_Zolotarev, good_theta2, result_type);
        r_Zolotarev = zol.pdf(myFloat(x), S1);
        error_Zolotarev = zol.error;
        result_type = zol.result_type;
        r_myFloat = std_dist.pdf(myFloat(x), log_flag, S1);
        r_abserr = std_dist.abserr;
        r_pdf.add_result(out, myFloat(alpha), myFloat(beta), myFloat(x),
                         r_myFloat, r_Zolotarev, g_theta2_error, r_abserr,
                         error_Zolotarev, good_theta2, result_type);
        r_Zolotarev = zol.ddx_pdf(myFloat(x), S1);
        error_Zolotarev = zol.error;
        result_type = zol.result_type;
        r_myFloat = std_dist.ddx_pdf(myFloat(x), S1);
        r_abserr = std_dist.abserr;
        r_ddx_pdf.add_result(out, myFloat(alpha), myFloat(beta), myFloat(x),
                         r_myFloat, r_Zolotarev, g_theta2_error, r_abserr,
                         error_Zolotarev, good_theta2, result_type);
      }
    }
  }
  out << endl << endl;
  r_cdf.print(out);
  out << endl;
  r_pdf.print(out);
  out << endl;
  r_ddx_pdf.print(out);
  out << endl;
}

int main(int argc, const char * argv[]) {

  mpfr::mpreal::set_default_prec(128);
  int n_gauss = 10;
  Kronrod<mpreal> g_k_big(n_gauss);
  mpfr::mpreal::set_default_prec(96);
  
  // Check the number of parameters
  if (argc > 2) {
    // Tell the user how to run the program
    show_usage(string(argv[0]));
    return 1;
  }
  
  int verbose = 0;
  if (argc > 0) {
    istringstream ss1((string(argv[4])));
    ss1 >> verbose;
  }
  
  bool noext = true;   //disable extrapolation
  mpreal epsabs = 0.;
  mpreal epsrel = 64*std::numeric_limits<mpreal>::epsilon();
  int subdivisions = 10000;
  
  IntegrationController<mpreal> controller(noext, g_k_big, epsabs, epsrel, subdivisions,
                                            0);
  IntegrationController<double> cntl_double(noext, g_k_big,
                                            static_cast<double>(epsabs),
                                            static_cast<double>(epsrel), subdivisions, 0);
  
  Controllers<double> ctls_double(cntl_double, cntl_double);
  
  string float_type="double";
  
  zolotarev_comparison<double>(float_type, ctls_double, verbose);
  
  Controllers<mpreal> ctls_mpreal(controller, cntl_double);
  
  float_type = "mpreal";
  
  zolotarev_comparison<mpreal>(float_type, ctls_mpreal, verbose);
  
  
 }
