/// \file dump_mpreal.cpp
/// Dumps high precision cdf, pdf, & ddx_pdf of standard stable distribution
/// \author Joseph Dunn
/// \copyright 2016, 2017 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::getline;

#define MPREAL
#include "stable_distribution.h"
using namespace stable_distribution;

#include <iomanip>
using std::setw;
using std::setprecision;
using std::right;
using std::scientific;
using std::fixed;

#include <string>
using std::string;

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <sstream>
using std::stringstream;

#include <vector>
using std::vector;

#include <array>
using std::array;

#include <deque>
using std::deque;

#include <thread>
using std::thread;
using std::this_thread::yield;

#include <mutex>
using std::mutex;
using std::unique_lock;

#include <condition_variable>
using std::condition_variable;

#include <chrono>
using std::chrono::high_resolution_clock;
using std::chrono::duration;

template<typename myFloat>
class result {
public:
  bool finished;
  myFloat alpha;
  myFloat beta;
  array<double, 407> xs;
  array<myFloat, 407> cdf;
  array<myFloat, 407> pdf;
  array<myFloat, 407> ddx_pdf;
  array<myFloat, 407> cdf_abserr;
  array<myFloat, 407> pdf_abserr;
  array<myFloat, 407> ddx_pdf_abserr;
  array<bool, 407> good_theta2;
  array<myFloat, 407> g_dd_theta2;
  result() : finished(false), alpha(std::numeric_limits<myFloat>::quiet_NaN()),
             beta(std::numeric_limits<myFloat>::quiet_NaN()) {
    xs.fill(std::numeric_limits<double>::quiet_NaN());
    cdf.fill(std::numeric_limits<myFloat>::quiet_NaN());
    pdf.fill(std::numeric_limits<myFloat>::quiet_NaN());
    ddx_pdf.fill(std::numeric_limits<myFloat>::quiet_NaN());
    cdf_abserr.fill(std::numeric_limits<myFloat>::quiet_NaN());
    pdf_abserr.fill(std::numeric_limits<myFloat>::quiet_NaN());
    ddx_pdf_abserr.fill(std::numeric_limits<myFloat>::quiet_NaN());
    good_theta2.fill(false);
    g_dd_theta2.fill(std::numeric_limits<myFloat>::quiet_NaN());
  }
  result(array<double, 407>& xs) : finished(false), alpha(std::numeric_limits<myFloat>::quiet_NaN()),
                                   beta(std::numeric_limits<myFloat>::quiet_NaN()), xs(xs){
    cdf.fill(std::numeric_limits<myFloat>::quiet_NaN());
    pdf.fill(std::numeric_limits<myFloat>::quiet_NaN());
    ddx_pdf.fill(std::numeric_limits<myFloat>::quiet_NaN());
    cdf_abserr.fill(std::numeric_limits<myFloat>::quiet_NaN());
    pdf_abserr.fill(std::numeric_limits<myFloat>::quiet_NaN());
    ddx_pdf_abserr.fill(std::numeric_limits<myFloat>::quiet_NaN());
    good_theta2.fill(false);
    g_dd_theta2.fill(std::numeric_limits<myFloat>::quiet_NaN());
  }
  void print(ostream& os) {
    for (int i=0; i<xs.size(); ++i) {
      os << setw(30) << setprecision(20) << scientific << alpha << " "
      << setw(30) << setprecision(20) << scientific << beta << " "
      << setw(30) << setprecision(20) << scientific << xs.at(i) << " "
      << setw(30) << setprecision(20) << scientific << cdf.at(i) << " "
      << setw(30) << setprecision(20) << scientific << pdf.at(i) << " "
      << setw(30) << setprecision(20) << scientific << ddx_pdf.at(i) << " "
      << setw(15) << setprecision(5) << scientific << cdf_abserr.at(i) << " "
      << setw(15) << setprecision(5) << scientific << pdf_abserr.at(i) << " "
      << setw(15) << setprecision(5) << scientific << ddx_pdf_abserr.at(i)
      << setw(10) << right << good_theta2.at(i) << " "
      << setw(15) << setprecision(5) << scientific << g_dd_theta2.at(i) << endl;
    }
  }
};
  

#define NBUFFERS 20

template<typename myFloat>
class result_buffer {
public:
  array<result<myFloat>, NBUFFERS> results;
  
  void initialize(array<double, 407>& xs) {
    front = 0;
    back = 0;
    number_in_use = 0;
    for (int i=0; i<NBUFFERS; ++i)
      results.at(i).xs = xs;
  }
  
  int reserve_buffer() {
    unique_lock<mutex> lock(buffer_mutex);
    not_full.wait(lock, [this]{return is_not_full();});
    int i_reserved = back;
    back = (back+1) % NBUFFERS;
    ++number_in_use;
    lock.unlock();
    return i_reserved;
  }
  
  void finish_buffer(int i) {
    unique_lock<mutex> lock(buffer_mutex);
    results.at(i).finished = true;
    if (i == front)
      next_available.notify_one();
  }
  
  void print_front(ostream& os) {
    unique_lock<mutex> lock(buffer_mutex);
    next_available.wait(lock, [this] {return next_ready();});
    results.at(front).print(os);
    results.at(front).finished = false;
    --number_in_use;
    front = (front+1) % NBUFFERS;
    lock.unlock();
    not_full.notify_one();
  }

private:
  int front;
  int back;

  size_t number_in_use;
  bool next_ready() const { return number_in_use > 0 && results.at(front).finished; }
  bool is_not_full() const { return number_in_use < NBUFFERS; }
  
  mutex buffer_mutex;
  condition_variable next_available;
  condition_variable not_full;
  
};

struct AB {
  double alpha;
  double beta;
};

class Jobs {
public:
  deque<AB> abs;
  mutex jobs_mutex;
};

mutex cout_mutex;

template<typename myFloat>
void reset_prec(int digits) {}

template<>
void reset_prec<mpreal>(int digits) {
  mpreal::set_default_prec(digits);
}

template<typename myFloat>
void free_cache() {}

template<>
void free_cache<mpreal>() {
  mpfr_free_cache();
}

template<typename myFloat, typename BigFloat>
void calculate_results(int thread_id, int noext, Kronrod<BigFloat> g_k_big,
                       myFloat epsabs, myFloat epsrel, int subdivisions, int verbose,
                       Jobs* jobs, result_buffer<myFloat>* res_buf, int digits) {
  reset_prec<myFloat>(digits);
  
  IntegrationController<myFloat> cntl(noext, g_k_big, epsabs, epsrel, subdivisions, verbose);
  IntegrationController<double> cntl_double(noext, g_k_big, static_cast<double>(epsabs),
                                            static_cast<double>(epsrel), subdivisions, verbose);
  Controllers<myFloat> ctls(cntl, cntl_double);
  {
    unique_lock<mutex> lock(cout_mutex);
    cout << "Starting thread " << thread_id << " with IntegraationController at " << &cntl << endl;
    cout << "Machine epsilon used: " << std::numeric_limits<myFloat>::epsilon() << endl;
    cout << "Digits10: " << int(digits * log(2)/log(10)) << endl;
    cout << "Input epsabs: " << epsabs << ", epsrel: " << epsrel << endl;

    cout << cntl << endl;
  }
  while(true) {
    int i;
    AB ab;
    int verbose = 0;
    int log_flag=false;
    int lower_tail = true;

    {
      unique_lock<mutex> lock(jobs->jobs_mutex);
      if (jobs->abs.size() == 0)
        return;
      else {
        ab = jobs->abs.front();
        jobs->abs.pop_front();
      }
    }
    i = res_buf->reserve_buffer();
    result<myFloat>& res = res_buf->results.at(i);
    myFloat alpha = static_cast<myFloat>(ab.alpha);
    myFloat beta = static_cast<myFloat>(ab.beta);
    StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctls, verbose);
    for (int j = 0; j<res.xs.size(); ++j) {
      myFloat x = static_cast<myFloat>(res.xs.at(j));
      res.alpha = alpha;
      res.beta = beta;
      res.cdf.at(j) = std_stable_dist.cdf(x, lower_tail, log_flag);
      res.cdf_abserr.at(j) = std_stable_dist.abserr;
      res.good_theta2.at(j) = std_stable_dist.good_theta2;
      res.g_dd_theta2.at(j) = std_stable_dist.g_dd_theta2;
      /*
       int p_termination_code = std_stable_dist.termination_code;
       */
      res.pdf.at(j) = std_stable_dist.pdf(x, log_flag);
      res.pdf_abserr.at(j) = std_stable_dist.abserr;
      /*
       int d_termination_code = std_stable_dist.termination_code;
       */
      res.ddx_pdf.at(j) = std_stable_dist.ddx_pdf(x);
      res.ddx_pdf_abserr.at(j) = std_stable_dist.abserr;
      /*
       int ddx_d_termination_code = std_stable_dist.termination_code;
       */
    }
    res_buf->finish_buffer(i);
  }
  free_cache<myFloat>();
}

template<typename myFloat, typename BigFloat>
int dump(Kronrod<BigFloat> g_k_big, int digits, string old_file_name, double alpha_low, double alpha_high) {
  high_resolution_clock::time_point start = high_resolution_clock::now();
  Eigen::initParallel();

  cout << "std::numeric_limits<myFloat>::epsilon() = " << std::numeric_limits<myFloat>::epsilon() << endl;
  StandardStableDistribution<myFloat>::initialize();
  int digits10 = digits * log(2)/log(10);
  cout.precision(digits10);
  cout << "zeta_tol:" << StandardStableDistribution<myFloat>::zeta_tol << endl;
  
  // Parameters for the integration controllers
  
  bool noext = true;   //disable extrapolation
  myFloat epsabs = 0;
  myFloat epsrel = 64*std::numeric_limits<myFloat>::epsilon();
  int subdivisions = 1000;
  int verbose = 0;
  
  // Set up iteration ranges for alpha, beta and x
  
  array<double, 24> alphas {{.01, .1, .2, .3, .4, .5, .6, .7, .8, .9, .99, 1.,
    1.01, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 1.99, 2.0}};
  
    array<double, 23> betas {{-1.0, -.99, -.9, -.8, -.7, -.6, -.5, -.4, -.3, -.2, -.1, 0,
      .1, .2, .3, .4, .5, .6, .7, .8, .9, .99, 1.}};
  
  string out_file = "../output/stable_mpreal.out";
  cout << "Writing output to " + out_file << endl;
  ofstream out(out_file) ;
  
  ifstream in(old_file_name);
  while (in) {
    string str_buf;
    getline(in, str_buf);
    stringstream ss(str_buf);
    double alpha; ss >> alpha;
    if (alpha < alpha_low) {
      out << str_buf << endl;
    } else {
      break;
    }
  }
  Jobs jobs;
  for (auto alpha : alphas) {
    if (alpha < alpha_low || alpha > alpha_high) continue;
    for (auto beta : betas)
      jobs.abs.push_back(AB{alpha, beta});
  }
  int number_of_jobs = static_cast<int>(jobs.abs.size());
    
  array<double, 407> xs;
  xs.at(0) = -std::numeric_limits<double>::infinity();
  xs.at(1) = -pow(10.,300);
  xs.at(2) = -pow(10.,200);
  xs.at(3) = -pow(10.,150);
  xs.at(4) = -pow(10.,120);
  for (int i = 5; i<104; ++i)
    xs.at(i) = -pow(10.,104-i);  //  -1e99, -1e98, ... -1e1
  for (int i = 104; i<=202; ++i)
    xs.at(i) = -static_cast<double>(203-i)/10; // -9.9, -9.8, ...,-.1
  xs.at(203) = 0;
  for (int i=204; i<=302; i++)
    xs.at(i) = static_cast<double>(i-203)/10; // .1, .2, ... 9.9
  for (int i=303; i<402; i++)
    xs.at(i) = pow(10.,i-302);   // 1e1, 1e2, ..., 1e99
  xs.at(402) = pow(10.,120);
  xs.at(403) = pow(10.,150);
  xs.at(404) = pow(10.,200);
  xs.at(405) = pow(10.,300);
  xs.at(406) = std::numeric_limits<double>::infinity();
  
  result_buffer<myFloat> res_buf;
  res_buf.initialize(xs);
  
  /*
  out << setw(30) << right << "alpha" << " "
  << setw(30) << right << "beta" << " "
  << setw(30) << right << "x" << " "
  << setw(30) << right << "cdf" << " "
  << setw(30) << right << "pdf" << " "
  << setw(30) << right << "ddx_pdf" << " "
  << setw(15) << right << "cdf.abserr" << " "
  << setw(15) << right << "pdf.abserr" << " "
  << setw(15) << right << "ddx_pdf.abserr" << " "
  << setw(10) << right <<  "good_th2" << " "
  << setw(15) << right << "g_dd_theta2" << endl;
   */
  
  thread t0(calculate_results<mpreal, mpreal>,
            1, noext, g_k_big, epsabs, epsrel, subdivisions, verbose, &jobs, &res_buf, digits);
  thread t1(calculate_results<mpreal, mpreal>,
            2, noext, g_k_big, epsabs, epsrel, subdivisions, verbose, &jobs, &res_buf, digits);
  thread t2(calculate_results<mpreal, mpreal>,
            3, noext, g_k_big, epsabs, epsrel, subdivisions, verbose, &jobs, &res_buf, digits);
  thread t3(calculate_results<mpreal, mpreal>,
            4, noext, g_k_big, epsabs, epsrel, subdivisions, verbose, &jobs, &res_buf, digits);
  thread t4(calculate_results<mpreal, mpreal>,
            5, noext, g_k_big, epsabs, epsrel, subdivisions, verbose, &jobs, &res_buf, digits);
  thread t5(calculate_results<mpreal, mpreal>,
            6, noext, g_k_big, epsabs, epsrel, subdivisions, verbose, &jobs, &res_buf, digits);
  thread t6(calculate_results<mpreal, mpreal>,
            7, noext, g_k_big, epsabs, epsrel, subdivisions, verbose, &jobs, &res_buf, digits);
  thread t7(calculate_results<mpreal, mpreal>,
            8, noext, g_k_big, epsabs, epsrel, subdivisions, verbose, &jobs, &res_buf, digits);
  
  for (int i=0; i<number_of_jobs; ++i) {
    res_buf.print_front(out);
  }
  
  t0.join();
  t1.join();
  t2.join();
  t3.join();
  t4.join();
  t5.join();
  t6.join();
  t7.join();
  
  while (in) {
    string str_buf;
    getline(in, str_buf);
    stringstream ss(str_buf);
    double alpha; ss >> alpha;
    if (alpha > alpha_high) {
      out << str_buf << endl;
    } else {
      continue;
    }
  }
  high_resolution_clock::time_point end = high_resolution_clock::now();
  duration<double> elapsed_time = end - start;
  cout << "Number of function calls = " << alphas.size() * betas.size() * xs.size() * 3 << endl;
  cout << "Elapsed time = " << fixed << setprecision(2) << elapsed_time.count() << " seconds" << endl;
  
  return 0;
}

static void show_usage (string name){
  cerr << "Usage: " << name << "[old_file_name alpha_low alpha_high]" << endl;
}

int main(int argc, char *argv[]) {
  
  // Check the number of parameters
  if (argc != 1 && argc != 4) {
    // Tell the user how to run the program
    show_usage(string(argv[0]));
    return 1;
  }
  string old_file_name;
  double alpha_low = 0;
  double alpha_high = 2;
  if (argc == 4) {
    stringstream ss1((string(argv[1]))), ss2((string(argv[2]))), ss3((string(argv[3])));
    ss1 >> old_file_name;
    ss2 >> alpha_low;
    ss3 >> alpha_high;
  }
  
  mpreal::set_default_prec(128);
  Kronrod<mpreal> g_k_big(10);
  mpreal::set_default_prec(96);
  return dump<mpreal, mpreal>(g_k_big, 96, old_file_name, alpha_low, alpha_high);

}



