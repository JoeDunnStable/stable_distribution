/// \file parallel_duality_check.cpp
/// Checks that Zolotarev duality theorem is satisfied
/// \author Joseph Dunn
/// \copyright 2017 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::getline;

#include "stable_config.h"

#define MPREAL
#define MPFR_FLOAT
#define CPP_BIN_FLOAT
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

#include <algorithm>
using std::sort;

#include <boost/filesystem.hpp>

template<typename myFloat>
myFloat epsdiff(myFloat r1, myFloat r2) {
  if (r1 == r2)
    return 0;
  else if (!boost::math::isfinite(r1) || !boost::math::isfinite(r2))
    return std::numeric_limits<myFloat>::infinity();
  else {
    myFloat one{1};
    return fabs(r1-r2)/(max(one, max(fabs(r1), fabs(r2)))*std::numeric_limits<myFloat>::epsilon());
  }
}

template<typename myFloat>
string fmt_eps(myFloat eps) {
  stringstream ss;
  if (eps < 1e4)
    ss << setw(10) << setprecision(1) << fixed << eps;
  else
    ss << setw(10) << right << "*******" ;
  return ss.str();
}

template<typename myFloat>
struct Result {
  double x_prime;
  myFloat alpha_prime;
  myFloat beta_star;
  myFloat ln_pdf1;
  myFloat x;
  double alpha;
  double beta;
  myFloat ln_pdf2;
  myFloat eps_diff;
  Result(double x_prime, myFloat alpha_prime, myFloat beta_star, myFloat ln_pdf1,
         myFloat x, double alpha, double beta, myFloat ln_pdf2, myFloat eps_diff) :
  x_prime(x_prime), alpha_prime(alpha_prime), beta_star(beta_star), ln_pdf1(ln_pdf1),
  x(x), alpha(alpha), beta(beta), ln_pdf2(ln_pdf2), eps_diff(eps_diff){}
};

void print_result_heading(std::ostream& os) {
  os << setw(15) << right << "x_prime"
  << setw(15) << right << "alpha_prime"
  << setw(15) << right << "beta_star"
  << setw(25) << right << "ln_pdf1"
  << setw(15) << right << "x"
  << setw(15) << right << "alpha"
  << setw(15) << right << "beta"
  << setw(25) << right << "ln_pdf2"
  << setw(10) << right << "eps_diff" << endl << endl;
}

template<typename myFloat>
std::ostream& operator<< (std::ostream& os, const Result<myFloat>& r) {
  os << setw(15) << setprecision(5) << scientific << r.x_prime
  << setw(15) << setprecision(10) << fixed << r.alpha_prime
  << setw(15) << setprecision(10) << fixed << r.beta_star
  << setw(25) << setprecision(16) << scientific << r.ln_pdf1
  << setw(15) << setprecision(5) << scientific << r.x
  << setw(15) << setprecision(10) << fixed << r.alpha
  << setw(15) << setprecision(10) << fixed << r.beta
  << setw(25) << setprecision(16) << scientific << r.ln_pdf2
  << setw(10) << fmt_eps(r.eps_diff) << endl;
  return os;
}

template<typename myFloat>
bool comp(Result<myFloat>& lhs, Result<myFloat>& rhs) {
  return lhs.eps_diff > rhs.eps_diff;
}

template<typename myFloat>
class Results;

template<typename myFloat>
ostream& operator<< (ostream&, Results<myFloat>&);

template<typename myFloat>
class Results {
public:
  static vector<double> probs;
  vector<Result<myFloat> > data;
  vector<myFloat> eps_all;
  double eps_best = -std::numeric_limits<double>::infinity();
  double eps_worst = std::numeric_limits<double>::infinity();
  double eps_sum = 0;
  int count = 0;
  int max_size;
  vector<myFloat> quantile() {
    long np = probs.size();
    vector<myFloat> qs(np);
    long n = eps_all.size();
    if (n > 0 && np > 0) {
      std::sort(eps_all.begin(), eps_all.end());
      for (int j=0; j<np; ++j) {
        double index = (n - 1) * probs.at(j);
        int lo = static_cast<int>(floor(index));
        int hi = static_cast<int>(ceil(index));
        double h = index - lo;
        qs.at(j) = (1-h) * eps_all.at(lo) + h * eps_all.at(hi);
      }
      return qs;
    } else {
      throw std::range_error("quantile: Both x and prob must be of length greater than 0");
    }
  }  //quantile
  
  Results(int max_size) : max_size(max_size) {}
  void add_result(Result<myFloat> res) {
    
    if (boost::math::isnan(res.ln_pdf1) || boost::math::isnan(res.ln_pdf2) )
      cout << res << endl;
    else {
      count++;
      eps_sum += static_cast<double>(res.eps_diff);
      eps_all.push_back(res.eps_diff);
      if (res.eps_diff > eps_best)
      {
        data.push_back(res);
        sort(data.begin(), data.end(), comp<myFloat>);
        if (data.size() > max_size) data.pop_back();
        eps_best = static_cast<double>(data.back().eps_diff);
        eps_worst = static_cast<double>(data.front().eps_diff);
      }
    }
  };
  
  bool pass() {
    return (eps_worst < 200) && (eps_sum/count < 10);
  }
  
  friend ostream& operator<< <myFloat>(ostream& os, Results<myFloat>& r);
};

template<typename myFloat>
vector<double> Results<myFloat>::probs = {.01, .05, .25, .5, .75, .95, .99};

template<typename myFloat>
ostream& operator<<(ostream& os, Results<myFloat>& r) {
  vector<myFloat> qs = r.quantile();
  os << "Table of worst " << r.data.size() << " out of " << r.count << " eps diffs" << endl << endl;
  print_result_heading(os);
  for (auto result : r.data) {
    os << result;
  }
  
  os << endl << setw(140) << "Average"
  << setw(10) << fmt_eps(r.eps_sum/r.count) << endl << endl;
  os << setw(140) << "Quantile" << endl << endl;
  for (int i =0; i<Results<myFloat>::probs.size(); ++i)
    os << setw(139) << fixed << setprecision(0) << Results<myFloat>::probs.at(i)*100 << "%"
       << setw(10) << fmt_eps(qs.at(i)) << endl;

  return os;
}

template<typename myFloat>
class JobOutput {
public:
  bool finished;
  vector<Result<myFloat> > results;

  JobOutput() : finished(false) {}
  void add_result(Result<myFloat> result) {results.push_back(result);}
  void clear() {results.clear();}
  void transfer(Results<myFloat>& final_results) {
    for (auto result : results)
      final_results.add_result(result);
    results.clear();
  }
};


#define NBUFFERS 20

template<typename myFloat>
class JobBuffer {
public:
  array<JobOutput<myFloat>, NBUFFERS> job_outputs;
  
  JobBuffer() {
    front = 0;
    back = 0;
    number_in_use = 0;
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
    job_outputs.at(i).finished = true;
    if (i == front)
      next_available.notify_one();
  }
  
  void transfer_front(Results<myFloat>& final_results) {
    unique_lock<mutex> lock(buffer_mutex);
    next_available.wait(lock, [this] {return next_ready();});
    job_outputs.at(front).transfer(final_results);
    job_outputs.at(front).finished = false;
    --number_in_use;
    front = (front+1) % NBUFFERS;
    lock.unlock();
    not_full.notify_one();
  }
  
private:
  int front;
  int back;
  
  size_t number_in_use;
  bool next_ready() const { return number_in_use > 0 && job_outputs.at(front).finished; }
  bool is_not_full() const { return number_in_use < NBUFFERS; }
  
  mutex buffer_mutex;
  condition_variable next_available;
  condition_variable not_full;
  
};

struct AB {
  double alpha_m_1;
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

template<typename myFloat>
void duality_check(double alpha_m_1, double beta, vector<double>& xs,
                   Controllers<myFloat> ctls, JobOutput<myFloat>& job_output) {
  myFloat pi2 = const_pi<myFloat>()/2;
  myFloat inf = std::numeric_limits<myFloat>::infinity();
  int log_flag = 1;
  // Zolotarev Theorem 2.3.2
  if (!(alpha_m_1 > 0 && -1 <= beta && beta <=1)) {
    throw std::range_error("duality check parameter error");
  }
  myFloat alpha = myFloat(alpha_m_1)+1;
  /*
  bool near = abs(alpha_m_1) < StandardStableDistribution<myFloat>::threshhold_1*abs(beta);
   */
  myFloat alpha_prime, alpha_prime_m_1, zeta, Q, beta_star, D;
  if (true /* near */) {
    alpha_prime_m_1 = expm1(-log1p(myFloat(alpha_m_1)));
    alpha_prime = 1 + alpha_prime_m_1;
    zeta = beta/tan(pi2*alpha_m_1);
    Q = (zeta<0) ? atan(1/abs(zeta))/pi2
    : 2-atan(1/abs(zeta))/pi2;
    beta_star = (alpha==2 || beta==-1)
    ? 1
    :-tan(pi2*alpha_prime_m_1)/tan(pi2*Q*alpha_prime);
    D = pow(1+zeta*zeta,1/(2*myFloat(alpha)))
        /sin(pi2*Q*alpha_prime);
  } else {
    alpha_prime = myFloat(1)/alpha;
    alpha_prime_m_1 = alpha_prime - 1;
    zeta = -beta*tan(pi2*alpha);
    Q = 1-atan(-zeta)/pi2;
    beta_star = (alpha==2 || beta==-1)
    ? 1
    :1/(tan(pi2*alpha_prime)*tan(pi2*Q*alpha_prime));
    D = pow(1+pow(beta*tan(pi2*alpha),2),1/(2*myFloat(alpha)))
    /sin(pi2*Q*alpha_prime);
  }
  myFloat min_ln = log(std::numeric_limits<myFloat>::min());
  StandardStableDistribution<myFloat> dist1(AlphaMinusOne<myFloat>(alpha_prime_m_1), beta_star, ctls,0);
  StandardStableDistribution<myFloat> dist2(AlphaMinusOne<myFloat>(alpha_m_1), beta, ctls,0);
  for (int i=0; i<xs.size(); ++i) {
    double x_prime = xs.at(i);
    myFloat x = D*pow(myFloat(x_prime),-alpha_prime);
    myFloat ln_pdf1 = dist1.pdf(x_prime, log_flag, S1);
    myFloat ln_pdf2 = log(D) + (-1-alpha_prime)*log(myFloat(x_prime))
    +dist2.pdf(x, log_flag, S1);
    if (ln_pdf1 < min_ln) ln_pdf1 = -inf;
    if (ln_pdf2 < min_ln) ln_pdf2 = -inf;
    myFloat eps_diff = epsdiff(ln_pdf1, ln_pdf2);
    Result<myFloat> result(x_prime, alpha_prime, beta_star, ln_pdf1, x,
                           static_cast<double>(alpha), beta, ln_pdf2, eps_diff);
    job_output.add_result(result);
  } // x
  return;
}

template<typename myFloat, typename BigFloat>
void calculate_job_outputs(int thread_id, int noext, Kronrod<BigFloat> g_k_big,
                       myFloat epsabs, myFloat epsrel, int subdivisions, int verbose,
                       vector<double> xs,
                       Jobs* jobs, JobBuffer<myFloat>* job_buffer, int digits) {
  reset_prec<myFloat>(digits);
  
  unique_lock<mutex> cout_lock(cout_mutex);
  IntegrationController<myFloat> cntl(noext, g_k_big, epsabs, epsrel, subdivisions, verbose);
  IntegrationController<double> cntl_double(noext, g_k_big, static_cast<double>(epsabs),
                                            static_cast<double>(epsrel), subdivisions, verbose);
  Controllers<myFloat> ctls(cntl, cntl_double);
  
  cout << "Starting thread " << thread_id << " with IntegraationController at " << &cntl << endl;
  cout << "Machine epsilon used: " << std::numeric_limits<myFloat>::epsilon() << endl;
  /*
   cout << "Digits10: " << int(digits * log(2)/log(10)) << endl;
   cout << "Input epsabs: " << epsabs << ", epsrel: " << epsrel << endl;
   
   cout << cntl << endl;
   */
  cout_lock.unlock();

  while(true) {
    int i;
    AB ab;
    
    {
      unique_lock<mutex> lock(jobs->jobs_mutex);
      if (jobs->abs.size() == 0)
        return;
      else {
        ab = jobs->abs.front();
        jobs->abs.pop_front();
      }
    }
    i = job_buffer->reserve_buffer();
    JobOutput<myFloat>& res = job_buffer->job_outputs.at(i);
    cout_lock.lock();
    duality_check(ab.alpha_m_1, ab.beta, xs, ctls, res);
    cout_lock.unlock();
    job_buffer->finish_buffer(i);
  }
  free_cache<myFloat>();
}

template<typename myFloat, typename BigFloat>
int create_jobs(Kronrod<BigFloat> g_k_big, ostream& os) {
  high_resolution_clock::time_point start = high_resolution_clock::now();
  Eigen::initParallel();
  
  Results<myFloat> final_results(100);
  
  cout << "std::numeric_limits<myFloat>::epsilon() = " << std::numeric_limits<myFloat>::epsilon() << endl;
  StandardStableDistribution<myFloat>::initialize();
  int digits = static_cast<int>(-log(std::numeric_limits<myFloat>::epsilon())/log(2.));
  int digits10 = static_cast<int>(-log(std::numeric_limits<myFloat>::epsilon())/log(10.));
  os.precision(digits10);
  
  // Parameters for the integration controllers
  
  bool noext = true;   //disable extrapolation
  myFloat epsabs = 0;
  myFloat epsrel = 64*std::numeric_limits<myFloat>::epsilon();
  int subdivisions = 1000;
  int verbose = 0;
  StandardStableDistribution<myFloat>::initialize();
  
  // Set up iteration ranges for alpha, beta and x
  
  vector<double> alphas_m_1;
  alphas_m_1.push_back(1./1024.);
  alphas_m_1.push_back(1./128.);  // 1/128
  for (int i=1; i<=12; ++i)
    alphas_m_1.push_back(i/12.);
  
  vector<double> betas;
  for (int i=-16; i<=16; ++i)
    betas.push_back(i/16.);
  Jobs jobs;
  for (auto alpha_m_1 : alphas_m_1) {
    for (auto beta : betas)
      jobs.abs.push_back(AB{alpha_m_1, beta});
  }
  int number_of_jobs = static_cast<int>(jobs.abs.size());
  
  vector<double> xs;
 
  for (int i=-32; i<=32; ++i)
    xs.push_back(pow(2.,i));

  JobBuffer<myFloat> job_buffer;
  
  thread t0(calculate_job_outputs<myFloat, BigFloat>,
            1, noext, g_k_big, epsabs, epsrel, subdivisions, verbose, xs, &jobs, &job_buffer, digits);
  thread t1(calculate_job_outputs<myFloat, BigFloat>,
            2, noext, g_k_big, epsabs, epsrel, subdivisions, verbose, xs, &jobs, &job_buffer, digits);
  thread t2(calculate_job_outputs<myFloat, BigFloat>,
            3, noext, g_k_big, epsabs, epsrel, subdivisions, verbose, xs, &jobs, &job_buffer, digits);
  thread t3(calculate_job_outputs<myFloat, BigFloat>,
            4, noext, g_k_big, epsabs, epsrel, subdivisions, verbose, xs, &jobs, &job_buffer, digits);
  thread t4(calculate_job_outputs<myFloat, BigFloat>,
            5, noext, g_k_big, epsabs, epsrel, subdivisions, verbose, xs, &jobs, &job_buffer, digits);
  thread t5(calculate_job_outputs<myFloat, BigFloat>,
            6, noext, g_k_big, epsabs, epsrel, subdivisions, verbose, xs, &jobs, &job_buffer, digits);
  thread t6(calculate_job_outputs<myFloat, BigFloat>,
            7, noext, g_k_big, epsabs, epsrel, subdivisions, verbose, xs, &jobs, &job_buffer, digits);
  thread t7(calculate_job_outputs<myFloat, BigFloat>,
            8, noext, g_k_big, epsabs, epsrel, subdivisions, verbose, xs, &jobs, &job_buffer, digits);
  
  for (int i=0; i<number_of_jobs; ++i) {
    job_buffer.transfer_front(final_results);
  }
  
  t0.join();
  t1.join();
  t2.join();
  t3.join();
  t4.join();
  t5.join();
  t6.join();
  t7.join();
  
  os << final_results;
  
  high_resolution_clock::time_point end = high_resolution_clock::now();
  duration<double> elapsed_time = end - start;
  os << "Number of function calls = " << alphas_m_1.size() * betas.size() * xs.size() * 2 << endl;
  os << "Elapsed time = " << fixed << setprecision(2) << elapsed_time.count() << " seconds" << endl;
  
  return !final_results.pass();
}

static void show_usage (string name){
  boost::filesystem::path p(name);
  cerr << "Usage: " << p.filename().string() << " test_name"<< endl << endl
  << "         where test_name is double"
#ifdef MPREAL
  << ", mpreal"
#endif
#ifdef CPP_BIN_FLOAT
  << ", cpp_bin_float"
#endif
#ifdef MPFR_FLOAT
  << ", mpfr_float"
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
  
  string out_dir = string("../output-") + 
	           string(PACKAGE_VERSION) + 
		   string("-") + 
		   string(PACKAGE_COMPILER);
  if (!boost::filesystem::is_directory(out_dir))
    boost::filesystem::create_directory(out_dir);
  string test_name{argv[1]};
  string out_file_name = out_dir + "/duality_check_" + test_name + ".out";
  string title = "Duality Check for " + test_name ;
  
#ifdef MPREAL
  // First create some high precision coeficients for
  mpreal::set_default_prec(128);
  Kronrod<mpreal> k_big(10);
  mpreal::set_default_prec(96);
#define BIGFLOAT mpreal
#elif MPFR_FLOAT
  Kronrod<BigMpfrFloat> k_big(10);
#define BIGFLOAT BigMpfrFloat
#elif CPP_BIN_FLOAT
  Kronrod<BigCppBinFloat> k_big(10);
#define BIGFLOAT BigCppBinFloat
#else
  Kronrod<double> k_big(10);
#define BIGFLOAT double
#endif

  if (test_name == "double") {
    cout << title << endl << endl;
    cout << "Writing to file " << out_file_name << endl << endl;
    ofstream os{out_file_name};
    os << title
       << ".  Machine epsilon = " << std::numeric_limits<double>::epsilon() << endl << endl;
    create_jobs<double, BIGFLOAT>(k_big, os);
  }
  
#ifdef MPREAL
  else if(test_name == "mpreal") {
    cout << title << endl << endl;
    cout << "Writing to file " << out_file_name << endl << endl;
    ofstream os{out_file_name};
    os << title
    << ".  Machine epsilon = " << std::numeric_limits<mpreal>::epsilon() << endl << endl;
    create_jobs<mpreal, mpreal>(k_big, os);
  }
#endif
  
#ifdef CPP_BIN_FLOAT
  else if (test_name == "cpp_bin_float") {
    Kronrod<BigCppBinFloat> k_big_cpp_bin_float(10);
    cout << title << endl << endl;
    cout << "Writing to file " << out_file_name << endl << endl;
    ofstream os{out_file_name};
    os << title
    << ".  Machine epsilon = " << std::numeric_limits<CppBinFloat>::epsilon() << endl << endl;
    create_jobs<CppBinFloat, BigCppBinFloat>(k_big_cpp_bin_float, os);
  }
#endif
  
#ifdef MPFR_FLOAT
  else if (test_name == "mpfr_float") {
    Kronrod<BigMpfrFloat> k_big_mpfr_float(10);
    cout << title << endl << endl;
    cout << "Writing to file " << out_file_name << endl << endl;
    ofstream os{out_file_name};
    os << title
    << ".  Machine epsilon = " << std::numeric_limits<MpfrFloat>::epsilon() << endl << endl;
    create_jobs<MpfrFloat, BigMpfrFloat>(k_big_mpfr_float, os);
  }
#endif
  else {
    show_usage(string(argv[0]));
    return 1;
    
  } //
  
  return 0;
}




