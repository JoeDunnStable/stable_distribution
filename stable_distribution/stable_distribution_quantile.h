//
/// \file stable_distribution_quantile.h
/// Implementation of quantile of standard stable distribution
/// Included in stable_distribution.h when LIBRARY is defined
/// \author Joseph Dunn
/// \copyright 2016, 2018 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include "q_guess.h"
#include "bracket_and_solve.h"

namespace stable_distribution {
  
using std::pair;
using boost::math::erf_inv;
using boost::math::erfc_inv;

// Functor passed to toms748 solve to find q for unit stable distribution
template<typename myFloat>
class p_solve {
private:
  myFloat p;
  StandardStableDistribution<myFloat>* std_stable_dist;
  int lower_tail;
  int log_p;
  Parameterization pm;
  
public:
  int neval;
  p_solve(myFloat p, StandardStableDistribution<myFloat>* std_stable_dist,
                 int lower_tail, int log_p,
          Parameterization pm) :
  p(p), std_stable_dist(std_stable_dist), lower_tail(lower_tail),
  log_p(log_p), pm(pm), neval(0) {}
  myFloat operator()(const myFloat q) {
    Fmt<myFloat> fmt;
    if (std_stable_dist->verbose)
      cout << "Calling cdf with parmerters" << endl
           << "q = "  << fmt << q << ", alpha = "  << fmt << std_stable_dist->alpha
           << ", beta = " << fmt << std_stable_dist->beta_input << endl
           << "targeting p = " << fmt << p << endl;
    
    myFloat r = std_stable_dist->cdf(q, lower_tail, log_p, pm)-p;
    if (std_stable_dist->verbose)
      cout << ", Resulting delta = " << fmt << r << endl;
    neval += std_stable_dist->neval;
    return r;
  }
};
  
template<typename myFloat>
double Controllers<myFloat>::q_guess(myFloat p, myFloat alpha,
                                myFloat beta, int lower_tail, int log_p) {
    Controllers<double> ctls_double(ctl_double, ctl_double);
    StandardStableDistribution<double> dist_double(static_cast<double>(alpha),
                                                   static_cast<double>(beta),
                                                   ctls_double, 0);
  double q_tol = 64*std::numeric_limits<double>::epsilon();
  return dist_double.quantile(static_cast<double>(p), lower_tail, log_p, q_tol);
  }
  
template<>
  double Controllers<double>::q_guess(double p, double alpha, double beta,
                                      int lower_tail, int log_p) {
    return stable_distribution::q_guess(p, alpha, beta, lower_tail, log_p);
  }
  
// Returns a quantile for the unit stable distribution.
template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::quantile(myFloat p, int lower_tail, int log_p, myFloat dbltol, Parameterization pm) {

  //Default values which will be reset if the integrator is used
  abserr = 0;
  neval = 0;
  termination_code = IntegrationController<myFloat>::TerminationCode::normal;
  last = 0;
  num_iter = 0;

  if (lower_tail) {
    if (fabs(log_p ? exp(p) : p) ==0)
      return NegInf;
    else if (fabs(1-(log_p ? exp(p) : p)) == 0)
      return PosInf;
  } else {
    if (fabs(log_p ? exp(p) : p) == 0)
      return PosInf;
    else if (fabs(1-(log_p ? exp(p) : p)) == 0)
      return NegInf;
  }
  if (alpha == 1 && beta_input == 0) {
    if (verbose)
      cout << "Returning inverse cdf of the Cauchy distribution. " << endl;
    neval = 0;
    abserr = 0;
    return lower_tail ? (log_p ? static_cast<myFloat>(tan((exp(p)-.5)*pi))
                               : tan((p-.5)*pi))
                      : (log_p ? static_cast<myFloat>(-tan((exp(p)-.5)*pi))
                               : -tan((p-.5)*pi)) ;
  }
  /*
  if (alpha == 2) {
    if (verbose)
      cout << "Returning inverse cdf of the normal distribution. " << endl;
    neval = 0;
    abserr = 0;
    return lower_tail ? (log_p ? static_cast<myFloat>(2*erf_inv(static_cast<myFloat>(2*exp(p)-1)))
                               : static_cast<myFloat>(2*erf_inv(2*p-1)))
                      : (log_p ? static_cast<myFloat>(-2*erf_inv(static_cast<myFloat>(2*exp(p)-1)))
                               : static_cast<myFloat>(-2*erf_inv(2*p-1)));
  }
  */
  p_solve<myFloat> p_s(p, this, lower_tail, log_p, pm);
  pair<myFloat,myFloat> r;
  RelativeComparisonTolerance<myFloat> tol(dbltol);
  myFloat guess = max(-1e300,min(1e300,controllers.q_guess(static_cast<double>(p),
                                              static_cast<double>(alpha),
                                              static_cast<double>(beta_input),
                                              lower_tail,log_p)));
  if (verbose)
    cout << "Guess for q " << fmt << guess << endl;
  myFloat factor = max<myFloat>(static_cast<myFloat>(1),static_cast<myFloat>(.1)*fabs(guess));
  bool rising = lower_tail;
  boost::uintmax_t maxiter = 1000;
  r=boost::math::tools::bracket_and_solve_root2(p_s,guess,factor,rising,tol,maxiter);
  if (verbose)
    cout << "r.first = " << fmt << r.first << ", r.second - " << fmt << r.second
         << " in " << maxiter << " iterations" << endl;
  neval = p_s.neval;
  num_iter = static_cast<int>(maxiter);
  return (r.first+r.second)/2;
}
} // namespace stable_distribution
