//
/// \file stable_distribution_quantile.h
/// \author Joseph Dunn
/// \copyright 2016 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include "q_guess.h"
#include <boost/math/tools/toms748_solve.hpp>

namespace boost { namespace math { namespace tools {
  
  /** modification of boost braket and solve root using -infinity as lower bound not 0 */
  template <class F, class T, class Tol, class Policy>
  std::pair<T, T> bracket_and_solve_root2(F& f, const T& guess, T factor, bool rising, Tol tol, boost::uintmax_t& max_iter, const Policy& pol)
  {
    BOOST_MATH_STD_USING
    static const char* function = "bracket_and_solve_root2<%1%>";
    //
    // Set up inital brackets:
    //
    T a = guess;
    T b = a;
    T fa = f(a);
    T fb = fa;
    //
    // Set up invocation count:
    //
    boost::uintmax_t count = max_iter - 1;
    
    int step = 32;
    
    if((fa < 0) == rising)
    {
      //
      // Zero is to the right of b, so walk upwards
      // until we find it:
      //
      while((boost::math::sign)(fb) == (boost::math::sign)(fa))
      {
        if(count == 0)
          return boost::math::detail::pair_from_single(policies::raise_evaluation_error(function, "Unable to bracket root, last nearest value was %1%", b, pol));
        //
        // Heuristic: normally it's best not to increase the step sizes as we'll just end up
        // with a really wide range to search for the root.  However, if the initial guess was *really*
        // bad then we need to speed up the search otherwise we'll take forever if we're orders of
        // magnitude out.  This happens most often if the guess is a small value (say 1) and the result
        // we're looking for is close to std::numeric_limits<T>::min().
        //
        if((max_iter - count) % step == 0)
        {
          factor *= 2;
          if(step > 1) step /= 2;
        }
        //
        // now go ahead and move our guess by "factor":
        //
        a = b;
        fa = fb;
        b += factor;
        fb = f(b);
        --count;
        BOOST_MATH_INSTRUMENT_CODE("a = " << a << " b = " << b << " fa = " << fa << " fb = " << fb << " count = " << count);
      }
    }
    else
    {
      //
      // Zero is to the left of a, so walk downwards
      // until we find it:
      //
      while((boost::math::sign)(fb) == (boost::math::sign)(fa))
      {
        if(count == 0)
          return boost::math::detail::pair_from_single(policies::raise_evaluation_error(function, "Unable to bracket root, last nearest value was %1%", a, pol));
        //
        // Heuristic: normally it's best not to increase the step sizes as we'll just end up
        // with a really wide range to search for the root.  However, if the initial guess was *really*
        // bad then we need to speed up the search otherwise we'll take forever if we're orders of
        // magnitude out.  This happens most often if the guess is a small value (say 1) and the result
        // we're looking for is close to numeric_limits<T>::min().
        //
        if((max_iter - count) % step == 0)
        {
          factor *= 2;
          if(step > 1) step /= 2;
        }
        //
        // n_gaussow go ahead and move are guess by "factor":
        //
        b = a;
        fb = fa;
        a -= factor;
        fa = f(a);
        --count;
        BOOST_MATH_INSTRUMENT_CODE("a = " << a << " b = " << b << " fa = " << fa << " fb = " << fb << " count = " << count);
      }
    }
    max_iter -= count;
    max_iter += 1;
    std::pair<T, T> r = toms748_solve(
                                 f, a, b, fa, fb,
                                 tol, count, pol);
    max_iter += count;
    BOOST_MATH_INSTRUMENT_CODE("max_iter = " << max_iter << " count = " << count);
    BOOST_MATH_LOG_COUNT(max_iter)
    return r;
  }
  
  template <class F, class T, class Tol>
  inline std::pair<T, T> bracket_and_solve_root2(F& f, const T& guess, const T& factor, bool rising, Tol tol, boost::uintmax_t& max_iter)
  {
    return bracket_and_solve_root2(f, guess, factor, rising, tol, max_iter, policies::policy<>());
  }
  
} // namespace tools
} // namespace math
} // namespace boost

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
  typename StandardStableDistribution<myFloat>::Parameterization pm;
  
public:
  int neval;
  p_solve(myFloat p, StandardStableDistribution<myFloat>* std_stable_dist,
                 int lower_tail, int log_p,
          typename StandardStableDistribution<myFloat>::Parameterization pm) :
  p(p), std_stable_dist(std_stable_dist), lower_tail(lower_tail),
  log_p(log_p), pm(pm), neval(0) {}
  myFloat operator()(const myFloat q) {
    if (std_stable_dist->verbose)
      cout << "Calling cdf with parmerters" << endl
           << "q = " << q << ", alpha = " << std_stable_dist->alpha
           << ", beta = " << std_stable_dist->beta_input << endl
           << "targeting p = " << p << endl;
    
    myFloat r = std_stable_dist->cdf(q, lower_tail, log_p, pm)-p;
    if (std_stable_dist->verbose)
      cout << ", Resulting delta = " << r << endl;
    neval += std_stable_dist->neval;
    return r;
  }
};
  
template<typename myFloat>
template<typename BigFloat>
AugmentedController<myFloat>::AugmentedController(bool noext, Kronrod<BigFloat> g_k_big,
                      myFloat epsabs,  myFloat epsrel,
                      int limit, int verbose)
  : IntegrationController<myFloat>(noext, g_k_big, epsabs, epsrel, limit, verbose),
  ctl_double(noext, g_k_big,
             min(64*std::numeric_limits<double>::epsilon(),static_cast<double>(epsabs)),
             min(64*std::numeric_limits<double>::epsilon(),static_cast<double>(epsrel)),
             limit, verbose) {}
  
template<>
template<typename BigFloat>
AugmentedController<double>::AugmentedController(bool noext, Kronrod<BigFloat> g_k_big,
                                                 double epsabs,  double epsrel,
                                                 int limit, int verbose)
  : IntegrationController<double>(noext, g_k_big, epsabs, epsrel, limit, verbose),
  ctl_double(noext, g_k_big, epsabs, epsrel, 0, verbose) {}
  
template<typename myFloat>
double AugmentedController<myFloat>::q_guess(myFloat p, myFloat alpha, myFloat beta,
                                              int lower_tail, int log_p) {
    StandardStableDistribution<double> dist(static_cast<double>(alpha),
                                            static_cast<double>(beta),
                                            ctl_double, 0);
    double q_tol = 64*std::numeric_limits<double>::epsilon();
    return dist.quantile(static_cast<double>(p), lower_tail, log_p, q_tol);
  }

template<>
double AugmentedController<double>::q_guess(double p, double alpha, double beta,
                                            int lower_tail, int log_p) {
    return q_guess(p, alpha, beta, lower_tail, log_p);
  }

// Returns a quantile for the unit stable distribution.
template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::quantile(myFloat p, int lower_tail, int log_p, myFloat dbltol, Parameterization pm) {

  //Default values which will be reset if the integrator is used
  abserr = 0;
  neval = 0;
  termination_code = IntegrationController<myFloat>::TerminationCode::normal;
  last = 0;

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
  myFloat guess = max(-1e300,min(1e300,q_guess(static_cast<double>(p),
                                              static_cast<double>(alpha),
                                              static_cast<double>(beta_input),
                                              lower_tail,log_p)));
  if (verbose)
    cout << "Guess for q " << guess << endl;
  myFloat factor = max<myFloat>(static_cast<myFloat>(1),static_cast<myFloat>(.1)*fabs(guess));
  bool rising = lower_tail;
  boost::uintmax_t maxiter = 1000;
  r=boost::math::tools::bracket_and_solve_root2(p_s,guess,factor,rising,tol,maxiter);
  if (verbose)
    cout << "r.first = " << r.first << ", r.second - " << r.second
         << " in " << maxiter << " iterations" << endl;
  neval = p_s.neval;
  num_iter = maxiter;
  return (r.first+r.second)/2;
}
} // namespace stable_distribution
