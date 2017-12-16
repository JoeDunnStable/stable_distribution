///
///  \file stable_distribution_mode.h
/// Implementation of mode of standard stable distribution
/// Included in stable_distribution.h when LIBRARY is defined
/// \author Joseph Dunn
/// \copyright 2016, 2017 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include <boost/math/tools/toms748_solve.hpp>
#include <boost/math/tools/minima.hpp>

namespace stable_distribution {
  
using std::max;
using std::min;

/// Functor passed to toms748_solve to find x such that ddx_pdf(x) == value
template<typename myFloat>
class DdxPdfSolve {
private:
  myFloat value;         ///< the target value
  StandardStableDistribution<myFloat>* std_stable_dist; ///< an instance of the StandardStableDistribution
  int verbose_mode;
  Parameterization pm;
  int x_flag;          ///< interpret input as log of x
public:
  /// constructor for functor
  DdxPdfSolve(myFloat value,                    ///< [in] the target value
                  StandardStableDistribution<myFloat>* std_stable_dist,  ///< [in] pointer to distribution
                  int verbose_mode,             ///< [in] indicator for verbose mode calculation
                  Parameterization pm,           ///< [in] the parameterization to use
                  int x_flag                  ///< [in] interpret the input as log of x or -x
                  )
  : value(value), std_stable_dist(std_stable_dist), verbose_mode(verbose_mode), pm(pm), x_flag(x_flag) {}
  /// return the value of the derivative wrt x os th pdf of the std. stable distribution
  myFloat operator()(const myFloat x_in) {
    myFloat ddx_min = std::numeric_limits<myFloat>::lowest();
    myFloat ddx_max = std::numeric_limits<myFloat>::max();
    myFloat x = x_in;
    switch (x_flag) {
      case 0:
        x=x_in;
        break;
      case 1:
        x=exp(x_in);
        break;
      case 2:
        x=-exp(x_in);
        break;
    }
    myFloat ret=max(ddx_min,min(ddx_max,std_stable_dist->ddx_pdf(x, pm)));
    if (verbose_mode) {
      myFloat zeta = -std_stable_dist->beta_input*tan(StandardStableDistribution<myFloat>::pi2*std_stable_dist->alpha);
      switch (pm) {
        case S0:
          cout << "x = " << x << ", x - zeta = " << x - zeta
               << ", ddx_pdf(x) = " << ret << endl;
          break;
        case S1:
          cout << "x = " << x + zeta << ", x - zeta = " << x
               << ", ddx_pdf(x) = " << ret << endl;
          break;
          
      }
    }
    return ret - value;
  }
};

/// Functor passed to brent_find_minima to find maximum of pdf
template<typename myFloat>
class Pdmaximum {
private:
  StandardStableDistribution<myFloat>* std_stable_dist; ///< an instance of the StandardStableDistribution
  int verbose_mode;
  Parameterization pm;
  int x_flag;          ///< interpret input as x, log of x or -x
public:
  /// constructor for functor
  Pdmaximum(StandardStableDistribution<myFloat>* std_stable_dist,  ///< [in] pointer to distribution
              int verbose_mode,             ///< [in] indicator for verbose mode calculation
              Parameterization pm,           ///< [in] the parameterization to use
              int x_flag                  ///< [in] interpret the input as log of x or -x
  )
  : std_stable_dist(std_stable_dist), verbose_mode(verbose_mode), pm(pm), x_flag(x_flag) {}
  /// return the - value of the pdf
  myFloat operator()(const myFloat x_in) {
    myFloat x = x_in;
    switch (x_flag) {
      case 0:
        x=x_in;
        break;
      case 1:
        x=exp(x_in);
        break;
      case 2:
        x=-exp(x_in);
        break;
    }
    int log_flag = 0;
    myFloat ret=-std_stable_dist->pdf(x, log_flag, pm);
    if (verbose_mode) {
      myFloat zeta = -std_stable_dist->beta_input*tan(StandardStableDistribution<myFloat>::pi2*std_stable_dist->alpha);
      switch (pm) {
        case S0:
          cout << "x = " << x << ", x - zeta = " << x - zeta
          << ", pdf(x) = " << -ret << endl;
          break;
        case S1:
          cout << "x = " << x + zeta << ", x - zeta = " << x
          << ", pdf(x) = " << -ret << endl;
          break;
          
      }
    }
    return ret;
  } // operator()
};
  
template<typename myFloat>
std::pair<myFloat,myFloat> StandardStableDistribution<myFloat>::mode(myFloat dbltol, int verbose_mode, Parameterization pm_requested)
{
  if(alpha * beta_input == 0){
    return std::pair<myFloat, myFloat>(0.,pdf(0.,false, pm_requested));
  }
  else  {
    myFloat alpha_proxy = 2 - 1e-5;
    myFloat upper, lower;
    myFloat outer, outer_low, inner_low;
    if (alpha < 1 && fabs(beta_input)==1) {
      myFloat eps=sqrt(std::numeric_limits<myFloat>::epsilon());
      myFloat beta_proxy = beta_input*(1-eps);
      if (verbose_mode)
        cout << "Using proxy beta = " << beta_proxy << endl;
      StandardStableDistribution<myFloat> proxy(alpha, beta_proxy, controllers, verbose);
      std::pair<myFloat, myFloat> mode = proxy.mode(dbltol, verbose_mode, pm_requested);
      mode.first = mode.first/(1-eps);
      return mode;
    } else if (alpha > alpha_proxy){
      if (verbose_mode)
        cout << "Using proxy alpha = " << alpha_proxy << endl;
      StandardStableDistribution<myFloat> proxy(alpha_proxy, beta_input, controllers, verbose);
      std::pair<myFloat, myFloat> mode = proxy.mode(dbltol, verbose_mode, pm_requested);
      mode.first = mode.first * (2-alpha)/(2-alpha_proxy);
      mode.second = pdf(mode.first, false, pm_requested);
      return mode;
    } else {
      Parameterization pm = (alpha < .5) ? S1 : S0;
      outer = .7;
      outer_low = max<myFloat>(log(std::numeric_limits<myFloat>::min())+5,-3*lgamma(1+1/alpha));
      inner_low = -lgamma(1+1/alpha);
      int x_flag = 0;
      if (tgamma(1+1/alpha) > std::numeric_limits<myFloat>::max())
          throw std::range_error("stable_mode: alpha is too small");
      
      
      if (beta_input>0){
        switch (pm) {
          case S0:
            lower=-outer;
            upper=0;
            break;
          case S1:
            lower = (alpha < .5) ? outer_low : -outer-zeta;
            upper = (alpha < .5) ? inner_low :  -zeta;
            x_flag = (alpha < .5) ? 1 : 0;
            break;
        }
      } else {
        switch (pm) {
          case S0:
            lower = 0;
            upper = outer;
            break;
          case S1:
            lower = (alpha < .5) ? outer_low : -zeta;
            upper = (alpha < .5) ? inner_low : outer - zeta;
            x_flag = (alpha < .5) ? 2 : 0;
            break;
        }
      }
      if (verbose_mode)
        cout << "mode: pm = " << pm << ", x_flag = " << x_flag << ", lower = " << lower << ", upper = " << upper << endl;
      boost::uintmax_t max_iter=10000;

/*
      DdxPdfSolve ddx_s(0, this, verbose_mode, pm, x_flag);
      RelativeComparisonTolerance tol(dbltol);
      std::pair<myFloat,myFloat> mode=boost::math::tools::toms748_solve(ddx_s,lower,upper,tol,max_iter);
      mode.first = .5 * (mode.first+mode.second);
      switch (x_flag) {
        case 1:
          mode.first = exp(mode.first);
          break;
        case 2:
          mode.first = -exp(mode.first);
          break;
        default:
          break;
      }
      mode.second = pdf(mode.first, false, pm);
 
 */
      int bits = std::numeric_limits<myFloat>::digits;
      Pdmaximum<myFloat> mode_s(this, verbose_mode, pm, x_flag);
      std::pair<myFloat, myFloat> mode=boost::math::tools::brent_find_minima(mode_s, lower, upper, bits, max_iter);
      mode.second = -mode.second;

      switch (x_flag) {
        case 1:
          mode.first = exp(mode.first);
          break;
        case 2:
          mode.first = -exp(mode.first);
          break;
        default:
          break;
      }

      if (pm != pm_requested) {
        switch (pm) {
          case S0:
            mode.first = mode.first-zeta;
            break;
          case S1:
            mode.first = mode.first + zeta;
            break;
        }
      }
      return mode;
    }
  }
}
  
} // namespace stable_distribution

