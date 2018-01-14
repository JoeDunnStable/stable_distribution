/// \file stable_distribution_cdf.h
/// Implementation of cdf for standard stable distribution
/// Included in stable_distribution.h when LIBRARY is defined
/// \author Joseph Dunn
/// \copyright 2016, 2017 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include <iostream>
#include <iomanip>
#include <boost/math/tools/toms748_solve.hpp>

namespace stable_distribution {
  
using std::string;
using std::vector;
using std::sort;
using std::setw;
using std::right;
using std::setprecision;
using std::pair;
using boost::math::tools::toms748_solve;

template<typename myFloat>
myFloat exp_m_x (myFloat x, StandardStableDistribution<myFloat>* std_stable_dist) {
  return (x<StandardStableDistribution<myFloat>::large_exp_arg)
            ? exp(-x)
            : static_cast<myFloat>(0);}

template<typename myFloat>
myFloat one_m_exp_m_x (myFloat x, StandardStableDistribution<myFloat>* std_stable_dist) {return -expm1(-x);}

template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::integrate_cdf() {
  // --- cdf()  for alpha < 2 ---
  // Returns the integral that preserves the most relative accuracy
  // cdf adjusts for whether this is the cdf or the complement of the cdf
  // or the difference between F.zeta and F.
  neval=0;
  if (verbose) cout << "integrate_cdf:" << endl;
  if(!good_theta2) {
    if (verbose)
      cout << "  Can't find theta2.  Proceeding without it." << endl;
  }
  myFloat r;

  // When x_m_zeta is large use the version that gets small as abs(x) becomes large
  // When x_m_zeta is small use the version that gets small as x approaches zeta
  bool use_one_m_exp_m_x = ((alpha<1) && !small_x_m_zet)
                            || ((alpha>1) && small_x_m_zet)
                            || ((alpha==1) && (beta>0));
  if (verbose)
    cout << "Integrand is " << (use_one_m_exp_m_x ? "1 - exp(-g)" : "exp(-g)") << endl;
  Integral_f_of_g<myFloat> int_g1 = use_one_m_exp_m_x
                                    ? Integral_f_of_g<myFloat>(&one_m_exp_m_x, this)
                                    : Integral_f_of_g<myFloat>(&exp_m_x, this);
  r=int_g1();
  myFloat c1;
  if (alpha != 1){
    c1 = 1/pi;
  } else{ // alpha == 1
    c1 = .5;
  }
  c_g_theta2_error = g_theta2_error;
  abserr=c1*int_g1.abserr;
  neval+=int_g1.neval;
  termination_code=int_g1.termination_code;
  last=int_g1.last;
  if (verbose) {
    cout << "  c1*sum(r)= " << c1 << " * " << r
         << " = " << c1  *(r) << endl
         << "  abs.err = " << c1*int_g1.abserr << endl
         << "  msg = " << int_g1.msg() << endl;
  }
  return c1*r;
} // StandardStableDistribution<myFloat>::integrate_cdf

// ------------------------------------------------------------------------------


template<typename myFloat>
myFloat retValue(myFloat F, int useF, int log_p) {
  return (useF) ? ((log_p) ? ((F==0)? StandardStableDistribution<myFloat>::NegInf : log(F)) : F)
  : ((log_p) ? ((F==1)? StandardStableDistribution<myFloat>::NegInf : boost::math::log1p(-F)) : 1 - F);
  
}

template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::cdf(myFloat x, int lower_tail, int log_p, Parameterization pm) {
  neval = 0;
  myFloat ret;
  set_x_m_zeta(x, pm);
  // default values which will be reset if the integrator is used
  c_g_theta2_error = 0;
  abserr = 0;
  neval = 0;
  termination_code = IntegrationController<myFloat>::TerminationCode::normal;
  last = 0;
  if (verbose)
    cout << "cdf: lower_tail = " << lower_tail << ", log_p = " << log_p << endl << *this;
  if(x==PosInf) {
    ret = retValue<myFloat>(1, lower_tail, log_p);
    if (verbose)
      cout << "  x is PosInf, returning " << ret << endl;
    return ret;
  }
  else if (x==NegInf) {
    ret = retValue<myFloat>(0, lower_tail, log_p);
    if (verbose)
      cout << "  x is NegInf, returing " << ret << endl;
    return ret;
  }
  switch (dist_type) {
    case Cauchy :
    {
      myFloat F;
      F = (fabs(x)<=1) ? static_cast<myFloat>(atan(-fabs(x))/pi + static_cast<myFloat>(.5))
                       : static_cast<myFloat>(atan(1/fabs(x))/pi);
      ret = retValue<myFloat>(F, (x>0)!=lower_tail, log_p);
      if (verbose)
        cout << "  Cauchy distribution, returning " << ret << endl;
      return ret;
    }
    case normal :
    {
      myFloat F;
      F = erfc(fabs(x/2))/2;
      ret = retValue<myFloat>(F, (x>0)!=lower_tail, log_p);
      if (verbose)
        cout << "  Normal distribution, returning " << ret << endl;
      return ret;
    }
    case fin_support :
      if(beta_input == 1 && x_m_zeta_input <=0) {
        ret = retValue<myFloat>(0, lower_tail, log_p);
      } else {
        ret = retValue<myFloat>(1, lower_tail, log_p);
      } 
      if (verbose)
        cout << "  x is outside of support, returning " << ret << endl;
      return ret;
    case other :
      if (verbose)
        cout << "cdf: General case:" << endl;
      myFloat F=0;
      if (alpha !=1 && small_x_m_zet) {
          myFloat F_zeta = (alpha<1 && fabs(beta)==1)
                           ? (lower_tail != (beta_input<1)) ? 0 : 1
                           :(lower_tail) ? static_cast<myFloat>(static_cast<myFloat>(.5) - theta0_x_gt_zeta/pi)
                                         : static_cast<myFloat>(static_cast<myFloat>(.5) + theta0_x_gt_zeta/pi);
          if (verbose)
            cout << "cdf: Using difference from F_zeta, " << F_zeta << endl;
          F = F_zeta;
        if (use_f_zeta) {
          if (verbose)
            cout << "cdf:: Using F_zeta" << endl;
          ret = (log_p) ? log(F) : F;
          return ret;
        } else {
          if (verbose)
            cout << "cdf: Using difference from F_zeta, " << F_zeta << endl;
        }
      }
      if (good_theta2) {
        myFloat integral = integrate_cdf();
        myFloat error = abserr;
        if ( (alpha_minus_one > 0 && beta == -1)
            || error < max(pow((alpha!=1)?x_m_zet:abs_x,-alpha),
                           4*controllers.controller.epsrel)* fabs(integral)
            || small_x_m_zet) {
          if (verbose)
            cout << "cdf: Integration error is small or x_m_zet is small.  Using integral" << endl;
          if (alpha != 1 && small_x_m_zet) {
            F -= ((lower_tail != (x_m_zeta_input > 0)) ? 1 : -1) * integral;
            ret = (log_p) ? log(F) : F;
            if (verbose)
              cout << "cdf returning " << ret << endl;
            return ret;
          } else {
            bool useF = !((x >= zeta && lower_tail) || (x < zeta && !lower_tail));
            F = min(static_cast<myFloat>(1),max(static_cast<myFloat>(0),integral));
            ret = retValue<myFloat>(F, useF, log_p);
            if (verbose)
              cout << "cdf: " << endl << "  Using tail integral, returning " << ret << endl;
            return ret;
          }
        }
      }
      if (alpha != 1 && small_x_m_zet) {
        if (verbose)
          cout << "cdf: Bad theta2 for small x_m_zet.  Using F_zeta" << endl;
        return (log_p) ? log(F) : F;
      } else {
        if (verbose)
          cout << "cdf: Bad theta2 or large integration error for large x_m_zet.  Using pPareto." << endl;
        // Note we're assuming here that a bad theta2 for alpha==1 implies a large x.
        return pPareto(x_m_zeta_input, alpha, beta_input, lower_tail, log_p);
      }
  } // switch on dist_type
} // StandardStableDistribution<myFloat>::cdf
  
template<typename myFloat>
myFloat pPareto(myFloat x, myFloat alpha, myFloat beta, bool lower_tail, bool log_p) {
  bool neg = x < 0;  // left tail
  beta = neg ? -beta : beta;
  if(log_p) {
    return (lower_tail != neg) ?
    log1p(-(1+beta)* C_stable_tail(alpha, false)* pow(fabs(x),-alpha)) :
    log1p(beta)+ C_stable_tail(alpha, true) - alpha*log(fabs(x));
  } else {
    myFloat iF = (1+beta)* C_stable_tail(alpha, false)* pow(fabs(x),-alpha);
    return (lower_tail != neg) ? 1-iF : iF;
  }
} // pPareto
  
} //namespace stable_distribution
