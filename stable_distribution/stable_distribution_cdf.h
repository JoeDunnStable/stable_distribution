/// \file stable_distribution_cdf.h
/// Implementation of cdf for standard stable distribution
/// Included in stable_distribution.h when LIBRARY is defined
/// \author Joseph Dunn
/// \copyright 2016, 2017, 2018 Joseph Dunn
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
      F = my_erfc(fabs(x/2))/2;
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
      if (use_series_small_x) {
        ret = series_small_x_cdf(x, lower_tail, pm);
        ret = (log_p) ? log(ret) : ret;
        abserr = fabs(error_series);
        if (verbose)
          cout << "  x is very small.  Using series = " << fmt << ret << endl;
        return ret;
      } else if (use_series_large_x) {
        ret = series_large_x_cdf(x, lower_tail, pm);
        ret = (log_p) ? log(ret) : ret;
        abserr = fabs(error_series);
        if (verbose)
          cout << "  x is very large.  Using series = " << fmt << ret << endl;
        return ret;
      }
      myFloat F=0;
      if (alpha !=1 && small_x_m_zet) {
          myFloat F_zeta = (alpha<1 && fabs(beta)==1)
                           ? (lower_tail != (beta_input<1)) ? 0 : 1
                           :(lower_tail) ? static_cast<myFloat>(static_cast<myFloat>(.5) - theta0_x_gt_zeta/pi)
                                         : static_cast<myFloat>(static_cast<myFloat>(.5) + theta0_x_gt_zeta/pi);
          if (verbose)
            cout << "cdf: Using difference from F_zeta, " << F_zeta << endl;
          F = F_zeta;
	 /*
        if (use_f_zeta) {
          if (verbose)
            cout << "cdf:: Using F_zeta" << endl;
          ret = (log_p) ? log(F) : F;
          return ret;
        } else {
          if (verbose)
            cout << "cdf: Using difference from F_zeta, " << F_zeta << endl;
        }
	*/
      }
      bool avoid_series = small_x_m_zet ? avoid_series_small_x : avoid_series_large_x;
      if (good_theta2 || avoid_series) {
        myFloat integral = integrate_cdf();
        myFloat error = abserr;
        if ( avoid_series || error < controllers.controller.epsrel* fabs(integral)) {
          if (verbose)
            cout << "cdf: " << (avoid_series?"Series is unreiable" : "Integration error is small") << ".  Using integral" << endl;
          if (alpha != 1 && small_x_m_zet) {
            F -= ((lower_tail != (x_m_zeta_input > 0)) ? 1 : -1) * integral;
            ret = (log_p) ? log(F) : F;
            if (verbose)
              cout << "cdf: " << endl
              << "integral is delta with F_zeta returning: " << ret << endl;
            return ret;
          } else {
            bool useF = !((x >= zeta && lower_tail) || (x < zeta && !lower_tail));
            F = min(static_cast<myFloat>(1),max(static_cast<myFloat>(0),integral));
            ret = retValue<myFloat>(F, useF, log_p);
            if (verbose)
              cout << "cdf: " << endl << "  Using tail integral, returning " << ret << endl;
            return ret;
          }
        } else { // We'll use the series if its better than  integration
          myFloat F_series = (alpha != 1 && small_x_m_zet) ? series_small_x_cdf(x, lower_tail, pm)
                                                         : series_large_x_cdf(x, lower_tail, pm);
          if (error_series < abserr) {
            abserr = fabs(error_series);
            if (verbose)
              cout << "cdf: Integration error is above threshhold and worse than series" << endl
              << "Using series for " << ((small_x_m_zet) ? "small x":"large x") << "." << endl;
            return (log_p) ? log(F_series) : F_series;
          } else {
            if (verbose)
              cout << "cdf:  Integration error is above threshhold but better that series.   Using integral" << endl;
            if (alpha != 1 && small_x_m_zet) {
              F -= ((lower_tail != (x_m_zeta_input > 0)) ? 1 : -1) * integral;
              ret = (log_p) ? log(F) : F;
              if (verbose)
                cout << "cdf: " << endl
                << "integral is delta with F_zeta returning: " << ret << endl;
              return ret;
            } else {
              bool useF = !((x >= zeta && lower_tail) || (x < zeta && !lower_tail));
              F = min(static_cast<myFloat>(1),max(static_cast<myFloat>(0),integral));
              ret = retValue<myFloat>(F, useF, log_p);
              if (verbose)
                cout << "cdf: " << endl << "  Using tail integral, returning " << ret << endl;
              return ret;
            }
          }  // error series is worse than integration error
        }
      } else {
        // bad_theta2 and we're not avoiding the series.
        if (verbose)
          cout << "cdf: Bad theta2.  Using series for " << ((alpha != 1 && small_x_m_zet) ? "small" : "large") << " x."  << endl;
        F = (alpha != 1 && small_x_m_zet) ? series_small_x_cdf(x, lower_tail, pm)
                                          : series_large_x_cdf(x, lower_tail, pm);
        abserr = fabs(error_series);
        return (log_p) ? log(F) : F;
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
  
/// cdf using zolotarev expansion for large x
template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::series_large_x_cdf(myFloat x0, int lower_tail, Parameterization pm) {
  if (verbose>0)
    cout << endl << "series_large_x_cdf: " << endl;
  myFloat Machine_eps = std::numeric_limits<myFloat>::epsilon();
  myFloat result_series;
  if (alpha < 1) {
    // Zolotarev formula 2.4.3 for upper tail, convergent for large x
    if (xB == 0) {
      result_series= std::numeric_limits<myFloat>::quiet_NaN();
      error_series = std::numeric_limits<myFloat>::infinity();
      n_series = 0;
    } else {
      result_series = 0;
      error_series = 0;
      myFloat abs_series = 0;
      myFloat term{std::numeric_limits<myFloat>::quiet_NaN()};
      myFloat abs_term{std::numeric_limits<myFloat>::quiet_NaN()};
      for (int n=1; n<=max_n; ++n) {
        abs_term = tgamma_ratio(n*alpha+1,myFloat(n+1)) * pow(xB,-n*alpha)/(n*alpha*pi);
        term = abs_term * (0 == n%2 ? -1 : 1) * sin(pi*n*rho*alpha);
        myFloat error = fabs(term) * ((n*alpha+1)*digamma(n*alpha+1) + n*alpha) * Machine_eps;
        if (verbose>1)
          cout << "n = " << n << ", term = " << fmt << term << endl;
        n_series = n;
        result_series += term;
        error_series += error;
        abs_series += fabs(term);
        if (!boost::math::isfinite(abs_term) || abs_term < 2*Machine_eps*fabs(result_series)) break;
      }
      result_series = (lower_tail == positive_x) ? 1-result_series : result_series;
      error_series += abs_term;
    } //xB !=0
    
  } else if (alpha == 1) {
    
    // Zolotarev 2.5.23 is asymptotic for large x
    // unlike the convergent formula this formula requires x > 0
    if (xB == 0) {
      result_series = std::numeric_limits<myFloat>::quiet_NaN();
      error_series = std::numeric_limits<myFloat>::max();
      n_series= 0;
    } else if (betaB == -1) {
      // Zolotarev formula 2.5.20
      myFloat xi = exp(xB-1);
      myFloat fac = exp(-xi)/sqrt(2*pi*xi);
      if (verbose > 1)
        cout << "xi = " << fmt << xi << endl
        << "fac = " << fmt << fac << endl;
      myFloat term = fac;
      if (verbose > 1)
        cout << "n = " << 0 << ", term = " << fmt << term << endl;
      result_series = term;
      myFloat old_term=term;
      myFloat alpha_xi_n = 1;
      n_series = 0;
      for (int n=1; n<Q_cdf.size(); ++n) {
        alpha_xi_n /= xi;
        term = fac * Q_cdf.at(n) * alpha_xi_n;
        if (fabs(term) > fabs(old_term)) {
          break;
        }
        n_series = n;
        if (verbose > 1)
          cout << "n = " << n << ", term = " << fmt << term << endl;
        result_series += term;
        if (fabs(term) <= eps*fabs(result_series) && fabs(alpha_xi_n) <= eps) break;
      }
      error_series = fabs(term)
/*      + fabs(result_series-fac)*exp(-pow(xi,.25)) */
      ;
      result_series = (lower_tail == positive_xB) ? 1-result_series : result_series;
    } else { // xB > 0 & betaB_ != -1
      myFloat log_x=log(xB);
      result_series = 0;
      myFloat term{std::numeric_limits<myFloat>::quiet_NaN()};
      myFloat old_term{std::numeric_limits<myFloat>::quiet_NaN()};
      int num_small_terms=0;
      for (int n=1; n<=max_n; ++n) {
        n_series=n;
        myFloat fac{0};
        for (int l=0; l<=n; ++l){
          myFloat r_l_n{0};
          for (int m=l; m<=n; ++m) {
            myFloat term0=binomial_coefficient<myFloat>(n,m)*binomial_coefficient<myFloat>(m,l);
            term0 *= (1-2*((m-l)%2)) * gamma_at_integers(m-l,n);
            term0 *=  pow(betaB, m) * pow(pi/2*(betaB_p_1),n-m) * sin(pi/2*(n-m));
            r_l_n += term0;
            if (verbose > 2){
              if (m == l)
                cout << "r(" << l << ", " << n << "} = " << fmt << term0;
              else
                cout << " + " << fmt << term0;
            }
          }
          if (verbose > 2) cout << endl;
          if (verbose>1)
            cout << "r(" << l << ", " << n << ") * log_x ^ l = "
            << fmt << r_l_n << " * " << fmt << log_x << " ^ " << l
            << " = " << fmt << r_l_n << " * " << fmt << pow(log_x,l)
            << " = " << fmt << r_l_n * pow(log_x,l) << endl;
          fac += r_l_n*pow(log_x,l);
          if (verbose>1)
            cout << "fac: " << fmt << fac << endl;
        }
        term=fac * pow(xB,-n) / (factorial<myFloat>(n)*pi);
        if (verbose>1)
          cout << "n = "<< n << ", term = " << fmt << term << endl;
        if (n>1 && fabs(term) > fabs(old_term)) break;
        result_series += term;
        if (verbose>0)
          cout << "n = " << n
          << ", result_series = " << fmt << result_series << endl;
        if (fabs(term) < fabs(result_series) * Machine_eps) {
          if (num_small_terms > 0) break;
          num_small_terms=num_small_terms+1;
        } else {
          num_small_terms=0;
          old_term=term;
        }
      }
      result_series = (lower_tail == positive_xB) ? 1-result_series : result_series;
      error_series=fabs(term) + fabs(result_series) * Machine_eps;
    } // xB !=0
    
  } else {
    // alpha > 1
    
    // Zolotarev Formula 2.5.4 is an asymptotic series for large x
    // for beta != -1
    if (xB == 0) {
      result_series = std::numeric_limits<myFloat>::quiet_NaN();
      error_series = std::numeric_limits<myFloat>::max();
      n_series = 0;
    } else if (beta == -1) {
      myFloat xi = fabs(alpha_m_1) * pow(xB/alpha, alpha/(alpha_m_1));
      myFloat fac = exp(-xi)/sqrt(2*pi*alpha*xi);
      if (verbose > 1)
        cout << "xi = " << fmt << xi << endl
        << "fac = " << fmt << fac << endl;
      myFloat term = fac;
      myFloat old_term = term;
      n_series=0;
      if (verbose > 1)
        cout << "n = " << 0 << ", term = " << fmt << term << endl;
      result_series = term;
      myFloat alpha_star = 1/alpha;
      myFloat alpha_xi_n = 1;
      for (int n = 1; n<Q_cdf.size(); ++n) {
        alpha_xi_n /= (alpha_star*xi);
        term = fac * Q_cdf.at(n) * alpha_xi_n;
        if (fabs(term) > fabs(old_term)) break;
        if (verbose > 1)
          cout << "n = " << n << ", term = " << fmt << term << endl;
        n_series = n;
        result_series += term;
        old_term = term;
        if (fabs(term) <= eps*fabs(result_series) && fabs(alpha_xi_n) <= eps) break;
      }
      error_series = fabs(term)
      /*        +fabs(result_series)*exp(-pow(xi,.25)) */
      ;
      result_series = (lower_tail==positive_x) ? 1 - result_series : result_series;
    } else {
      myFloat abs_term =(1/(pi))*tgamma(alpha)/tgamma(static_cast<myFloat>(2))*pow(xB,-alpha);
      myFloat term=abs_term*sin(pi/2*(2-alpha)*(betaB_p_1));
      result_series=term;
      myFloat old_abs_term=abs_term;
      n_series =1;
      if (verbose > 1)
        cout << "n = " << n_series
        << ", term = " << fmt << term << endl;
      for (int k=2; k<=max_n; ++k) {
        abs_term=(1/(pi))*tgamma_ratio(k*alpha,myFloat(k+1))* pow(xB,-alpha*k);
        term=abs_term*sin(pi*k/2*(2-alpha)*(betaB_p_1));
        if (abs_term > old_abs_term) break;
        n_series = k;
        if (verbose>1)
          cout <<"n = " << k << ", term = " << fmt << term << endl;
        result_series += term;
        if (abs_term < Machine_eps * fabs(result_series)) break;
        old_abs_term=abs_term;
      }
      error_series = abs_term + fabs(result_series)*Machine_eps;
      result_series = (lower_tail==positive_x) ? 1-result_series : result_series;
    } // xB != 0
    
  }  // alpha > 1
  if (verbose>0)
    cout <<"result_series = " << fmt << result_series
    << ", error_series = " << fmt << error_series
    << ", n_series = " << n_series << endl;

  return result_series;
} // series_large_x_cdf

/// cdf using zolotarev expansion for small x_m_zeta
template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::series_small_x_cdf(myFloat x0, int lower_tail, Parameterization pm) {
  if (verbose>0)
    cout << endl << "series_small_x_cdf: " << endl;
  myFloat result_series;
  myFloat Machine_eps = std::numeric_limits<myFloat>::epsilon();
  if (alpha < 1) {
    if (beta == 1) {
      // Zolotarev Theorem 2.5.3, asymptotic for small x
      myFloat xi = fabs(alpha_m_1) * pow(xB/alpha, alpha/(alpha_m_1));
      myFloat fac = exp(-xi)/sqrt(2*pi*alpha*xi);
      if (verbose > 1)
        cout << "xi = " << fmt << xi << endl
        << "fac = " << fmt << fac << endl;
      myFloat term = fac;
      if (verbose > 1)
        cout << "n = " << 0 << ", term = " << fmt << term << endl;
      myFloat old_term = term;
      result_series = term;
      myFloat alpha_star = alpha;
      myFloat alpha_xi_n = 1;
      for (int n = 1; n<Q_cdf.size(); ++n) {
        alpha_xi_n /= (alpha_star*xi);
        term = fac * Q_cdf.at(n) * alpha_xi_n;
        if (fabs(term) > fabs(old_term)) break;
        if (verbose > 1)
          cout << "n = " << n << ", term = " << fmt << term << endl;
        n_series = n;
        result_series += term;
        old_term = term;
        if (fabs(term) <= eps*fabs(result_series) && fabs(alpha_xi_n) <= eps) break;
      }
      error_series = fabs(term)
      /*        +fabs(result_series)*exp(-pow(xi,.25)) */
      ;
      result_series = (lower_tail == positive_x) ? result_series : 1 - result_series;
      
    } else { // beta != 1
      // Zolotarev Formula 2.5.3 asymptotic for x -> 0 adjusted to give upper tail
      n_series = 0;
      myFloat term = (beta==-1)?0:.5 * (1 + betaB);
      myFloat abs_term = 0;
      result_series = term;
      myFloat old_abs_term = fabs(term);
      myFloat abs_tail = fabs(term);
      if (verbose > 1)
        cout << "n = " << n_series << ", term = " << fmt << result_series << endl;
      if (xB != 0 && beta != -1) {
        for (int k=1; k<=max_n; ++k) {
          abs_term = tgamma_ratio(k/alpha,myFloat(k+1))* pow(xB,k)/(pi * alpha);
          term = -abs_term * sin(pi/2 * k * one_m_betaB);
          if (k > 1 && abs_term > old_abs_term) break;
          if (verbose>1)
            cout << "n = " << k << ", term = " << fmt << term << endl;
          n_series = k;
          result_series += term;
          abs_tail += fabs(term);
          if (abs_term < Machine_eps * result_series) break;
          old_abs_term=abs_term;
        }
      }
      result_series = (lower_tail==positive_x) ? 1 - result_series : result_series;
      error_series = ((betaB!=-1)?abs_term:fabs(term)) + abs_tail*Machine_eps;
    }
  } else if (alpha == 1) {
    n_series = 0;
    result_series = std::numeric_limits<myFloat>::quiet_NaN();
    error_series = std::numeric_limits<myFloat>::max();
    
  } else {
    // alpha > 1
    // Zolotarev 2.4.6
    
    //    Formula 2.4.3 is a convergent series that works well for x <=1, but
    //    but suffers from extreme rounding errors for larger x.
    n_series = 0;
    myFloat term{.5*(1+theta)};
    myFloat abs_term{fabs(term)};
    result_series = term;
    error_series = 0;
    if (verbose > 1)
      cout << "n = " << n_series << ", term = " << fmt << term << endl;
    myFloat abs_series = fabs(term);
    for (int n=1; n<=max_n; ++n) {
      abs_term = tgamma_ratio(n/alpha+1,myFloat(n+1))*pow(xB,n)/(n * pi);
      term = abs_term * (1==n%2 ? -1 : 1) * sin(pi*n*rho);
      myFloat error = fabs(term) * ((n/alpha+1)*digamma(n/alpha+1) + n)*Machine_eps;
      if (verbose>1)
        cout << "n = " << n << ", term = " << fmt << term << endl;
      n_series = n;
      result_series += term;
      error_series += error;
      abs_series += fabs(term);
      if (!boost::math::isfinite(abs_term) || abs_term < 2*Machine_eps*fabs(result_series)) break;
    }
    result_series = (lower_tail==positive_x) ? 1-result_series : result_series;
    error_series += abs_term;
    
  }  // alpha > 1
  if (verbose > 0)
    cout << "result_series = " << fmt << result_series
    << ", error_series = " << fmt << error_series
    << ", n_series = " << n_series << endl;
  return result_series;
} // series_small_x_cdf

} //namespace stable_distribution
