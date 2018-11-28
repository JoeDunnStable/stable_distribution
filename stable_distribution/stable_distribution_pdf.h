/// \file stable_distribution_pdf.h
/// Implementation of pdf of standard stable distribution
/// Included in stable_distribution.h when LIBRARY is defined
/// \author Joseph Dunn
/// \copyright 2016, 2017, 2018 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include <limits>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <boost/math/tools/toms748_solve.hpp>

namespace stable_distribution {

using std::string;
using std::vector;
using std::sort;
using std::setw;
using std::right;
using std::setprecision;
using std::ios;
using std::scientific;
using std::pair;
using std::max;
using std::min;
using boost::math::tools::toms748_solve;

/// x*exp(-x)  numerically stable, with correct limit 0 for x --> Inf
template<typename myFloat>
inline myFloat x_exp_m_x(myFloat x, StandardStableDistribution<myFloat>* ext) {
  myFloat r;
  if(boost::math::isnan(x))
    r = NAN;
  else if(x > StandardStableDistribution<myFloat>::large_exp_arg) // e.g. x == Inf
    r = 0;
  else
    r = x*exp(-x);
  return r;
}

template<typename myFloat>
myFloat dPareto(myFloat x, myFloat alpha, myFloat beta, bool log_flag) {
  if (x < 0) {
    x = -x;
    beta = -beta;
  }
  if (log_flag)
    return log(alpha) + log1p(beta) + C_stable_tail(alpha, true) -
      (1 + alpha) * log(x);
  else
    return alpha * (1 + beta) * C_stable_tail(alpha,log_flag) * pow(x,-(1 +
                  alpha));
}

template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::integrate_pdf(bool log_flag) {
  // --- pdf(x, alpha, beta, ..)  for alpha < 2 ---
  myFloat r;
  if (verbose) cout << "integrate_pdf:" << endl
                    << "Integrand is g * exp(-g)" << endl;
  
  Integral_f_of_g<myFloat> int_g1(&x_exp_m_x, this);
  r=int_g1();
  abserr=c2*int_g1.abserr;
  c_g_theta2_error = c2 * g_theta2_error;
  neval+=int_g1.neval;
  termination_code=int_g1.termination_code;
  last=int_g1.last;
  if (verbose) {
cout << "  c2*sum(r)= " << fmt << c2 << " * " << fmt << r << " = " << fmt << c2*(r) << endl
         << "  abs.err = " << fmt << c2*int_g1.abserr << endl
         << "  msg = " << fmt << int_g1.msg() << endl;
  }
  if (log_flag)
    return log(c2) + log(r);
  else ;
  return c2 * (r);
} //integrate_pdf

template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::pdf(myFloat x, int log_flag, Parameterization pm)
{
  //Default values which will be reset if the integrator is used
  abserr = 0;
  neval = 0;
  termination_code = IntegrationController<myFloat>::TerminationCode::normal;
  last = 0;
  
  set_x_m_zeta(x, pm);
  myFloat ret;
  if (verbose)
    cout << "pdf: log_flag = " << log_flag << endl << *this;
  if (!boost::math::isfinite(x)) {
    return log_flag ? NegInf: 0;
  }
  switch (dist_type) {
    case Cauchy :
      ret =log_flag ? static_cast<myFloat>(-log(pi) - log(1 + x*x))
                    : static_cast<myFloat>(1 / (pi * (1 + x*x)));
      if (verbose)
        cout << "  Using Cauchy Distribution = " << fmt << ret << endl;
      return ret;
    case normal :
      ret = log_flag ? static_cast<myFloat>(-x*x/4 -log(static_cast<myFloat>(2)) -log(pi)/2)
                     : exp(-x*x/4)/(2*sqrt(pi));
      if (verbose)
        cout << "  Using Normal Distribution = " << fmt << ret << endl;
      return ret;
    case fin_support :
      ret = log_flag ? NegInf : 0;
      if (verbose)
        cout << "  Outside of support.  Returning " << fmt << ret << endl;
      return ret;
    case other :
      myFloat ret;
      
      if (verbose)
        cout << "  General case" << endl;
      if (use_series_small_x) {
        ret = series_small_x_pdf(x, pm);
        ret = (log_flag) ? log(ret) : ret;
        abserr = fabs(error_series);
        if (verbose)
          cout << "  x is very small.  Using series = " << fmt << ret << endl;
        return ret;
      }
      if (use_series_large_x) {
        ret = series_large_x_pdf(x, pm);
        ret = (log_flag) ? log(ret) : ret;
        abserr = fabs(error_series);
        if (verbose)
          cout << "  x is very large.  Using series = " << fmt << ret << endl;
        return ret;
      }
/*
       if (alpha != 1) { // 0 < alpha < 2	&  |beta| <= 1 from above

        // For  x = zeta, have special case formula [Nolan(1997)];
        if (alpha < 1 && fabs(beta)==1)
          f_zeta = log_flag ? NegInf : 0;
        else
          f_zeta = (log_flag ? static_cast<myFloat>(lgamma(1 + 1/alpha) + log(cos(theta0)) - (log(pi) + log1p(pow(zeta,2))/(2 * alpha)))
                    : static_cast<myFloat>(tgamma(1 + 1/alpha) * cos(theta0)/(pi * pow(1 + pow(zeta,2),(1/(2 * alpha))))));

        // need to use it also for x ~= zeta
        if (use_f_zeta) {
          ret = f_zeta;
          if (verbose)
            cout << "  " << fmt << x_input
            << " ~= " << fmt << zeta
            << " Using f_zeta() = " << fmt << ret << endl;
          return ret;
        }
      } // alpha != 1
*/
      bool avoid_series = small_x_m_zet ? avoid_series_small_x : avoid_series_large_x;
      if (good_theta2 || avoid_series) {
        ret = integrate_pdf(log_flag);
        if (verbose) cout << "pdf:" << endl;
        myFloat error = max(c_g_theta2_error, abserr);
        if (avoid_series) {
          if (verbose)
            cout << "  Series is unreliable. Using integral = " << fmt << ret << endl;
        } else if (error < controllers.controller.epsrel * fabs(ret)) {
          if (verbose)
            cout << "  Integration error is below threshhold. Using integral = " << fmt << ret << endl;
        } else {
          myFloat result_series = (small_x_m_zet) ? series_small_x_pdf(x, pm) : series_large_x_pdf(x, pm);
          if (fabs(error_series) > abserr ) {
            if (verbose)
              cout << "  Integration error is above threshhold but better than series. Using integral = " << fmt << ret << endl;
          } else {
            ret = (log_flag) ? log(result_series) : result_series;
            abserr = fabs(error_series);
            if (verbose)
               cout << "  Integration error is above threshhold and worse than series.  Using series = " << fmt << ret << endl;
          }
        }
      } else { // bad theta2 and not avoiding series
        abserr = NAN;
        termination_code = IntegrationController<myFloat>::bad_integrand;
        last = 0;
        ret = (small_x_m_zet) ? series_small_x_pdf(x, pm) : series_large_x_pdf(x, pm);
        ret = (log_flag) ? log(ret) : ret;
        abserr = fabs(error_series);
        if (verbose)
          cout << "  Bad theta2.  Using series = " << fmt << ret << endl;
      } // bad theta2
      return ret;
  } // switch on dist_type

} // pdf
  
  using boost::math::binomial_coefficient;
  using boost::math::factorial;
  using std::endl;
  
/// pdf using zolotarev expansion for large x
template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::series_large_x_pdf(myFloat x0, Parameterization pm) {
  if (verbose)
    cout << endl << "series_large_x_pdf:" << endl;
  myFloat result_series;
  myFloat Machine_eps = std::numeric_limits<myFloat>::epsilon();
  if (alpha < 1) {
    // Zolotarev 2.4.8, a series that converges for large x
    
    if (rho == 0) {
      result_series = 0;
      error_series = 0;
      n_series = 0;
    } else if (xB==0) {
      result_series = std::numeric_limits<myFloat>::quiet_NaN();
      error_series = std::numeric_limits<myFloat>::infinity();
      n_series = 0;
    } else {
      myFloat term{0}, abs_term{0};
      result_series = 0;
      error_series = 0;
      for (int k=1; k<=max_n; ++k) {
        abs_term = tgamma_ratio(k*alpha+1,myFloat(k+1)) * pow(xB,-k*alpha-1)/(pi*gammaB);
        term = abs_term * (0 == k%2 ? -1 : 1) * sin(pi*k*rho*alpha);
        myFloat error = fabs(term) * ((k*alpha+1)*digamma(k*alpha+1) +k*alpha +1) * Machine_eps;
        if (k==1) { // check against dpareto
          myFloat dPareto = tgamma(alpha)/pi*sin(pi*alpha/2);
          dPareto *= alpha * (1+beta) * pow(x_m_zet,-alpha-1);
          if (verbose>0)
            cout << "1 - term1/dPareto = " << fmt << (1 -term/dPareto) << endl;
        }
        if (verbose>1)
          cout << "n = " << k << ", term = " << fmt << term << endl;
        result_series += term;
        error_series += error;
        n_series=k;
        if (!boost::math::isfinite(abs_term) || abs_term < 2*Machine_eps*fabs(result_series)) break;
      }
      error_series += abs_term;
    } // rho != 0
    
  } else if (alpha == 1) {
    
    // Zolotarev 2.5.23 is the asymptotic series for large x
    // For beta != -1,0
    // Unlike the convergent formula this one requires xB > 0
    
    if (xB == 0) {
      result_series = std::numeric_limits<myFloat>::quiet_NaN();
      error_series = std::numeric_limits<myFloat>::max();
      n_series = 0;
    } else if (betaB == -1) {
      // Zolotarev Theorem 2.5.2 asymptotic for x -> infinity
      myFloat xi = exp(xB-1);
      myFloat nu = 1;
      myFloat exp_m_xi = exp(-xi);
      myFloat fac = exp_m_xi !=0 ? nu*pow(xi,(2-alpha)/(2*alpha))*exp_m_xi/sqrt(2*pi*alpha)/gammaB
      : 0;
      if (verbose > 1)
        cout << "xi = " << fmt << xi << endl
        << "nu  = " << fmt << nu << endl
        << "fac = " << fmt << fac << endl;
      myFloat term = fac;
      if (verbose > 1)
        cout << "n = " << 0 << ", term = " << fmt << term << endl;
      myFloat old_term = term;
      result_series = term;
      myFloat alpha_star = 1/alpha;
      myFloat alpha_xi_n = 1;
      for (int n = 1; n<Q_pdf.size(); ++n) {
        alpha_xi_n /= (alpha_star*xi);
        term = fac * Q_pdf.at(n) * alpha_xi_n;
        if (fabs(term) > fabs(old_term)) break;
        if (verbose > 1)
          cout << "n = " << n << ", term = " << fmt << term << endl;
        n_series = n;
        result_series += term;
        old_term = term;
        if (fabs(term) <= eps*fabs(result_series) && fabs(alpha_xi_n) <= eps) break;
      }
      error_series = fabs(term)
/*      +fabs(result_series)*exp(-pow(xi,.25)) */
      ;
    } else {
      myFloat log_x=log(xB);
      result_series = 0;
      myFloat term{std::numeric_limits<myFloat>::quiet_NaN()};
      myFloat old_term{std::numeric_limits<myFloat>::quiet_NaN()};
      int num_small_terms=0;
      for (int n=1; n<=max_n; ++n) {
        myFloat fac{0};
        for (int l=0; l<=n; ++l){
          myFloat r_l_n{0};
          for (int m=l; m<=n; ++m) {
            myFloat term0=binomial_coefficient<myFloat>(n,m)*binomial_coefficient<myFloat>(m,l);
            term0 *= (1-2*((m-l)%2)) * gamma_at_integers(m-l,1+n);
            term0 *=  pow(betaB, m) * pow(pi/2*(betaB_p_1),n-m) * sin(pi/2*(n-m));
            r_l_n += term0;
            if (verbose > 2){
              if (m == l)
                cout << "r(" << l << ", " << n << "} = " << fmt << term0;
              else
                cout << " + " << fmt << term0;
            }
          }
          if (verbose > 2)
            cout << endl;
          if (verbose>1)
            cout << "r(" << l << ", " << n << ") * log_x ^ l = "
            << fmt << r_l_n << " * " << fmt << log_x << " ^ " << l
            << " = " << fmt << r_l_n << " * " << fmt << pow(log_x,l)
            << " = " << fmt << r_l_n * pow(log_x,l) << endl;
          fac += r_l_n*pow(log_x,l);
          if (verbose>1)
            cout << "fac: " << fmt << fac << endl;
        }
        term=fac * pow(xB,-n-1) / (factorial<myFloat>(n)*pi*gammaB);
        if (n==1) { // check against dpareto
          myFloat dPareto = tgamma(alpha)/pi*sin(pi*alpha/2);
          dPareto *= alpha * (betaB_p_1) * pow(fabs(x_m_zet),-alpha-1);
          if (verbose>0)
            cout << "1 - term1/dPareto = " << fmt << (1 -term/dPareto) << endl;
        }
        if (verbose>1)
          cout << "n = " << n << ": term = " << term << endl;
        if (n>1 && fabs(term) > fabs(old_term)) break;
        n_series = n;
        result_series += term;
        if (verbose>0)
          cout << "n = " << n
          << ": result_series = " << fmt << result_series << endl;
        if (fabs(term) < fabs(result_series) * Machine_eps) {
          if (num_small_terms > 0) break;
          num_small_terms=num_small_terms+1;
        } else {
          num_small_terms=0;
          old_term=term;
        }
      }
      error_series=fabs(term) + fabs(result_series) * Machine_eps;
    } // xB != 0
    
  } else {
    // alpha > 1
    
    if (xB == 0) {
      result_series = std::numeric_limits<myFloat>::quiet_NaN();
      error_series = std::numeric_limits<myFloat>::max();
      n_series = 0;
    } else if (beta == -1) {
      // Zolotarev Theorem 2.5.2 asymptotic for x -> infinity
      myFloat xi = fabs(1-alpha) * pow(xB/alpha, alpha/(alpha-1));
      myFloat nu = pow(abs(1-alpha),-1/alpha);
      myFloat exp_m_xi = exp(-xi);
      myFloat fac = exp_m_xi !=0 ? nu*pow(xi,(2-alpha)/(2*alpha))*exp(-xi)/sqrt(2*pi*alpha)/gammaB
      :0;
      if (verbose > 1)
        cout << "xi = " << fmt << xi << endl
        << "nu  = " << fmt << nu << endl
        << "fac = " << fmt << fac << endl;
      myFloat term = fac;
      myFloat old_term = term;
      if (verbose > 1)
        cout << "n = " << 0 << ", term = " << fmt << term << endl;
      result_series = term;
      myFloat alpha_star = 1/alpha;
      myFloat alpha_xi_n = 1;
      for (int n = 1; n<Q_pdf.size(); ++n) {
        alpha_xi_n /= (alpha_star*xi);
        term = fac * Q_pdf.at(n) * alpha_xi_n;
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
      
    } else {
      // Zolotarev 2.5.4, an asymptotic series that works well for large x
      n_series = 1;
      // myFloat alpha_star = 1/alpha;
      // myFloat beta_star = 1-(2-alpha)*(betaB_p_1);
      myFloat abs_term =(alpha/(gammaB*pi))*tgamma_ratio(alpha,myFloat(1))*pow(xB,-alpha-1);
      myFloat term=abs_term*sin(pi/2*(2-alpha)*(betaB_p_1));
      myFloat abs_tail{0};
      result_series = 0;
      // myFloat cos_ab = pow(cos((pi/2)*alpha_star*beta_star),1/alpha_star);
      myFloat cos_ab = 1;
      myFloat cos_ab_m_k = 1/cos_ab;
      error_series = std::numeric_limits<myFloat>::infinity();
      myFloat dPareto = tgamma(alpha)/pi*sin(pi*alpha/2);
      dPareto *= alpha * (1+beta) * pow(x_m_zet,-alpha-1);
      if (verbose>0) {
        cout << "1 - term1/dPareto = " << fmt << (1 -term/dPareto) << endl;
      }
      for (int k=2; k<=max_n; ++k) {
        abs_term=(alpha/(gammaB*pi))*tgamma_ratio(k*alpha,myFloat(k))* pow(xB,-alpha*k-1);
        myFloat new_term = abs_term*sin(pi*k/2*(2-alpha)*(betaB_p_1));
        cos_ab_m_k /= cos_ab;
        myFloat new_error = abs_term * cos_ab_m_k;
        if (new_error > error_series) break;
        n_series=k-1;
        result_series += term;
        abs_tail += fabs(term);
        error_series = new_error;
        if (verbose>0)
          cout << "n = " << k-1 << ", term = " << fmt << term 
		  << ", error = " << fmt << error_series << endl;
        term = new_term;
        if (abs_term < Machine_eps * result_series) break;
      }
      error_series +=  abs_tail*Machine_eps;
    } // xB != 0
  } // alpha > 1
  if (verbose>0)
    cout << "result_series = " << fmt << result_series
    << ", error_series = " << fmt << error_series << endl
    << ", n_series = " << n_series << endl;
  
  return result_series;
}
  
/// pdf using zolotarev expansion for small x
template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::series_small_x_pdf(myFloat x0, Parameterization pm) {
  if (verbose)
    cout << endl << "series_small_x_pdf:" << endl;
  myFloat result_series;
  myFloat Machine_eps = std::numeric_limits<myFloat>::epsilon();
  if (alpha < 1) {
    
    if (beta == 1) {
      // Zolotarev Theorem 2.5.2, asymptotic for small x
      myFloat xi = fabs(1-alpha) * pow(xB/alpha, alpha/(alpha-1));
      myFloat nu = pow(abs(1-alpha),-1/alpha);
      myFloat exp_m_xi = exp(-xi);
      myFloat fac = exp_m_xi !=0 ?nu*pow(xi,(2-alpha)/(2*alpha))*exp_m_xi/sqrt(2*pi*alpha)/gammaB
      :0;
      if (verbose > 1)
        cout << "xi = " << fmt << xi << endl
        << "nu  = " << fmt << nu << endl
        << "fac = " << fmt << fac << endl;
      myFloat term = fac;
      if (verbose > 1)
        cout << "n = " << 0 << ", term = " << fmt << term << endl;
      myFloat old_term = term;
      result_series = term;
      myFloat alpha_star = alpha;
      myFloat alpha_xi_n = 1;
      for (int n = 1; n<Q_pdf.size(); ++n) {
        alpha_xi_n /= (alpha_star*xi);
        term = fac * Q_pdf.at(n) * alpha_xi_n;
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
      
    } else { // beta != 1
      // Zolotarev formula 2.5.1, which is asmyptotic as x->0
      // for betaB != 1
      myFloat abs_term = tgamma(1/alpha)/(pi * alpha * gammaB);
      myFloat term = abs_term * ((beta!=-1)?sin(pi/2*one_m_betaB):0);
      myFloat abs_tail{fabs(term)};
      if (xB == 0 || beta== -1) {
        result_series = term;
        error_series = abs_tail * Machine_eps;
        n_series =0;
        if (verbose>0)
          cout << "n = " << n_series << ", term = " << term << endl;
      } else { // xB !=0
        result_series = 0;
        // myFloat cos_ab = pow(cos((pi/2)*alpha*betaB),1/alpha);
        myFloat cos_ab = 1;
        myFloat cos_ab_m_k = 1/cos_ab;
        error_series = std::numeric_limits<myFloat>::infinity();
        for (int k=1; k<=max_n; ++k) {
          // calculate the next term which is used to estimate the error
          abs_term = tgamma_ratio((k+1)/alpha,myFloat(k+1))*pow(xB,k)/(pi * alpha * gammaB);
          myFloat new_term = abs_term * sin(pi/2 * (k+1) * one_m_betaB);
          cos_ab_m_k /= cos_ab;
          myFloat new_error = abs_term * cos_ab_m_k;
          if (k>1 && (!boost::math::isfinite(new_error) || new_error > error_series)) break;
          n_series = k-1;
          result_series += term;
          abs_tail += fabs(term);
          error_series=new_error;
          if (verbose>0)
            cout << "n = " << k-1 << ", term = " << term  << ", error = " << error_series << endl;
          term = new_term;
          if (abs_term < Machine_eps * result_series) break;
        } // k loop
        error_series += abs_tail*Machine_eps;
      } // xB != 0 && beta != -1
    }
    
  } else if (alpha == 1) {
    // No useable series
    n_series = 0;
    result_series = std::numeric_limits<myFloat>::quiet_NaN();
    error_series = std::numeric_limits<myFloat>::max();

  } else {
    // alpha > 1
    // Zolotarev Formula 2.4.6 is a convergent series that works well for x <=1, but
    //    but suffers from extreme rounding errors for larger x.
    
    result_series = 0;
    error_series =0;
    myFloat term, abs_term{std::numeric_limits<myFloat>::quiet_NaN()};
    for (int n=1; n<=max_n; ++n) {
      abs_term = tgamma_ratio(n/alpha+1, myFloat(n+1))*pow(xB,n-1)/(pi*gammaB);
      term = abs_term * (0==n%2 ? -1 : 1) * sin(pi*n*rho);
      myFloat error = fabs(term)*((n/alpha+1)*digamma(n/alpha+1)+n-1) * Machine_eps;
      if (verbose>1)
        cout << "n = " << n << ", term = " << fmt << term << endl;
      n_series = n;
      result_series += term;
      error_series += error;
      if (!boost::math::isfinite(abs_term) || abs_term < 2*Machine_eps*fabs(result_series)) break;
    }
    error_series += abs_term;
    
  } // alpha > 1
  if (verbose>0)
    cout << "result_series = " << fmt << result_series
    << ", error_series = " << fmt << error_series << endl
    << ", n_series = " << n_series << endl;
  
  return result_series;
}

} // namespace stable_distribution
