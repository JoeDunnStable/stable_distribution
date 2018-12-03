/// \file stable_distribution_ddx_pdf.h
/// Implementation of derivative of pdf of standard stable distribution.
/// Included in stable_distribution.h when LIBRARY is defined
/// \author Joseph Dunn
/// \copyright 2016, 2017, 2018 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include <limits>
#include <iostream>
#include <iomanip>
#include <boost/math/tools/toms748_solve.hpp>

namespace stable_distribution {
  
using std::string;
using std::array;
using std::vector;
using std::sort;
using std::setw;
using std::right;
using std::setprecision;
using std::ios;
using std::pair;
using boost::math::tools::toms748_solve;

template<typename myFloat>
myFloat x_1_m_alpha_x_exp_m_x(myFloat x, StandardStableDistribution<myFloat>* std_stable_dist) {
  if (x > StandardStableDistribution<myFloat>::large_exp_arg)
    return 0;
  else
    return x * (1 - std_stable_dist->alpha*x) * exp(-x);
}

template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::integrate_ddx_pdf() {
  // --- pdf(x, alpha, beta, ..)  for alpha < 2 ---
  // For  x = zeta, have special case formula [n_gaussolan(1997)];
  // need to use it also for x ~= zeta, i.e., x.m.zet := |x - zeta| < delta :
  myFloat r;
  if (verbose)
    cout << "integrate_ddx_pdf:" << endl
         << "Integrand is g * (1 - alpha * g) * exp(-g)" << endl;
  // integral is zero near mode so need to use absolute tolerance near mode
  myFloat old_epsabs = controllers.controller.epsabs;
  controllers.controller.epsabs=(fabs(x_m_zeta_input) < 1) ? controllers.controller.epsrel/c_ddx : myFloat{0};
  Integral_f_of_g<myFloat> int_g1(x_1_m_alpha_x_exp_m_x,this, &controllers);
  controllers.controller.epsabs = old_epsabs;
  r=int_g1();
  abserr=fabs(c_ddx) * int_g1.abserr;
  c_g_theta2_error = fabs(c_ddx) * g_theta2_error;
  neval+=int_g1.neval;
  termination_code=int_g1.termination_code;
  last=int_g1.last;
  if (verbose) {
    if (verbose)
      cout << "  c_ddx*sum(r)= " << fmt << c_ddx << " * " << fmt << r << " = " << fmt << c_ddx*(r) << endl
      << "  abs.err = " << fmt << c_ddx*int_g1.abserr << endl
      << "  msg = " << fmt << int_g1.msg() << endl;
  }
  return c_ddx * (r);
} //integrate_ddx_pdf


template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::ddx_pdf(myFloat x, Parameterization pm)
{
  cout.precision(20);
  cout.setf(ios::scientific, ios::floatfield);
  // Default values which will be reset in the integrator is used
  abserr = 0;
  neval = 0;
  termination_code = IntegrationController<myFloat>::TerminationCode::normal;
  last = 0;
  
  myFloat ret;
  set_x_m_zeta(x, pm);
  if (verbose)
    cout << "ddx_pdf: pm = " << pm << endl << *this;
  if (!boost::math::isfinite(x)) {
    if (verbose)
      cout << "  x is infinite, returning " << 0 << endl;
     return 0;
  }
  switch (dist_type) {
    case Cauchy :
      ret = (-2*x / (pi *(1 + x*x)))*(1/(1+x*x));
      if (verbose)
        cout << "  Cauchy distribution, returning " << fmt << ret << endl;
      return  ret;
    case normal :
      ret = -x*exp(-x*x/4)/(4*sqrt(pi));
      if (verbose)
        cout << "  Normal distribution, returning " << fmt << ret << endl;
      return  ret;
    case fin_support :
      if (verbose)
        cout << "  Outside of support, returning " << 0 << endl;
      return 0;
    case other :
      // General Case
      if (verbose)
        cout << "  General Case"<< endl;
      if (use_series_small_x) {
        ret = series_small_x_ddx_pdf(x, pm);
        abserr = fabs(error_series);
        if (verbose)
          cout << "  x is very small.  Using series = " << fmt << ret << endl;
        return ret;
      } 
      if (use_series_large_x) {
        ret = series_large_x_ddx_pdf(x, pm);
        abserr = fabs(error_series);
        if (verbose)
          cout << "  x is very large.  Using series = " << fmt << ret << endl;
        return ret;
      }
/*      
      if (alpha != 1) { // 0 < alpha < 2	&  |beta| <= 1 from above

    // For  x = zeta, have special case formula [n_gaussolan(1997)];
        if (alpha < 1 && fabs(beta)==1)
          dfdx_zeta = 0;
        else if (theta0_x_gt_zeta==0)
          dfdx_zeta = 0;
        else
          dfdx_zeta = tgamma(1 + 2/alpha) * sin(2*theta0_x_gt_zeta)
          /(2*pi * pow(1 + pow(zeta,2),(1/(alpha))));

        // need to use it also for x ~= zeta
        if (use_f_zeta) {
          if (verbose)
            cout << "  " << fmt << x_input << " ~= " << fmt << zeta
            << " Using dfdx_zeta()" << endl;
          return dfdx_zeta;
        }
      } // alpha != 1
   */ 
      bool avoid_series = small_x_m_zet ? avoid_series_small_x : avoid_series_large_x;    
      if (good_theta2 || avoid_series) {
        ret = integrate_ddx_pdf();
        if (verbose) cout << "pdf:" << endl;
        myFloat error = max(c_g_theta2_error, abserr);
        if (avoid_series) {
          if (verbose)
            cout << "  Series is unreliable. Using integral = " << ret << endl;
        } else if (error < controllers.controller.epsrel * fabs(ret)) {
          if (verbose)
            cout << "  Integration error is below threshhold. Using integral = " << ret << endl;
        } else {
          myFloat result_series = (small_x_m_zet) ? series_small_x_ddx_pdf(x, pm) : series_large_x_ddx_pdf(x, pm);
          if (fabs(error_series) > abserr ) {
            if (verbose)
              cout << "  Integration error is above threshhold but better than series. Using integral = " << fmt << ret << endl;
          } else {
            ret = result_series;
            abserr = fabs(error_series);
            if (verbose)
               cout << "  Integration error is above threshhold and worse than series.  Using series = " << fmt << ret << endl;
          }
        }
      } else { // bad theta2 and not avoid series
        abserr = NAN;
        termination_code = IntegrationController<myFloat>::bad_integrand;
        last = 0;
        ret = (small_x_m_zet) ? series_small_x_ddx_pdf(x, pm) : series_large_x_ddx_pdf(x, pm);
        abserr = fabs(error_series);
        if (verbose)
          cout << "  Bad theta2.  Using series = " << fmt << ret << endl;
      } // bad theta2
      return ret;
  } // switch on dist_type
} //ddx_pdf

/// derivative of the pdf using zolotarev expansion for large x
template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::series_large_x_ddx_pdf(myFloat x0, Parameterization pm) {
  if (verbose)
    cout << endl << "series_large_x_ddx_pdf:" << endl;
  myFloat result_series;
  myFloat term{std::numeric_limits<myFloat>::quiet_NaN()};
  myFloat abs_term{std::numeric_limits<myFloat>::quiet_NaN()};
  myFloat Machine_eps = std::numeric_limits<myFloat>::epsilon();
  if (alpha < 1) {
    // Derivative of Zolotarev 2.4.8, a series that converges for large x
    
    myFloat abs_series =0;
    if (rho == 0) {
      result_series = 0;
      error_series = 0;
      n_series = 0;
    } else if (xB==0) {
      result_series = std::numeric_limits<myFloat>::quiet_NaN();
      error_series = std::numeric_limits<myFloat>::infinity();
      n_series = 0;
    } else {
      result_series = 0;
      error_series = 0;
      for (int k=1; k<=max_n; ++k) {
        abs_term = tgamma_ratio(k*alpha+2,myFloat(k+1)) * pow(xB,-k*alpha-2)/(pi*gammaB*gammaB);
        term = abs_term * (1 == k%2 ? -1 : 1) * sin(pi*k*rho*alpha);
        myFloat error = fabs(term) * ((k*alpha+2) * digamma(k*alpha+2) + k*alpha+2) * Machine_eps;
        if (k==1) { // check against dpareto
          myFloat ddx_dPareto = tgamma(alpha)/pi*sin(pi*alpha/2);
          ddx_dPareto *= -(alpha+1)*alpha * (1+beta) * pow(x_m_zet,-alpha-2);
          if (verbose>0)
            cout << "1 - term1/ddx_dPareto = " << fmt << (1 -term/ddx_dPareto) << endl;
        }
        if (verbose>1)
          cout << "n = " << k << ", term = " << fmt << term << ", error = " << fmt << error << endl;
        result_series += term;
        error_series += error;
        abs_series += fabs(term);
        n_series=k;
        if (!boost::math::isfinite(abs_term) || abs_term < 2*Machine_eps*fabs(result_series)) break;
      }
      result_series *= (positive_x) ? 1 : -1;
      error_series += abs_term;
    }
    
  } else if (alpha == 1) {

    // Unlike the convergent integrals, the asymptotic series requires xB>0
    if (xB == 0) {
      result_series = std::numeric_limits<myFloat>::quiet_NaN();
      error_series = std::numeric_limits<myFloat>::max();
      n_series = 0;
    } else if (betaB == -1) {
      // Derivative of Zolotarev formula 2.5.17 which is asymptotic for x large
      myFloat xi = exp(xB-1);
      if (exp(-xi) == 0) {
        result_series = 0;
        n_series = 0;
        error_series = 0;
      } else {
        myFloat nu = 1;
        myFloat fac = -nu*nu*pow(xi,(4-alpha)/(2*alpha))*exp(-xi)/sqrt(2*pi*alpha)/(gammaB*gammaB);
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
        for (int n = 1; n<Q_ddx_pdf.size(); ++n) {
          alpha_xi_n /= (alpha_star*xi);
          term = fac * Q_ddx_pdf.at(n) * alpha_xi_n;
          if (fabs(term) > fabs(old_term)) break;
          if (verbose > 1)
            cout << "n = " << n << ", term = " << fmt << term << endl;
          n_series = n;
          result_series += term;
          old_term = term;
          if (fabs(term) <= eps*fabs(result_series) && fabs(alpha_xi_n) <= eps) break;
        }
        result_series *= (positive_xB)? 1 : -1;
        error_series = fabs(term)
/*        +fabs(result_series)*exp(-xi-pow(xi,.25)) */
        ;
      }
    } else {
      // Zolotarev 2.5.23 is the series series for large x
      // For beta != -1,0
      myFloat log_x=log(xB);
      result_series = 0;
      myFloat term{std::numeric_limits<myFloat>::quiet_NaN()};
      myFloat old_term{std::numeric_limits<myFloat>::quiet_NaN()};
      int num_small_terms=0;
      for (int n=1; n<=max_n; ++n) {
        myFloat fac{0};
        myFloat fac_deriv{0};
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
          if (l>0)
            fac_deriv += r_l_n * l * pow(log_x, l-1) / (xB * gammaB);
          if (verbose>1)
            cout << "fac: " << fmt << fac << endl;
        }
        term= (-(n+1)*fac * pow(xB,-n-2)/gammaB + fac_deriv * pow(xB,-n-1)) / (factorial<myFloat>(n)*pi*gammaB);
        if (n==1) { // check against dpareto
          myFloat ddx_dPareto = tgamma(alpha)/pi*sin(pi*alpha/2);
          ddx_dPareto *= -(alpha+1)*alpha * (betaB_p_1) * pow(fabs(x_m_zet),-alpha-2);
          if (verbose>0)
            cout << "1 - term1/ddx_dPareto = " << fmt << (1 -term/ddx_dPareto) << endl;
        }
        if (verbose>1)
          cout << "n = " << n << ": term = " << fmt << term << endl;
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
      result_series *= (positive_xB) ? 1 : -1;
      error_series=fabs(term) + fabs(result_series) * Machine_eps;
    } // xB != 0
    
  } else {
    // alpha > 1
    
    if (beta == -1) {
      // Proof of Zolotarev Theorem 2.5.2 asymptotic for x -> infinity
      myFloat xi = fabs(alpha_m_1) * pow(xB/alpha, alpha/(alpha_m_1));
      if (xi < 1) {
        result_series = std::numeric_limits<myFloat>::quiet_NaN();
        error_series = std::numeric_limits<myFloat>::max();
        n_series = 0;
      } else {
        myFloat nu = pow(abs(alpha_m_1),-1/alpha);
        myFloat exp_m_xi = exp(-xi);
        myFloat fac = (exp_m_xi!=0) ? -nu*nu*pow(xi,(4-alpha)/(2*alpha))*exp_m_xi/sqrt(2*pi*alpha)/(gammaB*gammaB) : 0;
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
        for (int n = 1; n<Q_ddx_pdf.size(); ++n) {
          alpha_xi_n /= (alpha_star*xi);
          term = fac * Q_ddx_pdf.at(n) * alpha_xi_n;
          if (fabs(term) > fabs(old_term)) break;
          if (verbose > 1)
            cout << "n = " << n << ", term = " << fmt << term << endl;
          n_series = n;
          result_series += term;
          old_term = term;
          if (fabs(term) <= eps*fabs(result_series) && fabs(alpha_xi_n) <= eps) break;
        }
        result_series = (positive_x) ? result_series : -result_series;
        error_series = fabs(term)
        /*        +fabs(result_series)*exp(-pow(xi,.25)) */
        ;
      }
    } else if ( xB == 0) {
      result_series = std::numeric_limits<myFloat>::quiet_NaN();
      error_series = std::numeric_limits<myFloat>::max();
      n_series = 0;
    } else {
      // Derivative of Zolotarev 2.5.4, an asymptotic series that works well for large x
      // for betaB != -1
      n_series = 1;
      abs_term =(alpha/(gammaB*pi))*tgamma(alpha)/tgamma(myFloat(1))*(alpha+1)*pow(xB,-alpha-2)/gammaB;
      term=-abs_term*sin(pi/2*(2-alpha)*(betaB_p_1));
      myFloat ddx_dPareto = tgamma(alpha)/pi*sin(pi*alpha/2);
      ddx_dPareto *= -(alpha+1)*alpha * (1+beta) * pow(x_m_zet,-alpha-2);
      if (verbose>0) {
        cout << "1 - term1/ddx_dPareto = " << fmt << (1 -term/ddx_dPareto) << endl;
        cout << "n = "<< n_series << ", term = " << fmt << term << endl;
      }
      result_series = term;
      error_series = fabs(term) * ((alpha)*digamma(alpha) + alpha +2) * Machine_eps;
      myFloat old_abs_term=abs_term;
      for (int k=2; k<=max_n; ++k) {
        abs_term=(alpha/(gammaB*pi))*tgamma_ratio(k*alpha,myFloat(k))*
        (alpha*k+1)*pow(xB,-alpha*k-2)/gammaB;
        term=abs_term*-sin(pi*k/2*(2-alpha)*(betaB_p_1));
        myFloat error = fabs(term) * ((k*alpha)*digamma(k*alpha) + alpha*k +2) * Machine_eps;
        if (abs_term > old_abs_term) break;
        n_series=k;
        if (verbose>0)
          cout << "n = " << k << ", term = " << fmt << term << endl;
        result_series += term;
        error_series += error;
        if (abs_term < Machine_eps * result_series) break;
        old_abs_term=abs_term;
      }
      result_series *= (positive_x) ? 1 : -1;
      error_series += abs_term;
    } // xB != 0
    
  } // alpha > 1
  if (verbose>0)
    cout << "result_series = " << fmt << result_series
    << ", error_series = " << fmt << error_series << endl
    << ", n_series = " << n_series << endl;
  
  return result_series;
}

/// derivative of the pdf using zolotarev expansion for small x
template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::series_small_x_ddx_pdf(myFloat x0, Parameterization pm) {
  if (verbose)
    cout << endl << "series_small_x_ddx_pdf:" << endl;
  myFloat result_series;
  myFloat term{std::numeric_limits<myFloat>::quiet_NaN()};
  myFloat abs_term{std::numeric_limits<myFloat>::quiet_NaN()};
  myFloat Machine_eps = std::numeric_limits<myFloat>::epsilon();
  if (alpha < 1) {
    if (beta==1) {
      // Modification of the proof of Theorem 2.5.2, asymptotic for small x
      myFloat xi = fabs(alpha_m_1) * pow(xB/alpha, alpha/(alpha_m_1));
      if (xi < 1) {
        result_series = std::numeric_limits<myFloat>::quiet_NaN();
        error_series = std::numeric_limits<myFloat>::max();
        n_series = 0;
      } else {
        myFloat nu = pow(abs(alpha_m_1),-1/alpha);
        myFloat exp_m_xi = exp(-xi);
        myFloat fac = (exp_m_xi != 0) ? nu*nu*pow(xi,(4-alpha)/(2*alpha))*exp_m_xi/sqrt(2*pi*alpha)/(gammaB*gammaB) : 0;
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
        for (int n = 1; n<Q_ddx_pdf.size(); ++n) {
          alpha_xi_n /= (alpha_star*xi);
          term = fac * Q_ddx_pdf.at(n) * alpha_xi_n;
          if (fabs(term) > fabs(old_term)) break;
          if (verbose > 1)
            cout << "n = " << n << ", term = " << fmt << term << endl;
          n_series = n;
          result_series += term;
          old_term = term;
          if (fabs(term) <= eps*fabs(result_series) && fabs(alpha_xi_n) <= eps) break;
        }
        result_series *= (positive_x) ? 1 : -1;
        error_series = fabs(term)
/*        +fabs(result_series)*exp(-pow(xi,.25)) */
        ;
      }
    } else if (beta == 1 && xB >=10) {
      result_series=std::numeric_limits<myFloat>::quiet_NaN();
      error_series=std::numeric_limits<myFloat>::max();
      n_series=0;
    } else {
      // beta != 1
      myFloat abs_tail;
      if (xB == 0) {
        result_series = (betaB!=0 ? tgamma(2/alpha)/tgamma(myFloat(2))/(pi*alpha*gammaB*gammaB) * sin(pi * one_m_betaB) : 0);
        abs_tail = min(std::numeric_limits<myFloat>::max(),fabs(result_series));
        n_series = 0;
      } else {
        // Derivative of Zolotarev formula 2.5.1, which is asmyptotic as x->0
        // for betaB != 1
        abs_term = 0;
        term = 0;
        myFloat old_abs_term{abs_term};
        result_series = term;
        error_series = 0;
        abs_tail = fabs(term);
        n_series =0;
        if (verbose>0)
          cout << "n = " << n_series << ", term = " << fmt << term << endl;
        for (int k=1; k<=max_n; ++k) {
          abs_term = tgamma_ratio((k+1)/alpha,myFloat(k+1))*pow(xB,k-1)*k/(pi * alpha * gammaB*gammaB);
          term = abs_term * sin(pi/2 * (k+1) * one_m_betaB);
          if (k > 1 && abs_term > k*old_abs_term) break;
          if (verbose>0)
            cout << "n = " << k << ", term = " << fmt << term << endl;
          n_series = k;
          result_series += term;
          abs_tail += fabs(term);
          if (abs_term < Machine_eps * result_series) break;
          old_abs_term=abs_term;
        }
      }
      result_series *= (positive_x) ? 1 : -1;
      error_series = ((xB!=0 && beta!=-1)?abs_term:0) + abs_tail*Machine_eps;
    }
    
  } else if (alpha == 1) {
    // Derivatie of Zolotarev Formula 2.2.9 using sin as Im(e^ix)
    // given by infinite integral
    n_series = 0;
    result_series=std::numeric_limits<myFloat>::quiet_NaN();
    error_series=std::numeric_limits<myFloat>::max();
    
  } else {
    // alpha > 1
    // Derivative of Zolotarev Formula 2.4.6 is a convergent series that works well for x <=1, but
    //    but suffers from extreme rounding errors for larger x.
    
    result_series = 0;
    error_series = 0;
    myFloat abs_series =0;
    myFloat term{std::numeric_limits<myFloat>::quiet_NaN()};
    myFloat abs_term{std::numeric_limits<myFloat>::quiet_NaN()};
    for (int n=2; n<=max_n; ++n) {
      abs_term = tgamma_ratio(n/alpha+1,myFloat(n+1)) * pow(xB,n-2)*(n-1)/(pi*gammaB*gammaB);
      term = abs_term * (0==n%2 ? -1 : 1) * ((n*rho)!=int(n*rho) ? sin(pi*n*rho) : 0);
      myFloat error = fabs(term) * ((n/alpha+1) * digamma(n/alpha+1) + (n-1))* Machine_eps;
      if (verbose>1)
        cout << "n = " << n << ", term = " << fmt << term << ", error = " << fmt << error << endl;
      n_series = n;
      result_series += term;
      error_series += error;
      abs_series += fabs(term);
      if (!boost::math::isfinite(abs_term) || abs_term < 2*Machine_eps*fabs(result_series)) break;
    }
    result_series *= (positive_x) ? 1 : -1;
    error_series += abs_term;
    
  } // alpha > 1
  if (verbose>0)
    cout << "result_series = " << fmt << result_series
    << ", error_series = " << fmt << error_series << endl
    << ", n_series = " << n_series << endl;
  
  return result_series;
}
} // namespace stable_distribution
