/// \file zolotarev_pdf.h
/// Implementation of pdf of standard stable distribution per Zolotarev
/// Included in zolotarev.h when LIBRARY is defined
/// \author Joseph Dunn
/// \copyright 2017 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include <boost/math/distributions/cauchy.hpp>
// This will be included at the end of zolotarev.h if LIBRARY is defined

namespace stable_distribution {

using boost::math::binomial_coefficient;
using boost::math::factorial;
using std::endl;

/// pdf using zolotarev expansion
template<typename myFloat>
myFloat Zolotarev<myFloat>::pdf(myFloat x0, Parameterization pm) {
  switch (pm) {
    case S0:
      set_x_m_zeta(x0-zeta);
      break;
    case S1:
      set_x_m_zeta(x0);
  }
  myFloat Machine_eps = std::numeric_limits<myFloat>::epsilon();
  if (alpha < 1) {
    // Zolotarev 2.4.8, a series that converges for large x
    
    if (rho == 0) {
      result_convergent = 0;
      error_convergent = 0;
      n_convergent = 0;
    } else if (xB==0) {
      result_convergent = std::numeric_limits<myFloat>::quiet_NaN();
      error_convergent = std::numeric_limits<myFloat>::infinity();
      n_convergent = 0;
    } else {
      myFloat term{0}, abs_term{0};
      result_convergent = 0;
      error_convergent = 0;
      for (int k=1; k<=max_n_conv; ++k) {
        abs_term = tgamma_ratio(k*alpha+1,myFloat(k+1)) * pow(xB,-k*alpha-1)/(pi*gammaB);
        term = abs_term * (0 == k%2 ? -1 : 1) * sin(pi*k*rho*alpha);
        myFloat error = fabs(term) * ((k*alpha+1)*digamma(k*alpha+1) +k*alpha +1) * Machine_eps;
        if (k==1) { // check against dpareto
          myFloat dPareto = tgamma(alpha)/pi*sin(pi*alpha/2);
          dPareto *= alpha * (1+beta) * pow(x_m_zet,-alpha-1);
          if (verbose>0)
            cout << "1 - term1/dPareto = " << (1 -term/dPareto) << endl;
        }
        if (verbose>1)
          cout << "n = " << k << ", term = " << term << endl;
        result_convergent += term;
        error_convergent += error;
        n_convergent=k;
        if (!boost::math::isfinite(abs_term) || abs_term < 2*Machine_eps*fabs(result_convergent)) break;
      }
      error_convergent += abs_term;
    } // rho != 0
    if (verbose>0)
      cout << "result convergent = " << result_convergent
      << ", error_convergent = " << error_convergent << endl;
    
    if (beta == 1) {
      // Zolotarev Theorem 2.5.2, asymptotic for small x
      myFloat psi = fabs(1-alpha) * pow(xB/alpha, alpha/(alpha-1));
      myFloat nu = pow(abs(1-alpha),-1/alpha);
      myFloat exp_m_psi = exp(-psi);
      myFloat fac = exp_m_psi !=0 ?nu*pow(psi,(2-alpha)/(2*alpha))*exp_m_psi/sqrt(2*pi*alpha)/gammaB
                                  :0;
      if (verbose > 1)
        cout << "psi = " << psi << endl
        << "nu  = " << nu << endl
        << "fac = " << fac << endl;
      myFloat term = fac;
      if (verbose > 1)
        cout << "n = " << 0 << ", term = " << term << endl;
      myFloat old_term = term;
      result_asymptotic = term;
      myFloat alpha_star = alpha;
      myFloat alpha_psi_n = 1;
      for (int n = 1; n<Q_pdf.size(); ++n) {
        alpha_psi_n /= (alpha_star*psi);
        term = fac * Q_pdf.at(n) * alpha_psi_n;
        if (fabs(term) > fabs(old_term)) break;
        if (verbose > 1)
          cout << "n = " << n << ", term = " << term << endl;
        n_asymptotic = n;
        result_asymptotic += term;
        old_term = term;
      }
      error_asymptotic = fabs(term)+fabs(result_asymptotic)*exp(-pow(psi,.25));
      
    } else { // beta != 1
      // Zolotarev formula 2.5.1, which is asmyptotic as x->0
      // for betaB != 1
      myFloat abs_term = tgamma(1/alpha)/(pi * alpha * gammaB);
      myFloat term = abs_term * ((beta!=-1)?sin(pi/2*(1-betaB)):0);
      myFloat abs_tail{fabs(term)};
      if (xB == 0 || beta== -1) {
        result_asymptotic = term;
        error_asymptotic = abs_tail * Machine_eps;
        n_asymptotic =0;
        if (verbose>0)
          cout << "n = " << n_asymptotic << ", term = " << term << endl;
      } else { // xB !=0
        result_asymptotic = 0;
        // myFloat cos_ab = pow(cos((pi/2)*alpha*betaB),1/alpha);
        myFloat cos_ab = 1;
        myFloat cos_ab_m_k = 1/cos_ab;
        error_asymptotic = std::numeric_limits<myFloat>::infinity();
        for (int k=1; k<=max_n_asymp; ++k) {
          abs_term = tgamma_ratio((k+1)/alpha,myFloat(k+1))*pow(xB,k)/(pi * alpha * gammaB);
          myFloat new_term = abs_term * sin(pi/2 * (k+1) * (1-betaB));
          cos_ab_m_k /= cos_ab;
          myFloat new_error = abs_term * cos_ab_m_k;
          if (new_error > error_asymptotic) break;
          n_asymptotic = k-1;
          result_asymptotic += term;
          abs_tail += fabs(term);
          error_asymptotic=new_error;
          if (verbose>0)
            cout << "n = " << k-1 << ", term = " << term  << ", error = " << error_asymptotic << endl;
          term = new_term;
          if (abs_term < Machine_eps * result_asymptotic) break;
        } // k loop
        error_asymptotic += abs_tail*Machine_eps;
      } // xB != 0 && beta != -1
    }
    if (verbose > 0)
      cout << "result asymptotic = " << result_asymptotic
      << ", error_asymptotic = " << error_asymptotic << endl
      << "result_asymptotic - result_convergent = " << result_asymptotic-result_convergent << endl;
  
  } else if (alpha == 1) {
    // Zolotarev Formula 2.2.9 using sin as Im(e^ix)
    // given by infinite integral
    n_convergent = 0;
    if (betaB == 0) {
      using namespace boost::math;
      cauchy_distribution<myFloat> dist_cauchy(0,1);
      result_convergent = boost::math::pdf(dist_cauchy, x_m_zet);
      error_convergent=5*std::numeric_limits<myFloat>::epsilon()*result_convergent;
    } else if (!boost::math::isfinite(x0)) {
      result_convergent = 0;
      error_convergent = 0;
    } else {
      result_convergent = pdf_alpha_1();
      error_convergent = pdf_alpha_1.abserr;
    }
    if (verbose>0)
      cout << "result convergent = " << result_convergent
      << ", error_convergent = " << error_convergent << endl;
    
    // Zolotarev 2.5.23 is the asymptotic series for large x
    // For beta != -1,0
    // Unlike the convergent formula this one requires xB > 0

    if (xB == 0) {
      result_asymptotic = std::numeric_limits<myFloat>::quiet_NaN();
      error_asymptotic = std::numeric_limits<myFloat>::max();
      n_asymptotic = 0;
    } else if (betaB == -1) {
      // Zolotarev Theorem 2.5.2 asymptotic for x -> infinity
      myFloat psi = exp(xB-1);
      myFloat nu = 1;
      myFloat exp_m_psi = exp(-psi);
      myFloat fac = exp_m_psi !=0 ? nu*pow(psi,(2-alpha)/(2*alpha))*exp_m_psi/sqrt(2*pi*alpha)/gammaB
                                  : 0;
      if (verbose > 1)
        cout << "psi = " << psi << endl
        << "nu  = " << nu << endl
        << "fac = " << fac << endl;
      myFloat term = fac;
      if (verbose > 1)
        cout << "n = " << 0 << ", term = " << term << endl;
      myFloat old_term = term;
      result_asymptotic = term;
      myFloat alpha_star = 1/alpha;
      myFloat alpha_psi_n = 1;
      for (int n = 1; n<Q_pdf.size(); ++n) {
        alpha_psi_n /= (alpha_star*psi);
        term = fac * Q_pdf.at(n) * alpha_psi_n;
        if (fabs(term) > fabs(old_term)) break;
        if (verbose > 1)
          cout << "n = " << n << ", term = " << term << endl;
        n_asymptotic = n;
        result_asymptotic += term;
        old_term = term;
      }
      error_asymptotic = fabs(term)+fabs(result_asymptotic)*exp(-pow(psi,.25));
    } else {
      myFloat log_x=log(xB);
      result_asymptotic = 0;
      myFloat term{std::numeric_limits<myFloat>::quiet_NaN()};
      myFloat old_term{std::numeric_limits<myFloat>::quiet_NaN()};
      int num_small_terms=0;
      for (int n=1; n<=max_n_asymp; ++n) {
        myFloat fac{0};
        for (int l=0; l<=n; ++l){
          myFloat r_l_n{0};
          for (int m=l; m<=n; ++m) {
            myFloat term0=binomial_coefficient<myFloat>(n,m)*binomial_coefficient<myFloat>(m,l);
            term0 *= (1-2*((m-l)%2)) * gamma_at_integers(m-l,1+n);
            term0 *=  pow(betaB, m) * pow(pi/2*(1+betaB),n-m) * sin(pi/2*(n-m));
            r_l_n += term0;
            if (verbose > 2){
              if (m == l)
                cout << "r(" << l << ", " << n << "} = " << term0;
              else
                cout << " + " << term0;
            }
          }
          if (verbose > 2)
            cout << endl;
          if (verbose>1)
            cout << "r(" << l << ", " << n << ") * log_x ^ l = "
            << r_l_n << " * " << log_x << " ^ " << l
            << " = " << r_l_n << " * " << pow(log_x,l)
            << " = " << r_l_n * pow(log_x,l) << endl;
          fac += r_l_n*pow(log_x,l);
          if (verbose>1)
            cout << "fac: " << fac << endl;
        }
        term=fac * pow(xB,-n-1) / (factorial<myFloat>(n)*pi*gammaB);
        if (n==1) { // check against dpareto
          myFloat dPareto = tgamma(alpha)/pi*sin(pi*alpha/2);
          dPareto *= alpha * (1+betaB) * pow(fabs(x_m_zet),-alpha-1);
          if (verbose>0)
            cout << "1 - term1/dPareto = " << (1 -term/dPareto) << endl;
        }
        if (verbose>1)
          cout << "n = " << n << ": term = " << term << endl;
        if (n>1 && fabs(term) > fabs(old_term)) break;
        n_asymptotic = n;
        result_asymptotic += term;
        if (verbose>0)
          cout << "n = " << n
          << ": result_asymptotic = " << result_asymptotic << endl;
        if (fabs(term) < fabs(result_asymptotic) * Machine_eps) {
          if (num_small_terms > 0) break;
          num_small_terms=num_small_terms+1;
        } else {
          num_small_terms=0;
          old_term=term;
        }
      }
      error_asymptotic=fabs(term) + fabs(result_asymptotic) * Machine_eps;
    } // xB != 0
    if (verbose>0)
      cout << "result_asymptotic = " << result_asymptotic
      << ", error_asymptotic = " << error_asymptotic << endl
      << " result_asymptotic - result_convergent = " << result_asymptotic - result_convergent << endl;
    
    if (boost::math::isnan(result_convergent) || !boost::math::isfinite(result_convergent)
        ||  (!boost::math::isnan(result_asymptotic)
             && error_convergent > max(100*Machine_eps*fabs(result_convergent),error_asymptotic))) {
      result_type = asymptotic;
      n = n_asymptotic;
      error = error_asymptotic;
      result = result_asymptotic;
    } else {
      result_type = convergent;
      n = n_convergent;
      error = error_convergent;
      result = result_convergent;
    }
    return result;
    
  } else {
    // alpha > 1
    // Zolotarev Formula 2.4.6 is a convergent series that works well for x <=1, but
    //    but suffers from extreme rounding errors for larger x.
    
    result_convergent = 0;
    error_convergent =0;
    myFloat term, abs_term{std::numeric_limits<myFloat>::quiet_NaN()};
    for (n=1; n<=max_n_conv; ++n) {
      abs_term = tgamma_ratio(n/alpha+1, myFloat(n+1))*pow(xB,n-1)/(pi*gammaB);
      term = abs_term * (0==n%2 ? -1 : 1) * sin(pi*n*rho);
      myFloat error = fabs(term)*((n/alpha+1)*digamma(n/alpha+1)+n-1) * Machine_eps;
      if (verbose>1)
        cout << "n = " << n << ", term = " << term << endl;
      n_convergent = n;
      result_convergent += term;
      error_convergent += error;
      if (!boost::math::isfinite(abs_term) || abs_term < 2*Machine_eps*fabs(result_convergent)) break;
    }
    error_convergent += abs_term;
    if (verbose>0)
      cout << "result_convergent = " << result_convergent
      << ", error_convergent = " << error_convergent << endl;
    
    if (xB == 0) {
      result_asymptotic = std::numeric_limits<myFloat>::quiet_NaN();
      error_asymptotic = std::numeric_limits<myFloat>::max();
      n_asymptotic = 0;
    } else if (beta == -1) {
      // Zolotarev Theorem 2.5.2 asymptotic for x -> infinity
      myFloat psi = fabs(1-alpha) * pow(xB/alpha, alpha/(alpha-1));
      myFloat nu = pow(abs(1-alpha),-1/alpha);
      myFloat exp_m_psi = exp(-psi);
      myFloat fac = exp_m_psi !=0 ? nu*pow(psi,(2-alpha)/(2*alpha))*exp(-psi)/sqrt(2*pi*alpha)/gammaB
                                  :0;
      if (verbose > 1)
        cout << "psi = " << psi << endl
        << "nu  = " << nu << endl
        << "fac = " << fac << endl;
      myFloat term = fac;
      myFloat old_term = term;
      if (verbose > 1)
        cout << "n = " << 0 << ", term = " << term << endl;
      result_asymptotic = term;
      myFloat alpha_star = 1/alpha;
      myFloat alpha_psi_n = 1;
      for (int n = 1; n<Q_pdf.size(); ++n) {
        alpha_psi_n /= (alpha_star*psi);
        term = fac * Q_pdf.at(n) * alpha_psi_n;
        if (fabs(term) > fabs(old_term)) break;
        if (verbose > 1)
          cout << "n = " << n << ", term = " << term << endl;
        n_asymptotic = n;
        result_asymptotic += term;
        old_term = term;
      }
      error_asymptotic = fabs(term)+result_asymptotic*exp(-pow(psi,.25));
      
    } else {
      // Zolotarev 2.5.4, an asymptotic series that works well for large x
      n_asymptotic = 1;
      // myFloat alpha_star = 1/alpha;
      // myFloat beta_star = 1-(2-alpha)*(1+betaB);
      myFloat abs_term =(alpha/(gammaB*pi))*tgamma_ratio(alpha,myFloat(1))*pow(xB,-alpha-1);
      myFloat term=abs_term*sin(pi/2*(2-alpha)*(1+betaB));
      myFloat abs_tail{0};
      result_asymptotic = 0;
      // myFloat cos_ab = pow(cos((pi/2)*alpha_star*beta_star),1/alpha_star);
      myFloat cos_ab = 1;
      myFloat cos_ab_m_k = 1/cos_ab;
      error_asymptotic = std::numeric_limits<myFloat>::infinity();
      myFloat dPareto = tgamma(alpha)/pi*sin(pi*alpha/2);
      dPareto *= alpha * (1+beta) * pow(x_m_zet,-alpha-1);
      if (verbose>0) {
        cout << "1 - term1/dPareto = " << (1 -term/dPareto) << endl;
      }
      for (int k=2; k<=max_n_asymp; ++k) {
        abs_term=(alpha/(gammaB*pi))*tgamma_ratio(k*alpha,myFloat(k))* pow(xB,-alpha*k-1);
        myFloat new_term = abs_term*sin(pi*k/2*(2-alpha)*(1+betaB));
        cos_ab_m_k /= cos_ab;
        myFloat new_error = abs_term * cos_ab_m_k;
        if (new_error > error_asymptotic) break;
        n_asymptotic=k-1;
        result_asymptotic += term;
        abs_tail += fabs(term);
        error_asymptotic = new_error;
        if (verbose>0)
          cout << "n = " << k-1 << ", term = " << term << ", error = " << error_asymptotic << endl;
        term = new_term;
        if (abs_term < Machine_eps * result_asymptotic) break;
      }
      error_asymptotic +=  abs_tail*Machine_eps;
    } // xB != 0
    if (verbose>0)
      cout << "result_asymptotic = " << result_asymptotic
      << ", error_asymptotic = " << error_asymptotic << endl
      << "result_asymptotic - result_convergent = " << result_asymptotic-result_convergent << endl;
    
  } // alpha > 1
  if (boost::math::isnan(result_convergent)
      || (error_convergent > error_asymptotic)) {
    result_type = asymptotic;
    n = n_asymptotic;
    error = error_asymptotic;
    result = result_asymptotic;
  } else {
    result_type = convergent;
    n = n_convergent;
    error = error_convergent;
    result = result_convergent;
  }
  return result;
}
}  // namespace stable_distribution

