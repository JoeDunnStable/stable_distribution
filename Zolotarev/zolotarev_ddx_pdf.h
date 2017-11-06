/// \file zolotarev_ddx_pdf.h
/// \author Joseph Dunn
/// \copyright 2017 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include <iostream>
// This will be included at the end of zolotarev.h if LIBRARY is defined

namespace stable_distribution {

using std::cout;
using std::endl;

/// derivative of the pdf using zolotarev expansion
template<typename myFloat>
myFloat Zolotarev<myFloat>::ddx_pdf(myFloat x0, Parameterization pm) {
  switch (pm) {
    case S0:
      set_x_m_zeta(x0-zeta);
      break;
    case S1:
      set_x_m_zeta(x0);
  }
  myFloat term, abs_term;
  myFloat Machine_eps = std::numeric_limits<myFloat>::epsilon();
  if (alpha < 1) {
    // Derivative of Zolotarev 2.4.8, a series that converges for large x
    
    myFloat abs_series =0;
    if (rho == 0) {
      result_convergent = 0;
      error_convergent = 0;
      n_convergent = 0;
    } else if (xB==0) {
      result_convergent = std::numeric_limits<myFloat>::quiet_NaN();
      error_convergent = std::numeric_limits<myFloat>::infinity();
      n_convergent = 0;
    } else {
      result_convergent = 0;
      error_convergent = 0;
      for (int k=1; k<=max_n_conv; ++k) {
        abs_term = tgamma_ratio(k*alpha+2,myFloat(k+1)) * pow(xB,-k*alpha-2)/(pi*gammaB*gammaB);
        term = abs_term * (1 == k%2 ? -1 : 1) * sin(pi*k*rho*alpha);
        myFloat error = fabs(term) * ((k*alpha+2) * digamma(k*alpha+2) + k*alpha+2) * Machine_eps;
        if (k==1) { // check against dpareto
          myFloat ddx_dPareto = tgamma(alpha)/pi*sin(pi*alpha/2);
          ddx_dPareto *= -(alpha+1)*alpha * (1+beta) * pow(x_m_zet,-alpha-2);
          if (verbose>0)
            cout << "1 - term1/ddx_dPareto = " << (1 -term/ddx_dPareto) << endl;
        }
        if (verbose>1)
          cout << "n = " << k << ", term = " << term << ", error = " << error << endl;
        result_convergent += term;
        error_convergent += error;
        abs_series += fabs(term);
        n_convergent=k;
        if (!boost::math::isfinite(abs_term) || abs_term < 2*Machine_eps*fabs(result_convergent)) break;
      }
      result_convergent *= (positive_x) ? 1 : -1;
      error_convergent += abs_term;
    }
    if (verbose>0)
      cout << "result convergent = " << result_convergent
      << ", error_convergent = " << error_convergent << endl;
    
    if (beta==1) {
      // Modification of the proof of Theorem 2.5.2, asymptotic for small x
      myFloat psi = fabs(1-alpha) * pow(xB/alpha, alpha/(alpha-1));
      if (psi < 1) {
        result_asymptotic = std::numeric_limits<myFloat>::quiet_NaN();
        error_asymptotic = std::numeric_limits<myFloat>::max();
        n_asymptotic = 0;
      } else {
        myFloat nu = pow(abs(1-alpha),-1/alpha);
        myFloat exp_m_psi = exp(-psi);
        myFloat fac = (exp_m_psi != 0) ? nu*nu*pow(psi,(4-alpha)/(2*alpha))*exp_m_psi/sqrt(2*pi*alpha)/(gammaB*gammaB) : 0;
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
        for (int n = 1; n<Q_ddx_pdf.size(); ++n) {
          alpha_psi_n /= (alpha_star*psi);
          term = fac * Q_ddx_pdf.at(n) * alpha_psi_n;
          if (fabs(term) > fabs(old_term)) break;
          if (verbose > 1)
            cout << "n = " << n << ", term = " << term << endl;
          n_asymptotic = n;
          result_asymptotic += term;
          old_term = term;
        }
        result_asymptotic *= (positive_x) ? 1 : -1;
        error_asymptotic = fabs(term)+fabs(result_asymptotic)*exp(-pow(psi,.25));
      }
    } else if (beta == 1 && xB >=10) {
      result_asymptotic=std::numeric_limits<myFloat>::quiet_NaN();
      error_asymptotic=std::numeric_limits<myFloat>::max();
      n_asymptotic=0;
    } else {
      // beta != 1
      myFloat abs_tail;
      if (xB == 0) {
        result_asymptotic = (betaB!=0 ? tgamma(2/alpha)/tgamma(myFloat(2))/(pi*alpha*gammaB*gammaB) * sin(pi * (1-betaB)) : 0);
        abs_tail = min(std::numeric_limits<myFloat>::max(),fabs(result_asymptotic));
        n_asymptotic = 0;
      } else {
        // Derivative of Zolotarev formula 2.5.1, which is asmyptotic as x->0
        // for betaB != 1
        abs_term = 0;
        term = 0;
        myFloat old_abs_term{abs_term};
        result_asymptotic = term;
        error_asymptotic = 0;
        abs_tail = fabs(term);
        n_asymptotic =0;
        if (verbose>0)
          cout << "n = " << n_asymptotic << ", term = " << term << endl;
        for (int k=1; k<=max_n_asymp; ++k) {
          abs_term = tgamma_ratio((k+1)/alpha,myFloat(k+1))*pow(xB,k-1)*k/(pi * alpha * gammaB*gammaB);
          term = abs_term * sin(pi/2 * (k+1) * (1-betaB));
          if (k > 1 && abs_term > k*old_abs_term) break;
          if (verbose>0)
            cout << "n = " << k << ", term = " << term << endl;
          n_asymptotic = k;
          result_asymptotic += term;
          abs_tail += fabs(term);
          if (abs_term < Machine_eps * result_asymptotic) break;
          old_abs_term=abs_term;
        }
      }
      result_asymptotic *= (positive_x) ? 1 : -1;
      error_asymptotic = ((xB!=0 && beta!=-1)?abs_term:0) + abs_tail*Machine_eps;
    }
    if (verbose > 0)
      cout << "result asymptotic = " << result_asymptotic
      << ", error_asymptotic = " << error_asymptotic << endl
      << "result_asymptotic - result_convergent = " << result_asymptotic-result_convergent << endl;
    
  } else if (alpha == 1) {
    // Derivatie of Zolotarev Formula 2.2.9 using sin as Im(e^ix)
    // given by infinite integral
    n_convergent = 0;
    if (betaB == 0) { // Derivative of Cauchy Distribution
      result_convergent=-2*x_m_zet/(pi*pow(1+x_m_zet*x_m_zet,2));
      result_convergent *= positive_x ? 1 : -1;
      error_convergent=5*std::numeric_limits<myFloat>::epsilon()*result_convergent;
    } else {
      result_convergent = ((positive_x)? 1 : -1) * ddx_pdf_alpha_1();
      error_convergent = ddx_pdf_alpha_1.abserr;
    }
    if (verbose>0)
      cout << "result convergent = " << result_convergent
      << ", error_convergent = " << error_convergent << endl;
    // Unlike the convergent integrals, the asymptotic series requires x>0
    myFloat xB_ = fabs(xB) + log(gammaB);
    myFloat betaB_ =(xB>=0) ? betaB : -betaB;
    myFloat beta_ = (xB>=0) ? beta : -beta;
    bool positive_x_ = (xB < 0) != positive_x;
    if (xB == 0) {
      result_asymptotic = std::numeric_limits<myFloat>::quiet_NaN();
      error_asymptotic = std::numeric_limits<myFloat>::max();
      n_asymptotic = 0;
    } else if (beta_ == -1) {
      // Derivative of Zolotarev formula 2.5.17 which is asymptotic for x large
      myFloat psi = exp(xB_-1);
      if (exp(-psi) == 0) {
        result_asymptotic = 0;
        n_asymptotic = 0;
        error_asymptotic = 0;
      } else {
        myFloat nu = 1;
        myFloat fac = -nu*nu*pow(psi,(4-alpha)/(2*alpha))*exp(-psi)/sqrt(2*pi*alpha)/(gammaB*gammaB);
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
        for (int n = 1; n<Q_ddx_pdf.size(); ++n) {
          alpha_psi_n /= (alpha_star*psi);
          term = fac * Q_ddx_pdf.at(n) * alpha_psi_n;
          if (fabs(term) > fabs(old_term)) break;
          if (verbose > 1)
            cout << "n = " << n << ", term = " << term << endl;
          n_asymptotic = n;
          result_asymptotic += term;
          old_term = term;
        }
        result_asymptotic *= (positive_x_)? 1 : -1;
        error_asymptotic = fabs(term)+fabs(result_asymptotic)*exp(-psi-pow(psi,.25));
      }
    } else {
      // Zolotarev 2.5.23 is the asymptotic series for large x
      // For beta != -1,0
      myFloat log_x=log(xB_);
      result_asymptotic = 0;
      myFloat term, old_term;
      int num_small_terms=0;
      for (int n=1; n<=max_n_asymp; ++n) {
        myFloat fac{0};
        myFloat fac_deriv{0};
        for (int l=0; l<=n; ++l){
          myFloat r_l_n{0};
          for (int m=l; m<=n; ++m) {
            myFloat term0=binomial_coefficient<myFloat>(n,m)*binomial_coefficient<myFloat>(m,l);
            term0 *= (1-2*((m-l)%2)) * gamma_at_integers(m-l,1+n);
            term0 *=  pow(betaB_, m) * pow(pi/2*(1+betaB_),n-m) * sin(pi/2*(n-m));
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
          if (l>0)
            fac_deriv += r_l_n * l * pow(log_x, l-1) / (xB_ * gammaB);
          if (verbose>1)
            cout << "fac: " << fac << endl;
        }
        term= (-(n+1)*fac * pow(xB_,-n-2)/gammaB + fac_deriv * pow(xB_,-n-1)) / (factorial<myFloat>(n)*pi*gammaB);
        if (n==1) { // check against dpareto
          myFloat ddx_dPareto = tgamma(alpha)/pi*sin(pi*alpha/2);
          ddx_dPareto *= -(alpha+1)*alpha * (1+beta) * pow(x_m_zet,-alpha-2);
          if (verbose>0)
            cout << "1 - term1/ddx_dPareto = " << (1 -term/ddx_dPareto) << endl;
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
      result_asymptotic *= (positive_x_) ? 1 : -1;
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
    // Derivative of Zolotarev Formula 2.4.6 is a convergent series that works well for x <=1, but
    //    but suffers from extreme rounding errors for larger x.
    
    result_convergent = 0;
    error_convergent = 0;
    myFloat abs_series =0;
    myFloat term, abs_term;
    for (n=2; n<=max_n_conv; ++n) {
      abs_term = tgamma_ratio(n/alpha+1,myFloat(n+1)) * pow(xB,n-2)*(n-1)/(pi*gammaB*gammaB);
      term = abs_term * (0==n%2 ? -1 : 1) * ((n*rho)!=int(n*rho) ? sin(pi*n*rho) : 0);
      myFloat error = fabs(term) * ((n/alpha+1) * digamma(n/alpha+1) + (n-1))* Machine_eps;
      if (verbose>1)
        cout << "n = " << n << ", term = " << term << ", error = " << error << endl;
      n_convergent = n;
      result_convergent += term;
      error_convergent += error;
      abs_series += fabs(term);
      if (!boost::math::isfinite(abs_term) || abs_term < 2*Machine_eps*fabs(result_convergent)) break;
    }
    result_convergent *= (positive_x) ? 1 : -1;
    error_convergent += abs_term;
    if (verbose>0)
      cout << "result_convergent = " << result_convergent
      << ", error_convergent = " << error_convergent << endl;
    
    if (beta == -1) {
      // Proof of Zolotarev Theorem 2.5.2 asymptotic for x -> infinity
      myFloat psi = fabs(1-alpha) * pow(xB/alpha, alpha/(alpha-1));
      if (psi < 1) {
        result_asymptotic = std::numeric_limits<myFloat>::quiet_NaN();
        error_asymptotic = std::numeric_limits<myFloat>::max();
        n_asymptotic = 0;
      } else {
        myFloat nu = pow(abs(1-alpha),-1/alpha);
        myFloat exp_m_psi = exp(-psi);
        myFloat fac = (exp_m_psi!=0) ? -nu*nu*pow(psi,(4-alpha)/(2*alpha))*exp_m_psi/sqrt(2*pi*alpha)/(gammaB*gammaB) : 0;
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
        for (int n = 1; n<Q_ddx_pdf.size(); ++n) {
          alpha_psi_n /= (alpha_star*psi);
          term = fac * Q_ddx_pdf.at(n) * alpha_psi_n;
          if (fabs(term) > fabs(old_term)) break;
          if (verbose > 1)
            cout << "n = " << n << ", term = " << term << endl;
          n_asymptotic = n;
          result_asymptotic += term;
          old_term = term;
        }
        result_asymptotic = (positive_x) ? result_asymptotic : -result_asymptotic;
        error_asymptotic = fabs(term)+fabs(result_asymptotic)*exp(-pow(psi,.25));
      }
    } else if ( xB == 0) {
      result_asymptotic = std::numeric_limits<myFloat>::quiet_NaN();
      error_asymptotic = std::numeric_limits<myFloat>::max();
      n_asymptotic = 0;
    } else {
      // Derivative of Zolotarev 2.5.4, an asymptotic series that works well for large x
      // for betaB != -1
      n_asymptotic = 1;
      abs_term =(alpha/(gammaB*pi))*tgamma(alpha)/tgamma(myFloat(1))*(alpha+1)*pow(xB,-alpha-2)/gammaB;
      term=-abs_term*sin(pi/2*(2-alpha)*(1+betaB));
      myFloat ddx_dPareto = tgamma(alpha)/pi*sin(pi*alpha/2);
      ddx_dPareto *= -(alpha+1)*alpha * (1+beta) * pow(x_m_zet,-alpha-2);
      if (verbose>0) {
        cout << "1 - term1/ddx_dPareto = " << (1 -term/ddx_dPareto) << endl;
        cout << "n = "<< n_asymptotic << ", term = " << term << endl;
      }
      result_asymptotic = term;
      myFloat old_abs_term=abs_term;
      for (int k=2; k<=max_n_asymp; ++k) {
        abs_term=(alpha/(gammaB*pi))*tgamma_ratio(k*alpha,myFloat(k))*
        (alpha*k+1)*pow(xB,-alpha*k-2)/gammaB;
        term=abs_term*-sin(pi*k/2*(2-alpha)*(1+betaB));
        myFloat error = fabs(term) * ((k*alpha)*digamma(k*alpha) + alpha*k +2) * Machine_eps;
        if (abs_term > old_abs_term) break;
        n_asymptotic=k;
        if (verbose>0)
          cout << "n = " << k << ", term = " << term << endl;
        result_asymptotic += term;
        error_asymptotic += error;
        if (abs_term < Machine_eps * result_asymptotic) break;
        old_abs_term=abs_term;
      }
      result_asymptotic *= (positive_x) ? 1 : -1;
      error_asymptotic += abs_term;
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

