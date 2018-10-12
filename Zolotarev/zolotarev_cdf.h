/// \file zolotarev_cdf.h
/// Implementation of cdf for standard stable distribution per Zolotarev
/// Included in zolotarev.h when LIBRARY is defined
/// \author Joseph Dunn
/// \copyright 2017 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include <boost/math/distributions/cauchy.hpp>
// Thie will be included at the end of zolotarev.h if LIBRARY is defined

namespace stable_distribution {
  
using std::endl;

/// cdf using zolotarev expansion
template<typename myFloat>
myFloat Zolotarev<myFloat>::cdf(myFloat x0, int lower_tail, Parameterization pm) {
  myFloat x_m_zeta_in;
  switch (pm) {
    case S0:
      x_m_zeta_in = x0-zeta;
      break;
    case S1:
      x_m_zeta_in = x0;
  }
  set_x_m_zeta(x_m_zeta_in);
  myFloat Machine_eps = std::numeric_limits<myFloat>::epsilon();
  if (alpha < 1) {
    // Zolotarev formula 2.4.3 for upper tail, convergent for large x
    if (xB == 0) {
      result_convergent = std::numeric_limits<myFloat>::quiet_NaN();
      error_convergent = std::numeric_limits<myFloat>::infinity();
      n_convergent = 0;
    } else {
      result_convergent = 0;
      error_convergent = 0;
      myFloat abs_series = 0;
      myFloat term{std::numeric_limits<myFloat>::quiet_NaN()};
      myFloat abs_term{std::numeric_limits<myFloat>::quiet_NaN()};
      for (int n=1; n<=max_n_conv; ++n) {
        abs_term = tgamma_ratio(n*alpha+1,myFloat(n+1)) * pow(xB,-n*alpha)/(n*alpha*pi);
        term = abs_term * (0 == n%2 ? -1 : 1) * sin(pi*n*rho*alpha);
        myFloat error = fabs(term) * ((n*alpha+1)*digamma(n*alpha+1) + n*alpha) * Machine_eps;
        if (verbose>1)
          cout << n << "   " << setprecision(digits10) << scientific<< term << endl;
        n_convergent = n;
        result_convergent += term;
        error_convergent += error;
        abs_series += fabs(term);
        if (!boost::math::isfinite(abs_term) || abs_term < 2*Machine_eps*fabs(result_convergent)) break;
      }
      result_convergent = (lower_tail == positive_x) ? 1-result_convergent : result_convergent;
      error_convergent += abs_term;
    } //xB !=0
    if (verbose>0)
      cout <<"result_convergent = " << setprecision(digits10) << scientific<< result_convergent
      << ", error_convergent = " << setprecision(digits10) << scientific<< error_convergent << endl;
    
    if (beta == 1) {
      // Zolotarev Theorem 2.5.3, asymptotic for small x
      myFloat psi = fabs(1-alpha) * pow(xB/alpha, alpha/(alpha-1));
      myFloat fac = exp(-psi)/sqrt(2*pi*alpha*psi);
      if (verbose > 1)
        cout << "psi = " << psi << endl
        << "fac = " << fac << endl;
      myFloat term = fac;
      if (verbose > 1)
        cout << "n = " << 0 << ", term = " << setprecision(digits10) << scientific<< term << endl;
      myFloat old_term = term;
      result_asymptotic = term;
      myFloat alpha_star = alpha;
      myFloat alpha_psi_n = 1;
      for (int n = 1; n<Q_cdf.size(); ++n) {
        alpha_psi_n /= (alpha_star*psi);
        term = fac * Q_cdf.at(n) * alpha_psi_n;
        if (fabs(term) > fabs(old_term)) break;
        if (verbose > 1)
          cout << "n = " << n << ", term = " << setprecision(digits10) << scientific<< term << endl;
        n_asymptotic = n;
        result_asymptotic += term;
        old_term = term;
      }
      error_asymptotic = fabs(term) + fabs(result_asymptotic-fac)*exp(-pow(psi,.25));
      result_asymptotic = (lower_tail == positive_x) ? result_asymptotic : 1 - result_asymptotic;

    } else { // beta != 1
      // Zolotarev Formula 2.5.3 asymptotic for x -> 0 adjusted to give upper tail
      n_asymptotic = 0;
      myFloat term = (beta==-1)?0:.5 * (1 + betaB);
      myFloat abs_term = 0;
      result_asymptotic = term;
      myFloat old_abs_term = fabs(term);
      myFloat abs_tail = fabs(term);
      if (verbose > 0)
        cout << "n = " << n_asymptotic
             << ", term = " << setprecision(digits10) << scientific<< result_asymptotic << endl;
      if (xB != 0 && beta != -1) {
        for (int k=1; k<=max_n_asymp; ++k) {
          abs_term = tgamma_ratio(k/alpha,myFloat(k+1))* pow(xB,k)/(pi * alpha);
          term = -abs_term * sin(pi/2 * k * (1-betaB));
          if (k > 1 && abs_term > old_abs_term) break;
          if (verbose>0)
            cout << "n = " << k
            << ", term = " << setprecision(digits10) << scientific << term << endl;
          n_asymptotic = k;
          result_asymptotic += term;
          abs_tail += fabs(term);
          if (abs_term < Machine_eps * result_asymptotic) break;
          old_abs_term=abs_term;
        }
      }
      result_asymptotic = (lower_tail==positive_x) ? 1 - result_asymptotic : result_asymptotic;
      error_asymptotic = ((betaB!=-1)?abs_term:fabs(term)) + abs_tail*Machine_eps;
    }
    if (verbose > 0)
      cout << "result_asmyptotic = " << setprecision(digits10) << scientific << result_asymptotic
      << ", error_asymptotic = " << setprecision(digits10) << scientific<< error_asymptotic << endl
      << "result_asymptotic - result_convergent = " << result_asymptotic-result_convergent << endl;
  } else if (alpha == 1) {
    n_convergent = 0;
    if (betaB == 0) {
      using namespace boost::math;
      cauchy_distribution<myFloat> dist_cauchy(0,1);
      result_convergent = (lower_tail==positive_x) ? boost::math::cdf(dist_cauchy, x_m_zet)
                                      : boost::math::cdf(complement(dist_cauchy, x_m_zet));
      error_convergent=std::numeric_limits<myFloat>::epsilon()*result_convergent;
    } else if (x_m_zeta_in == std::numeric_limits<myFloat>::infinity()) {
      result_convergent = (lower_tail)?0:1;
      error_convergent = 0;
    } else if (x_m_zeta_in == -std::numeric_limits<myFloat>::infinity()) {
      result_convergent = (lower_tail)?1:0;
      error_convergent =0;
    } else {
      // Formula from Zolotarev's Corollary 1 to Theorem 2.2.1
      // requires integrals over an infinite range
      myFloat tmp = cdf_alpha_1();
      // Needed for consistency with aymptotic formulae
      result_convergent = (lower_tail == positive_x) ? 1-tmp : tmp;
      error_convergent = cdf_alpha_1.abserr;
    }
    if (verbose>0)
      cout << "result_convergent = " << setprecision(digits10) << scientific<< result_convergent
      << ", error_convergent = " << setprecision(digits10) << scientific<< error_convergent << endl;
        
    // Zolotarev 2.5.23 is asymptotic for large x
    // unlike the convergent formula this formula requires x > 0
    if (xB == 0) {
      result_asymptotic = std::numeric_limits<myFloat>::quiet_NaN();
      error_asymptotic = std::numeric_limits<myFloat>::max();
      n_asymptotic = 0;
    } else if (betaB == -1) {
      // Zolotarev formula 2.5.20
      myFloat psi = exp(xB-1);
      myFloat fac = exp(-psi)/sqrt(2*pi*psi);
      if (verbose > 1)
        cout << "psi = " << setprecision(digits10) << scientific<< psi << endl
        << "fac = " << setprecision(digits10) << scientific<< fac << endl;
      myFloat term = fac;
      if (verbose > 1)
        cout << "n = " << 0 << ", term = " << setprecision(digits10) << scientific << term << endl;
      result_asymptotic = term;
      myFloat old_term=term;
      myFloat psi_n = 1;
      n_asymptotic = 0;
      for (n=1; n<Q_cdf.size(); ++n) {
        psi_n /= psi;
        term = fac * Q_cdf.at(n) * psi_n;
        if (fabs(term) > fabs(old_term)) {
          break;
        }
        if (verbose > 1)
          cout << "n = " << n << ", term = " << setprecision(digits10) << scientific << term << endl;
       n_asymptotic = n;
        result_asymptotic += term;
      }
      error_asymptotic = fabs(term) + fabs(result_asymptotic-fac)*exp(-pow(psi,.25));
      result_asymptotic = (lower_tail == positive_xB) ? 1-result_asymptotic : result_asymptotic;
    } else { // xB > 0 & betaB_ != -1
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
            term0 *= (1-2*((m-l)%2)) * gamma_at_integers(m-l,n);
            term0 *=  pow(betaB, m) * pow(pi/2*(1+betaB),n-m) * sin(pi/2*(n-m));
            r_l_n += term0;
            if (verbose > 2){
              if (m == l)
                cout << "r(" << l << ", " << n << "} = " << setprecision(digits10) << scientific << term0;
              else
                cout << " + " << setprecision(digits10) << scientific<< term0;
            }
          }
          if (verbose > 2) cout << endl;
          if (verbose>1)
            cout << "r(" << l << ", " << n << ") * log_x ^ l = "
            << setprecision(digits10) << scientific << r_l_n << " * "
            << setprecision(digits10) << scientific << log_x << " ^ " << l
            << " = " << setprecision(digits10) << scientific << r_l_n << " * "
            << setprecision(digits10) << scientific << pow(log_x,l)
            << " = " << setprecision(digits10) << scientific << r_l_n * pow(log_x,l) << endl;
          fac += r_l_n*pow(log_x,l);
          if (verbose>1)
            cout << "fac: " << fac << endl;
        }
        term=fac * pow(xB,-n) / (factorial<myFloat>(n)*pi);
        if (verbose>1)
          cout << "n = "<< n << ", term = " << setprecision(digits10) << scientific << term << endl;
        if (n>1 && fabs(term) > fabs(old_term)) break;
        n_asymptotic = n;
        result_asymptotic += term;
        if (verbose>0)
          cout << "n = " << n
          << ", result_asymptotic = " << setprecision(digits10) << scientific << result_asymptotic << endl;
        if (fabs(term) < fabs(result_asymptotic) * Machine_eps) {
          if (num_small_terms > 0) break;
          num_small_terms=num_small_terms+1;
        } else {
          num_small_terms=0;
          old_term=term;
        }
      }
      result_asymptotic = (lower_tail == positive_xB) ? 1-result_asymptotic : result_asymptotic;
      error_asymptotic=fabs(term) + fabs(result_asymptotic) * Machine_eps;
    } // xB !=0
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
    if (verbose>0)
      cout << "result_asymptotic = " << setprecision(digits10) << scientific << result_asymptotic
      << ", error_asymptotic = " << setprecision(digits10) << scientific << error_asymptotic << endl
      << " result_asymptotic - result_convergent = " << result_asymptotic - result_convergent << endl;
    return result;
    
  } else {
    // alpha > 1
    // Zolotarev 2.4.6
    
    //    Formula 2.4.3 is a convergent series that works well for x <=1, but
    //    but suffers from extreme rounding errors for larger x.
    n_convergent = 0;
    myFloat term{.5*(1+theta)};
    myFloat abs_term{fabs(term)};
    result_convergent = term;
    error_convergent = 0;
    if (verbose > 1)
      cout << "n = " << n_convergent
      << ", term = " << setprecision(digits10) << scientific << term << endl;
    myFloat abs_series = fabs(term);
    for (int n=1; n<=max_n_conv; ++n) {
      abs_term = tgamma_ratio(n/alpha+1,myFloat(n+1))*pow(xB,n)/(n * pi);
      term = abs_term * (1==n%2 ? -1 : 1) * sin(pi*n*rho);
      myFloat error = fabs(term) * ((n/alpha+1)*digamma(n/alpha+1) + n)*Machine_eps;
      if (verbose>1)
        cout << "n = " << n
        << ", term = " << setprecision(digits10) << scientific << term << endl;
      n_convergent = n;
      result_convergent += term;
      error_convergent += error;
      abs_series += fabs(term);
      if (!boost::math::isfinite(abs_term) || abs_term < 2*Machine_eps*fabs(result_convergent)) break;
    }
    result_convergent = (lower_tail==positive_x) ? 1-result_convergent : result_convergent;
    error_convergent += abs_term;
    if (verbose>0)
      cout << "result_convergent = " << setprecision(digits10) << scientific << result_convergent
      << ", error_convergent = " << setprecision(digits10) << scientific << error_convergent << endl;
    
    // Zolotarev Formula 2.5.4 is an asymptotic series for large x
    // for beta != -1
    if (xB == 0) {
      result_asymptotic = std::numeric_limits<myFloat>::quiet_NaN();
      error_asymptotic = std::numeric_limits<myFloat>::max();
      n_asymptotic = 0;
    } else if (beta == -1) {
      myFloat psi = fabs(1-alpha) * pow(xB/alpha, alpha/(alpha-1));
      myFloat fac = exp(-psi)/sqrt(2*pi*alpha*psi);
      if (verbose > 1)
        cout << "psi = " << setprecision(digits10) << scientific << psi << endl
        << "fac = " << setprecision(digits10) << scientific << fac << endl;
      myFloat term = fac;
      myFloat old_term = term;
      if (verbose > 1)
        cout << "n = " << 0 << ", term = " << term << endl;
      result_asymptotic = term;
      myFloat alpha_star = 1/alpha;
      myFloat alpha_psi_n = 1;
      for (int n = 1; n<Q_cdf.size(); ++n) {
        alpha_psi_n /= (alpha_star*psi);
        term = fac * Q_cdf.at(n) * alpha_psi_n;
        if (fabs(term) > fabs(old_term)) break;
        if (verbose > 1)
          cout << "n = " << n
          << ", term = " << setprecision(digits10) << scientific << term << endl;
       n_asymptotic = n;
        result_asymptotic += term;
        old_term = term;
      }
      error_asymptotic = fabs(term)+fabs(result_asymptotic-fac)*exp(-pow(psi,.25));
      result_asymptotic = (lower_tail==positive_x) ? 1 - result_asymptotic : result_asymptotic;
    } else {
      abs_term =(1/(pi))*tgamma(alpha)/tgamma<myFloat>(2)*pow(xB,-alpha);
      term=abs_term*sin(pi/2*(2-alpha)*(1+betaB));
      result_asymptotic=term;
      myFloat old_abs_term=abs_term;
      n_asymptotic=1;
      if (verbose)
        cout << "n = " << n_asymptotic
             << ", term = " << setprecision(digits10) << scientific << term << endl;
      for (int k=2; k<=max_n_asymp; ++k) {
        abs_term=(1/(pi))*tgamma_ratio(k*alpha,myFloat(k+1))* pow(xB,-alpha*k);
        term=abs_term*sin(pi*k/2*(2-alpha)*(1+betaB));
        if (abs_term > old_abs_term) break;
        n_asymptotic = k;
        if (verbose>0)
          cout <<"n = " << k
          << ", term = " << setprecision(digits10) << scientific << term << endl;
        result_asymptotic += term;
        if (abs_term < Machine_eps * fabs(result_asymptotic)) break;
        old_abs_term=abs_term;
      }
      error_asymptotic = abs_term + fabs(result_asymptotic)*Machine_eps;
      result_asymptotic = (lower_tail==positive_x) ? 1-result_asymptotic : result_asymptotic;
    } // xB != 0
    if (verbose>0)
      cout << "result_asymptotic = " << setprecision(digits10) << scientific << result_asymptotic
      << ", error_asymptotic = " << setprecision(digits10) << scientific << error_asymptotic << endl
      << "result_asymptotic - result_convergent = " << result_asymptotic-result_convergent << endl;
    
  }  // alpha > 1
  
  if (boost::math::isnan(result_convergent)
      || error_convergent > error_asymptotic) {
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
} //cdf
} // namespace stable_distribution
