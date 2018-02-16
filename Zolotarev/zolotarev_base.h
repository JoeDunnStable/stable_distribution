/// \file zolotarev_base.h
/// Implemenation of ommon routines for standard stable distribution per Zolotarev
/// Included in zolotarev.h when LIBRARY is defined
/// \author Joseph Dunn
/// \copyright 2017 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include <iomanip>
#include <boost/math/special_functions/bernoulli.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/tools/polynomial.hpp>
#include "gamma_derivative_at_integers.h"
// This will be included at the end of zolotarev.h if LIBRARY is defined
namespace stable_distribution {

using std::endl;
using std::setw;
using std::setprecision;
using std::scientific;
using boost::math::factorial;
using boost::math::double_factorial;
using boost::math::binomial_coefficient;
using boost::math::bernoulli_b2n;
using boost::math::tools::polynomial;

template<typename myFloat> bool Zolotarev<myFloat>::initialized=false;
template<typename myFloat> myFloat Zolotarev<myFloat>::pi;
template<typename myFloat>
Array<myFloat, Dynamic, Dynamic> Zolotarev<myFloat>::gamma_at_integers;
template<typename myFloat> int Zolotarev<myFloat>::max_n_conv;
template<typename myFloat> int Zolotarev<myFloat>::max_n_asymp;
template<typename myFloat> vector<myFloat> Zolotarev<myFloat>::points;

/// constructor for zolotarev
template<typename myFloat>
Zolotarev<myFloat>::Zolotarev(myFloat alpha, myFloat beta_input,
                    IntegrationController<myFloat>* cntl,
                     int verbose,
                     int verbose_integration) :
  alpha(alpha), beta_input(beta_input), verbose(verbose),
  cdf_alpha_1(Q_integrand<myFloat>(*this), points, *cntl, verbose_integration ),
  pdf_alpha_1(q_integrand<myFloat>(*this), points, *cntl, verbose_integration),
  ddx_pdf_alpha_1(ddx_q_integrand<myFloat>(*this), points, *cntl, verbose_integration){
  if (!initialized) {
    pi = const_pi<myFloat>();
    points = {0, pi/2};
    max_n_conv = 1000;
    max_n_asymp = 50;
    gamma_at_integers = gamma_derivative_at_integers<myFloat>(max_n_asymp);
    initialized=true;
  }
  if (alpha!=1) {
    zeta=-beta_input*tan(pi*alpha/2);
    theta0_x_gt_zeta = atan(-zeta)/alpha;
  } else {
    zeta=0;
  }
  myFloat alpha_star = (alpha < 1) ? alpha : 1/alpha;
  std::vector<myFloat> bernoulli_2n;
  
  if (fabs(beta_input)==1) {
    // Stuff we may need if alpha < 1, beta=1 and x -> 0, or
    // alpha >=1 beta=-1 and x -> infinity
    // See Zolotarev's Theorems 2.5.2 and 2.5.3
    
    // The Bernoulli for 2n =  0 to max_n_asymp
    bernoulli_b2n<myFloat>(0, max_n_asymp+1, std::back_inserter(bernoulli_2n)); // Fill vector with even Bernoulli numbers.
    // Zolotarev formula 2.5.13
    std::vector<myFloat> a_n;
    a_n.push_back(0);
    for (int n = 1; n<=max_n_asymp; ++n) {
      myFloat tmp =pow(myFloat(2),2*n)*fabs(bernoulli_2n.at(n))/(2*n*factorial<myFloat>(2*n));
      if (alpha_star==1)
        tmp *= 2*n+1;
      else
        tmp *= (alpha_star*(1-pow(alpha_star,2*n))/(1-alpha_star) + 1 -pow(1-alpha_star, 2*n));
      a_n.push_back(tmp);
    }
    // Zolotarev formula 2.5.14
    std::vector<myFloat> C_n;   // the complete exponential Bell polynomial
    std::vector<myFloat> b_n;
    C_n.push_back(1);
    b_n.push_back(1);
    for (int n = 1; n<=max_n_asymp; ++n) {
      myFloat tmp{0};
      for ( int i = 0; i<n; ++i)
        tmp += binomial_coefficient<myFloat>(n-1,i) * C_n.at(n-1-i) * (factorial<myFloat>(i+1)*a_n.at(i+1));
      C_n.push_back(tmp);
      b_n.push_back(tmp/factorial<myFloat>(n));
    }
    // Zolotarev formula 2.5.16
    std::vector<polynomial<myFloat> > d_n;
    myFloat c0[] = {0};
    myFloat c1[] = {1};
    myFloat x1[] = {0,1};
    d_n.push_back(polynomial<myFloat>{c0,0});
    polynomial<myFloat> x{x1, 1};
    polynomial<myFloat> x_2n_p_2 = x * x;
    for (int n = 1; n<=max_n_asymp-1; ++n) {
      x_2n_p_2 *= x * x;
      polynomial<myFloat> tmp = (-b_n.at(n+1)/alpha_star) * x_2n_p_2;
      d_n.push_back(tmp);
    }
    std::vector<polynomial<myFloat> > C_n_poly;
    std::vector<polynomial<myFloat> > q_cdf_n;
    C_n_poly.push_back(polynomial<myFloat>{c1, 0});
    q_cdf_n.push_back(polynomial<myFloat>{c1, 0});
    for (int n = 1; n<=max_n_asymp-1; ++n) {
      polynomial<myFloat> tmp{0};
      // Recurrence relation for complete exponential Bell polynomials
      for ( int i = 0; i<n; ++i)
        tmp += binomial_coefficient<myFloat>(n-1,i) * C_n_poly.at(n-1-i) * (factorial<myFloat>(i+1)*d_n.at(i+1));
      C_n_poly.push_back(tmp);
      q_cdf_n.push_back((1/factorial<myFloat>(n))*tmp);
    }
    Q_cdf.push_back(1);
    Q_pdf.push_back(1);
    Q_ddx_pdf.push_back(1);
    if (verbose > 1)
      cout << setw(10) << "n"
      << setw(25) << "Q_cdf"
      << setw(25) << "Q_pdf"
      << setw(25) << "Q_ddx_pdf" << endl << endl;
    if (verbose > 1)
      cout << setw(10) << 0
      << setw(25) << setprecision(15) << scientific << Q_cdf.back()
      << setw(25) << setprecision(15) << scientific << Q_pdf.back()
      << setw(25) << setprecision(15) << scientific << Q_ddx_pdf.back()
      << endl;

    for (int n = 1; n<=max_n_asymp-1; ++n ) {
      myFloat tmp{0};
      for (int i = 2; i<=4*n; i+=2)
        // Zolotarev formula 2.5.8
        // the 2 * i moment of normal distriubtion is (i-1)!!
        tmp += double_factorial<myFloat>(i-1)*q_cdf_n.at(n)[i];
      Q_cdf.push_back(tmp);
      Q_pdf.push_back((alpha_star/2)*(2*n+1)*Q_cdf.at(n-1)+Q_cdf.at(n));
      Q_ddx_pdf.push_back(((alpha_star/2)*(2*n+1)-1)*Q_pdf.at(n-1) + Q_pdf.at(n));
      if (verbose > 1)
        cout << setw(10) << n
             << setw(25) << setprecision(15) << scientific << Q_cdf.back()
             << setw(25) << setprecision(15) << scientific << Q_pdf.back()
             << setw(25) << setprecision(15) << scientific << Q_ddx_pdf.back()
             << endl;
    }
  }
}

template<typename myFloat>
void Zolotarev<myFloat>::set_x_m_zeta(myFloat x_m_zeta_in) {
  if (alpha!=1){
    x_m_zet=fabs(x_m_zeta_in);
    if (x_m_zeta_in>=0){
      beta=beta_input;
      theta0 =theta0_x_gt_zeta;
      positive_x = true;
    } else {
      beta=-beta_input;
      theta0 = -theta0_x_gt_zeta;
      positive_x = false;
    }
    theta = (beta == -1 && alpha < 1) ? -1 : theta0/(pi/2);
    rho = (1+theta)/2;
    myFloat k_alpha = (alpha>1) ? alpha - 2 : alpha;
    betaB =alpha * theta/ k_alpha;
    if (verbose>0) {
      cout << "theta = " << theta << endl;
      cout << "rho = " << rho << endl;
      cout << "betaB = " << betaB << endl;
    }
    gammaB = pow(cos(alpha*theta0),-1/alpha); // = Zolotarev lambda^(1/alpha)
    xB = x_m_zet/gammaB;
  } else { // alpha = 1
    // This following is used by the integrals
    x_m_zet=fabs(x_m_zeta_in);
    if (beta_input > 0 || (beta_input == 0 && x_m_zeta_in >=0)) {
      beta = beta_input;
      x_m_zet=x_m_zeta_in;
      positive_x = true;
    } else {
      beta = -beta_input;
      x_m_zet=-x_m_zeta_in;
      positive_x = false;
    }
    vector<myFloat> grid;
    for (int i=1; i<=1000; ++i)
      grid.push_back(myFloat(i)/10);
    points.clear();
    points.push_back(0);
    for (auto u : grid)
      points.push_back((u/fabs(max<myFloat>(1,x_m_zet)))/(1+u/fabs(max<myFloat>(1,x_m_zet))));
         
    points.push_back(1);
    // This is used by the asymptotic series
    // Zolotarev 1.1.8
    gammaB=2/pi;
    xB = (x_m_zeta_in - beta_input*gammaB*log(gammaB))/gammaB;
    if (xB >= 0) {
      betaB = beta_input;
      positive_xB = true;
    } else {
      betaB = -beta_input;
      positive_xB=false;
      xB = -xB;
    }
    if (verbose>0)
      cout << "betaB = " << betaB << endl;
  }
  if (verbose>0)
    cout << "gammaB = " << gammaB << endl
    << "xB = " << xB << endl;
  
}
  
} // namespace stable_distribution

