/// @file stable_distribution_base.h
/// Implementation of common routines for standard stable distribution
/// Included in stable_distribution.h when LIBRARY is defined
/// \author Joseph Dunn
/// \copyright 2017, 2018 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include <limits>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <boost/math/tools/toms748_solve.hpp>
#include <boost/math/special_functions/bernoulli.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/tools/polynomial.hpp>
#include <cassert>
#include "gamma_derivative_at_integers.h"

namespace stable_distribution {

using std::string;
using std::vector;
using std::sort;
using std::setw;
using std::right;
using std::setprecision;
using std::ios;
using std::scientific;
using std::ostringstream;
using std::pair;
using std::max;
using std::min;
using boost::math::tools::toms748_solve;
using boost::math::factorial;
using boost::math::double_factorial;
using boost::math::binomial_coefficient;
using boost::math::bernoulli_b2n;
using boost::math::tools::polynomial;

/* This is an ancillary functions used by pdf for small alpha */
template<typename myFloat>
myFloat pdf_smallA(myFloat x, myFloat alpha, myFloat beta, bool log_flag=false);

template<typename myFloat> mutex StandardStableDistribution<myFloat>::stable_mutex;
template<typename myFloat> myFloat StandardStableDistribution<myFloat>::pi;
template<typename myFloat> myFloat StandardStableDistribution<myFloat>::pi2;
template<typename myFloat> myFloat StandardStableDistribution<myFloat>::eps;
template<typename myFloat> myFloat StandardStableDistribution<myFloat>::xmin;
template<typename myFloat> myFloat StandardStableDistribution<myFloat>::large_exp_arg;
template<typename myFloat> myFloat StandardStableDistribution<myFloat>::PosInf;
template<typename myFloat> myFloat StandardStableDistribution<myFloat>::NegInf;
template<typename myFloat> double StandardStableDistribution<myFloat>::threshhold_1;
template<typename myFloat> bool StandardStableDistribution<myFloat>::initialized = false;
template<typename myFloat>
Array<myFloat, Dynamic, Dynamic> StandardStableDistribution<myFloat>::gamma_at_integers;
template<typename myFloat> int StandardStableDistribution<myFloat>::max_n;

template<typename myFloat>
void StandardStableDistribution<myFloat>::initialize(){
  unique_lock<mutex> stable_lock(stable_mutex);
  if (!initialized) {
    pi = const_pi<myFloat>();
    pi2 = pi/2;
    eps = std::numeric_limits<myFloat>::epsilon();
    xmin = std::numeric_limits<myFloat>::min();
    large_exp_arg = log(std::numeric_limits<myFloat>::max());
    PosInf = std::numeric_limits<myFloat>::infinity();
    NegInf = -PosInf;
    threshhold_1 = .5;
    max_n = 50;
    gamma_at_integers = gamma_derivative_at_integers<myFloat>(max_n);
    initialized = true;
  }
}

/// consturctor
template<typename myFloat>StandardStableDistribution<myFloat>::StandardStableDistribution(
                             myFloat alpha,                  ///< [in] the structural parameter of the distribution
                             myFloat beta,                       ///< [in] the skewness parameter of the distribution
                             Controllers<myFloat> ctls,        ///< [in] reference to integration controllers
                             int verbose                        ///< [in] indicator for verbose output
  ) : alpha(alpha), alpha_m_1(alpha-1), beta_input(beta), x_input(NAN),
  x_m_zeta_input(NAN), controllers(ctls), verbose(verbose) {
  if (!initialized) initialize();
  if (fabs(alpha_m_1)>64*std::numeric_limits<myFloat>::epsilon()){
    if (fabs(alpha_m_1) > threshhold_1 * fabs(beta)) {
      zeta = -beta_input*tan(alpha*pi2);
      theta0_x_gt_zeta = atan(beta_input*tan(alpha*pi2))/alpha;
    } else {
      zeta = beta_input/tan(pi2*alpha_m_1);
      theta0_x_gt_zeta = ((zeta<0)? 1 : -1) *
      (pi2 - (pi2*alpha_m_1 + atan(1/fabs(zeta)))/alpha);
    }
    theta0_x_gt_zeta = min(pi2,max<myFloat>(-pi2,theta0_x_gt_zeta));
    cat0=1/sqrt(1+zeta*zeta);
  } else {
    zeta=0;
    this->alpha=1;
    this->alpha_m_1=0;
    c2=pi2*fabs(1/(2*beta_input));
    c_ddx=-c2*pi2/beta_input;
  }
  if (fabs(beta_input) == 1)
    Q_initializer();
}
  
template<typename myFloat>
StandardStableDistribution<myFloat>::StandardStableDistribution(
                             AlphaMinusOne<myFloat> a_m_1,      ///< [in] the structural parameter minus one
                             myFloat beta,                       ///< [in] the skewness parameter of the distribution
                             Controllers<myFloat> ctls,        ///< [in] reference to integration controller
                             int verbose                        ///< [in] indicator for verbose output
  ) : alpha(a_m_1.alpha_m_1 + 1), alpha_m_1(a_m_1.alpha_m_1), beta_input(beta), x_input(NAN),
  x_m_zeta_input(NAN), controllers(ctls), verbose(verbose) {
  if (!initialized) initialize();
  if (fabs(alpha_m_1)> 64 * std::numeric_limits<myFloat>::min()){
    if (fabs(alpha_m_1) > threshhold_1 * fabs(beta)) {
      zeta = -beta_input*tan(alpha*pi2);
      theta0_x_gt_zeta = atan(beta_input*tan(alpha*pi2))/alpha;
    } else {
      zeta = beta_input/tan(pi2*alpha_m_1);
      theta0_x_gt_zeta = ((zeta<0)? 1 : -1) *
      (pi2 - (pi2*alpha_m_1 + atan(1/fabs(zeta)))/alpha);
    }
    theta0_x_gt_zeta = min(pi2,max<myFloat>(-pi2,theta0_x_gt_zeta));
    cat0=1/sqrt(1+zeta*zeta);
  } else {
    zeta=0;
    this->alpha=1;
    this->alpha_m_1=0;
    c2=pi2*fabs(1/(2*beta_input));
    c_ddx=-c2*pi2/beta_input;
  }
  if (fabs(beta_input) == 1)
    Q_initializer();
  
}

template<typename myFloat>
void StandardStableDistribution<myFloat>::Q_initializer() {
  myFloat alpha_star = (alpha < 1) ? alpha : 1/alpha;
  std::vector<myFloat> bernoulli_2n;
  
  // Stuff we may need if alpha < 1, beta=1 and x -> 0, or
  // alpha >=1 beta=-1 and x -> infinity
  // See Zolotarev's Theorems 2.5.2 and 2.5.3
  
  // The Bernoulli for 2n =  0 to max_n_asymp
  bernoulli_b2n<myFloat>(0, max_n+1, std::back_inserter(bernoulli_2n)); // Fill vector with even Bernoulli numbers.
  // Zolotarev formula 2.5.13
  std::vector<myFloat> a_n;
  a_n.push_back(0);
  for (int n = 1; n<=max_n; ++n) {
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
  for (int n = 1; n<=max_n; ++n) {
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
  for (int n = 1; n<=max_n-1; ++n) {
    x_2n_p_2 *= x * x;
    polynomial<myFloat> tmp = (-b_n.at(n+1)/alpha_star) * x_2n_p_2;
    d_n.push_back(tmp);
  }
  std::vector<polynomial<myFloat> > C_n_poly;
  std::vector<polynomial<myFloat> > q_cdf_n;
  C_n_poly.push_back(polynomial<myFloat>{c1, 0});
  q_cdf_n.push_back(polynomial<myFloat>{c1, 0});
  for (int n = 1; n<=max_n-1; ++n) {
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
    << setw(fmt.width) << "Q_cdf"
    << setw(fmt.width) << "Q_pdf"
    << setw(fmt.width) << "Q_ddx_pdf" << endl << endl;
  if (verbose > 1)
    cout << setw(10) << 0
    << fmt << Q_cdf.back()
    << fmt << Q_pdf.back()
    << fmt<< Q_ddx_pdf.back()
    << endl;
  
  for (int n = 1; n<=max_n-1; ++n ) {
    myFloat tmp{0};
    for (int i = 2; i<=4*n; i+=2)
      // Zolotarev formula 2.5.8
      // the 2 * i moment of normal distriubtion is (i-1)!!
      tmp += double_factorial<myFloat>(i-1)*q_cdf_n.at(n)[i];
    Q_cdf.push_back(tmp);
    Q_pdf.push_back((alpha_star/2)*(2*n-1)*Q_cdf.at(n-1)+Q_cdf.at(n));
    Q_ddx_pdf.push_back(((alpha_star/2)*(2*n-1-2/alpha))*Q_pdf.at(n-1) + Q_pdf.at(n));
/*
    Q_ddx_pdf.push_back(((alpha_star/2)*(2*n-1)-1)*Q_pdf.at(n-1) + Q_pdf.at(n));
*/
    if (verbose > 1)
      cout << setw(10) << n
      << fmt << Q_cdf.back()
      << fmt << Q_pdf.back()
      << fmt << Q_ddx_pdf.back()
      << endl;
  }
} // Q_initializer
  
template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::g_l(myFloat th_l) const{
  // Similar to g except the input variable is th_l = th+theta0
  if ((alpha < 1 && th_l<4*xmin) || (alpha > 1 && th_l>th_span*(1-4*eps))) {
    myFloat g0 = 0;
    if ((alpha<1 && beta == 1) || (alpha > 1 && beta == -1))
      g0 = pow(cat0*pow(x_m_zet/alpha,alpha),(1/(alpha_m_1)))*fabs(alpha_m_1);
    return g0;
  }
  else if ((alpha < 1 && th_l > th_span*(1-4*eps)) || (alpha >1 && th_l<4*xmin))
    return PosInf;
  else {
    myFloat cos_sin_att;
    myFloat cos_costh;
    if (add_l != 0) {
      myFloat costh = max<myFloat>(static_cast<myFloat>(0),sin(th_l-add_l));
      cos_sin_att = costh/sin(alpha*th_l);
      cos_costh = max<myFloat>(static_cast<myFloat>(0),-sin((alpha_m_1)*th_l+add_l))/costh;
    } else {
      myFloat sincth = sinc_pi(th_l);
      cos_sin_att = sincth/(alpha*sinc_pi(alpha*th_l));
      cos_costh = -(alpha_m_1)*sinc_pi((alpha_m_1)*th_l)/sincth;
    }
    myFloat pow2;
    
    if (fabs(zeta) < 1 || fabs(x_m_zeta_input+zeta) > .1 * fabs(zeta)) {
      myFloat x_cos_sin = x_m_zet*(cos_sin_att);
      myFloat pow1 = pow(x_cos_sin,alpha);
      pow2 = pow(cat0*pow1,(1/(alpha_m_1)));
    } else {
      myFloat ln_pow1 = alpha * (log(fabs(zeta))+log1p(-x_m_zeta_input/zeta-1)+log(cos_sin_att));
      myFloat ln_pow2 = (log(cat0) + ln_pow1)/(alpha_m_1);
      pow2 = exp(ln_pow2);
    }
    
    myFloat ret = pow2*cos_costh;
    return ret;
    
  }
}

template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::g_r(myFloat th_r) const{
  // Similar to g except the input variable is th_r = pi/2 - th
  if ((alpha>1 && th_r > th_span*(1-4*eps)) || (alpha<1 && th_r<4*xmin) )
    return PosInf;
  else if ((alpha>1 && th_r<4*xmin) || (alpha <1 && th_r > th_span*(1-4*eps))){
    myFloat g0 = 0;
    if ((alpha<1 && beta == 1) || (alpha > 1 && beta == -1))
      g0=pow(cat0*pow(x_m_zet/alpha,alpha),(1/(alpha_m_1)))*fabs(alpha_m_1);
    return g0;
  }
  else {
    myFloat th_l = max<myFloat>(0, th_span - th_r);
    // rounding errors cause problems when th_r is near th_span.  th_l works better.
    if (th_l< min(256*eps, th_span/2))
      return g_l(th_l);
    myFloat cos_sin_att;
    myFloat cos_costh;
    if (add_r != 0) {
      myFloat costh = max<myFloat>(0,sin(th_r));
      cos_sin_att = costh/sin(alpha*th_r+add_r);
      cos_costh = max(static_cast<myFloat>(0.),static_cast<myFloat>(sin((alpha_m_1)*th_r+add_r)))/costh;
    } else {
      myFloat sincth = sinc_pi(th_r);
      cos_sin_att = sincth/(alpha*sinc_pi(alpha*th_r));
      cos_costh = (alpha_m_1)*sinc_pi((alpha_m_1)*th_r)/sincth;
    }
    myFloat pow2;
    if (x_m_zet < 1e100) {
      myFloat pow1 = pow(x_m_zet*(cos_sin_att),alpha);
      pow2 = pow(cat0 * pow1,1/(alpha_m_1));
    } else {
      myFloat ln_pow1 = alpha*(log(x_m_zet) + log(cos_sin_att));
      pow2 = exp((log(cat0) + ln_pow1)/(alpha - 1));
    }
    myFloat ret = pow2*cos_costh;
    return ret;
    
  }
}

template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::ga1_r(myFloat u_r) const{

  myFloat h = p2b+pi2-u_r*pi2;
  myFloat h2b = h/p2b;
  myFloat tanth;
  myFloat h_tanth;
  if(u_r==2){
    tanth = NegInf;
    if (beta == 1)
      h_tanth = -1;
    else
      h_tanth = NegInf;
  } else if (u_r==0) {
    tanth = PosInf;
    if (beta == -1)
      h_tanth = -1;
    else
      h_tanth = PosInf;
  } else {
    tanth = 1/tan(pi*u_r/2);
    h_tanth = h*tanth;
  }
  myFloat exp_ea_p_h_tan_th = exp(ea+h_tanth);
  myFloat costh = sin(pi2*u_r);
  if (u_r==2) {
    if (beta>0) {
      if (beta==1)
        return exp_ea_p_h_tan_th/pi2;
      else
        return 0;
    } else {
      return PosInf;
    }
  } else if (u_r==0){
    if (beta<0){
      if (beta==-1){
        return exp_ea_p_h_tan_th/pi2;
      } else
        return 0;
    } else
      return PosInf;
  } else if (exp_ea_p_h_tan_th ==0)
    return 0;
  else {
    return h2b*exp_ea_p_h_tan_th/costh;
  }
}

template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::g_r_c(myFloat th_c) const{
  if (th_c< .5*th_min || th_c > .5 * th_max) {
    // We don't need the accuracy of th_c away from the peak
    // and dealing with the endpoints is a pain.
    return g_r(th_c-th_min);
  }
  myFloat ln_g = ln_g_th2
  +(1/alpha_m_1)*log1p(-2*pow(sin(th_c/2),2)+sin(th_c)*cot_th2_1)
  -(alpha/alpha_m_1)*log1p(-2*pow(sin(alpha*th_c/2),2)+sin(alpha*th_c)*cot_th2_2)
  +log1p(-2*pow(sin(alpha_m_1*th_c/2),2)+sin(alpha_m_1*th_c)*cot_th2_3);
  return exp(ln_g);
}

template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::g_l_c(myFloat th_c) const {
  if (th_c< .5*th_min || th_c > .5 * th_max) {
    // We don't need the accuracy of th_c away from the peak
    // and dealing with the endpoints is a pain.
    return g_l(th_c-th_min);
  }
  myFloat ln_g = ln_g_th2
  +(1/alpha_m_1)*log1p(-2*pow(sin(th_c/2),2)+sin(th_c)*cot_th2_1)
  -(alpha/alpha_m_1)*log1p(-2*pow(sin(alpha*th_c/2),2)+sin(alpha*th_c)*cot_th2_2)
  +log1p(-2*pow(sin(alpha_m_1*th_c/2),2)+sin(alpha_m_1*th_c)*cot_th2_3);
  return exp(ln_g);
}
  
template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::ga1_c(myFloat u_c) const{
  
  if (fabs(u_c) > .01) {
    // We don't need the high resolution away for u_c==0
    return ga1_r(u_c-th_min);
  }
  myFloat h = h_u2-u_c*pi2;
  myFloat del_ln_h2b = log1p(-u_c/(h_u2/pi2));
  myFloat sin_u_c = sin(u_c * pi2);
  myFloat cos_u_c = cos(u_c * pi2);
  myFloat del_h_tanth;
  
  myFloat tan_u_c;
  myFloat ret;
  if(u_c==th_max){
    if (beta == 1)
      del_h_tanth = -1-h_u2*tanth_u2;
    else
      del_h_tanth = NegInf;
  } else if (u_c == th_min) {
    if (beta == -1)
      del_h_tanth = -1-h_u2*tanth_u2;
    else
      del_h_tanth = PosInf;
  } else {
    tan_u_c = sin_u_c/cos_u_c;
    if (costh_u2*costh_u2 !=0)
      del_h_tanth = -h*tan_u_c/pow(costh_u2,2)/(1+tanth_u2*tan_u_c)-pi2*u_c*tanth_u2;
    else {
      if (h*tan_u_c == 0)
        del_h_tanth = 0;
      else if (tan_u_c < 1)
        del_h_tanth = -h*tan_u_c/pow(costh_u2,2)/(1+tanth_u2*tan_u_c)-pi2*u_c*tanth_u2;
      else
        del_h_tanth = -h/pow(costh_u2,2)/(1/tan_u_c+tanth_u2)-pi2*u_c*tanth_u2;
    }
  }
  myFloat del_ln_costh = log1p(tanth_u2*sin_u_c-2*pow(sin(pi*u_c/4),2));
  if (u_c==th_max) {
    if (beta>0) {
      if (beta==1)
        ret = exp(ln_g_u2-log(h_u2/p2b)+log(costh_u2)+del_h_tanth);
      else
        ret = 0;
    } else {
      ret = PosInf;
    }
  } else if (u_c==th_min){
    if (beta<0){
      if (beta==-1){
        ret = exp(ln_g_u2-log(h_u2/p2b)+log(costh_u2)+del_h_tanth);
      } else
        ret = 0;
    } else
      ret = PosInf;
  }
  else {
    ret = exp(ln_g_u2 + del_ln_h2b + del_h_tanth - del_ln_costh);
  }
  if (boost::math::isnan(ret)) {
    cerr << "ga1_c(alpha = " << alpha << ", beta = " << beta_input
         << ", x = " << x_m_zeta_input+zeta << ")" << endl
         << "u_c = " << u_c << endl;
    throw(std::range_error("ga1_c: Returning nan"));
  }
  return ret;
}

template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::g_a_near_1_r(myFloat th_r) const {
  myFloat tol = 4 * eps;
  if ((alpha_m_1 < 0 && fabs(th_r - th_span) < tol)
      ||(alpha_m_1 > 0 && th_r == 0)) {
    myFloat g0 = 0;
    if ((alpha_m_1 < 0  && beta == 1) || (alpha_m_1 >0 && beta == -1))
      g0=pow(cat0*pow(x_m_zet/alpha,alpha),(1/(alpha_m_1)))*fabs(alpha_m_1);
    return g0;
  } else if ((alpha_m_1 < 0 & th_r == 0)
             || (alpha_m_1 > 0 & fabs(th_r - th_span) < tol)) {
    return PosInf;
  } else {
    myFloat ea = (alpha/alpha_m_1) * log1p(-x_input/zeta);
    myFloat h_tanh = -(alpha/alpha_m_1)*log1p(max<myFloat>(-1,-2 * pow(sin((alpha_m_1*th_r+add_r)/2),2) +
                                                            1/tan(th_r)*sin(alpha_m_1*th_r + add_r)));
    myFloat ln_h2b = log(abs(zeta)*sin(alpha_m_1*th_r+add_r));
    myFloat ln_costh = log(sin(th_r));
    myFloat residual = -1/(2*alpha_m_1)*log1p(pow(zeta,-2));
    return exp(ea+h_tanh+ln_h2b-ln_costh+residual);
  }
} // g_a_near_1_r
  
template<typename myFloat>
  myFloat StandardStableDistribution<myFloat>::g_a_near_1_l(myFloat th_l) const{
  myFloat tol = 4*eps;
  if ((alpha_m_1 < 0 && th_l == 0)
      ||(alpha_m_1 > 0 && fabs(th_l - th_span) < tol)) {
    myFloat g0 = 0;
    if ((alpha_m_1 < 0  && beta == 1) || (alpha_m_1 >0 && beta == -1))
      g0=pow(cat0*pow(x_m_zet/alpha,alpha),(1/(alpha_m_1)))*fabs(alpha_m_1);
    return g0;
  } else if ((alpha_m_1 < 0 && fabs(th_l - th_span) < tol)
             ||(alpha_m_1 > 0 & th_l == 0)) {
    return PosInf;
  } else {
    myFloat ea = (alpha/alpha_m_1) * log1p(-x_input/zeta);
    myFloat h_tanh = -(alpha/alpha_m_1)*log1p(max<myFloat>(-1,-2 * pow(sin((alpha_m_1*th_l+add_l)/2),2) +
                                                                1/tan(th_l-add_l)*sin(alpha_m_1*th_l + add_l)));
    myFloat ln_h2b = log(fabs(zeta)*-sin(alpha_m_1*th_l+add_l));
    myFloat ln_costh = log(sin(th_l-add_l));
    myFloat residual = -1/(2*alpha_m_1)*log1p(pow(zeta,-2));
    return exp(ea+h_tanh+ln_h2b-ln_costh+residual);
  }
} // g_a_near_1_l

template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::dlng_dth_r(myFloat th_r){
  myFloat t1 = + cos(th_r)/(sin(th_r)*(alpha_m_1));
  myFloat t2 = - (alpha*alpha)/(alpha_m_1)*cos(alpha*th_r+add_r)/sin(alpha*th_r+add_r);
  myFloat t3 = + (alpha_m_1) * cos((alpha_m_1)*th_r+add_r)/sin((alpha_m_1)*th_r+add_r);
  return t1+t2+t3;
}

template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::dlng_dth_l(myFloat th_l) {
  myFloat t1 = cos(th_l-add_l)/sin(th_l-add_l);
  myFloat t2 = -(alpha*alpha)/(alpha_m_1)*cos(alpha*th_l)/sin(alpha*th_l);
  myFloat t3 = (-alpha_m_1) * cos((-alpha_m_1)*th_l - add_l)/sin((-alpha_m_1)*th_l - add_l);
  return t1+t2+t3;
}

template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::dlnga1_du_r(myFloat u_r) {
  myFloat t1 = -1/(1/beta + 1 - u_r);
  myFloat t2 = - pi * cos(pi2 * u_r)/sin(pi2 * u_r);
  myFloat  t3 = - pi2*pi2* (1/beta + 1 - u_r) * pow(sin(pi2 * u_r),-2);
  return t1+t2+t3;
}

template<typename myFloat>
ostream& operator<<(ostream& os, const StandardStableDistribution<myFloat>& dist) {
const int w = 21;
  os << "StandardStableDistribution: " << endl
     << setw(w) << " alpha = " << dist.fmt << dist.alpha << endl
     << setw(w) << " beta_input = " << dist.fmt << dist.beta_input << endl
     << setw(w) << " x_input = " << dist.fmt << dist.x_input << endl
     << setw(w) << " x_m_zeta_input = " << dist.fmt << dist.x_m_zeta_input << endl;
  if (!boost::math::isfinite(dist.x_input)) return os;
  os << setw(w) << " betaB = " << dist.fmt << dist.betaB << endl
     << setw(w) << " betaB_p_1 = " << dist.fmt << dist.betaB_p_1 << endl;
  if (dist.alpha != 1) {
    os << setw(w) << " theta = " << dist.fmt << dist.theta << endl
       << setw(w) << " rho = " << dist.fmt << dist.rho << endl
       << setw(w) << " one_m_betaB = " << dist.fmt << dist.one_m_betaB << endl;
  }
  os << setw(w) << " gammaB = " << dist.fmt << dist.gammaB << endl
     << setw(w) << " xB = " << dist.fmt << dist.xB << endl
     << setw(w) << " beta = " << dist.fmt << dist.beta << endl;
  switch (dist.dist_type) {
    case StandardStableDistribution<myFloat>::Cauchy :
      os << "Cauchy distribution." << endl;
      return os;
    case StandardStableDistribution<myFloat>::normal :
      os << "Normal distribution with std. dev. = sqrt(2)." << endl;
      return os;
    case StandardStableDistribution<myFloat>::fin_support :
      os << "Outside support of distribution." << endl;
      return os;
    case StandardStableDistribution<myFloat>::other :
      if (dist.use_series_small_x) return os << "Use_series_small x = true"<< endl;
      if (dist.use_series_large_x) return os << "Use_series_large x = true"<< endl;
      if (dist.alpha!=1) {
        os << setw(w) << " zeta = " << dist.fmt << dist.zeta << endl
           << setw(w) << " theta0_x_gt_zeta = " << dist.fmt << dist.theta0_x_gt_zeta << endl
           << setw(w) << " cos(alpha*theta0) = " << dist.fmt << dist.cat0 << endl
           << setw(w) << " theta0 = " << dist.fmt << dist.theta0 << endl
           << setw(w) << " x_m_zet = " << dist.fmt << dist.x_m_zet << endl
           << setw(w) << " small_x_m_zet = " << dist.fmt << dist.small_x_m_zet << endl
           << setw(w) << " add_l = " << dist.fmt << dist.add_l << endl
           << setw(w) << " add_r = " << dist.fmt << dist.add_r << endl;
      }
      /*
        if (!dist.use_f_zeta)
      */
      os << setw(w) << " th_span = " << dist.fmt << dist.th_span << endl
         << setw(w) << " th_min = " << dist.fmt << dist.th_min << endl
         << setw(w) << " th_max = " << dist.fmt << dist.th_max << endl         //Both th_l and th_r range from 0 to th_max
         << setw(w) << " c2 = " << dist.fmt << dist.c2 << endl
         << setw(w) << " c_ddx = " << dist.fmt << dist.c_ddx << endl
         << setw(w) << " fun_type = " << dist.fmt << dist.fun_type << endl;
      if (dist.fun_type == StandardStableDistribution<myFloat>::fun_g_l_c
          || dist.fun_type == StandardStableDistribution<myFloat>::fun_g_r_c) {
        os << setw(w) << " ln_g_th2 = " << dist.fmt << dist.ln_g_th2 << endl
           << setw(w) << " cot_th2_1 = " << dist.fmt << dist.cot_th2_1 << endl
           << setw(w) << " cot_th2_2 = " << dist.fmt << dist.cot_th2_2 << endl
           << setw(w) << " cot_th2_3 = " << dist.fmt << dist.cot_th2_3 << endl;
      }
      
      if (dist.alpha==1) {  // variable used when alpha = 1
        os << setw(w) << " abs_x = " << dist.fmt << dist.abs_x << endl
           << setw(w) << " i2b = " << dist.fmt << dist.i2b << endl
           << setw(w) << " p2b = " << dist.fmt << dist.p2b << endl
           << setw(w) << " ea = " << dist.fmt << dist.ea << endl;
        if (dist.fun_type==StandardStableDistribution<myFloat>::fun_ga1_c) {
          os << setw(w) << " ln_g_u2 - " << dist.fmt << dist.ln_g_u2 << endl
             << setw(w) << " costh_u2 = " << dist.fmt << dist.costh_u2 << endl
             << setw(w) << " tanth_u2 = " << dist.fmt << dist.tanth_u2 << endl
             << setw(w) << " h_u2 = " << dist.fmt << dist.h_u2 << endl;
        }
      }
      os << setw(w) << " good_theta2 = " << dist.fmt << dist.good_theta2 << endl
         << setw(w) << " g(theta2) error = " << dist.fmt <<  dist.g_theta2_error << endl
         << setw(w) << " g_dd_theta2 = " << dist.fmt << dist.g_dd_theta2 << endl;
        
      os << "g_map: " << endl << dist.fmt << "theta" << dist.fmt << "g(theta)" <<endl;
      for (typename vector<myFloat>::const_iterator ppoint=dist.points.begin(); ppoint<dist.points.end(); ppoint++) {
        os << dist.fmt << *ppoint
           << dist.fmt << dist.g(*ppoint) << ((*ppoint==dist.theta2)? " *":"") << endl;
      }
      return os;
  }   //switch on dist.dist_type
}

template<typename myFloat>
void f_of_g (myFloat *th, int n, void *ext) {
  Integral_f_of_g<myFloat> * int_f_of_g = (Integral_f_of_g<myFloat> *) ext;
  for (int i=0; i<n; i++) {
    th[i]= int_f_of_g->f_of_g(th[i]);
  }
}

template<typename myFloat>
myFloat Integral_f_of_g<myFloat>::operator() () {
  Fmt<myFloat> fmt; 
  int verbose = controllers->controller.get_verbose();
  if (verbose==4){
    cout << endl
         << "IntegrationController::integrate(f_of_g,..)" << endl;
  }
  controllers->controller.integrate(*this, std_stable_dist->points,
                        result, abserr, neval, termination_code, last);
  
  if (verbose>=3){
    myFloat rsum=0, esum=0;
    for (int i=0; i<last; i++) {
      rsum += controllers->controller.subs.at(i).r;
      esum += controllers->controller.subs.at(i).e;
    }
    
    if (termination_code > 0)
      cout << msg() << ":" << endl;
    cout << "Integral of f_of_g from theta = " << fmt << std_stable_dist->points.front()
         << " to theta = " << fmt << std_stable_dist->points.back()
         << " = " << fmt << result
         << ", with absolute error = "<< fmt  << abserr
         << ", subintervals = " << last << endl
         << "rsum = " << fmt << rsum << ", esum = " << fmt << esum << endl;
    print_subs_summary(cout, controllers->controller.subs, last, std_stable_dist->points);
  }
  if (verbose>=4){
    print_subs(cout, controllers->controller.subs, last, std_stable_dist->points);
  }
  return result;
}

template<typename myFloat>
void StandardStableDistribution<myFloat>::set_x_m_zeta(myFloat x, Parameterization pm) {
  myFloat x_m_zeta_in;
  switch (pm) {
    case S0:
      x_input=x;
      x_m_zeta_in = x-zeta;
      break;
    case S1:
      x_input = x + zeta;
      x_m_zeta_in = x;
  }

  if (x_m_zeta_in != x_m_zeta_input) {
    // Default values which will be redetermined if g_map is called
    good_theta2 = true;
    g_dd_theta2=0;
    c_g_theta2_error = 0;

    x_m_zeta_input=x_m_zeta_in;
    if (!boost::math::isfinite(x_m_zeta_in)) {
      if (verbose)
        cout << "set_x_m_zeta: x_m_zeta is not finite." << endl;
      return;
    }
    if (alpha == 1 && beta_input == 0) {
      if (verbose) {
        cout << "set_x_m_zeta: Cauchy distribution" <<endl;
      }
      dist_type=Cauchy;
      return;
    }
    if (alpha == 2) {
      if (verbose) {
        cout << "set_x_m_zeta: normal distribution with std dev of sqrt(2)" <<endl;
      }
      dist_type=normal;
      return;
    }
    if (alpha < 1 && ((beta_input==1 && x_m_zeta_in<=0) || (beta_input==-1 && x_m_zeta_in>=0))) {
      // Not used but it's nice to have them
      beta=-beta_input;
      theta0 = -theta0_x_gt_zeta;
      
      if (verbose) {
        cout << "set_x_m_zeta: outside of distribution support" << endl;
      }
      dist_type=fin_support;
      return;
    }
    dist_type=other;
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
      if (fabs(betaB + 1) > .1)
        betaB_p_1 = betaB + 1;
      else
        betaB_p_1 = atan((1+beta)*tan(alpha*pi/2)/(1-beta*pow(tan(alpha*pi/2),2)))
                    /(k_alpha*pi/2);
      if (fabs(1-betaB) > .1)
        one_m_betaB = 1-betaB;
      else
        one_m_betaB = atan((1-beta)*tan(alpha*pi/2)/(1+beta*pow(tan(alpha*pi/2),2)))
                      /(k_alpha*pi/2);
      gammaB = pow(cat0,-1/alpha); // = Zolotarev lambda^(1/alpha)
      xB = x_m_zet/gammaB;
    } else {
      gammaB = 2/pi;
      xB = (x_m_zeta_in - beta_input*gammaB*log(gammaB))/gammaB;
      if (xB >= 0) {
        betaB = beta_input;
        positive_xB = true;
      } else {
        betaB = -beta_input;
        positive_xB=false;
        xB = -xB;
      }
      betaB_p_1=betaB+1;
      one_m_betaB=1-betaB;  // Not used
      if (x_m_zeta_in >= 0) {
        beta = beta_input;
      } else {
        beta = -beta_input;
      }
    }  //alpha == 1
    myFloat eps_term_avoid = .01;
    myFloat eps_term_use = pow(eps,.2);
    if ((alpha<1 && beta == 1) || (alpha>=1 && beta == -1)) {
      myFloat xbound_avoid = alpha==1 ? -log(eps_term_avoid)+1
                           : alpha*pow(1/(eps_term_avoid*alpha*fabs(alpha_m_1)),1-1/alpha);
      avoid_series_small_x = (xB!=0 && ((alpha == 1) || (alpha <= .1))) || ((alpha<1 && xB > xbound_avoid)  || (alpha>1 && xB > eps_term_avoid));
      avoid_series_large_x = ((alpha>=1 && xB < xbound_avoid) || (alpha<1 && pow(xB,-alpha) > eps_term_avoid));
      myFloat xbound_use = alpha==1 ? -log(eps_term_use)+1
                        : alpha*pow(1/(eps_term_use*alpha*fabs(alpha_m_1)),1-1/alpha);
      use_series_small_x = (!avoid_series_small_x) && ((alpha<1 && xB < xbound_use) || (alpha>1 && xB < eps_term_use)) ;
      use_series_large_x = (!avoid_series_large_x) && ((alpha>=1 && xB > xbound_use) || (alpha<1 && pow(xB,-alpha) < eps_term_use));
    } else {
      avoid_series_small_x = (xB!=0 && ((alpha == 1) || (alpha <= .1))) || (xB > eps_term_avoid);
      avoid_series_large_x = pow(xB,-alpha) > eps_term_avoid;
      use_series_small_x = (!avoid_series_small_x) && xB < eps_term_use;
      use_series_large_x = (!avoid_series_large_x) && pow(xB,-alpha) < eps_term_use;
    }

    if (use_series_small_x || use_series_large_x) return; 
    
    // Set up the items needed for the integration
    if (alpha !=1) {
      th_min=0;
      if (fabs(alpha_m_1)> threshhold_1 * fabs(beta)) {
        th_span=(pi2)+theta0;
        add_l = theta0-(pi2);
        add_r = max<myFloat>(static_cast<myFloat>(0),pi-alpha*(th_span));
      } else {
        myFloat zeta_adj = (x_m_zeta_in>=0) ? zeta : -zeta;
        th_span = (zeta_adj<0)
                    ? pi-(pi2*alpha_m_1+atan(1/abs(zeta_adj)))/alpha
                    : (pi2*alpha_m_1+atan(1/abs(zeta_adj)))/alpha;
        add_l = (zeta_adj<0)
                    ? -(pi2*alpha_m_1 + atan(1/abs(zeta_adj)))/alpha
                    : (pi2*alpha_m_1 + atan(1/abs(zeta_adj)))/alpha-pi;
        add_r = (zeta_adj<0)
                    ? -alpha_m_1*pi+(pi2*alpha_m_1+atan(1/abs(zeta_adj)))
                    : pi-(pi2*alpha_m_1+atan(1/abs(zeta_adj)));
      }
      th_max=th_span;
      if (fabs(alpha_m_1)>= 0 || fabs(x_input)>.5*fabs(zeta)) {
        small_x_m_zet = x_m_zet<max<myFloat>(static_cast<myFloat>(1),fabs(zeta));
        if (small_x_m_zet)
          fun_type=fun_g_l;
        else
          fun_type=fun_g_r;
      } else {
        small_x_m_zet = (x_input>0 && zeta>0) || (x_input<0 && zeta<0);
        if (small_x_m_zet)
          fun_type=fun_g_a_near_1_l;
        else
          fun_type=fun_g_a_near_1_r;
      }
      if (alpha_m_1 < 0) {
        if (beta==1)
          add_l=0;
        else if (beta==-1)
          add_l=pi;
      } else if (alpha_m_1>0) {
        if (beta==-1)
          add_r=0;
      }
      if (fun_type==fun_g_l || fun_type==fun_g_r) {
        myFloat g_hi = g(th_max);
        myFloat g_lo = g(th_min);
        if (g_hi <=1 || g_lo<=1) {
          gSolve<myFloat> g_s(0., this, true);
          myFloat lower = th_min;
          myFloat upper = th_max;
          myFloat g_upper = g_hi;
          myFloat g_lower = g_lo;
          th_guess(1,lower, g_lower, upper, g_upper);
          boost::uintmax_t max_iter = 1000;
          RelativeComparisonTolerance<myFloat> rel_tol(8*eps);
          
          pair<myFloat,myFloat> th1_pair = toms748_solve(g_s, lower, upper, rel_tol, max_iter);
          myFloat ln_g_first=log(g(th1_pair.first));
          myFloat ln_g_second=log(g(th1_pair.second));
          if (verbose) {
            cout << setw(15) << "th_first = " << fmt << th1_pair.first << endl
                 << setw(15) << "ln_g_first = " << fmt << ln_g_first << endl
                 << setw(15) << "th_second = " << fmt << th1_pair.second << endl
                 << setw(15) << "ln_g_second = " << fmt << ln_g_second << endl;
          }
          myFloat th2;
          bool center_on_th2;
          if (boost::math::isfinite(ln_g_first) && boost::math::isfinite(ln_g_second)) {
            th2 = (th1_pair.first+th1_pair.second)/2;
            ln_g_th2 = log(g(th2));
            myFloat d = dlng_dth(th2);
            center_on_th2 = isnan(d) || !boost::math::isfinite(d) || fabs(th2*d) > controllers.controller.epsrel/8;
          } else if (boost::math::isfinite(ln_g_first)) {
            th2 = th1_pair.first;
            ln_g_th2 = ln_g_first;
            center_on_th2 = true;
          } else if (boost::math::isfinite(ln_g_second)){
            th2 = th1_pair.second;
            ln_g_th2 = ln_g_second;
            center_on_th2 = true;
          } else {
            th2 = (th1_pair.first+th1_pair.second)/2;
            ln_g_th2 = 0;
            center_on_th2 = true;
          }
          if (center_on_th2){
            th_min -= th2;
            th_max -= th2;
            if (fun_type==fun_g_l) {
              fun_type=fun_g_l_c;
              cot_th2_1 = 1/tan(th2-add_l);
              cot_th2_2 = 1/tan(alpha*th2);
              cot_th2_3 = 1/tan(alpha_m_1*th2+add_l);
            } else {
              fun_type=fun_g_r_c;
              cot_th2_1 = 1/tan(th2);
              cot_th2_2 = 1/tan(alpha*th2+add_r);
              cot_th2_3 = 1/tan(alpha_m_1*th2 + add_r);
            }
          }
        }
      }
      
      c2 = (alpha/(pi * fabs(alpha_m_1) * x_m_zet));
      c_ddx = c2/((alpha_m_1)*(x_m_zeta_in));
    } else { // alpha = 1
      abs_x=fabs(x_m_zeta_in);
      small_x_m_zet = abs_x < 10;
      i2b=1/(2*beta);
      p2b=pi*i2b;
      ea = -p2b*abs_x;
      th_span = 2;
      th_min=0;
      th_max=2;
      fun_type=fun_ga1_r;
      if (abs_x > 0){
        myFloat g_hi = g(th_max);
        myFloat g_lo = g(th_min);
        if (g_hi <=1 || g_lo<=1) {
          gSolve<myFloat> g_s(0., this, true);
          myFloat lower = th_min;
          myFloat upper = th_max;
          myFloat g_upper = g_hi;
          myFloat g_lower = g_lo;
          th_guess(1,lower, g_lower, upper, g_upper);
          boost::uintmax_t max_iter = 1000;
          RelativeComparisonTolerance<myFloat> rel_tol(8*eps);
          
          pair<myFloat,myFloat> ur1_pair = toms748_solve(g_s, lower, upper, rel_tol, max_iter);
          myFloat ln_g_first=log(g(ur1_pair.first));
          myFloat ln_g_second=log(g(ur1_pair.second));
          myFloat u2;
          bool center_on_u2;
          if (boost::math::isfinite(ln_g_first) && boost::math::isfinite(ln_g_second)) {
            u2 = (ur1_pair.first+ur1_pair.second)/2;
            ln_g_u2 = log(g(u2));
            myFloat d = dlng_dth(u2);
            center_on_u2 = isnan(d) || !boost::math::isfinite(d) || fabs(u2*d) > controllers.controller.epsrel/8;
          } else if (boost::math::isfinite(ln_g_first)) {
            u2 = ur1_pair.first;
            ln_g_u2 = ln_g_first;
            center_on_u2 = true;
          } else if (boost::math::isfinite(ln_g_second)){
            u2 = ur1_pair.second;
            ln_g_u2 = ln_g_second;
            center_on_u2 = true;
          } else {
            u2 = (ur1_pair.first+ur1_pair.second)/2;
            ln_g_u2 = 0;
            center_on_u2 = true;
          }
          if (center_on_u2){
            th_min -= u2;
            th_max -= u2;
            costh_u2 = sin(u2 * pi2);
            tanth_u2 = cos(u2 * pi2)/sin(u2 * pi2);
            h_u2 = p2b+pi2-u2*pi2;
            fun_type = fun_ga1_c;
          }
        }
        
      }
    }
    map_g();
  }
}

template<typename myFloat>
void StandardStableDistribution<myFloat>::map_g() {
  if (verbose > 1) cout << "map_g:" << endl;
  boost::uintmax_t max_iter;
  RelativeComparisonTolerance<myFloat> rel_tol_th2(8*eps);
  RelativeComparisonTolerance<myFloat> rel_tol(controllers.controller.epsrel);
  
  //' g() is strictly monotone -- Nolan(1997) ["3. Numerical Considerations"]
  //'     alpha >= 1  <==>  g() is falling, ie. from Inf --> 0;  otherwise growing from 0 to +Inf
  
  points.resize(0);
  points.push_back(th_min);
  points.push_back(th_max);
  if (th_min==th_max) return;
  bool do_hi=true, do_lo=true;
  myFloat g_hi = g(th_max);
  myFloat g_lo = g(th_min);
  // log(0) is broken on CPP_BIN_FLOAT
  myFloat log_g_hi = (g_hi == PosInf) ? PosInf : (g_hi==0) ? NegInf : log(g_hi);
  myFloat log_g_lo = (g_lo == PosInf) ? PosInf : (g_lo==0) ? NegInf : log(g_lo);
  myFloat ln_g_theta2{0};
  neval+=2;
  
  if (verbose>1)
    cout << "  g_hi = " << fmt << g_hi << ", g_lo = " << fmt << g_lo << endl;
  gSolve<myFloat> g_s(0., this, true);
  if (g_hi>=g_lo && g_lo>1){
    if (verbose>1)
      cout << "  Theta2 is at th_min" << endl;
    good_theta2=true;
    theta2=th_min;
    g_theta2=g_lo;
    ln_g_theta2 = log_g_lo;
    g_theta2_error = 0;
    max_iter=0;
    do_lo=false;
  } else if (g_lo>g_hi && g_hi>1){
    if (verbose>1)
      cout << "  Theta2 is at th_max" << endl;
    good_theta2=true;
    theta2=th_max;
    ln_g_theta2 = log_g_hi;
    g_theta2_error = 0;
    max_iter=0;
    do_hi=false;
  } else {
    if (verbose>1)
      cout << "  Theta2 is in the interior" << endl;
    myFloat upper =th_max;
    myFloat lower = th_min;
    myFloat g_upper = g_hi;
    myFloat g_lower = g_lo;
    th_guess(1,lower, g_lower, upper, g_upper);
    
    if (verbose>1)
      cout << "  theta 2 range from th_guess = " << fmt << lower << " to " << fmt << upper << endl;
    // toms748_solve passes lower, upper and max_iter by reference and changes them.
    max_iter = 1000;
    
    pair<myFloat,myFloat> ur1_pair = toms748_solve(g_s, lower, upper, rel_tol_th2, max_iter);
    theta2 = (ur1_pair.first+ur1_pair.second)/2;
    g_theta2 = g(theta2);
    ln_g_theta2 = log(g_theta2);
    neval+=max_iter+1;
    if (theta2!=th_min && theta2!=th_max) points.push_back(theta2);
    g_theta2_error = fabs(g(ur1_pair.second)-g(ur1_pair.first));
    myFloat cap_g_theta2_error = .01;
    if ( boost::math::isnan(g_theta2_error) || g_theta2_error > cap_g_theta2_error ){
      if (verbose>1) {
        cout << endl << "  theta2 is not good, g_theta_error: "
        << fmt << g_theta2_error;
        if (boost::math::isnan(g_theta2_error))
          cout << endl;
        else
          cout << " > " << fmt << cap_g_theta2_error << endl;
      }
      good_theta2=false;
//      sort(points.begin(),points.end());
//      return;
    } else
      good_theta2=true;
  }
  
  g_dd_theta2=dlng_dth(theta2);
  myFloat th;
  if (verbose>=2)
    cout << endl << "    theta2 = " << fmt << theta2
    << ", log(g(theta2)) = " << fmt << ln_g_theta2 << ", iterations = " << max_iter << endl
    << "    ddx_lng(theta2) = " << fmt << g_dd_theta2 <<endl;
  vector<double> ln_g_lo={-.1,-.2,-.3,-.4,-.5,-.75,
    -1,-2,-3,-4,-6,-8,-10,-12,-14,-16,
    -18,-20,-24,-28,-32,-36,-40,-50, -60,-70,-80, -100, -120, -140};
  vector<double> ln_g_hi={.1, .5, 1, 2, 3, 4, 5, 6};
  th=theta2;
  if (do_lo) {
    vector<double> *ln_g = (boost::math::isfinite(g_lo)) ? &ln_g_lo : &ln_g_hi;
    for (int j=0; j<(*ln_g).size() && (!boost::math::isfinite(g_lo) || (*ln_g)[j] > log_g_lo); j++) {
      myFloat target = ln_g_theta2+(*ln_g)[j];
      if (target >= log(g_s.g_max)) break;
      if (verbose >=2) {
        cout << "    target log(g) = " << fmt << target << endl;
      }
      g_s.set_value(target);
      pair<myFloat,myFloat> th_pair;
      max_iter=1000;
      myFloat lower = th_min;
      myFloat upper = th;
      myFloat ln_g_th = log(g(th));
      if (!boost::math::isfinite(ln_g_th)) break;
      bool bracketed = (boost::math::isfinite(g_lo)) ?(log_g_lo <= target) && (target < ln_g_th)
      :(ln_g_th <= target) && (target <= log_g_lo);
      if (!bracketed) continue;
      th_pair = toms748_solve(g_s, lower, upper, rel_tol, max_iter);
      neval+=max_iter;
      if (max_iter==1000) break;
      myFloat th_new=th_pair.second;
      if (th_new==th_min) break;
      if (rel_tol(th,th_new)) continue;
      points.push_back(th_new);
      th=th_new;
      if (verbose>=2){
        cout << "    theta = " << fmt << th
        << ", log(g(theta)) = " << fmt << log(g(th)) << ", Iterations = " << max_iter << endl;
      }
    }
    if (points.size() == 3) {
      // no intervals were added on the low side.  Add one if there's room
      if ( theta2 > 100 * std::numeric_limits<myFloat>::min()) {
        points.push_back((th_min+theta2)/2);
      }
    }
  }
  long npts_low = points.size();
  th=theta2;
  if (do_hi){
    vector<double> *ln_g = (boost::math::isfinite(g_hi)) ? &ln_g_lo : &ln_g_hi;
    for (int j=0; j<(*ln_g).size() && (!boost::math::isfinite(g_hi) || g_hi == 0 || (*ln_g)[j] > log(g_hi)); j++) {
      myFloat target = ln_g_theta2+(*ln_g)[j];
      if (target >= log(g_s.g_max)) break;
      if (verbose >=2) {
        cout << "    target log(g) = " << target << endl;
      }
      g_s.set_value(target);
      pair<myFloat,myFloat> th_pair;
      max_iter=1000;
      myFloat lower = th;
      myFloat upper = th_max;
      myFloat ln_g_th = log(g(th));
      if (!boost::math::isfinite(ln_g_th)) break;
      bool bracketed = (boost::math::isfinite(g_hi)) ? (log_g_hi <= target) && (target < ln_g_th)
      : (ln_g_th <= target) && (target <= log_g_hi);
      if (!bracketed) continue;
      th_pair = toms748_solve(g_s, lower, upper, rel_tol, max_iter);
      neval+=max_iter;
      if (max_iter==1000) break;
      myFloat th_new=th_pair.first;
      if (fabs(th_new-th_max)< 200*eps*th_max) {
        myFloat del_th{2/controllers.controller.edge_spacing()};
        while (fabs(th-th_max) >= 2*del_th*eps*th_max) {
          th_new = th_max*(1-del_th*eps);
          points.push_back(th_new);
          del_th = 2*del_th;
        }
        break;
      }
      if (rel_tol(th,th_new)) continue;
      points.push_back(th_new);
      th=th_new;
      if (verbose>=2){
        cout << "    theta = " << fmt << th
        << ", log(g(theta)) = " << fmt << log(g(th)) << ", Iterations = " << max_iter << endl;
      }
    }
    if (points.size() == npts_low) {
      // no intervals were added on the high side.  Add one if there's room
      if ( (th_max-theta2) > 100 * eps * th_max) {
        points.push_back(theta2 + (th_max - theta2)/2);
      }
    }
  }
  sort(points.begin(),points.end());
  return;
}

template<typename myFloat>
void StandardStableDistribution<myFloat>::th_guess(const myFloat &value,
                                  myFloat &lower, myFloat &g_lower,
                                  myFloat &upper, myFloat &g_upper) {
  myFloat theta_guess;
  switch(fun_type){
    case fun_g_l :
    case fun_g_a_near_1_l:
      theta_guess= x_m_zet*pow(value,(-alpha_m_1)/alpha)*pow(cat0,1/alpha)*sin(-add_l)/alpha;
      break;
    case fun_g_r :
    case fun_g_a_near_1_r:
      theta_guess=sin(add_r)/cat0*pow(value,alpha_m_1)/pow(x_m_zet,alpha);
      break;
    case fun_ga1_r:
      theta_guess = 2 * (1 + beta) / abs_x;
      break;
    case fun_g_l_c:
    case fun_g_r_c:
    case fun_ga1_c:
      theta_guess = 0;
  }
  if (lower < theta_guess && theta_guess < upper) {
    myFloat g_theta_guess = g(theta_guess);
    if (g_theta_guess >= value) {
      if (g_upper==PosInf) {
        upper = theta_guess;
        g_upper = g_theta_guess;
      } else {
        lower = theta_guess;
        g_lower = g_theta_guess;
      }
    } else {
      if (g_upper==PosInf) {
        lower = theta_guess;
        g_lower = g_theta_guess;
      } else {
        upper = theta_guess;
        g_upper = g_theta_guess;
      }
    }
  }
}

// ------------------------------------------------------------------------------

template<typename myFloat>
myFloat C_stable_tail(myFloat alpha, bool log_flag) {
    if (!(0 <= alpha && alpha <= 2)) {
        throw std::range_error("C_stable_tail: alpha is not between 0 and 2 inclusive");
    }
  myFloat r = alpha;
  if (alpha == 0)
    r = (log_flag) ? -log(static_cast<myFloat>(2)) : static_cast<myFloat>(0.5);
  else if (alpha == 2)
      r = (log_flag) ? StandardStableDistribution<myFloat>::NegInf : static_cast<myFloat>(0);
  else
      r = (log_flag) ? static_cast<myFloat>(lgamma(alpha) - log(StandardStableDistribution<myFloat>::pi) + log(sin(alpha * StandardStableDistribution<myFloat>::pi2)))
                     : static_cast<myFloat>(tgamma(alpha)/StandardStableDistribution<myFloat>::pi * sin(alpha * StandardStableDistribution<myFloat>::pi2));
  return r;
}

 
} // namespace stable_distribution
