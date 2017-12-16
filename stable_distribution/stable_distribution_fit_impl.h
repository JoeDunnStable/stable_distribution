/// \file stable_distribution_fit_impl.h
/// Implementation of routines to fit stable distribution
/// Included in stable_distribution_fit.h when LIBRARY is defined
/// \author Joseph Dunn
/// \copyright 2016, 2017 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include "stable_distribution_fit.h"

#include "stable_distribution_Vec.h"
#include "parameter_check.h"

#include <chrono>

#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>

//#define BOOST_MATH_INSTRUMENT
#include <boost/math/tools/toms748_solve.hpp>

#include "neldermeadsolver.h"
#include "stable_distribution.h"

namespace boost { namespace math {
  
  template <class RealType, class Policy>
  inline RealType mycdf(const students_t_distribution<RealType, Policy>& dist, const RealType& x)
  {
    if((boost::math::isinf)(x))
    {
      if(x < 0) return 0; // -infinity
      return 1; // + infinity
    }
    return cdf(dist, x);
  }
  
  template <class RealType, class Policy>
  inline RealType mycdf(const complemented2_type<students_t_distribution<RealType, Policy>, RealType>& c)
  {
    return mycdf(c.dist, -c.param);
  }
  
} } // namespace boost::math

namespace stable_distribution {
  
using std::setw;
using std::setprecision;
using std::fixed;
using std::right;
using std::scientific;

using std::stringstream;
using std::ofstream;

using std::sort;
using std::string;

using std::chrono::high_resolution_clock;
using std::chrono::duration;
using std::min;


using boost::math::tools::toms748_solve;

using cppoptlib::Problem;
using cppoptlib::NelderMeadSolver;

using boost::math::mycdf;

template<typename myFloat>
DstableQuick<myFloat>::DstableQuick(StandardStableDistribution<myFloat> *std_stable_dist): std_stable_dist(std_stable_dist), dist_t(std_stable_dist->alpha){
  // the splines don't work well near the mode so we'll use pdf for everything near the mode
  myFloat x_mode=std_stable_dist->mode(1e-9).first;
  myFloat p_break_0 = 1e-6;
  x_break_0 = std_stable_dist->quantile(p_break_0, true, false, 1e-9);
  myFloat x_break_1 = max(x_break_0,x_mode-10);
  myFloat x_break_2 = max(x_break_1,x_mode-1);
  x_break_3 = max(x_break_2, x_mode-.1);
  //        myFloat p_break_7 = min(p_high, 1e-6);
  myFloat p_break_7 = 1e-6;
  //        myFloat x_break_7 = (p_break_7 < p_high) ? std_stable_dist.quantile(p_break_7, false, false, 1e-9)
  //                                                : x_high;
  x_break_7 = std_stable_dist->quantile(p_break_7, false, false, 1e-9);
  myFloat x_break_6 = min(x_break_7, x_mode+10);
  myFloat x_break_5 = min(x_break_6, x_mode+1);
  x_break_4 = min(x_break_5, x_mode+.1);
  
  myFloat pt_break_0 = mycdf(dist_t, x_break_0);
  myFloat pt_break_1 = mycdf(dist_t,x_break_1);
  myFloat pt_break_2 = mycdf(dist_t,x_break_2);
  myFloat pt_break_3 = mycdf(dist_t,x_break_3);
  unsigned int n_knots_1 = (pt_break_0 < pt_break_1) ? 10 : 0;
  unsigned int n_knots_2 = (pt_break_1 < pt_break_2) ? 50 : 0;
  unsigned int n_knots_3 = (pt_break_2 < pt_break_3) ? 30 : 0;
  unsigned int n_knots = n_knots_1 + n_knots_2 + n_knots_3 + 1;
  Vec pt_knot(n_knots);
  int j=0;
  myFloat del_pt = (n_knots_1 > 0) ? (pt_break_1-pt_break_0)/n_knots_1 : 0;
  pt_knot(j) = pt_break_0;
  j = 1;
  for (int i=0; i<n_knots_1; i++) {
    pt_knot(j) = pt_knot(j-1)+del_pt;
    j++;
  }
  pt_knot(j-1) = pt_break_1; // Otherwise it's slightly off because of rounding
  myFloat del_x = (n_knots_2 > 0) ? (x_break_2-x_break_1)/n_knots_2: 0.;
  myFloat x_last=x_break_1;
  for (int i=0; i<n_knots_2; i++) {
    x_last += del_x;
    pt_knot(j) = mycdf(dist_t,x_last);
    j++;
  }
  pt_knot(j-1) = pt_break_2;  //It's otherwise slightly off because of rounding
  del_x = (n_knots_3 > 0) ? (x_break_3 - x_break_2)/n_knots_3 : 0;
  x_last=x_break_2;
  for (int i=0; i<n_knots_3; i++) {
    x_last += del_x;
    pt_knot(j) = mycdf(dist_t, x_last);
    j++;
  }
  pt_knot(n_knots-1) = pt_break_3;  //It's otherwise slightly off because of rounding
  Vec y_knot(n_knots);
  j=0;
  for (unsigned int i=0; i<n_knots; ++i) {
    if (pt_knot(i) != 0 && pt_knot(i) != 1) {
      myFloat x_tmp = quantile(dist_t,pt_knot(i));
      myFloat y_tmp = std_stable_dist->pdf(x_tmp,true);
      if (boost::math::isfinite(y_tmp)) {
        pt_knot(j)=pt_knot(i);
        y_knot(j) = y_tmp - log(pdf(dist_t, x_tmp));
        ++j;
      }
    }
  }
  pt_knot.conservativeResize(j);
  y_knot.conservativeResize(j);
  if (j>=2) {
    x_break_0 = quantile(dist_t,pt_knot(0));
    x_break_3 = quantile(dist_t,pt_knot(j-1));
  } else {
    x_break_0 = x_break_3;  // Do not use splinelow
  }
  Vec der(2);  der << 0,0;
  spline_low = CubicSpline<myFloat>(pt_knot, y_knot, false, der);
  
  myFloat pt_break_7 = mycdf(complement(dist_t, x_break_7));
  myFloat pt_break_6 = mycdf(complement(dist_t,x_break_6));
  myFloat pt_break_5 = mycdf(complement(dist_t,x_break_5));
  myFloat pt_break_4 = mycdf(complement(dist_t,x_break_4));
  n_knots_1 = (pt_break_6 > pt_break_7) ? 10 : 0;
  n_knots_2 = (pt_break_5 > pt_break_6) ? 50 : 0;
  n_knots_3 = (pt_break_4 > pt_break_5) ? 30 : 0;
  n_knots = n_knots_1 + n_knots_2 + n_knots_3 + 1;
  pt_knot.resize(n_knots);
  y_knot.resize(n_knots);
  j=0;
  del_pt = (n_knots_1 > 0) ? (pt_break_6-pt_break_7)/n_knots_1 : 0;
  pt_knot(j) = pt_break_7;
  j = 1;
  for (int i=0; i<n_knots_1; i++) {
    pt_knot(j) = pt_knot(j-1)+del_pt;
    j++;
  }
  pt_knot(j-1) = pt_break_6; // Otherwise it's slightly off because of rounding
  del_x = (n_knots_2 > 0) ? (x_break_5-x_break_6)/n_knots_2: 0.;
  x_last=x_break_6;
  for (int i=0; i<n_knots_2; i++) {
    x_last += del_x;
    pt_knot(j) = mycdf(complement(dist_t,x_last));
    j++;
  }
  pt_knot(j-1) = pt_break_5;  //It's otherwise slightly off because of rounding
  del_x = (n_knots_3 > 0) ? (x_break_4 - x_break_5)/n_knots_3 : 0;
  x_last=x_break_5;
  for (int i=0; i<n_knots_3; i++) {
    x_last += del_x;
    pt_knot(j) = mycdf(complement(dist_t, x_last));
    j++;
  }
  pt_knot(n_knots-1) = pt_break_4;  //It's otherwise slightly off because of rounding
  j=0;
  for (unsigned int i=0; i<n_knots; ++i) {
    if (pt_knot(i) != 0 && pt_knot(i) != 1) {
      myFloat x_tmp = quantile(complement(dist_t,pt_knot(i)));
      myFloat y_tmp = std_stable_dist->pdf(x_tmp,true);
      if (boost::math::isfinite(y_tmp)) {
        pt_knot(j)=pt_knot(i);
        y_knot(j) = y_tmp - log(pdf(dist_t,x_tmp));
        ++j;
      };
    }
  }
  pt_knot.conservativeResize(j);
  y_knot.conservativeResize(j);
  if (j>=2) {
    x_break_7 = quantile(complement(dist_t,pt_knot(0)));
    x_break_4 = quantile(complement(dist_t,pt_knot(j-1)));
  } else {
    x_break_7 = x_break_4;
  }
  spline_high = CubicSpline<myFloat>(pt_knot, y_knot, false, der);
}

// This function returns an approximation to the log likelihood of the observations x
template<typename myFloat>
Vec DstableQuick<myFloat>::operator() (const Vec& x)
{
  unsigned int n = static_cast<unsigned int>(x.size());
  Vec ret(n);
  if (std_stable_dist->alpha < 2 - std::numeric_limits<myFloat>::epsilon()*10){
    pMat idx(sort_indexes(x));
    Vec xx = idx.inverse() * x;  //xx contains the values from x sorted in ascending order
    unsigned int n0;
    for (n0=0; xx(n0)<x_break_0; ++n0)
    {
      ret(n0) = std_stable_dist->pdf(xx(n0), true);
    }
    unsigned int n3;
    for (n3=n; n3 > 1 && xx(n3-1) > x_break_7; n3--){
      ret(n3-1) = std_stable_dist->pdf(xx(n3-1),true);
    }
    unsigned int n1;
    for (n1=n0; n1<n-1 && xx(n1)<x_break_3; ++n1);
    Vec pt_x_low(n1-n0);
    Vec ln_dt_x_low(n1-n0);
    for (unsigned int i = n0; i<n1; ++i) {
      pt_x_low(i-n0) = mycdf(dist_t,xx(i));
      ln_dt_x_low(i-n0) = log(pdf(dist_t, xx(i)));
    }
    ret.segment(n0,n1-n0) = spline_low(pt_x_low) + ln_dt_x_low;
    unsigned int n2;
    for (n2=n1; n2 < n3 && xx(n2)<x_break_4; ++n2) {
      ret(n2) = std_stable_dist->pdf(xx(n2),true);
    }
    Vec pt_x_high(n3-n2);
    Vec ln_dt_x_high(n3-n2);
    for (unsigned int i = n2; i<n3; ++i) {
      pt_x_high(i-n2) = mycdf(complement(dist_t,xx(i)));
      ln_dt_x_high(i-n2) = log(pdf(dist_t, xx(i)));
    }
    ret.segment(n2,n3-n2)= spline_high(pt_x_high) + ln_dt_x_high;
    return idx * ret;
  } else {
    for (int i=0; i<n; ++i) {
      ret(i) = -x(i)*x(i)/4 -log(2) -log(StandardStableDistribution<myFloat>::pi)/2;
    }
    return ret;
  }
}

template<typename myFloat>
myFloat DstableQuick<myFloat>::operator() (const myFloat& x)
{
  if (std_stable_dist->alpha < 2 - std::numeric_limits<myFloat>::epsilon() * 10){
    if (x<x_break_0) {
      return std_stable_dist->pdf(x, true);
    } else if (x >= x_break_7){
      return std_stable_dist->pdf(x,true);
    } else if (x >= x_break_0 && x < x_break_3) {
      Vec pt_x_low(1);
      myFloat ln_dt_x_low;
      pt_x_low(0) = mycdf(dist_t,x);
      ln_dt_x_low =log(pdf(dist_t,x));
      return spline_low(pt_x_low)(0) + ln_dt_x_low;
    } else if (x >= x_break_3 && x < x_break_4) {
      return std_stable_dist->pdf(x,true);
    } else {
      Vec pt_x_high(1);
      pt_x_high(0) = mycdf(complement(dist_t,x));
      myFloat ln_dt_x_high = log(pdf(dist_t,x));
      return spline_high(pt_x_high)(0) + ln_dt_x_high;
    }
  } else {
    return -x*x/4 -log(2) -log(StandardStableDistribution<myFloat>::pi)/2;
  }
}

template<typename myFloat>
Vec DstableQuick<myFloat>::pt(const Vec& x) {
  Vec ret(x.size());
  for (int i=0; i<x.size(); i++) {
    if (x(i) >=x_break_0 && x(i) < x_break_3) ret(i) = mycdf(dist_t,x(i));
    else if (x(i)>=x_break_4 && x(i) < x_break_7) ret(i) = mycdf(complement(dist_t,x(i)));
    else ret(i) = NAN;
  }
  return ret;
}

template<typename myFloat>
Vec pdf_quick(const Vec& x, const myFloat alpha, const myFloat beta,
              const Vec& gamma, const Vec& delta, const int pm,
              int log_flag, Controllers<myFloat> ctls, const int verbose){
  parameter_check("pdf", x, alpha, beta, gamma, delta, pm);
  
  // Switch to S0 parameterization
  Vec gamma0(gamma.size());
  Vec delta0(delta.size());
  switch_to_S0(alpha, beta, gamma, delta, pm, ctls, verbose,
                   gamma0, delta0);
  
  // Shift and Scale:
  Vec x_std = (x - delta0).array()/gamma0.array();
  StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctls, verbose);
  DstableQuick<myFloat> std_stable_dist_quick(&std_stable_dist);
  Vec ll = std_stable_dist_quick(x_std);
  if (log_flag)
    return ll.array() - gamma0.array().log();
  else
    return exp(ll.array())/gamma0.array();
}

template<typename myFloat>
myFloat capped_pdf(const Vec& y, const myFloat alpha, const myFloat beta, const myFloat gamma, const myFloat delta,
                        const bool quick, Controllers<myFloat> ctls, const int verbose) {
  StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctls, verbose);
  Vec yy = (y.array() - delta)/gamma;
  myFloat log_gamma = log(gamma);
  myFloat out =0;
  myFloat p_low = std_stable_dist.cdf(-1e300, true, true);
  myFloat p_high = std_stable_dist.cdf(1e300, false, true);
  if (quick) {
    DstableQuick<myFloat> std_stable_dist_quick(&std_stable_dist);
    int j = 0;
    for (int i=0; i<yy.size(); i++) {
      if (yy(i) > 1e300)
        out += p_high;
      else if (yy(i) < -1e300)
        out += p_low;
      else {
        yy(j) = yy(i);
        j++;
      }
    }
    yy.conservativeResize(j);
    out += std_stable_dist_quick(yy).sum() - j*log_gamma;
  } else {
    for (int i=0; i<yy.size(); i++) {
      if (yy(i) > 1e300)
        out += p_high;
      else if (yy(i) < -1e300)
        out += p_low;
      else
        out += std_stable_dist.pdf(yy(i),true) - log_gamma;
    }
  }
  out = -2 * out / y.size();
  return out;
}

template<typename myFloat>
class McCullochFit : public Problem<myFloat> {
public:
  myFloat q_kurt;
  myFloat q_skew;
  myFloat alpha_min;
  myFloat alpha_max;
  ostream *trace;
  myFloat dbltol;
  Controllers<myFloat> ctls;
  int verbose;
  myFloat value(const Vec &par) {
    myFloat alpha = par_to_alpha(par(0));
    myFloat beta = atan(par(1))/StandardStableDistribution<myFloat>::pi2;
    *trace << setw(18) << setprecision(8) << alpha
    << setw(18) << setprecision(8) << beta;
    (*trace).flush();
    StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctls, verbose);
    Vec probs(5); probs << .05, .25, .5, .75, .95;
    Vec gamma; gamma.setOnes(probs.size());
    Vec delta; delta.setZero(probs.size());
    Parameterization pm = S0;
    bool lower_tail = true;
    bool log_p = false;
    Vec qs = quantile(probs, alpha, beta, gamma, delta, pm, lower_tail, log_p, dbltol, ctls, verbose);
    myFloat qs_kurt = (qs(4) - qs(0))/(qs(3)-qs(1));
    myFloat qs_skew = (qs(4) -2 * qs(2) + qs(0))/(qs(4)-qs(0));
    *trace << setw(18) << setprecision(8) << qs_kurt
    << setw(18) << setprecision(8) << qs_skew;
    myFloat out = pow(q_kurt/qs_kurt-1,2)+100*pow(qs_skew-q_skew,2);
    *trace << setw(18) << setprecision(8) << out << endl;
    return min(out,static_cast<myFloat>(1e100));  //Optim doesn't like infinite numbers
  }
  myFloat alpha_to_par(const myFloat alpha) {
    myFloat a = min(alpha_max,max(alpha_min,alpha));
    return tan((alpha_max-a)/(alpha_max-alpha_min)*(-StandardStableDistribution<myFloat>::pi2)+(a-alpha_min)/(alpha_max-alpha_min)*StandardStableDistribution<myFloat>::pi2);
    
  }
  myFloat par_to_alpha(const myFloat par) {
    return alpha_min + (alpha_max-alpha_min)*(.5+atan(par)/StandardStableDistribution<myFloat>::pi);
  }
  myFloat norm(const Vec& par1, const Vec& par2) {
    Vec del(2);
    del(0) = par_to_alpha(par1(0)) - par_to_alpha(par2(0));
    del(1) = (atan(par1(1))-atan(par2(1)))/StandardStableDistribution<myFloat>::pi2;
    return del.array().abs().maxCoeff();
  }
  
  McCullochFit(const myFloat q_kurt, const myFloat q_skew, const myFloat alpha_min, const myFloat alpha_max, ostream *trace, const myFloat dbltol, Controllers<myFloat> ctls, const int verbose )
  : q_kurt(q_kurt), q_skew(q_skew), alpha_min(alpha_min), alpha_max(alpha_max), trace(trace), dbltol(dbltol), ctls(ctls), verbose(verbose) {};
};

template<typename myFloat>
class DunnFit : public Problem<myFloat> {
public:
  myFloat q_kurt;
  myFloat q_mode;
  myFloat skew;
  myFloat alpha_min;
  myFloat alpha_max;
  ostream *trace;
  myFloat dbltol;
  Controllers<myFloat> ctls;
  int verbose;
  myFloat value(const Vec &par) {
    myFloat alpha = par_to_alpha(par(0));
    myFloat beta = atan(par(1))/StandardStableDistribution<myFloat>::pi2;
    *trace << setw(18) << setprecision(8) << alpha
    << setw(18) << setprecision(8) << beta;
    (*trace).flush();
    StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctls, verbose);
    Vec probs(5); probs << .05, .25, .5, .75, .95;
    Vec gamma; gamma.setOnes(probs.size());
    Vec delta; delta.setZero(probs.size());
    Parameterization pm = S0;
    bool lower_tail = true;
    bool log_p = false;
    Vec qs = quantile(probs, alpha, beta, gamma, delta, pm, lower_tail, log_p, dbltol, ctls, verbose);
    myFloat mode = std_stable_dist.mode(dbltol).first;
    myFloat q_delta = .5*(qs(4)-qs(0));
    myFloat p_minus = std_stable_dist.cdf(mode-q_delta,true,false);
    myFloat p_plus = std_stable_dist.cdf(mode + q_delta, false, false);
    myFloat s_skew = (p_plus - p_minus)/(p_plus + p_minus);
    myFloat qs_kurt = (qs(4) - qs(0))/(qs(3)-qs(1));
    *trace << setw(18) << setprecision(8) << qs_kurt
    << setw(18) << setprecision(8) << s_skew;
    myFloat out = pow(q_kurt/qs_kurt-1,2)+pow(s_skew-skew,2);
    *trace << setw(18) << setprecision(8) << out << endl;
    return min(out,static_cast<myFloat>(1e100));  //Optim doesn't like infinite numbers
  }
  myFloat alpha_to_par(const myFloat alpha) {
    myFloat a = min(alpha_max,max(alpha_min,alpha));
    return tan((alpha_max-a)/(alpha_max-alpha_min)*(-StandardStableDistribution<myFloat>::pi2)+(a-alpha_min)/(alpha_max-alpha_min)*StandardStableDistribution<myFloat>::pi2);
    
  }
  myFloat par_to_alpha(const myFloat par) {
    return alpha_min + (alpha_max-alpha_min)*(.5+atan(par)/StandardStableDistribution<myFloat>::pi);
  }
  myFloat norm(const Vec& par1, const Vec& par2) {
    Vec del(2);
    del(0) = par_to_alpha(par1(0)) - par_to_alpha(par2(0));
    del(1) = (atan(par1(1))-atan(par2(1)))/StandardStableDistribution<myFloat>::pi2;
    return del.array().abs().maxCoeff();
  }
  
  DunnFit(const myFloat q_kurt, const myFloat q_mode, const myFloat skew, const myFloat alpha_min, const myFloat alpha_max,
          ostream *trace, const myFloat dbltol, Controllers<myFloat> ctls, const int verbose )
  : q_kurt(q_kurt), q_mode(q_mode), skew(skew), alpha_min(alpha_min), alpha_max(alpha_max), trace(trace), dbltol(dbltol), ctls(ctls), verbose(verbose) {};
};

template<typename myFloat>
class MLEFit : public Problem<myFloat> {
public:
  Vec y;
  myFloat alpha_min;
  myFloat alpha_max;
  bool quick;
  ostream *trace;
  Controllers<myFloat> ctls;
  int verbose;
  myFloat value(const Vec &par) {
    myFloat alpha = par_to_alpha(par(0));
    myFloat beta = atan(par(1))/StandardStableDistribution<myFloat>::pi2;
    myFloat gamma = exp(par(2));
    myFloat delta = par(3);
    *trace << setw(18) << setprecision(8) << alpha
    << setw(18) << setprecision(8) << beta
    << setw(18) << setprecision(8) << gamma
    << setw(18) << setprecision(8) << delta;
    (*trace).flush();
    myFloat out = 0;
    Vec yy = (y.array()-delta)/gamma;
    StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctls, verbose);
    myFloat p_high = std_stable_dist.cdf(1e300, false, true);
    myFloat p_low = std_stable_dist.cdf(-1e300, true, true);
    myFloat log_gamma = log(gamma);
    if (quick) {
      DstableQuick<myFloat> std_stable_dist_quick(&std_stable_dist);
      int j = 0;
      for (int i=0; i<yy.size(); i++) {
        if (yy(i) > static_cast<myFloat>(1e300))
          out += p_high;
        else if (yy(i) < static_cast<myFloat>(-1e300))
          out += p_low;
        else {
          yy(j) = yy(i);
          j++;
        }
      }
      yy.conservativeResize(j);
      out += std_stable_dist_quick(yy).sum() - j*log_gamma;
    } else {
      for (int i=0; i<yy.size(); i++) {
        if (yy(i) > 1e300)
          out += p_high;
        else if (yy(i) < -1e300)
          out += p_low;
        else
          out += std_stable_dist.pdf(yy(i),true) - log_gamma;
      }
    }
    out = -2 * out / y.size();
    *trace << setw(18) << setprecision(8) << out << endl;
    return max(min(out,static_cast<myFloat>(1e100)),static_cast<myFloat>(-1e+100));  //Optim doesn't like infinite numbers
  }
  myFloat alpha_to_par(const myFloat alpha) {
    myFloat a = min(alpha_max,max(alpha_min,alpha));
    return tan((alpha_max-a)/(alpha_max-alpha_min)*(-StandardStableDistribution<myFloat>::pi2)+(a-alpha_min)/(alpha_max-alpha_min)*StandardStableDistribution<myFloat>::pi2);
    
  }
  myFloat par_to_alpha(const myFloat par) {
    return alpha_min + (alpha_max-alpha_min)*(.5+atan(par)/StandardStableDistribution<myFloat>::pi);
  }
  myFloat norm(const Vec& par1, const Vec& par2) {
    Vec del(4);
    del(0) = par_to_alpha(par1(0)) - par_to_alpha(par2(0));
    del(1) = (atan(par1(1))-atan(par2(1)))/StandardStableDistribution<myFloat>::pi2;
    del(2) = exp(par1(2))-exp(par2(2));
    del(3) = par1(3)-par2(3);
    return del.array().abs().maxCoeff();
  }
  
  MLEFit(const Vec &y, const myFloat alpha_min, const myFloat alpha_max, bool quick, ostream *trace,
         Controllers<myFloat> ctls, const int verbose )
  : y(y), alpha_min(alpha_min), alpha_max(alpha_max), quick(quick), trace(trace), ctls(ctls), verbose(verbose) {};
};


template<typename myFloat>
class QMLEFit : public Problem<myFloat> {
public:
  Vec y;
  StandardStableDistribution<myFloat> std_stable_dist;
  bool quick;
  ostream *trace;
  DstableQuick<myFloat> std_stable_dist_quick;
  
  myFloat value(const Vec &par) {
    myFloat gamma = exp(par(0));
    myFloat delta = par(1);
    *trace << setw(18) << setprecision(8) << std_stable_dist.alpha
    << setw(18) << setprecision(8) << std_stable_dist.beta_input
    << setw(18) << setprecision(8) << gamma
    << setw(18) << setprecision(8) << delta;
    myFloat out = 0;
    Vec yy = (y.array()-delta)/gamma;
    myFloat p_high = std_stable_dist.cdf(1e300, false, true);
    myFloat p_low = std_stable_dist.cdf(-1e300, true, true);
    myFloat log_gamma = log(gamma);
    if (quick) {
      int j = 0;
      for (int i=0; i<yy.size(); i++) {
        if (yy(i) > 1e300)
          out += p_high;
        else if (yy(i) < -1e300)
          out += p_low;
        else {
          yy(j) = yy(i);
          j++;
        }
      }
      yy.conservativeResize(j);
      out += std_stable_dist_quick(yy).sum() - j*log_gamma;
    } else {
      for (int i=0; i<y.size(); i++) {
        if (yy(i) > 1e300)
          out += p_high;
        else if (yy(i) < -1e300)
          out += p_low;
        else
          out += std_stable_dist.pdf(yy(i),true) - log_gamma;
      }
    }
    out = -2 * out / y.size();
    *trace << setw(18) << setprecision(8) << out << endl;
    return max(min(out,static_cast<myFloat>(1e100)),static_cast<myFloat>(-1e+100));  //Optim doesn't like infinite numbers
  }
  myFloat norm(const Vec& par1, const Vec& par2) {
    Vec del(2);
    del(0) = exp(par1(0))-exp(par2(0));
    del(1) = par1(1)-par2(1);
    return del.array().abs().maxCoeff();
  }
  
  QMLEFit(const Vec &y, const myFloat &alpha, const myFloat &beta, bool quick, ostream *trace, Controllers<myFloat> ctls,
          const int verbose )
  : y(y), std_stable_dist(alpha, beta, ctls, verbose), quick(quick), trace(trace), std_stable_dist_quick(&std_stable_dist) {};
  
};

template<typename myFloat>
inline myFloat to_prob(myFloat p){
  return max(static_cast<myFloat>(0.),min(static_cast<myFloat>(1.),p));
}

template<typename myFloat>
Vec quantile(Vec &x, const Vec& probs)
{
  myFloat eps = 100 * std::numeric_limits<myFloat>::epsilon();
  long np = probs.size();
  if ((probs.array() < -eps || probs.array() > 1 + eps).any())
    throw std::range_error("quantile: 'probs' outside [0,1]");
  probs.array().unaryExpr(std::ptr_fun(to_prob<myFloat>));
  Vec qs(np);
  long n = x.size();
  if (n > 0 && np > 0) {
    std::sort(x.data(), x.data()+n);
    for (int j=0; j<np; ++j) {
      myFloat index = (n - 1) * probs(j);
      int lo = static_cast<int>(floor(index));
      int hi = static_cast<int>(ceil(index));
      myFloat h = index - lo;
      qs(j) = (1-h) * x(lo) + h * x(hi);
    }
    return qs;
  } else {
    throw std::range_error("quantile: Both x and prob must be of length greater than 0");
  }
}

template<typename myFloat>
myFloat p_sample(Vec x_sample, myFloat x ) {
  long n = x_sample.size();
  long n_lower = 0;
  for (int i=0; i<n; ++i) {
    if (x_sample(i) <= x) {
      ++n_lower;
    }
  }
  return static_cast<myFloat>(n_lower)/n;
}

template<typename myFloat>
myFloat mode_sample(Vec x_sample) {
  int n = static_cast<int>(x_sample.size());
  std::sort(x_sample.data(),x_sample.data()+n);
  int m = n/100;
  int i_mode = -1;
  myFloat x_delta = std::numeric_limits<myFloat>::infinity();
  for (int i=0; i+m<n; ++i) {
    if (x_sample(i+m) - x_sample(i) < x_delta) {
      i_mode = i;
      x_delta = x_sample(i+m) - x_sample(i);
    }
  }
  return x_sample(i_mode+m/2);
}

template<typename myFloat>
ostream& operator<< (ostream &os, const FitResult<myFloat> &fr) {
  os << setw(10) << fr.method
  << setw(14) << setprecision(5) << scientific << fr.alpha
  << setw(14) << setprecision(5) << fr.beta
  << setw(14) << setprecision(5) << fr.gamma
  << setw(14) << setprecision(5) << fr.delta
  << setw(14) << setprecision(5) << fr.two_ll_n
  << setw(7) << fr.n
  << setw(14) << setprecision(5) << fr.q_kurt
  << setw(14) << setprecision(5) << fr.q_skew
  << setw(14) << setprecision(5) << fr.q_scale
  << setw(14) << setprecision(5) << fr.q_location
  << setw(6) << fr.convergence
  << setw(6)  << fr.iterations
  << setw(6) << setprecision(1) << fixed << fr.cpu_time << endl;
  return os;
}

void result_heading(ostream &os) {
  os << setw(10) << right << "method"
  << setw(14) << right << "alpha"
  << setw(14) << right << "beta"
  << setw(14) << right << "gamma"
  << setw(14) << right << "delta"
  << setw(14) << right << "two_ll_n"
  << setw(7) << right << "n"
  << setw(14) << right << "q_kurt"
  << setw(14) << right << "q_skew"
  << setw(14) << right << "q_scale"
  << setw(14) << right << "q_location"
  << setw(6) << right << "conv?"
  << setw(6)  << right << "iter"
  << setw(6) << right << "time" << endl;
}

template<typename myFloat>
std::vector<FitResult<myFloat> > stable_fit(const Vec& yy, Controllers<myFloat> ctls, const myFloat dbltol,
                                  const string type,const bool quick, const int verbose) {
  int n=static_cast<int>(yy.size());
  ofstream trace("../output/stable_fit_trace.txt");
  
  // First McCulloch's method
  
  Vec probs(5); probs << .05, .25, .5, .75, .95;
  Vec y(yy);
  Vec q = quantile(y,probs);  // y is sorted as a side effect
  myFloat alpha;
  myFloat beta;
  string convergence;
  duration<double> phase1_times;
  int iterations;
  myFloat q_kurt=(q(4)-q(0))/(q(3)-q(1));
  myFloat q_skew=(q(4)+q(0)-2*q(2))/(q(4)-q(0));
  myFloat q_mode{std::numeric_limits<myFloat>::quiet_NaN()}, skew;
  NelderMeadSolver<myFloat> solver;
  typename NelderMeadSolver<myFloat>::Info &myctrl = solver.ctrl();
  myctrl.iterations = 1000;
  myctrl.obj_spread = 1e-10;
  myctrl.x_spread = 1e10;   // Effectively not used.
  const typename NelderMeadSolver<myFloat>::Info &myinfo = solver.info();
  if (q_kurt < 1000) {
    high_resolution_clock::time_point t0 = high_resolution_clock::now();
    Vec par_McCulloch(2); par_McCulloch << 1., q_skew;
    McCullochFit<myFloat> mcculloch_fit(q_kurt, q_skew, .1, 2., &trace, dbltol, ctls, verbose);
    solver.minimize(mcculloch_fit, par_McCulloch);
    iterations = static_cast<int>(myinfo.iterations);
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    phase1_times = t1 - t0;
    if (myinfo.obj_spread <= myctrl.obj_spread && myinfo.x_spread <= myctrl.x_spread)
      convergence = "True";
    else
      convergence = "False";
    alpha = mcculloch_fit.par_to_alpha(par_McCulloch(0));
    beta = atan(par_McCulloch(1))/StandardStableDistribution<myFloat>::pi2;
  } else {
    q_mode = mode_sample(y);
    myFloat p_upper = 1-p_sample(y,q_mode+.5*(q(4)-q(0)));
    myFloat p_lower = p_sample(y,q_mode-.5*((q(4)-q(0))));
    skew = (p_upper-p_lower)/(p_upper + p_lower);
    high_resolution_clock::time_point t0 = high_resolution_clock::now();
    Vec par_Dunn(2); par_Dunn << 1., skew;
    DunnFit<myFloat> dunn_fit(q_kurt, q_mode, skew, .1, 2., &trace, dbltol, ctls, verbose);
    solver.minimize(dunn_fit, par_Dunn);
    iterations = static_cast<int>(myinfo.iterations);
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    phase1_times = t1 - t0;
    if (myinfo.obj_spread <= myctrl.obj_spread && myinfo.x_spread <= myctrl.x_spread)
      convergence = "True";
    else
      convergence = "False";
    alpha = dunn_fit.par_to_alpha(par_Dunn(0));
    beta = atan(par_Dunn(1))/StandardStableDistribution<myFloat>::pi2;
  }
  Vec ps(5);
  ps << .05, .25, .5, .75, .95;
  Vec gamma0; gamma0.setOnes(ps.size());
  Vec delta0; delta0.setZero(ps.size());
  Parameterization pm = S0;
  Vec qs = quantile(ps, alpha, beta, gamma0, delta0, pm, true, false, dbltol, ctls, verbose);
  myFloat gamma = (q(3)-q(1))/(qs(3)-qs(1));
  myFloat delta;
  if (q_kurt < 1000) {
    delta = q(2) - gamma * qs(2);
  } else {
    StandardStableDistribution<myFloat> McCulloch_dist(alpha, beta, ctls, verbose);
    delta = q_mode - gamma * McCulloch_dist.mode(dbltol).first;
  }
  myFloat alpha_min = .1, alpha_max = 2.;
  MLEFit<myFloat> mle_fit(y, alpha_min, alpha_max, quick, &trace, ctls, verbose);
  Vec par_mle(4);
  par_mle << mle_fit.alpha_to_par(alpha), tan(StandardStableDistribution<myFloat>::pi2*beta), log(gamma), delta;
  
  QMLEFit<myFloat> q_mle_fit(y, alpha, beta, quick, &trace, ctls, verbose);
  Vec par_q_mle(2);
  par_q_mle << log(gamma), delta;
  
  qs= gamma * qs.array() + delta;
  FitResult<myFloat> fr_phase1((q_kurt < 1000) ? "McCulloch" : "Dunn", alpha, beta, gamma, delta, -mle_fit.value(par_mle), n, qs, convergence, iterations,
                      phase1_times.count());
  
  std::vector<FitResult<myFloat> > results;
  results.push_back(fr_phase1);
  
  iterations = 0;
  if (type=="mle") {
    duration<double> mle_times;
    high_resolution_clock::time_point t0 = high_resolution_clock::now();
    myctrl.iterations = 1000;
    solver.minimize(mle_fit, par_mle);
    iterations += myinfo.iterations;
    if (myinfo.obj_spread <= myctrl.obj_spread && myinfo.x_spread <= myctrl.x_spread)
      convergence = "True";
    else
      convergence = "False";
    alpha = mle_fit.par_to_alpha(par_mle(0));
    beta = atan(par_mle(1))/StandardStableDistribution<myFloat>::pi2;
    gamma = exp(par_mle(2));
    delta = par_mle(3);
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    mle_times = t1 - t0;
    qs=gamma*quantile(ps, alpha, beta, gamma0, delta0, pm, true, false, dbltol, ctls, verbose).array() + delta;
    FitResult<myFloat> fr_mle("mle", alpha, beta, gamma, delta,-mle_fit.value(par_mle),n,qs,convergence,iterations,
                     mle_times.count());
    results.push_back(fr_mle);
  }
  else if (type=="q_mle") {
    // alpha and beta are from McCulloch's method
    duration<double> q_mle_times;
    high_resolution_clock::time_point t0 = high_resolution_clock::now();
    myctrl.iterations=1000;
    solver.minimize(q_mle_fit, par_q_mle);
    iterations += myinfo.iterations;
    string convergence;
    if (myinfo.obj_spread <= myctrl.obj_spread && myinfo.x_spread <= myctrl.x_spread)
      convergence = "True";
    else
      convergence = "False";
    par_mle(2)=par_q_mle(0);
    gamma = exp(par_q_mle(0));
    par_mle(3) = delta = par_q_mle(1);
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    q_mle_times = t1-t0;
    qs=gamma*quantile(ps, alpha, beta, gamma0, delta0, pm, true, false, dbltol, ctls, verbose).array()+delta;
    FitResult<myFloat> fr_q_mle("q_mle", alpha, beta, gamma, delta, -mle_fit(par_mle), n, qs, convergence, iterations,
                       q_mle_times.count());
    results.push_back(fr_q_mle);
  }
  return results;
}

} //namespace stable_distribution


