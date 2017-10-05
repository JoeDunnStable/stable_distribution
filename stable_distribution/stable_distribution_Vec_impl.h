//
/// \file stable_distribution_Vec.cpp
/// \author Joseph Dunn
/// \copyright 2016 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include "parameter_check.h"

namespace stable_distribution {
  
using namespace adaptive_integration;
using Eigen::Matrix;
using Eigen::Dynamic;

template<typename myFloat>
Vec std_pdf(const Vec& x, const myFloat alpha, const myFloat beta, int log_flag,
            Controllers<myFloat> ctls, const int verbose) {
  StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctls, verbose);
  Vec ret(x.size());
  for (int i=0; i<x.size(); i++)
    ret(i) = std_stable_dist.pdf(x(i), log_flag, S1);
  return ret;
}
  
template<typename myFloat>
Vec pdf(const Vec& x, const myFloat alpha, const myFloat beta,
        const Vec& gamma, const Vec& delta, const int pm, const int log_flag,
             Controllers<myFloat> ctls, const int verbose){

  parameter_check("pdf", x, alpha, beta, gamma, delta, pm);

  // Switch to S1 parameterization
  Vec gamma1(gamma.size());
  Vec delta1(delta.size());
  switch_to_S1(alpha, beta, gamma, delta, pm, ctls, verbose,
                   gamma1, delta1);
  
// Shift and Scale:
  Vec x_std = (x - delta1).array()/gamma1.array();
  
  Vec ret = std_pdf(x_std, alpha, beta, log_flag, ctls, verbose);
  if (log_flag)
    return ret.array() - gamma1.array().log();
  else
    return ret.array()/gamma1.array();
}

template<typename myFloat>
Vec std_cdf(const Vec& x, const myFloat alpha, const myFloat beta, int lower_tail, int log_flag,
            Controllers<myFloat> ctls, const int verbose) {
  StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctls, verbose);
  Vec ret(x.size());
  for (int i=0; i<x.size(); i++)
    ret(i) = std_stable_dist.cdf(x(i), lower_tail, log_flag, S1);
  return ret;
}
  
template<typename myFloat>
Vec cdf(const Vec& x, const myFloat alpha, const myFloat beta,
        const Vec& gamma, const Vec& delta, const int pm, const int lower_tail, const int log_p,
             Controllers<myFloat> ctls, const int verbose) {

  parameter_check("pdf", x, alpha, beta, gamma, delta, pm);
  
  // Switch to S0 parameterization
  Vec gamma1(gamma.size());
  Vec delta1(delta.size());
  switch_to_S1(alpha, beta, gamma, delta, pm, ctls, verbose,
                   gamma1, delta1);
  
  // Shift and Scale:
  Vec x_std = (x - delta1).array()/gamma1.array();
  return std_cdf(x_std, alpha, beta, lower_tail, log_p, ctls, verbose);
}

template<typename myFloat>
Vec std_quantile(const Vec& p, const myFloat alpha, const myFloat beta,
             const int lower_tail, const int log_p,
             const myFloat dbltol, Controllers<myFloat> ctls,
             const int verbose) {
  StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctls, verbose);
  Vec ret(p.size());
  for (int i=0; i<p.size(); i++)
    ret(i)=std_stable_dist.quantile(p(i), lower_tail, log_p, dbltol, S1);
  return ret;
}

template<typename myFloat>
Vec quantile(const Vec& p, const myFloat alpha, const myFloat beta,
              const Vec& gamma, const Vec& delta, const int pm, const int lower_tail, const int log_p,
             const myFloat dbltol, Controllers<myFloat> ctls,
             const int verbose) {

  parameter_check("quantile", p, alpha, beta, gamma, delta, pm);
  
  // Switch to S1 parameterization
  Vec gamma1(gamma.size());
  Vec delta1(delta.size());
  switch_to_S1(alpha, beta, gamma, delta, pm, ctls, verbose,
                   gamma1, delta1);
  Vec ret = std_quantile(p, alpha, beta, lower_tail, log_p, dbltol, ctls, verbose);
  return gamma1.array()*ret.array() + delta1.array();
}

template<typename myFloat>
Vec std_ddx_pdf(const Vec& x, const myFloat alpha, const myFloat beta,
            Controllers<myFloat> ctls, const int verbose){
  Vec ret(x.size());
  StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctls, verbose);
  for (int i=0; i<x.size(); i++)
    ret(i)=std_stable_dist.ddx_pdf(x(i), S1);
  return ret;
}

template<typename myFloat>
Vec ddx_pdf(const Vec& x, const myFloat alpha, const myFloat beta,
            const Vec& gamma, const Vec& delta, const int pm,
                 Controllers<myFloat> ctls, const int verbose){
  parameter_check("ddx_pdf", x, alpha, beta, gamma, delta, pm);
  
  // Switch to S1 parameterization
  Vec gamma1(gamma.size());
  Vec delta1(delta.size());
  switch_to_S1(alpha, beta, gamma, delta, pm, ctls, verbose,
                   gamma1, delta1);
  
  // Shift and Scale:
  Vec x_std = (x - delta1).array()/gamma1.array();
  Vec ret = std_ddx_pdf(x_std, alpha, beta, ctls, verbose);
  return ret.array()/gamma1.array();
}

template<typename myFloat>
Vec std_random_stable(const myFloat alpha, const myFloat beta,
                  const Vec &u1, const Vec &u2) {
  size_t n = u1.size();
  Vec ret(n);
  for (size_t i = 0; i<n; ++i)
    ret(i)=random_stable(alpha, beta, u1(i), u2(i), S1);
  return ret;
}
  
template<typename myFloat>
Vec random_stable(const myFloat alpha, const myFloat beta,
                  const Vec& gamma, const Vec& delta, const int pm,
                  Controllers<myFloat> ctls, const int verbose,
                  const Vec &u1, const Vec &u2) {
  parameter_check("random_stable", u1, alpha, beta, gamma, delta, pm);
  if (u1.size() != u2.size()) {
    throw std::out_of_range("random_stable: u1.size() != u2.size()");
  }

  // Switch to S1 parameterization
  Vec gamma1(gamma.size());
  Vec delta1(delta.size());
  
  // We need a controller to calculate the mode
  switch_to_S1(alpha, beta, gamma, delta, pm, ctls, verbose,
                   gamma1, delta1);
  return gamma1.array()* std_random_stable(alpha, beta, u1, u2).array() + delta1.array();
  
}
  
} //namespace stable_distribution

