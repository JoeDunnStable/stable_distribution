//
///  \file parameter_check.h
///
///  \author Joseph Dunn on 7/27/16.
///  \copyright 2016 Joseph Dunn.
//

#ifndef parameter_check_h
#define parameter_check_h

#include <Eigen/Dense>
#include "stable_distribution.h"

namespace stable_distribution {
  
using Eigen::Matrix;
  using Eigen::Dynamic;
#define Vec Matrix<myFloat, Dynamic, 1>
  
template<typename myFloat>
void inline parameter_check(const string fcn, const Vec& x, const myFloat alpha, const myFloat beta,
                            const Vec& gamma, const Vec& delta, const int pm) {
  if (alpha <= 0 || alpha > 2) {
    throw std::out_of_range(fcn + ": alpha < 0 or alpha > 2");
  }
  if (beta < -1 || beta > 1) {
    throw std::out_of_range(fcn + ": beta < -1 or beta > 1");
  }
  if ((gamma.array() < 0).any()) {
    throw std::out_of_range(fcn + ": gamma < 0");
  }
  if (pm < 0 || pm > 2) {
    throw std::out_of_range(fcn + ": pm is not 0, 1 or 2");
  }
  if (x.size() != delta.size() || x.size() != delta.size()) {
    throw std::out_of_range(fcn + ": x, gamma and delta must be the same length");
  }
}

template<typename myFloat>
void inline switch_to_S1(const myFloat alpha, const myFloat beta,
                      const Vec& gamma, const Vec& delta, const int pm,
                      Controllers<myFloat> ctls, const int verbose,
                      Vec& gamma1, Vec& delta1) {
  if (!StandardStableDistribution<myFloat>::initialized)
    StandardStableDistribution<myFloat>::initialize();
  myFloat& pi2 = StandardStableDistribution<myFloat>::pi2;
  if (pm == 0) {
    gamma1=gamma;
    if (alpha != round(alpha))
      delta1 = delta - beta * gamma * tan(pi2*alpha);
    else if (alpha == 1)
      delta1 = delta.array() - beta * gamma.array() * gamma.array().log()/pi2;
    else
      delta1 = delta;
  } else if (pm == 1) {
    gamma1 = gamma;
    delta1 = delta;
  } else { // pm=2
    gamma1 = pow(alpha, -1/alpha) * gamma;
    StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctls, verbose);
    delta1 = delta - gamma1 * std_stable_dist.mode(sqrt(std::numeric_limits<myFloat>::epsilon()),
                                                   S1).first;
  }
}

template<typename myFloat>
void inline switch_to_S0(const myFloat alpha, const myFloat beta,
                         const Vec& gamma, const Vec& delta, const int pm,
                         Controllers<myFloat> ctls, const int verbose,
                         Vec& gamma0, Vec& delta0) {
  if (!StandardStableDistribution<myFloat>::initialized)
    StandardStableDistribution<myFloat>::initialize();
  myFloat& pi2 = StandardStableDistribution<myFloat>::pi2;
  if (pm == 0) {
    gamma0 = gamma;
    delta0 = delta;
  } else if (pm == 1) {
    gamma0=gamma;
    if (alpha != round(alpha))
      delta0 = delta + beta * gamma * tan(pi2*alpha);
    else if (alpha == 1)
      delta0 = delta.array() + beta * gamma.array() * gamma.array().log()/pi2;
    else
      delta0 = delta;
  } else { // pm=2
    gamma0 = pow(alpha, -1/alpha) * gamma;
    StandardStableDistribution<myFloat> std_stable_dist(alpha, beta, ctls, verbose);
    delta0 = delta - gamma0 * std_stable_dist.mode(sqrt(std::numeric_limits<myFloat>::epsilon()),
                                                   S0).first;
  }
}
} //namespace stable_distribution

#endif /* parameter_check_h */
