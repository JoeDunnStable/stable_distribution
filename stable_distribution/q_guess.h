//
/// \file q_guess.h
/// Initial guess for quantile of standard stable distribution
/// \author Joseph Dunn
/// \copyright 2016 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#ifndef q_guess_h
#define q_guess_h

namespace stable_distribution {

/** Generate an initial guess of the quantile q of a stable distribution given value of cdf p
 * @returns the quantile
 */
double q_guess(double p,       /**< the target value of the cdf */
               double alpha,   /**< the shape parameter of the stable distribtuon */
               double beta,    /**< the skewness parameter of the stable distribution */
               int lower_tail, /**< true => return lower tail, else return upper tail */
               int log_p       /**< true if p should be interpreted as log of the cdf */
               );

} //namespace stable_distribution
#endif /* q_guess_h */
