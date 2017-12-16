/// \file  gamma_derivative_at_integers.h
/// Derivatives of gamma function at integers
/// \author Joseph Dunn
/// \copyright 2017 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#ifndef gamma_derivative_at_integers_h
#define gamma_derivative_at_integers_h

#include <Eigen/Dense>

using Eigen::Array;
using Eigen::Dynamic;

#include "myFloat.h"

template<typename myFloat>
Array<myFloat, Dynamic, Dynamic> gamma_derivative_at_integers(const int nn) {
  // First do x = 1
  // http://dlmf.nist.gov/5.4.E12 http://dlmf.nist.gov/5.15.E2
  Array<myFloat, Dynamic, Dynamic> polygamma(nn+1, nn+2);
  polygamma(0,1) = -const_euler<myFloat>();
  for (int n=1; n<=nn; ++n)
    polygamma(n,1) = (2*(n % 2)-1) * tgamma(myFloat(n+1)) * zeta(myFloat(n+1));
  // http://dlmf.nist.gov/5.15.E5
  for (int i=2; i<=nn+1; ++i) {
    for (int n = 0; n<=nn; ++n)
      polygamma(n,i) = polygamma(n,i-1)+(1-2*(n%2))*tgamma(myFloat(n+1))*pow(myFloat(i-1),-n-1);
  }
  //  n<-c(mpfr(0,166),n+1)
  Array<myFloat, Dynamic, Dynamic> lgamma_taylor(nn+2,nn+2);
  for (int i = 1; i<=nn+1; ++i) {
    lgamma_taylor(0,i) = lgamma(myFloat(i));
    for (int n=1; n<=nn+1; ++n)
      lgamma_taylor(n, i) = polygamma(n-1,i)/tgamma(myFloat(n+1));
  }
  Array<myFloat, Dynamic, Dynamic> gamma_taylor(nn+2, nn+2);
  gamma_taylor.setZero();
  gamma_taylor.row(0) = exp(lgamma_taylor.row(0));
  // j is power of x in lgamma_taylor its coefficients in row j of lgamma_taylor
  for (int j =1; j<=nn+1; ++j) {
    Array<myFloat, Dynamic, Dynamic> tmp = gamma_taylor; // pick up the zero order term(1) in exp
    // k i spower of x in gamma_taylor its coefficents in row k of gamma_taylor
    for (int k=1; k<=nn+1; ++k) {
      for (int i=1; i<=nn+1; ++i ){
        int kk = k - i*j;
        if (kk < 0) break;
        tmp.row(k) += gamma_taylor.row(kk) * pow(lgamma_taylor.row(j),i)/tgamma(myFloat(i+1));
      }
    }
    gamma_taylor = tmp;
  }
  for (int n=0; n<=nn+1; ++n)
    gamma_taylor.row(n) *= tgamma(myFloat(n+1));
  gamma_taylor.col(0).setConstant(NAN);
  return gamma_taylor;
}


#endif /* gamma_derivative_at_integers_h */
