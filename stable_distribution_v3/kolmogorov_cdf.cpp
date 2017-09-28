/// \file kolomogorov.cpp
/// \author Joseph Dunn
/// \copyright 2017 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include "kolmogorov.h"
#include <Eigen/Dense>
using Matrix = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>;

/// Returns complementary cummulative probability for D * sqrt(n)
double kolmogorov_asymptotic_cdf(double x) {
  double ret{0};
  double eps{50*std::numeric_limits<double>::epsilon()};
  if (x<=0)
    ret=1;
  else {
    for (int i=1; i<100; ++i) {
      double term = 2 * pow(-1, i-1) * exp(-2*pow(i*x,2));
      ret += term;
      if (fabs(term) < eps * fabs(ret)) break;
    }
  }
  return ret;
}

// The follow routine is from Evaluating Kolmogorov’s Distribution
//  George Marsaglia
// The Florida State University
// Wai Wan Tsang†
//  Jingbo Wang
// The University of Hong Kong
// http://www.jstatsoft.org/v08/i18/paper
// Modified by J. Dunn to use C++ and Eigen

class Matrix_Exponent {
public:
  Matrix mat;   // Should be square
  double e;
  Matrix_Exponent() : mat(), e(0) {};
};

void mPower(const Matrix_Exponent& A, Matrix_Exponent& V,int n)
{
  if (n==1) {
    V=A;
    return;
  }
  mPower(A, V, n /2);
  if(n%2==0){
    V.mat = V.mat*V.mat;
    V.e = 2 * V.e;
  } else {
    V.mat = A.mat * V.mat * V.mat;
    V.e = A.e + 2*V.e;
  }
  size_t m = A.mat.rows();
  if(V.mat(m/2, m/2)>1e140) {
    V.mat = V.mat*1e-140;
    V.e+=140;
  }
}

double kolmogorov_cdf(int n,double d)
{ int k,m,i,j,g;
  double h,s;
  //OMIT NEXT LINE IF YOU REQUIRE >7 DIGIT ACCURACY IN THE RIGHT TAIL
  // s=d*d*n; if(s>7.24||(s>3.76&&n>99)) return 1-2*exp(-(2.000071+.331/sqrt(n)+1.409/n)*s);
  k=(int)(n*d)+1;
  m=2*k-1;
  h=k-n*d;
  Matrix_Exponent H;
  H.mat.resize(m,m);
  H.e=0;
  Matrix_Exponent Q;
  for(i=0;i<m;i++)
    for(j=0;j<m;j++)
     if(i-j+1<0)
       H.mat(i,j)=0;
     else
       H.mat(i,j)=1;
  for(i=0;i<m;i++) {
    H.mat(i,0)-=pow(h,i+1);
    H.mat(m-1,i)-=pow(h,(m-i));
  }
  H.mat(m-1, 0)+=(2*h-1>0?pow(2*h-1,m):0);
  for(i=0;i<m;i++)
    for(j=0;j<m;j++)
      if(i-j+1>0)
        for(g=1;g<=i-j+1;g++)
          H.mat(i,j)/=g;
  mPower(H, Q, n);
  s=Q.mat(k-1,k-1);
  for(i=1;i<=n;i++) {
    s=s*i/n;
    if(s<1e-140) {
      s*=1e140;
      Q.e-=140;
    }
  }
  s*=pow(10.,Q.e);
  return s;
}
