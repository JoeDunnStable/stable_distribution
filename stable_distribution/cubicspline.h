/// \file cubicspline.h
/// General purpose cubic splines
/// \author Joseph Dunn
/// \copyright 2016, 2017, 2018 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#ifndef CUBICSPLINE_H
#define CUBICSPLINE_H

#include "myFloat.h"
#include <Eigen/Dense>
using Eigen::Sequential;
#include <algorithm>
using std::sort;

#define Mat Eigen::Matrix<myFloat, Eigen::Dynamic, Eigen::Dynamic>
#define Vec Eigen::Matrix<myFloat, Eigen::Dynamic, 1>
typedef Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> uVec;
typedef Eigen::Matrix<int, Eigen::Dynamic, 1> iVec;
typedef Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> pMat;

/** compare indices using values in another vector */
template<typename myFloat>
class CompareIndicesByAnotherVectorValues {
    const Vec& _values;     ///< the input values driving the ordering
public:
  /// Constructor which simply copies the values
  CompareIndicesByAnotherVectorValues(
                                      const Vec &values ///< [in] reference to input values
                                      )
     : _values(values) {}
    /// Given two indices, compare the values
    bool operator() (
                     const unsigned int& a,  ///< [in] reference to index for lhs
                     const unsigned int& b   ///< [in] reference to index for rhs
                     ) const { return (_values)(a) < (_values)(b); }
};

/** given vector v returns a vector of sort indices ordering v in ascending order */
template<typename myFloat>
inline iVec sort_indexes(const Vec &v /**< the vector of values to be sorted*/) {
    
    CompareIndicesByAnotherVectorValues<myFloat> comp(v);
    
    // initialize original index locations
    int n = static_cast<int>(v.size());
    iVec idx = iVec::LinSpaced(n,0,n-1);
    
    // sort indexes based on comparing values in v
    sort(idx.data(), idx.data()+idx.size(), comp);
    
    return idx;
}

/** construct and evaluate a cubic spline with specified values at given knots */
template<typename myFloat>
class CubicSpline {
private:
  Vec knots;
  Mat coefs;
public:
    /** default constructor doing nothing */
    CubicSpline() {};
  /** constructor taking knots, values and possibly derivative at end */
  CubicSpline(const Vec& x,        /**< [in] the vector of knots */
              const Vec& y,        /**< [in] the value at the corresponding knot */
              const bool endp2nd, /**< [in] flag indicating endpoints are clamped */
              const Vec& der       /**< [in] a 2 vector with the value of the derivative at endpoints */
             )
  {
    unsigned int n = static_cast<unsigned int>(x.size());
    if (n<2) {
        return;
    }
    Vec h = x.tail(n-1)-x.head(n-1);
    if ((h.array()<=0).any())
        throw std::range_error("cubic_spline; the knots must be distinct and in ascending order.\n");
    knots=x;
    if (n==2) {
      coefs.setZero(1,4);
      coefs(0,0)=y(0);
      coefs(0,1)=(x(1)!=x(0)) ? (y(1)-y(0))/(x(1)-x(0)) : 0;
      return;
    }
    Vec e(n);
    Mat A;
    A.setZero(n,n);
    Vec d(n-1);
    e(0)= 2 *h(0);
    e.segment(1,n-2) = 2*(h.segment(0,n-2) + h.segment(1,n-2));
    e(n-1)= 2*h(n - 2);
    A.diagonal()=e;
    A.diagonal(-1)=h;
    A.diagonal(1)=h;
    d = (y.tail(n-1)-y.head(n-1)).array()/h.array();

    Vec rhs(n);
    rhs.segment(1,n-2) = 3 * (d.segment(1,n-2) - d.segment(0,n - 2));
    myFloat der0 = endp2nd ? der(0) : 0;
    myFloat dern = endp2nd ? der(1) : 0;
    if (endp2nd) {
      A(0,0) = 2 * h[0];
      A(0,1) =  h(0);
      A(n-1, n-1) = 2 * h(n - 2);
      A(n-2,n-3) = h(n-2);
      rhs(0) = 3 * (d(0) - der0);
      rhs(n-1) = 3 * (dern - d(n-2));
    }
    else {
      A.topRows(1).setZero();
      A(0,0) = 1;
      A.bottomRows(1).setZero();
      A(n-1, n-1) = 1;
      rhs(0) = der0;
      rhs(n-1)= dern;
    }

    coefs.setZero(n, 4);
      coefs.col(2) = A.colPivHouseholderQr().solve(rhs);
    unsigned int m;
    for (m=0; m<n-1; m++) {
      coefs(m, 3) = (coefs(m+1,2) - coefs(m, 2))/3/h(m);
      coefs(m,1) = d(m) - h(m)/3 * (coefs(m + 1, 2) + 2 * coefs(m, 2));
      coefs(m,0) = y(m);
    }
    coefs.conservativeResize(n-1,4);
//      cout << "knots:" << endl << knots << endl;
//      cout << "y_knots" << endl << y << endl;
//      cout << "coefs: " << endl << coefs << endl;
  } //constructor

    /** evaluate the spline at points x. */
    Vec operator() (const Vec& x /**< [in] the points at which to evaluate the spline */){
    Vec ret(x.size());
    iVec s_index = sort_indexes(x);
    uVec indices;
    indices.setZero(x.size());
    iVec isgood;
    isgood.setOnes(x.size());
    unsigned int i=0, j=0;
    for (j=0; j<x.size() && (x(s_index(j))< knots(0)) ; j++)
      isgood(s_index(j)) = 0;
    while (j<x.size()){
      if (x(s_index(j))==knots(i)){
        indices(s_index(j)) = std::min(i,static_cast<unsigned int>(knots.size())-2);
        j++;
      } else if (x(s_index(j)) > knots(i)){
        i++;
        if (i>knots.size()-1){
          for (;j<x.size();j++)
            isgood(s_index(j))=0;
          break;
        }
      } else {
        indices(s_index(j)) = std::min(i-1,static_cast<unsigned int>(knots.size())-2);
        j++;
      }
    }
    for (i=0; i<x.size(); i++) {
      if (isgood(i)) {
        j=indices(i);
        myFloat dx=x(i)-knots(j);
        ret(i) = ((coefs(j,3)*dx + coefs(j,2))*dx + coefs(j,1))*dx + coefs(j,0);
      } else
        ret(i) = NAN;
    }
    return ret;
  }
  /** returns the knots used for the spline */
  Vec get_knots() {
     return knots;
  }
  /** returns the coefficient matrix of the splines */
  Mat get_coefs() {
      return coefs;
  }
  
  /** returns the number of knots in the spline */
  unsigned int get_n_knots() {
      return static_cast<unsigned int>(knots.size());
  }
};

#undef Mat
#undef Vec
#endif //CUBICSPLINE_H
