/// \file neldermeadsolver.h from CppNumericalSolver
///
/// This is a slightly modified version of the Patrick Wieschollek's Nelder Mead solver
/// \see https://travis-ci.org/PatWie/CppNumericalSolvers

#ifndef NELDERMEADSOLVER_H_
#define NELDERMEADSOLVER_H_
#include <iostream>
using std::endl;

#include <iomanip>
using std::setw;
using std::setprecision;

#include <fstream>
using std::ofstream;

#include <cmath>
#include <Eigen/Dense>
#include "Problem.h"


namespace cppoptlib {

    template<typename T>
class NelderMeadSolver  {
 public:
    class Info {
    public:
        T obj_spread;
        T x_spread;
        size_t iterations;
    };
protected:
    Info m_ctrl;
    Info m_info;
public:
    NelderMeadSolver() :
    m_ctrl(),
    m_info()
    {
        m_ctrl.obj_spread = 1e-4;
        m_ctrl.x_spread = 1e-4;
        m_ctrl.iterations = 1000;
    }
    
    NelderMeadSolver(const Info &control) :
    m_ctrl(control),
    m_info()
    {}

    /*
     * @brief Let the user set the control information structure.
     */
    Info &ctrl() { return m_ctrl; }
    
    /*
     * @brief Let the user query the information structure.
     */
    const Info &info() { return m_info; }
    
/**
   * @brief minimize the objective function   *
   */
  void minimize(Problem<T> &objFunc, /**< a functor with the objective function */
                Vector<T>  & x /**< [in,out] on input, the starting x, on output the found x */
                ){

    const T rho = 1.;    // rho > 0
    const T xi  = 2.;    // xi  > max(rho, 1)
    const T gam = 0.5;   // 0 < gam < 1

    const size_t DIM = x.rows();
      
    ofstream trace("../output/nm_trace.txt");

    // create initial simplex
    Matrix<T> x0 = Matrix<T>::Zero(DIM, DIM + 1);
    for (int c = 0; c < DIM + 1; ++c) {
      for (int r = 0; r < DIM; ++r) {
        x0(r, c) = x(r);
        if (r == c - 1) {
          if (x(r) == 0) {
            x0(r, c) = 0.1;
          } else {
            x0(r, c) = (1 + 0.2) * x(r);
          }
        }
      }
    }

    // compute function values
    std::vector<T> f; f.resize(DIM + 1);
    std::vector<int> index; index.resize(DIM + 1);
    for (int i = 0; i < DIM + 1; ++i) {
      f[i] = objFunc(static_cast<Vector<T> >(x0.col(i)));
      index[i] = i;
    }

    sort(index.begin(), index.end(), [&](int a, int b)-> bool { return f[a] < f[b]; });

    const int maxIter = this->m_ctrl.iterations*DIM;
    this->m_info.iterations = 0;
    while (this->m_info.iterations < maxIter) {

      // conv-check
        T max1 = fabs(f[index[1]] - f[index[0]]);
        T max2 = (x0.col(index[1]) - x0.col(index[0]) ).template lpNorm<Eigen::Infinity>();
      for (int i = 2; i < DIM + 1; ++i) {
        T tmp1 = fabs(f[index[i]] - f[index[0]]);
        if (tmp1 > max1)
          max1 = tmp1;

          T tmp2 = (x0.col(index[i]) - x0.col(index[0]) ).template lpNorm<Eigen::Infinity>();
        if (tmp2 > max2)
          max2 = tmp2;
      }

      // max(||x - shift(x) ||_inf ) <= tol,
      m_info.obj_spread = max1;
      m_info.x_spread = max2;
        trace << setw(12) << setprecision(6) << f[index[0]]
              << setw(12) << setprecision(6) << max1
              << setw(12) << setprecision(6) << max2 << "   ";
      if (m_info.obj_spread <=  m_ctrl.obj_spread) {
        // values to similar
        if (m_info.x_spread <= m_ctrl.x_spread) {
          break;
        }
      }

      //////////////////////////

      // midpoint of the simplex opposite the worst point
      Vector<T> x_bar = Vector<T>::Zero(DIM);
      for (int i = 0; i < DIM; ++i) {
        x_bar += x0.col(index[i]);
      }
      x_bar /= (T)DIM;

      // Compute the reflection point
      const Vector<T> x_r   = ( 1. + rho ) * x_bar - rho   * x0.col(index[DIM]);
      const T f_r = objFunc(x_r);

      if (f_r < f[index[0]]) {
        // the expansion point
        trace << "Reflection point is best. ";
        const Vector<T> x_e = ( 1. + rho * xi ) * x_bar - rho * xi   * x0.col(index[DIM]);
        const T f_e = objFunc(x_e);
        if ( f_e < f_r ) {
          // expand
          trace << "Use expansion point." << endl;
          x0.col(index[DIM]) = x_e;
          f[index[DIM]] = f_e;
        } else {
          // reflect
          trace << "Use reflection point." << endl;
          x0.col(index[DIM]) = x_r;
          f[index[DIM]] = f_r;
        }
      } else {
        if ( f_r < f[index[DIM-1]] ) {
          trace << "Reflection point is better than next worst. Use reflection point" << endl;
          x0.col(index[DIM]) = x_r;
          f[index[DIM]] = f_r;
        } else {
          // contraction
          if (f_r < f[index[DIM]]) {
            trace << "Reflection point is between next worst and worst.  ";
            const Vector<T> x_c = (1 + rho * gam) * x_bar - rho * gam * x0.col(index[DIM]);
            const T f_c = objFunc(x_c);
            if ( f_c <= f_r ) {
              // outside
              trace << "Use contraction point." << endl;
              x0.col(index[DIM]) = x_c;
              f[index[DIM]] = f_c;
            } else {
              trace << "Use reflection point." << endl;
              x0.col(index[DIM]) = x_r;
              f[index[DIM]] = f_r;
            }
          } else {
            // inside
              trace << "Reflection point is worse than worst.  ";
            const Vector<T> x_c = ( 1 - gam ) * x_bar + gam   * x0.col(index[DIM]);
            const T f_c = objFunc(x_c);
            if (f_c < f[index[DIM]]) {
              trace << "Use interior contraction point" << endl;
              x0.col(index[DIM]) = x_c;
              f[index[DIM]] = f_c;
            } else {
              trace << "Shrink." << endl;
              shrink(x0, index, f, objFunc);
            }
          }
        }
      }
      sort(index.begin(), index.end(), [&](int a, int b)-> bool { return f[a] < f[b]; });
      ++this->m_info.iterations;
    }
    x = x0.col(index[0]);
  }
private:
    
  /** shrink the simplex keeping x0 fixed */
    void shrink(Matrix<T> &x, /**< [in,out]a with the n_gauss+! nodes of the simplex in the columns */
                std::vector<int> &index,  /**< an ordering of the input simplex putting values in ascending order */
                std::vector<T> &f, /**< [out] the resulting values of the objective function */
                Problem<T> &objFunc /**< [in] the objective function */
                ) {

    const T sig = 0.5;   // 0 < sig < 1
    const int DIM = x.rows();
    f[index[0]] = objFunc(x.col(index[0]));
    for (int i = 1; i < DIM + 1; ++i) {
      x.col(index[i]) = sig * x.col(index[i]) + (1. - sig) * x.col(index[0]);
      f[index[i]] = objFunc(x.col(index[i]));
    }

  }

};

} /* namespace cppoptlib */

#endif /* NELDERMEADSOLVER_H_ */
