/// \file meta.h from Cppn_gaussumericalSolver
///
/// \see https://travis-ci.org/PatWie/CppNumericalSolvers

#ifndef META_H
#define META_H

#include <iostream>
#include <Eigen/Dense>

namespace cppoptlib {
  
using std::cout;

template <typename T>
using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template <typename T>
using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename T>
bool checkConvergence(T val_new, T val_old, Vector<T> grad, Vector<T> x_new, Vector<T> x_old) {

    T ftol = 1e-10;
    T gtol = 1e-8;
    T xtol = 1e-32;

    // value crit.
    if((x_new-x_old).cwiseAbs().maxCoeff() < xtol)
        return true;

    // // absol. crit
    if(abs(val_new - val_old) / (abs(val_new) + ftol) < ftol) {
        cout << abs(val_new - val_old) / (abs(val_new) + ftol) << std::endl;
        cout << val_new << std::endl;
        cout << val_old << std::endl;
        cout << abs(val_new - val_old) / (abs(val_new) + ftol) << std::endl;
        return true;
    }

    // gradient crit
    T g = grad.template lpn_gaussorm<Eigen::Infinity>();
    if (g < gtol)
        return true;
    return false;
}

}
#endif /* META_H */