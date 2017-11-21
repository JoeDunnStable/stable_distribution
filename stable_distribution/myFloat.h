//
///  \file myFloat.h
/// Declares four different types of floats to be used:
/// 1. double, which is the default
/// 2. mpreal, when the preprocessor variable MPREAL is defined
/// 3. mprf_float, the boost multiprecision wrapper around mpfr when MPFR_FLOAT_50 is defined
/// 4. cpp_bin_float, the boost multiprecion cpp_bin_float when CPP_BIN_FLOAT is defined
///
/// \author Joseph Dunn
/// \copyright 2016, 2014 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#ifndef myFloat_h
#define myFloat_h

#undef BOOST_MATH_OVERFLOW_ERROR_POLICY
#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error

#include <Eigen/Eigenvalues>

#include <boost/multiprecision/cpp_bin_float.hpp>
using CppBinFloat = boost::multiprecision::cpp_bin_float_quad;
using BigCppBinFloat = boost::multiprecision::cpp_bin_float_50;

#include <boost/multiprecision/mpfr.hpp>

using MpfrFloat = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<50>, boost::multiprecision::et_off>;
using BigMpfrFloat = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<60>, boost::multiprecision::et_off>;

#include <mpreal.h>
using mpfr::mpreal;
using mpfr::digamma;
using mpfr::const_pi;
using mpfr::const_euler;

// the boost version of tgamma_ratio is broken for variable precision mpreal
inline mpreal tgamma_ratio(mpreal num, mpreal denom) {
    return tgamma(num)/tgamma(denom);
}

#include <boost/math/tools/real_cast.hpp>
namespace mpfr {
  template <class Policy>
  inline long long lltrunc(mpfr::mpreal const& x, const Policy& pol)
  {
    return boost::math::tools::real_cast<long long>(boost::math::trunc(x, pol));
  }

}


// These are for double

#include <cmath>
#include <boost/math/special_functions/digamma.hpp>
using boost::math::digamma;
#include <boost/math/special_functions/zeta.hpp>
using boost::math::zeta;

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/zeta.hpp>
#include <boost/math/special_functions/erf.hpp>
//using fmax = boost::multiprecision::max;
//using fmin = boost::multiprecision::min;

using boost::math::isfinite;
using boost::math::isnan;
using boost::math::isinf;
using boost::math::expm1;
using boost::math::log1p;
using boost::math::erf;
using boost::math::erfc;
using boost::math::tgamma;
using boost::math::lgamma;
using boost::math::digamma;
using boost::math::tgamma_ratio;
using boost::math::binomial_coefficient;
using boost::math::factorial;
using boost::math::zeta;

using boost::math::policies::policy;
using boost::math::policies::max_root_iterations;
typedef policy<max_root_iterations<1000> > my_erf_inv_policy;

template<typename myFloat>
myFloat pow(myFloat x, myFloat y) {
  if (x<0)
    return NAN;
  else if (x == 0 && y == 0)
    return NAN;
  else if (x == 0 && y < 0)
    return std::numeric_limits<myFloat>::infinity();
  else if (x == 0 && y>0)
    return 0;
  else
    return std::pow(x, y);
}

template<> double pow<double>(double x, double y){
  return pow(x,y);
}

template<> mpreal pow<mpreal>(mpreal x, mpreal y) {
  return pow(x,y);
}

#include <boost/math/constants/constants.hpp>
template<typename myFloat>
inline myFloat const_pi() {
  return boost::math::constants::pi<myFloat>();
}

template<>
inline mpreal const_pi<mpreal>() {return const_pi();}

template<typename myFloat>
inline myFloat const_euler() {
  return boost::math::constants::euler<myFloat>();
}

template<>
inline mpreal const_euler<mpreal>() {return const_euler();}

template<typename myFloat>
inline myFloat erf_inv(myFloat p) {
  return boost::math::erf_inv(p);
}

template<typename myFloat>
inline myFloat erfc_inv(myFloat p) {
  return boost::math::erfc_inv(p);
}

#include <boost/math/tools/roots.hpp>
using boost::math::tools::newton_raphson_iterate;

class ErfSolve {
  mpreal c;
  mpreal target;
  bool lower_tail;
public:
  ErfSolve(mpreal target, bool lower_tail)
     : target(target), lower_tail(lower_tail), c(2/sqrt(const_pi<mpreal>())) {}
  std::pair<mpreal,mpreal> operator() (mpreal x) {
    return (lower_tail) ? std::make_pair(erf(x)-target,c * exp(-x*x))
    :std::make_pair(erfc(x)-target, -c * exp(-x*x));
  }
};

template<>
inline mpreal erf_inv<mpreal> (mpreal p) {
  ErfSolve erf_s(p, true);
  mpreal guess = erf_inv(static_cast<double>(p));
  mpreal min_x = 0;
  mpreal max_x = std::numeric_limits<mpreal>::max();
  int digits = static_cast<int>(mpreal::get_default_prec()) - 4;
  boost::uintmax_t max_iter=100;
  mpreal ret= newton_raphson_iterate(erf_s, guess, min_x, max_x, digits, max_iter);
  return ret;
}

template<>
inline mpreal erfc_inv<mpreal> (mpreal p) {
  ErfSolve erfc_s(p, false);
  mpreal guess = erfc_inv(static_cast<double>(p));
  mpreal min_x = 0;
  mpreal max_x = std::numeric_limits<mpreal>::max();
  int digits = static_cast<int>(mpreal::get_default_prec()) - 4;
  boost::uintmax_t max_iter=100;
  mpreal ret = newton_raphson_iterate(erfc_s, guess, min_x, max_x, digits, max_iter);
  return ret;
}


#include <string>
using std::string;

namespace Eigen {
  
  template<> struct NumTraits<CppBinFloat> : GenericNumTraits<CppBinFloat>
  {
    static inline Real dummy_precision() { return epsilon()*static_cast<CppBinFloat>(1024); }
  };
  
  template<> struct NumTraits<MpfrFloat> : GenericNumTraits<MpfrFloat>
  {
    static inline Real dummy_precision() { return epsilon()*static_cast<MpfrFloat>(1024); }
  };
  
  template<> struct NumTraits<mpreal> : GenericNumTraits<mpreal>
  {
    static inline Real dummy_precision() { return epsilon()*static_cast<mpreal>(1024); }
  };
  
} // namespace eigen







#endif /* myFloat_h */

