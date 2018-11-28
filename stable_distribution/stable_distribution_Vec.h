//
/// \file  stable_distribution_Vec.h
/// Declaration of routines calculating 4 parameter stable distributions
/// x, gamma and delta can be vectors
/// \author Joseph Dunn
/// \copyright 2016, 2017, 2018 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#ifndef stable_distribution_Vec_h
#define stable_distribution_Vec_h

#include <vector>
#include <Eigen/Dense>
#include "stable_distribution.h"

namespace stable_distribution {

using std::vector;
using adaptive_integration::IntegrationController;
  
/// an Eigen vector containing myFloats
using Eigen::Matrix;
using Eigen::Dynamic;
#define Vec Matrix<myFloat, Dynamic, 1>

/// calculate the pdf of the standard stable distribution for vector of x's
template<typename myFloat>
Vec std_pdf(const Vec& x,                     ///< [in] the vector of points at which to calculate the pdf
            const myFloat alpha,              ///< [in] the structural parameter of the distribution
            const myFloat beta,               ///< [in] the skewness parameter of the distribution
            const Parameterization pm_std,    ///< [in] the parameterization
            const int log_flag,               ///< [in] return log of the pdf, otherwise the pdf
            Controllers<myFloat> ctls,        ///< [in,out] reference to the integration controller
            const int verbose                 ///< [in] indicator for verbose output
);

/// calculate the pdf of the stable distribution for vector of x's
template<typename myFloat>
Vec pdf(const Vec& x,                     ///< [in] the vector of points at which to calculate the pdf
             const myFloat alpha,         ///< [in] the structural parameter of the distribution
             const myFloat beta,          ///< [in] the skewness parameter of the distribution
             const Vec& gamma,            ///< [in] a vector of scale parameters for distribution
             const Vec& delta,            ///< [in] a vector of location parameters for distribution
             const int pm,                ///< [in] the parameterization, 0, 1 or 2
             const int log_flag,          ///< [in] return log of the pdf, otherwise the pdf
             Controllers<myFloat> ctls,   ///< [in,out] reference to the integration controller
             const int verbose            ///< [in] indicator for verbose output
             );

/// return the cdf of the standard stable distribution for a vector of points z
template<typename myFloat>
Vec std_cdf(const Vec& z,                     ///< [in] a vector of points at which to calculate the distribution
            const myFloat alpha,              ///< [in] the structural parameter of the distribution
            const myFloat beta,               ///< [in] the skewness parameter of the distribution
            const Parameterization pm_std,    ///< [in] the parameterization
            const int lower_tail,             ///< [in] return the lower tail, otherwise the upper tail
            const int log_p,                  ///< [in] return the log(pdf), otherwise the pdf
            Controllers<myFloat> ctls,        ///< [in,out] reference to the integration controller
            const int verbose                 ///< [in] indicator for verbose output
);

/// return the cdf of the stable distribution for a vector of points z
template<typename myFloat>
Vec cdf(const Vec& z,                     ///< [in] a vector of points at which to calculate the distribution
             const myFloat alpha,         ///< [in] the structural parameter of the distribution
             const myFloat beta,          ///< [in] the skewness parameter of the distribution
             const Vec& gamma,            ///< [in] a vector of scale parameters for distribution
             const Vec& delta,            ///< [in] a vector of location parameters for distribution
             const int pm,                ///< [in] the parameterization, 0, 1 or 2
             const int lower_tail,        ///< [in] return the lower tail, otherwise the upper tail
             const int log_p,             ///< [in] return the log(pdf), otherwise the pdf 
             Controllers<myFloat> ctls,   ///< [in,out] reference to the integration controller
             const int verbose            ///< [in] indicator for verbose output
             );

/// return the inverse of the cdf of the std. stable distribution in S1 for a vector of probabilities p
template<typename myFloat>
Vec std_quantile(const Vec& p,                     ///< [in] vector of probabilities at which to calculate the inverse
                 myFloat alpha,                    ///< [in] the structural parameter of the distribution
                 myFloat beta,                     ///< [in] the skewness parameter of the distribution
                 const Parameterization pm_std,    ///< [in] the parameterization
                 int lower_tail,                   ///< [in] return the lower tail, otherwise the upper tail
                 int log_p,                        ///< [in] p is the log(pdf), otherwise the pdf
                 const myFloat dbltol,             ///< [in] the tolerance for the inverse
                 Controllers<myFloat> ctls,        ///< [in,out] reference to the integration controller
                 const int verbose                 ///< [in] indicator for verbose output
);

/// return the inverse of the cdf of the stable distribution for a vector of probabilities p
template<typename myFloat>
Vec quantile(const Vec& p,               ///< [in] vector of probabilities at which to calculate the inverse
             const myFloat alpha,         ///< [in] the structural parameter of the distribution
             const myFloat beta,          ///< [in] the skewness parameter of the distribution
             const Vec& gamma,            ///< [in] a vector of scale parameters for distribution
             const Vec& delta,            ///< [in] a vector of location parameters for distribution
             const int pm,                ///< [in] the parameterization, 0, 1 or 2
             const int lower_tail,        ///< [in] return the lower tail, otherwise the upper tail
             const int log_p,             ///< [in] p is the log(pdf), otherwise the pdf
             const myFloat dbltol,        ///< [in] the tolerance for the inverse
             Controllers<myFloat> ctls,   ///< [in,out] reference to the integration controller
             const int verbose            ///< [in] indicator for verbose output
             );

/// return the derivative wrt x of the pdf of the standard stable distribution for vector of x
template<typename myFloat>
Vec std_ddx_pdf(const Vec& x,                     ///< [in] the vector of points at which to calculate the pdf
                const myFloat alpha,              ///< [in] the structural parameter of the distribution
                const myFloat beta,               ///< [in] the skewness parameter of the distribution
                const Parameterization pm_std,    ///< [in] the parameterization
                Controllers<myFloat> ctls,        ///< [in,out] reference to the integration controller
                const int verbose                 ///< [in] indicator for verbose output
);

/// return the derivative wrt x of the pdf of the stable distribution for vector of x
template<typename myFloat>
Vec ddx_pdf(const Vec& x,                     ///< [in] the vector of points at which to calculate the pdf
                 const myFloat alpha,         ///< [in] the structural parameter of the distribution
                 const myFloat beta,          ///< [in] the skewness parameter of the distribution
                 const Vec& gamma,            ///< [in] a vector of scale parameters for distribution
                 const Vec& delta,            ///< [in] a vector of location parameters for distribution
                 const int pm,                ///< [in] the parameterization, 0, 1 or 2
                 Controllers<myFloat> ctls,   ///< [in,out] reference to the integration controller
                 const int verbose            ///< [in] indicator for verbose output
                 );

/// return a vector of deviates from the standard stable distribution
template<typename myFloat>
Vec std_random_stable(const myFloat alpha,              ///< [in] the structural parameter of the distribution
                      const myFloat beta,               ///< [in] the skewness parameter of the distribution
                      const Parameterization pm_std,    ///< [in] the parameterization
                      const Vec &u1,                    ///< [in] a vector of the first std uniform deviates used
                      const Vec &u2                     ///< [in] a vector of the second std uniform deviates used
);

/// return a vector of deviates from the stable distribution
template<typename myFloat>
Vec random_stable(const myFloat alpha,    ///< [in] the structural parameter of the distribution
             const myFloat beta,          ///< [in] the skewness parameter of the distribution
             const Vec& gamma,            ///< [in] a vector of scale parameters for distribution
             const Vec& delta,            ///< [in] a vector of location parameters for distribution
             const int pm,                ///< [in] the parameterization, 0, 1 or 2
             Controllers<myFloat> ctls,   ///< [in,out] reference to the integration controller
             const int verbose,           ///< [in] indicator for verbose output
             const Vec &u1,               ///< [in] a vector of the first std uniform deviates used
             const Vec &u2                ///< [in] a vector of the second std uniform deviates used
             );
  
} //namespace stable_distribution
  
#define VEC_TEMPLATES(EXT, EXP, T) \
EXT template EXP Matrix<T,Dynamic,1> std_pdf(const Matrix<T,Dynamic,1>&, const T, const T, const Parameterization, const int, Controllers<T>, const int); \
EXT template EXP Matrix<T,Dynamic,1> pdf(const Matrix<T,Dynamic,1>&, const T, const T, const Matrix<T,Dynamic,1>&, const Matrix<T,Dynamic,1>&, const int, const int, Controllers<T>, const int); \
EXT template EXP Matrix<T,Dynamic,1> std_cdf(const Matrix<T,Dynamic,1>&, const T, const T, const Parameterization, const int, const int, Controllers<T>, const int); \
EXT template EXP Matrix<T,Dynamic,1> cdf(const Matrix<T,Dynamic,1>&, const T, const T, const Matrix<T,Dynamic,1>&, const Matrix<T,Dynamic,1>&, const int, const int, const int, Controllers<T>, const int); \
EXT template EXP Matrix<T,Dynamic,1> std_quantile(const Matrix<T,Dynamic,1>&, const T, const T, const Parameterization, const int, const int, const T, Controllers<T>, const int); \
EXT template EXP Matrix<T,Dynamic,1> quantile(const Matrix<T,Dynamic,1>&, const T, const T, const Matrix<T,Dynamic,1>&, const Matrix<T,Dynamic,1>&, const int, const int, const int, const T, Controllers<T>, const int); \
EXT template EXP Matrix<T,Dynamic,1> std_ddx_pdf(const Matrix<T,Dynamic,1>&, const T, const T, const Parameterization, Controllers<T>, const int); \
EXT template EXP Matrix<T,Dynamic,1> ddx_pdf(const Matrix<T,Dynamic,1>&, const T, const T, const Matrix<T,Dynamic,1>&, const Matrix<T,Dynamic,1>&, const int, Controllers<T>, const int); \
EXT template EXP Matrix<T,Dynamic,1> std_random_stable(const T, const T, const Parameterization, const Matrix<T,Dynamic,1>&, const Matrix<T,Dynamic,1>&); \
EXT template EXP Matrix<T,Dynamic,1> random_stable(const T, const T, const Matrix<T,Dynamic,1>&, const Matrix<T,Dynamic,1>&, const int, Controllers<T>, const int, const Matrix<T,Dynamic,1>&, const Matrix<T,Dynamic,1>&);

#ifdef LIBRARY

#include "stable_distribution_Vec_impl.h"

#endif

namespace stable_distribution {
VEC_TEMPLATES(STABLE_EXT, STABLE_EXP, double)
#ifdef CPP_BIN_FLOAT
  VEC_TEMPLATES(STABLE_EXT, STABLE_EXP, CppBinFloat)
#endif
#ifdef MPFR_FLOAT
  VEC_TEMPLATES(STABLE_EXT, STABLE_EXP, MpfrFloat)
#endif
#ifdef MPREAL
  VEC_TEMPLATES(STABLE_EXT, STABLE_EXP, mpreal)
#endif
}

#undef Vec

#endif // stable_distribution_Vec_h
