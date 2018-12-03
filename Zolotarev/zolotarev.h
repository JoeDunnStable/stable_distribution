/// \file zolotarev.h
/// cdf, pdf and ddx_pdf of standard stable distribution per Zolotarev
/// \author Joseph Dunn
/// \copyright 2017, 2018 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3


#ifndef zolotarev_h
#define zolotarev_h

#include "myFloat.h"
#include "stable_distribution.h"
#include <Eigen/Dense>
#include <vector>
#include <string>

namespace stable_distribution {
  
using adaptive_integration::Integral;
using adaptive_integration::IntegrationController;
using Eigen::Array;
using Eigen::Dynamic;
using std::vector;
using std::string;
  
template<typename myFloat> class Zolotarev;
  
template<typename myFloat>
class Q_integrand {
  Zolotarev<myFloat>& zol;
public:
  Q_integrand(Zolotarev<myFloat>& zol) : zol(zol) {}
  void operator() (vector<myFloat>& u);
};

/// integrand for  use in calculating pdf for alpha=1 and small x
template<typename myFloat>
class q_integrand{
  Zolotarev<myFloat>& zol;
public:
  q_integrand(Zolotarev<myFloat>& zol) : zol(zol) {}
  void operator() (vector<myFloat>& u);
};

/// integrand for  use in calculating ddx_pdf for alpha=1 and small x
template<typename myFloat>
class ddx_q_integrand {
  Zolotarev<myFloat>& zol;
public:
  ddx_q_integrand(Zolotarev<myFloat>& zol) : zol(zol) {}
  void operator() (vector<myFloat>& u);
};
  
enum ResultType {asymptotic, convergent};  ///< enum for the result

  /// The data and functions needed to calculated Zolotarev's series for
/// the stable distrivution
template<typename myFloat>
class Zolotarev {
public:
  static Array<myFloat, Dynamic, Dynamic> gamma_at_integers;  ///<gamma function and derivatives at integers
  static bool initialized;   ///< Are the static elements initialized
  static myFloat pi;         ///< Cached pi to requisite precision
  static int digits10;           ///< Precision at time of initialization
  static int max_n_conv;     ///< the maximum number of terms for convergent series
  static int max_n_asymp;    ///< the maximum number of terms for asymptotic series
  static vector<myFloat> points;     ///< the initial subdivision of [0, pi/2] for integration
  IntegrationController<myFloat>* cntl; ///< pointer to integration controller to use
  myFloat alpha;             ///< shape parameter for stable distribution
  myFloat alpha_m_1;         ///< alpah - 1
  myFloat beta_input;        ///< beta, the skewness parameter as input
  myFloat beta;              ///< beta actually use.  +/- beta_input
  bool positive_x;           ///< indicates whether x_m_zet has been negated for calculation
  bool positive_xB;          ///< indicatew wheter xB has been negated for calculation
  myFloat zeta;              ///< beta * tan( pi * alpha / 2)
  myFloat theta0_x_gt_zeta;  ///< the lower limit of integration when x  > zeta
  std::vector<myFloat> Q_cdf; ///< the coefficients in the asumptotic expansion for cdf
                              ///< when abs(beta) == 1
  std::vector<myFloat> Q_pdf;///< the coefficents in the asymptotic expansion for pdf
                             ///< when abs(beta) == 1
  std::vector<myFloat> Q_ddx_pdf;    ///< the coefficents in the asymptotic expansion for ddx_pdf
                             ///< when abs(beta) == 1

  myFloat theta0;            ///< the lower limit of integration
  myFloat theta;             ///< for alpha!=1, theta0/(pi/2)
  myFloat rho;               ///< for alpah!=1, (1+theta)/2
  myFloat gammaB;            ///< Zolotarev lambda^(1/alpha)
  myFloat betaB;             ///< for alpha==1, beta for Zolotarev's B representation
  myFloat x_m_zet;           ///< x - zeta
  myFloat xB;                ///< x in Zolotarev's B representation = x_m_zeta/gammaB
  int n_convergent;          ///< the number of terms used in the convergent series
  int n_asymptotic;          ///< the number of terms used in the asymptotic series
  int n;                     ///< the number of terms used in the best series
  myFloat error_convergent;  ///< the estimated error in the convergent series
  myFloat error_asymptotic;  ///< the estimated error in the asymptotic series
  myFloat error;             ///< estimated error in the best series
  myFloat result_convergent; ///< the result of the convergent series
  myFloat result_asymptotic; ///< the result of the asymptotic series
  myFloat result;            ///< the result of the best series
  
  ResultType result_type;    ///< indicator giving which series was used
  int verbose;               ///< verbosity of output
  
  /// constructor for zolotarev
  Zolotarev(myFloat alpha, myFloat beta_input, IntegrationController<myFloat>* cntl,
            int verbose=0, int verbose_integration=0);
  
  /// common logic to set x - zeta
  void set_x_m_zeta(myFloat x_m_zeta_in);
  
  /// pdf using zolotarev expansion
  myFloat pdf(myFloat x0, Parameterization pm=S0);
  
  /// ddx_pdf using zolotarev expansion
  myFloat ddx_pdf(myFloat x0, Parameterization pm=S0);
  
  /// cdf using zolotarev expansion
  myFloat cdf(myFloat x0, int lower_tail, Parameterization pm=S0);
  
  /// Functor to evaluate cdf when alpha=1 and x is small
  Integral<myFloat, Q_integrand<myFloat> > cdf_alpha_1;

  /// Functor to evaluate pdf when alpha=1 and x is small
  Integral<myFloat, q_integrand<myFloat> > pdf_alpha_1;
  
  /// Functor to evaluate ddx_pdf when alpha=1 and x is small
  Integral<myFloat, ddx_q_integrand<myFloat> > ddx_pdf_alpha_1;
  
};  //Zolotarev

}  //namespace stable_distribution

#define ZOLOTAREV_TEMPLATES(EXT, T) \
EXT template class Zolotarev<T>; \
EXT template class Q_integrand<T>; \
EXT template class q_integrand<T>; \
EXT template class ddx_q_integrand<T>;

#ifdef LIBRARY
#include "zolotarev_base.h"
#include "zolotarev_cdf.h"
#include "zolotarev_pdf.h"
#include "zolotarev_ddx_pdf.h"
#include "zolotarev_integral.h"

namespace stable_distribution {
ZOLOTAREV_TEMPLATES(,double)
#ifdef CPP_BIN_FLOAT
  ZOLOTAREV_TEMPLATES(,CppBinFloat)
#endif
#ifdef MPFR_FLOAT
  ZOLOTAREV_TEMPLATES(,MpfrFloat)
#endif
#ifdef MPREAL
ZOLOTAREV_TEMPLATES(,mpreal)
#endif
}

#else

namespace stable_distribution {
ZOLOTAREV_TEMPLATES(extern, double)
#ifdef CPP_BIN_FLOAT
  ZOLOTAREV_TEMPLATES(extern, CppBinFloat)
#endif
#ifdef MPFR_FLOAT
  ZOLOTAREV_TEMPLATES(extern, MpfrFloat)
#endif
#ifdef MPREAL
ZOLOTAREV_TEMPLATES(extern, mpreal)
#endif
}

#endif

#endif /* zolotarev_h */
