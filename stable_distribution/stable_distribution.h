/// \file stable_distribution.h
/// Class for standard stable distriubution
/// \author Joseph Dunn
/// \copyright 2016, 2017, 2018 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#ifndef stable_distribution_H
#define stable_distribution_H

// Generic definitions for shared library support
#ifdef LIBRARY
  #define STABLE_EXT
  #if defined _MSC_VER
    #define STABLE_EXP __declspec(dllexport)
  #else
    #define STABLE_EXP __attribute__ ((visibility("default")))
  #endif
#else
  #define STABLE_EXT extern
  #if defined _MSC_VER
    #define STABLE_EXP __declspec(dllimport)
  #else
    #define STABLE_EXP
  #endif
#endif

#include <string>
#include <vector>
#include "adaptive_integration.h"
#include "myFloat.h"
#include <Eigen/Dense>
#include <mutex>
using std::mutex;
using std::unique_lock;
using Eigen::Array;
using Eigen::Dynamic;

namespace stable_distribution {

using std::endl;
using std::ostream;

using std::string;
using std::vector;
using std::max;
using std::min;
using adaptive_integration::Kronrod;
using adaptive_integration::IntegrationController;
using adaptive_integration::Subinterval;
  
template<typename myFloat> class StandardStableDistribution;
template<typename myFloat> ostream& operator<< (ostream& os,
                        const StandardStableDistribution<myFloat>& dist);

/// A pair of integration controllers the 
template<typename myFloat>
  class Controllers  {
public:
  IntegrationController<myFloat>& controller;
  IntegrationController<double>& ctl_double;
  Controllers(IntegrationController<myFloat>& controller,
              IntegrationController<double>& ctl_double) :
    controller(controller), ctl_double(ctl_double){}
  double q_guess(myFloat p, myFloat alpha, myFloat beta, int lower_tail, int log_p);
};

/// the types of parameterizations supported by StandardStableDistribution
enum Parameterization {S0=0, S1=1};
  
/// a class allowing the use of alpha-1 to initialize StandardStableDistribution
template<typename myFloat>
struct AlphaMinusOne {
  myFloat alpha_minus_one;
  AlphaMinusOne(myFloat alpha_minus_one) : alpha_minus_one(alpha_minus_one) {}
};
  
/// the data and functions to calculate the standard stable distribution
template<typename myFloat>
class STABLE_EXP StandardStableDistribution {
public:
  static mutex stable_mutex;                ///< protect duing initialization of static variables
  static myFloat pi;                        ///< pi
  static myFloat pi2;                       ///< pi/2
  static myFloat eps;                       ///< machine epsilon
  static myFloat xmin;                      ///< the minimum positive number
  static myFloat large_exp_arg;             ///< the largest argument for exp()
  static myFloat PosInf;                    ///< positive inifinity 
  static myFloat NegInf;                    ///< negative inifinity
  static double threshhold_1;               ///< used to define near=abs(alpha-1) < threshhold_1 * abs(beta)
  static Array<myFloat, Dynamic, Dynamic> gamma_at_integers;  ///<gamma function and derivatives at integers
  static int max_n;     ///< the maximum number of terms for series expansions
  static bool initialized;                  ///< have the static member been initialized
  myFloat alpha;                            ///< the structural parameter to stable distribution
  myFloat alpha_minus_one;                  ///< alpha minus one.  Used near 1
  myFloat beta_input;                       ///< beta as input,  sign might be flipped in the computation 
private:
  Fmt<myFloat> fmt;                         /// a manipulator to format myFloat to maximum precision
  myFloat x_input;                          ///< the input x in the S0 parameterization, sign might be flipped in computation
  myFloat x_m_zeta_input;                   ///< x - zeta as input, sign might be flipped in computation

  // variables independent of x used when alpha != 1
  myFloat zeta;                             ///< beta * tan( pi * alpha / 2 
  myFloat theta0_x_gt_zeta;                 ///< the lower limit of integration when x  > zeta 
  myFloat cat0;                             ///< cos(alpha * theta0) 
  myFloat gammaB;                           ///< Zolotarev lambda^(1/alpha)
  std::vector<myFloat> Q_cdf;               ///< the coefficients in the asumptotic expansion for cdf
                                            ///< when abs(beta) == 1
  std::vector<myFloat> Q_pdf;               ///< the coefficents in the asymptotic expansion for pdf
                                            ///< when abs(beta) == 1
  std::vector<myFloat> Q_ddx_pdf;           ///< the coefficents in the asymptotic expansion for ddx_pdf
                                            ///< when abs(beta) == 1

  // variable dependent on x used when alpha != 1
  myFloat beta;                             ///< the skewness parameter actually used in computation 
  myFloat betaB;                            ///< for alpha==1, beta for Zolotarev's B representation
  myFloat betaB_p_1;                        ///< betaB + 1 via numerically stable method
  myFloat one_m_betaB;                      ///< 1 - betaB vis numerically stable method
  myFloat theta0;                           ///< the lower limit of integration used
  bool positive_x;                          ///< Indicator the sign x has been flipped
  bool positive_xB;                         ///< Indicator the sign of xB has been flipped
  myFloat theta;                            ///< for alpha!=1, theta0/(pi/2)
  myFloat rho;                              ///< for alpah!=1, (1+theta)/2
  myFloat x_m_zet;                          ///< x - zeta actually used
  myFloat xB;                               ///< x in Zolotarev's B representation = x_m_zeta/gammaB
  bool avoid_series_small_x;                ///< avoid the use of the small x series even if integral is inaccurate
  bool use_series_small_x;                  ///< use the Zolotarev small x serries
  bool avoid_series_large_x;                ///< avoid the use of the large x series even if integral is inaccurate
  bool use_series_large_x;                  ///< use the Zolotarev large x series
  myFloat add_l;                            ///< adjustment parameter used when lower limit is moved to zero
  myFloat add_r;                            ///< adjustment parameter used when upper limit is moved to zero
  myFloat th_span;                          ///< the length of the interval of integration
  myFloat th_min;                           ///< the lower limit of integration
  myFloat th_max;                           ///< the upper limit of integration
  myFloat c2;                               ///< scaling parameter for pdf 
  myFloat c_ddx;                            ///< scaling parameter for ddx_pdf

  /// special cases 
  enum dist {Cauchy=1, normal=2, fin_support=3, other=4};
  dist dist_type;                          ///< the type indicating special case 

  /// flag used to determine which calculation to use
  bool small_x_m_zet;
  
  /// types of g function
  enum Fun {fun_g_l=1,        ///< g uses th_l for alpha != 1
            fun_g_r=2,        ///< g uses th_r for alpha !=1
            fun_ga1_r=3,      ///< g use u_r for alpha == 1
            fun_g_l_c=4,      ///< g w th_l primary and th_c secondary for alpha ! 1
            fun_g_r_c=5,      ///< g w th_r primary and th_c secondary for alpha != 1
            fun_ga1_c=6,      ///< g use u_c for alpha == 1
            fun_g_a_near_1_l=7, ///< g uses th_c for small x_m_zet and alpha ~ 1
            fun_g_a_near_1_r=8  ///< g uses th_c for large x_m_zet and alpha ~ 1
           };
  Fun fun_type;                            ///< function type 

  // variable used when alpha = 1
  myFloat abs_x;                            ///< abs(x) 
  myFloat i2b;                              ///< 1 / (2*beta) 
  myFloat p2b;                              ///< pi / (2*beta) 
  myFloat ea;                               ///< -p2b*abs_x
  // variables used when using th_c when alpha != 1
  myFloat ln_g_th2;                         ///< ln of g at th2
  myFloat cot_th2_1;                        ///< first cotan in formuala for g(th_c)
  myFloat cot_th2_2;                        ///< second cotan in formula for g(th_c)
  myFloat cot_th2_3;                        ///< third cotan in formula for g(th_c)
  // variables only used when alpha == 1 and abs(x) is large
  myFloat ln_g_u2;                          ///< the residual ln(g) when u==u2
  myFloat costh_u2;                         ///< cos(th) when u_r==u2
  myFloat tanth_u2;                         ///< tan(th) when u_r == u2
  myFloat h_u2;                             ///< h when u_r == u2
  
    
  // variables controling the integrator
public:
  Controllers<myFloat> controllers;       ///< pointer to integration controllers
    int verbose;                           ///< indicator to request detailed output

  // the results of the map_g process. Used to establish the initial subintervals.
  bool good_theta2;                         ///< was map_g able to find theta 2, point where g==1
  myFloat theta2;                           ///< point where g == 1
  myFloat g_theta2_error;                   ///< uncertainty in g at theta2
  myFloat c_g_theta2_error;                ///< uncertainty in result because of g
  myFloat g_theta2;                         ///< actual value of g at theta 2
  myFloat g_dd_theta2;                      ///< derivative of g wrt theta at theta2 
  vector<myFloat> points;                   ///< list of endpoints of initial subdivision 
    
  int last;                                ///< The last subinterval actually used 
  // additional output from integrator
  myFloat abserr;                           ///< the estimated absolute error of calculation 
  int neval;                               ///< the number of function evaluations 
  // additional output from series
  int n_series;                             ///< the number of terms used in the series
  myFloat error_series;                     ///< the estimated error in the series expansion
  typename IntegrationController<myFloat>::TerminationCode termination_code;                                 ///< the error code from the integrator
  int num_iter;                             ///< the number of iterations required by quantile
  /// initialize the constants pi and pi2
  static void initialize();
  /// consturctor
  StandardStableDistribution(myFloat alpha,                  ///< [in] the structural parameter of the distribution
                             myFloat beta,                       ///< [in] the skewness parameter of the distribution
                             Controllers<myFloat> ctls,        ///< [in] reference to integration controller
                             int verbose                        ///< [in] indicator for verbose output
  );
  /// consturctor using AlphaMinusOne
  StandardStableDistribution(AlphaMinusOne<myFloat> a_m_1,      ///< [in] the structural parameter minus one
                             myFloat beta,                       ///< [in] the skewness parameter of the distribution
                             Controllers<myFloat> ctls,        ///< [in] reference to integration controller
                             int verbose                        ///< [in] indicator for verbose output
  );
  
  /// set or reset x - zeta
  void set_x_m_zeta(myFloat x,  ///< [in] x - zeta whether positive or negative
                    Parameterization pm  ///< [in] the type of parameterization
                   );

  /// returns the value of the integrand at th
  myFloat g(myFloat th ///< [in] the variable of integration
          ) const
  {
    switch(fun_type){
      case fun_g_l:
        return g_l(th);
      case fun_g_r:
        return g_r(th);
      case fun_ga1_r:
        return ga1_r(th);
      case fun_g_l_c:
        return g_l_c(th);
      case fun_g_r_c:
        return g_r_c(th);
      case fun_ga1_c:
        return ga1_c(th);
      case fun_g_a_near_1_l:
        return g_a_near_1_l(th);
      case fun_g_a_near_1_r:
        return g_a_near_1_r(th);
    }
  }
  bool use_series() {
    return (boost::math::isfinite(x_input) && dist_type==other && (fabs(xB) <= 1e-3 || fabs(xB) >= 1e3));
  }
private:
  /// Q initializer used by series when |beta| == 1
  void Q_initializer() ;
  
  /// g when th_l is theta 0 
  myFloat g_l(myFloat th_l ///< [in] point of integration when theta0 is moved to zero
             ) const;

  /// g when th_r is pi/2 - th 
  myFloat g_r(myFloat th_r ///< [in]the distance down from pi/2
             ) const;

  /// g when alpha=1 and u_r = 1-u = 1 - th/pi2 
  myFloat ga1_r(myFloat u_r ///< [in]the distance down from 1
               ) const;
  
  /// g when alpha!=1, th_l is primary and th_c is secondary
  myFloat g_l_c(myFloat th_c ///< [in] the th_l-theta2
               ) const;
    
  /// g when alpha!=1, th_r is primary and th_c is secondary
  myFloat g_r_c(myFloat th_c ///< [in] the th_r-theta2
               ) const;
  
  /// g when alpha=1 & x large with u_c = u_r - u_r_2 where g(u_r_2) = 1
  myFloat ga1_c(myFloat u_c ///< [in]the distance up from u_r_2
  ) const;
  
  /// g when alpha near 1 & sign(x) = sign(zeta)
  myFloat g_a_near_1_l(myFloat th_l  ///< [in] th_l = th - theta0
                      ) const;
  
  /// g when alpha near 1 & sign(x) != sign(zeta)
  myFloat g_a_near_1_r(myFloat th_r  ///< [in] th_r = pi/2 -th
                      )const;
  
  /// return the derivative of ln(g) wrt to th
  myFloat dlng_dth(myFloat th ///< the point of integration
                    ) {
    switch(fun_type){
      case fun_g_l:
      case fun_g_a_near_1_l :
        return dlng_dth_l(th);
      case fun_g_r:
      case fun_g_a_near_1_r :
        return dlng_dth_r(th);
      case fun_ga1_r:
        return dlnga1_du_r(th);
      case fun_g_l_c:
        return dlng_dth_l(th-th_min);
      case fun_g_r_c:
        return dlng_dth_r(th-th_min);
      case fun_ga1_c:
        return dlnga1_du_r(th-th_min);
    }
  }
 
  myFloat dlng_dth_l(myFloat th_l);
  myFloat dlng_dth_r(myFloat th_r);
  myFloat dlnga1_du_r(myFloat u_r);

  /// determine the initial subdivisions for integration 
  void map_g();

    /// guess a range for theta given value of g(theta) 
  void th_guess(const myFloat &value, myFloat &lower, myFloat &g_lower, myFloat &upper, myFloat &g_upper);
public:
  /// calculate the integral used for pdf 
  myFloat integrate_pdf(bool log_flag ///< [in] return the log of the answer
                           );

  /// calculate the integral used for ddx_pdf 
  myFloat integrate_ddx_pdf();

  /// calculate the integral used for cdf 
  myFloat integrate_cdf();
  
  /// calculate the pdf of the standard stable distribution
  myFloat pdf(myFloat x,              ///< [in] the point at which to evaluate the distribution
              int log_flag,           ///< [in] return the log of the result
              Parameterization pm=S0  ///< [in] the parameterization
              );

  /// calculate the cdf of the standard stable distribution 
  myFloat cdf(myFloat z,              ///< [in] the point at which to evaluate the distributioin
              int lower_tail,         ///< [in] return the lower tail, otherwise the upper tail
              int log_p,              ///< [in] return the log of p
              Parameterization pm=S0  ///< [in] the parameterization
              );

  /// caluclate the inverse of the cdf of the standard stable distribution
  myFloat quantile(myFloat p,              ///< [in] the target probability or its log
                   int lower_tail,         ///< [in] return lower tail, other return upper tail
                   int log_p,              ///< [in] interpret p a the log of the probability
                   myFloat dbltol,         ///< [in] the tolerance used in determining the inverse
                   Parameterization pm=S0  ///< [in] the parameterization
                   );

  /// calculate the derivative wrt x of the pdf of the standard stable distribution 
  myFloat ddx_pdf(myFloat x,    ///< [in] the point at which to evaluate the distribution
                  Parameterization pm=S0  ///< [in] the parameterization
                  );
  /// return the location and size of the mode of the standard stable distribution
  
  /// pdf using zolotarev expansion for very large x
  myFloat series_large_x_pdf(myFloat x0, Parameterization pm=S0);
  
  /// ddx_pdf using zolotarev expansion for very large x
  myFloat series_large_x_ddx_pdf(myFloat x0, Parameterization pm=S0);
  
  /// cdf using zolotarev expansion for very large x
  myFloat series_large_x_cdf(myFloat x0, int lower_tail, Parameterization pm=S0);
  
  /// pdf using zolotarev expansion for very small x
  myFloat series_small_x_pdf(myFloat x0, Parameterization pm=S0);
  
  /// ddx_pdf using zolotarev expansion for very small x
  myFloat series_small_x_ddx_pdf(myFloat x0, Parameterization pm=S0);
  
  /// cdf using zolotarev expansion for very small x
  myFloat series_small_x_cdf(myFloat x0, int lower_tail, Parameterization pm=S0);
  

  std::pair<myFloat, myFloat> mode(myFloat dbltol,         ///< [in] the tolerance needed for the solution
               int verbose=0,          ///< [in] verbose indicator for mode calculation
               Parameterization pm=S0  ///< [in] the parameterization
              );
  
  /// send an instance of StandardStableDistribution to an output stream
  friend ostream& operator<< <>(ostream& os,                  ///< [in,out] the output stream
                             const StandardStableDistribution<myFloat>& gc  ///< [in] the StableDistriubion to print
                             );
}; // StandardStableDistribution

/// generate a sample from the stable distribution
template<typename myFloat>
myFloat random_stable(myFloat alpha,     ///< [in] the structure parameter of the distribution
                 myFloat beta,      ///< [in]the skewness parameter of the distributuon
                 myFloat u1,        ///< [in] a random number uniform of 0,1
                 myFloat u2,         ///< [in] a second random number uniform on 0,1
                 Parameterization pm = S0  ///< the parameterization used
                 );

/// everything needed to integrate f(g) except the endpoints
/// and starting subintervals of the integration range.
template<typename myFloat>
class Integral_f_of_g {
private:
  StandardStableDistribution<myFloat> *std_stable_dist;                ///< a pointer to an instance of StandardStableDistribution
  myFloat (*f)(myFloat, StandardStableDistribution<myFloat>*); ///< a pointer to function which will be applied to g
  Controllers<myFloat>* controllers;
  
public:
  /// return f of g 
  void operator() (vector<myFloat>& th ///< [in, out] point at which to evaluate
                ) {
                  for (int i=0; i<th.size(); ++i)
                  th.at(i) = (*f)(std_stable_dist->g(th.at(i)), std_stable_dist);
                  
                }

  /// construct the functor
  Integral_f_of_g(myFloat (*f)(myFloat, StandardStableDistribution<myFloat>*),  ///< [in]pointer to function to be evaluated
                  StandardStableDistribution<myFloat>* std_stable_dist,                 ///< [in] pointer to an instance of StandardStableDistribution
                  Controllers<myFloat>* ctls = nullptr ///< [in] pointer to the integration controller
                 ) : std_stable_dist(std_stable_dist), f(f),
                     controllers((ctls==nullptr) ? &std_stable_dist->controllers : ctls) {}
  
  /// return a human readable error message
  string msg() {
    string msgs[]={"OK","Maximum subdivisions reached","Roundoff error detected",
      "Bad integrand behavior","Roundoff error in the extrapolation table",
      "Integral probably divergent","Input is invalid"};
    return msgs[termination_code];};
  
  /// return a pointer to the instance of StandardStableDistribution 
  StandardStableDistribution<myFloat>* get_std_stable_dist() {return std_stable_dist;}
  myFloat result;        ///< the approximation of integral of f(g(x)) 
  myFloat abserr;        ///< the estimated absolute error of the approximation 
  int neval;            ///< the number of function evaluations used 
  typename IntegrationController<myFloat>::TerminationCode termination_code; ///< the error indicator returned from IntegrationController::integrate
  int last;             ///< the number of subintervals required 
  
  /// return the approximation to the integral of f(g(x)) 
  myFloat operator() ();
};

/// Functor passed to toms748_solve to find x such that g(x) == value
template<typename myFloat>
class gSolve {
private:
  StandardStableDistribution<myFloat>* std_stable_dist;   ///< pointer to an instance of StandardStableDistribution
  myFloat value;            ///< the target value 
  bool log_flag;           ///< interpret the value as the log of the target 
  myFloat g_min;            ///< the floor on the allowable g 
public:
  myFloat g_max;            ///< the cap on the allowable g 
  
  /// constructor
  gSolve(myFloat value,              ///< [in] the value or log of the value to solve for
         StandardStableDistribution<myFloat>* std_stable_dist, ///< [in] pointer to an instance of StandardStableDistribution
         bool log_flag=false        ///< [in] interpret value as the log of g
  ) :
  std_stable_dist(std_stable_dist), value(value), log_flag(log_flag), g_min(std::numeric_limits<myFloat>::min()),
  g_max(std::numeric_limits<myFloat>::max()) {}
  
  /// return the difference between g(th) and the target value 
  myFloat operator()(const myFloat th ///< [in] the point at which to evaluate g
                    ) {
    myFloat g_ = std_stable_dist->g(th);
    g_=max(g_min,min(g_,g_max));
    if (std_stable_dist->verbose >=3)
      cout << "      theta = " << th
      << ", " << (log_flag?"ln_":"") << "g(theta) = " << (log_flag?log(g_):g_) << endl;
    if (log_flag)
      return log(g_)-value;
    else
      return g_ - value;
  }
  /// set the target value for g(th) 
  void set_value(myFloat value_in) {value=value_in;}
  
}; //gSolve

/// calculate the constant used for the asymptotoic behavior of the tail
template<typename myFloat>
myFloat C_stable_tail(myFloat alpha,  ///< [in] the structural parameter
                      bool log_flag   ///< [in] return the log
);
  
/// return the pdf of the Pareto distribution
template<typename myFloat>
myFloat dPareto(myFloat x,      ///< [in]the point at which to evaluate the distribution 
               myFloat alpha,  ///< [in] the structural parameter of the stable distribution 
               myFloat beta,   ///< [in] the skewness parameter of the stable distribution 
               bool log_flag  ///< [in] return the low of the distribution 
               );

/// return the cdf of the Pareto distribution 
template<typename myFloat>
myFloat pPareto(myFloat x,        ///< [in] the point at which to evaluate the distribution
               myFloat alpha,    ///< [in] the structural parameter of the stable distribution 
               myFloat beta,     ///< [in] the skewness parameter of the stable idstribution 
               bool lower_tail, ///< [in] return the lower tail, otherwise the upper tail 
               bool log_p       ///< [in] return the log of the pdf, otherwise the pdf 
               );

/// class to determine relative tolerance
template<typename myFloat>
class RelativeComparisonTolerance
{
public:
  /// constructor 
  RelativeComparisonTolerance(myFloat eps ///< the specified tolerance
                             ) : eps(eps) {};
  /// return true if a and b are within the tolerance 
  inline bool operator()(const myFloat& a, const myFloat& b)
  {
    return fabs(a - b) <= (eps * (std::min)(fabs(a), fabs(b)));
  }
private:
  myFloat eps;
};

/// class to determine absolute tolerance 
template<typename myFloat>
class AbsComparisonTolerance
{
public:
  /// constructor 
  AbsComparisonTolerance(myFloat eps_in ///< the specified tolerance
                        ) : eps(eps_in){};
  /// return true if a and b are within the tolerance
  inline bool operator()(const myFloat& a, const myFloat& b)
  {
    bool r=fabs(a - b) <= eps;
    return r;
  }
private:
  myFloat eps;
};
} // namespace stable_distribution

#define STABLE_TEMPLATES(EXT,EXP,T) \
EXT template class EXP adaptive_integration::IntegrationController<T>; \
EXT template class EXP stable_distribution::StandardStableDistribution<T>; \
EXT template EXP T stable_distribution::random_stable<T>(T, T, T, T, Parameterization); \
EXT template EXP T stable_distribution::C_stable_tail<T>(T, bool); \
EXT template EXP T stable_distribution::dPareto<T>(T, T, T, bool); \
EXT template EXP T stable_distribution::pPareto<T>(T, T, T, bool, bool);


#ifdef LIBRARY

#include "stable_distribution_base.h"
#include "stable_distribution_cdf.h"
#include "stable_distribution_pdf.h"
#include "stable_distribution_ddx_pdf.h"
#include "stable_distribution_quantile.h"
#include "stable_distribution_random.h"
#include "stable_distribution_mode.h"

#endif

STABLE_TEMPLATES(STABLE_EXT, STABLE_EXP ,double)
#ifdef CPP_BIN_FLOAT
STABLE_TEMPLATES(STABLE_EXT, STABLE_EXP ,CppBinFloat)
#endif
#ifdef MPFR_FLOAT
STABLE_TEMPLATES(STABLE_EXT, STABLE_EXP ,MpfrFloat)
#endif
#ifdef MPREAL
STABLE_TEMPLATES(STABLE_EXT, STABLE_EXP ,mpreal)
#endif
  

#endif // ndef stable_distribution.H
