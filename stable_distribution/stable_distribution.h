/// \file stable_distribution.h
/// \author Joseph Dunn
/// \copyright 2016 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#ifndef stable_distribution_H
#define stable_distribution_H

#include <iostream>
#include <string>
#include <vector>
#include "adaptive_integration.h"
#include "myFloat.h"


namespace stable_distribution {

using std::cout;
using std::cerr;
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
  
/// the data and functions to calculate the standard stable distribution
template<typename myFloat>
class StandardStableDistribution {
public:
  static myFloat pi;                        ///< pi
  static myFloat pi2;                       ///< pi/2
  static myFloat eps;                       ///< machine epsilon
  static myFloat large_exp_arg;             ///< the largest argument for exp()
  static myFloat PosInf;                    ///< positive inifinity 
  static myFloat NegInf;                    ///< negative inifinity
  static myFloat zeta_tol;                  ///< tolerance for the use of f_zeta
  static bool initialized;                  ///< have the static member been initialized
  myFloat alpha;                            ///< the structural parameter to stable distribution
  myFloat beta_input;                       ///< beta as input,  sign might be flipped in the computation 
private:
  myFloat x_input;                          ///< the input x in the S0 parameterization, sign might be flipped in computation
  myFloat x_m_zeta_input;                   ///< x - zeta as input, sign might be flipped in computation

  // variables independent of x used when alpha != 1
  myFloat zeta;                             ///< beta * tan( pi * alpha / 2 
  myFloat theta0_x_gt_zeta;                 ///< the lower limit of integration when x  > zeta 
  myFloat cat0;                             ///< cos(alpha * theta0) 

  // variable dependent on x used when alpha != 1
  myFloat beta;                             ///< the skewness parameter actually used in computation 
  myFloat theta0;                           ///< the lower limit of integration used 
  myFloat x_m_zet;                          ///< x - zeta actually used 
  myFloat add_l;                            ///< adjustment parameter used when lower limit is moved to zero 
  myFloat add_r;                            ///< adjustment parameter used when upper limit is moved to zero
  myFloat th_min;                           ///< the lower limit of integration, usually 0
  myFloat th_max;                           ///< the limit of integration, both th_l and th_r range from th_min to th_max
  myFloat c2;                               ///< scaling parameter for pdf 
  myFloat c_ddx;                            ///< scaling parameter for ddx_pdf
  bool use_f_zeta;                          ///< use result for zeta w/o integrating

  /// special cases 
  enum dist {Cauchy=1, normal=2, fin_support=3, other=4};
  dist dist_type;                          ///< the type indicating special case 

  /// types of g function 
  enum Fun {fun_g_l=1,   ///< g uses th_l for alpha != 1 
            fun_g_r=2,   ///< g uses th_r for alpha !=1 
            fun_ga1_r=3, ///< g use u_r for alpha == 1
            fun_ga1_c=4  ///< g use u_c for alpha == 1
           };
  Fun fun_type;                            ///< function type 

  // variable used when alpha = 1
  myFloat abs_x;                            ///< abs(x) 
  myFloat i2b;                              ///< 1 / (2*beta) 
  myFloat p2b;                              ///< pi / (2*beta) 
  myFloat ea;                               ///< -p2b*abs_x 
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
  typename IntegrationController<myFloat>::TerminationCode termination_code;                                 ///< the error code from the integrator
  int num_iter;                             ///< the number of iterations required by quantile
  /// initialize the constants pi and pi2
  static void initialize();
  /// consturctor
  StandardStableDistribution(myFloat alpha,                  ///< [in] the structural parameter of the distribution
                             myFloat beta,                       ///< [in] the skewness parameter of the distribution
                             Controllers<myFloat> ctls,        ///< [in] reference to integration controller
                             int verbose                        ///< [in] indicator for verbose output
  ) : alpha(alpha), beta_input(beta), x_input(NAN),
  x_m_zeta_input(NAN), controllers(ctls), verbose(verbose) {
    if (!initialized) initialize();
    if (fabs(alpha-1)>100*std::numeric_limits<myFloat>::epsilon()){
      zeta = -beta_input*tan(alpha*pi2);
      theta0_x_gt_zeta = atan(beta_input*tan(alpha*pi2))/alpha;
      theta0_x_gt_zeta = min(pi2,max<myFloat>(-pi2,theta0_x_gt_zeta));
      cat0=1/sqrt(1+zeta*zeta);
    } else {
      zeta=0;
      this->alpha=1;
      c2=pi2*fabs(1/(2*beta_input));
      c_ddx=-c2*pi2/beta_input;
    }

  };

  /// set or reset x - zeta
  void set_x_m_zeta(myFloat x,  ///< [in] x - zeta whether positive or negative
                    Parameterization pm  ///< [in] the type of parameterization
                   );

  /// returns the value of the integrand at th
  myFloat g(myFloat th ///< [in] the variable of integration
          ) const{
    switch(fun_type){
    case fun_g_l:
      return g_l(th);
    case fun_g_r:
      return g_r(th);
    case fun_ga1_r:
      return ga1_r(th);
    case fun_ga1_c:
      return ga1_c(th);
    }
  }
private:
  /// g when th_l is theta 0 
  myFloat g_l(myFloat th_l ///> [in] point of integration when theta0 is moved to zero
             ) const;

  /// g when th_r is pi/2 - th 
  myFloat g_r(myFloat th_r ///< [in]the distance down from pi/2
             ) const;

  /// g when alpha=1 and u_r = 1-u = 1 - th/pi2 
  myFloat ga1_r(myFloat u_r ///< [in]the distance down from 1
               ) const;
    
  /// g when alpha=1 & x large with u_c = u_r - u_r_2 where g(u_r_2) = 1
  myFloat ga1_c(myFloat u_c ///< [in]the distance up from u_r_2
  ) const;
  
  /// return the derivative of ln(g) wrt to th
  myFloat dlng_dth(myFloat th ///< the point of integration
                    ) {
    switch(fun_type){
      case fun_g_l:
        return dlng_dth_l(th);
      case fun_g_r:
        return dlng_dth_r(th);
      case fun_ga1_r:
        return dlnga1_du_r(th);
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
  myFloat f_of_g(myFloat th ///< [in] point at which to evaluate
                ) {return (*f)(std_stable_dist->g(th), std_stable_dist);}

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
  myFloat g_max;            ///< the cap on the allowable g 
  
public:
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
      << ", g(theta) = " << g_ << endl;
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

#define STABLE_TEMPLATES(EXT,T) \
EXT template class adaptive_integration::IntegrationController<T>; \
EXT template class stable_distribution::StandardStableDistribution<T>; \
EXT template T stable_distribution::random_stable<T>(T, T, T, T, Parameterization); \
EXT template T stable_distribution::C_stable_tail<T>(T, bool); \
EXT template T stable_distribution::dPareto<T>(T, T, T, bool); \
EXT template T stable_distribution::pPareto<T>(T, T, T, bool, bool);


#ifdef LIBRARY

#include "stable_distribution_base.h"
#include "stable_distribution_cdf.h"
#include "stable_distribution_pdf.h"
#include "stable_distribution_ddx_pdf.h"
#include "stable_distribution_quantile.h"
#include "stable_distribution_random.h"
#include "stable_distribution_mode.h"

STABLE_TEMPLATES(,double)
STABLE_TEMPLATES(,mpreal)
  
#else

STABLE_TEMPLATES(extern,double)
STABLE_TEMPLATES(extern,mpreal)

#endif

#endif // ndef stable_distribution.H
