/// \file stable_distribution_ddx_pdf.h
/// Implementation of derivative of pdf of standard stable distribution.
/// Included in stable_distribution.h when LIBRARY is defined
/// \author Joseph Dunn
/// \copyright 2016, 2017 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include <limits>
#include <iostream>
#include <iomanip>
#include <boost/math/tools/toms748_solve.hpp>

namespace stable_distribution {
  
using std::string;
using std::array;
using std::vector;
using std::sort;
using std::setw;
using std::right;
using std::setprecision;
using std::ios;
using std::pair;
using boost::math::tools::toms748_solve;

template<typename myFloat>
myFloat x_1_m_alpha_x_exp_m_x(myFloat x, StandardStableDistribution<myFloat>* std_stable_dist) {
  if (x > StandardStableDistribution<myFloat>::large_exp_arg)
    return 0;
  else
    return x * (1 - std_stable_dist->alpha*x) * exp(-x);
}

template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::integrate_ddx_pdf() {
  // --- pdf(x, alpha, beta, ..)  for alpha < 2 ---
  // For  x = zeta, have special case formula [n_gaussolan(1997)];
  // need to use it also for x ~= zeta, i.e., x.m.zet := |x - zeta| < delta :
  myFloat r;
  if (verbose)
    cout << "integrate_ddx_pdf:" << endl
         << "Integrand is g * (1 - alpha * g) * exp(-g)" << endl;
  // integral is zero near mode so need to use absolute tolerance near mode
  myFloat old_epsabs = controllers.controller.epsabs;
  controllers.controller.epsabs=(fabs(x_m_zeta_input) < 1) ? controllers.controller.epsrel/c_ddx : myFloat{0};
  Integral_f_of_g<myFloat> int_g1(x_1_m_alpha_x_exp_m_x,this, &controllers);
  controllers.controller.epsabs = old_epsabs;
  r=int_g1();
  abserr=fabs(c_ddx) * int_g1.abserr;
  c_g_theta2_error = fabs(c_ddx) * g_theta2_error;
  neval+=int_g1.neval;
  termination_code=int_g1.termination_code;
  last=int_g1.last;
  if (verbose) {
    if (verbose)
      cout << "  c_ddx*sum(r)= " << c_ddx << " * " << r << " = " << c_ddx*(r) << endl
      << "  abs.err = " << c_ddx*int_g1.abserr << endl
      << "  msg = " << int_g1.msg() << endl;
  }
  return c_ddx * (r);
} //integrate_ddx_pdf


template<typename myFloat>
myFloat StandardStableDistribution<myFloat>::ddx_pdf(myFloat x, Parameterization pm)
{
  cout.precision(20);
  cout.setf(ios::scientific, ios::floatfield);
  // Default values which will be reset in the integrator is used
  abserr = 0;
  neval = 0;
  termination_code = IntegrationController<myFloat>::TerminationCode::normal;
  last = 0;
  
  myFloat ret;
  set_x_m_zeta(x, pm);
  if (verbose)
    cout << "ddx_pdf: pm = " << pm << endl << *this;
  if (!boost::math::isfinite(x)) {
    if (verbose)
      cout << "  x is infinite, returning " << 0 << endl;
     return 0;
  }
  switch (dist_type) {
    case Cauchy :
      ret = (-2*x / (pi *(1 + x*x)))*(1/(1+x*x));
      if (verbose)
        cout << "  Cauchy distribution, returning " << ret << endl;
      return  ret;
    case normal :
      ret = -x*exp(-x*x/4)/(4*sqrt(pi));
      if (verbose)
        cout << "  Normal distribution, returning " << ret << endl;
      return  ret;
    case fin_support :
      if (verbose)
        cout << "  Outside of support, returning " << 0 << endl;
      return 0;
    case other :
      myFloat dfdx_zeta{0};
      
      if (alpha != 1) { // 0 < alpha < 2	&  |beta| <= 1 from above
        // General Case
        if (verbose)
          cout << "  General Case"<< endl;
        // For  x = zeta, have special case formula [n_gaussolan(1997)];
        if (alpha < 1 && fabs(beta)==1)
          dfdx_zeta = 0;
        else if (theta0_x_gt_zeta==0)
          dfdx_zeta = 0;
        else
          dfdx_zeta = tgamma(1 + 2/alpha) * sin(2*theta0_x_gt_zeta)
          /(2*pi * pow(1 + pow(zeta,2),(1/(alpha))));

        // need to use it also for x ~= zeta
        if (use_f_zeta) {
          if (verbose)
            cout << "  " << x_input << " ~= " << zeta
            << " Using dfdx_zeta()" << endl;
          return dfdx_zeta;
        }
        
      } // alpha != 1
      if (good_theta2) {
        ret = integrate_ddx_pdf();
        if (abserr > 1e-6 * fabs(ret) && fun_type>1) {
          ret = -(1+alpha)*dPareto(x_m_zeta_input+zeta, alpha, beta_input, false)/x;
          if (verbose)
            cout << "ddx_pdf:" << endl << "  Integral has large error and x is large. using ddx_dPareto = " << ret << endl;
        } else {
          if (verbose)
            cout << "ddx_pdf:" << endl << "  Using integral = " << ret << endl;
        }
      } else { // bad theta2
        abserr = NAN;
        termination_code = IntegrationController<myFloat>::bad_integrand;
        last = 0;
        if (fun_type>1){
          ret = -(1+alpha)*dPareto(x_m_zeta_input+zeta, alpha, beta_input, false)/x;
          if (verbose){
            cout<< "ddx_pdf:" << endl << "  Theta2 is bad & x is large so using ddx_dPareto = " << ret << endl;
          }
        } else {
          ret = dfdx_zeta;
          if (verbose)
            cout<< "ddx_pdf:" << endl << "  Theta2 is bad & x is small, so using dfdx_zeta = " << ret << endl;
        }
      } // bad theta2
      return ret;
  } // switch on dist_type
} //ddx_sdstaple1

} // namespace stable_distribution
