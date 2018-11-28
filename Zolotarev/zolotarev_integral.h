/// \file zolotarev_integral.h
/// Implemenation of convergent integral used for alpha = 1 & small x
/// Included in zolotarev.h when LIBRARY is defined
/// \author Joseph Dunn
/// \copyright 2017, 2018 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

// Allows the calculation of pdf and cdf for alpha = 1 via integrals

// Theis will be included at the end of zolotarev.h if LIBRARY is defined
#include "adaptive_integration.h"

namespace stable_distribution {

template<typename myFloat>
  void Q_integrand<myFloat>::operator() (vector<myFloat>& u){
  myFloat eps0 = 50*std::numeric_limits<myFloat>::min();
  myFloat eps = 50*std::numeric_limits<myFloat>::epsilon();
  myFloat& x{zol.x_m_zet};
  myFloat& beta{zol.beta};
  myFloat& pi{Zolotarev<myFloat>::pi};
  size_t n = u.size();
  for (int i=0; i<n; ++i) {
    myFloat uu{u.at(i)/(1-u.at(i))};
    if (u.at(i) > eps0 && fabs(u.at(i)-1) > eps)
      // Formula from Zolotarev's Corollary 1 to Theorem 2.2.1
      // with slight modification to remap to finite interval
      u.at(i) = exp(-x*uu-beta*(2/pi)*uu*log(uu)) * sin((1+beta)*uu)
            /(uu*pi*(1-u.at(i))*(1-u.at(i)));
    else if (fabs(u.at(i)-1)<eps)
      u.at(i) = 0;
    else
      u.at(i) = (1+beta)/pi;
  }
}

template<typename myFloat>
void q_integrand<myFloat>::operator() (vector<myFloat>& u) {
  myFloat eps0 = 50*std::numeric_limits<myFloat>::min();
  myFloat eps = 50*std::numeric_limits<myFloat>::epsilon();
  myFloat& x{zol.x_m_zet};
  myFloat& beta{zol.beta};
  myFloat& pi{Zolotarev<myFloat>::pi};
  size_t n = u.size();
  for (int i=0; i<n; ++i) {
    myFloat uu{u.at(i)/(1-u.at(i))};
    if (u.at(i) > eps0 && fabs(u.at(i)-1) > eps)
      // Zolotarev Formula 2.2.9 using sin as Im(e^ix)
      // and adjusting to S1 parameterization from Zolotarev's B
      u.at(i) = exp(-x*uu-beta*(2/pi)*uu*log(uu)) * sin((1+beta)*uu)
             /(pi*(1-u.at(i))*(1-u.at(i)));
    else
      u.at(i) = 0;
  }
}
  
template<typename myFloat>
  void ddx_q_integrand<myFloat>::operator() (vector<myFloat>& u) {
  myFloat eps0 = 50*std::numeric_limits<myFloat>::min();
  myFloat eps = 50*std::numeric_limits<myFloat>::epsilon();
  myFloat& x{zol.x_m_zet};
  myFloat& beta{zol.beta};
  myFloat& pi{Zolotarev<myFloat>::pi};
    size_t n = u.size();
  for (int i=0; i<n; ++i) {
    myFloat uu{u.at(i)/(1-u.at(i))};
    if (u.at(i) > eps0 && fabs(u.at(i)-1) > eps)
      // Zolotarev Formula 2.2.9 using sin as Im(e^ix)
      // and adjusting to S1 parameterization from Zolotarev's B
      u.at(i) = -uu*exp(-x*uu-beta*(2/pi)*uu*log(uu)) * sin((1+beta)*uu)
      /(pi*(1-u.at(i))*(1-u.at(i)));
    else
      u.at(i) = 0;
  }
}

} // namespace stable_distribution
