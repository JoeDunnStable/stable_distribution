// \file zolotarev_integral.cpp
/// \author Joseph Dunn
/// \copyright 2017 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

// Allows the calculation of pdf and cdf for alpha = 1 via integrals

// Theis will be included at the end of zolotarev.h if LIBRARY is defined
#include "adaptive_integration.h"

namespace stable_distribution {

template<typename myFloat>
void Q_integrand(myFloat* u, int n, void* ex){
  Zolotarev<myFloat>* zp = (Zolotarev<myFloat>*)ex;
  myFloat eps0 = 50*std::numeric_limits<myFloat>::min();
  myFloat eps = 50*std::numeric_limits<myFloat>::epsilon();
  myFloat& x{zp->x_m_zet};
  myFloat& beta{zp->beta};
  myFloat& pi{Zolotarev<myFloat>::pi};
  for (int i=0; i<n; ++i) {
    myFloat uu{u[i]/(1-u[i])};
    if (u[i] > eps0 && fabs(u[i]-1) > eps)
      // Formula from Zolotarev's Corollary 1 to Theorem 2.2.1
      // with slight modification to remap to finite interval
      u[i] = exp(-x*uu-beta*(2/pi)*uu*log(uu)) * sin((1+beta)*uu)
            /(uu*pi*(1-u[i])*(1-u[i]));
    else if (fabs(u[i]-1)<eps)
      u[i] = 0;
    else
      u[i] = (1+beta)/pi;
  }
}

template<typename myFloat>
void q_integrand(myFloat* u, int n, void* ex) {
  Zolotarev<myFloat>* zp = (Zolotarev<myFloat>*)ex;
  myFloat eps0 = 50*std::numeric_limits<myFloat>::min();
  myFloat eps = 50*std::numeric_limits<myFloat>::epsilon();
  myFloat& x{zp->x_m_zet};
  myFloat& beta{zp->beta};
  myFloat& pi{Zolotarev<myFloat>::pi};
  for (int i=0; i<n; ++i) {
    myFloat uu{u[i]/(1-u[i])};
    if (u[i] > eps0 && fabs(u[i]-1) > eps)
      // Zolotarev Formula 2.2.9 using sin as Im(e^ix)
      // and adjusting to S1 parameterization from Zolotarev's B
      u[i] = exp(-x*uu-beta*(2/pi)*uu*log(uu)) * sin((1+beta)*uu)
             /(pi*(1-u[i])*(1-u[i]));
    else
      u[i] = 0;
  }
}
  
template<typename myFloat>
void ddx_q_integrand(myFloat* u, int n, void* ex) {
  Zolotarev<myFloat>* zp = (Zolotarev<myFloat>*)ex;
  myFloat eps0 = 50*std::numeric_limits<myFloat>::min();
  myFloat eps = 50*std::numeric_limits<myFloat>::epsilon();
  myFloat& x{zp->x_m_zet};
  myFloat& beta{zp->beta};
  myFloat& pi{Zolotarev<myFloat>::pi};
  for (int i=0; i<n; ++i) {
    myFloat uu{u[i]/(1-u[i])};
    if (u[i] > eps0 && fabs(u[i]-1) > eps)
      // Zolotarev Formula 2.2.9 using sin as Im(e^ix)
      // and adjusting to S1 parameterization from Zolotarev's B
      u[i] = -uu*exp(-x*uu-beta*(2/pi)*uu*log(uu)) * sin((1+beta)*uu)
      /(pi*(1-u[i])*(1-u[i]));
    else
      u[i] = 0;
  }
}

} // namespace stable_distribution
