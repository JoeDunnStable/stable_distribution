///
/// @file stable_distribution_random.h
/// Implementation of random generator for standard stable distribution
/// Included in stable_distribution.h when LIBRARY is defined
/// \author Joseph Dunn
/// \copyright 2016, 2017, 2018 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

namespace stable_distribution {
  
template<typename myFloat>
myFloat random_stable(myFloat alpha, myFloat beta,
                        myFloat u1, myFloat u2,
                        Parameterization pm)
  {
    // Description:
    //	 Returns one random variates for standard stable DF
    if (!StandardStableDistribution<myFloat>::initialized)
      StandardStableDistribution<myFloat>::initialize();
    auto &pi =  StandardStableDistribution<myFloat>::pi;
    auto &pi2 =  StandardStableDistribution<myFloat>::pi2;
    
    // Calculate uniform and exponential distributed random numbers:
    myFloat theta = pi * (u1-1./2.);
    myFloat w = -log(u2);
    
    myFloat result;
    if (alpha == 1 ) {
      result = (1/pi2)*((pi2+beta*theta)*tan(theta)-beta*log(pi2*w*cos(theta)/(pi2+beta*theta)));
    } else if (alpha==0) {
      if (w > 1)
        return 0;
      else if (theta + beta * pi2 > 0)
        return std::numeric_limits<myFloat>::infinity();
      else
        return -std::numeric_limits<myFloat>::infinity();
    } else if (fabs(alpha-1)< 1./128.) {
      myFloat a_m_1 = alpha - 1;
      myFloat zeta = beta/tan(pi2*a_m_1);
      myFloat cos_a_th0 = pow(1+zeta*zeta,-.5);
      myFloat sin_a_th0 = ((a_m_1>0) ? -1 : 1) * sqrt(1-cos_a_th0*cos_a_th0);
/*
      cout.precision(16);
      cout << "a_m_1 = " << a_m_1 << endl
      << "zeta = " << zeta << endl
      << "theta0 = " << theta0 << endl
      << "cos_a_th0 = " << cos_a_th0 << endl
      << "sin_a_th0 = " << sin_a_th0 << endl
      << "pm = " << pm << endl;
 */
      myFloat cos_th_ath =cos(a_m_1*theta)*cos_a_th0-sin(a_m_1*theta)*sin_a_th0;
      myFloat eps1 = expm1(a_m_1/alpha*(log(cos_a_th0)+log(cos(theta))-log(cos_th_ath/w)));
      myFloat term1 =(sin(a_m_1*theta)+cos(a_m_1*theta)*tan(theta))+zeta*sin(a_m_1*theta)*tan(theta);
      term1 *= (1+eps1);
      myFloat term2 = -zeta*(cos(a_m_1*theta))*eps1;
      myFloat term3;
      if (pm==S0)
        term3 = zeta*2*pow(sin(a_m_1*theta/2),2);
      else
        term3 = -zeta*cos(a_m_1*theta);
      result = term1 + term2 +term3;
/*
      cout << "cos_th_ath = " << cos_th_ath << endl
      << " eps1 = " << eps1 << endl
      << "term1 = " << term1 << endl
      << "term2 = " << term2 << endl
      << "term3 = " << term3 << endl
      << "result = " << result << endl;
 */
    } else {
      myFloat b_tan_pa = beta*tan(pi2*alpha);
      myFloat theta0 = min(max<myFloat>(-pi2, atan(b_tan_pa) / alpha), pi2);
      myFloat c = pow(1+pow(b_tan_pa,2),1/(2*alpha));
      myFloat a_tht = alpha*(theta+theta0);
      myFloat r = c*sin(a_tht)*exp(-log(cos(theta))/alpha+((1-alpha)/alpha)*log(cos(theta-a_tht)/w));
      if (pm == S0)
        result = r - b_tan_pa;
      else
        result = r;
    }
    return result;
  }
} // namespace stable_distribution

