/// \file kolmogorov.h
/// The cdf of the Kolmogorov Smirnov statistic
/// \author Joseph Dunn
/// \copyright 2017 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#ifndef kolmogorov_h
#define kolmogorov_h

#include <vector>
using std::vector;

/// Return the complementary cummulative probability for D*sqrt(n) for large n
double kolmogorov_asymptotic_cdf(double x);

/// Returns the exact distribution of Kolmogorov Smirnov for all n
double kolmogorov_cdf(int n, double d);

#endif
