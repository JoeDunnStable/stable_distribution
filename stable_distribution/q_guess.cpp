///
/// \file q_guess.cpp
/// Implementation of rough guess for quantile of stable distribution
/// \author Joseph Dunn
/// \copyright 2016, 2017, 2018 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include "q_guess.h"
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/tools/toms748_solve.hpp>
#include <boost/math/constants/constants.hpp>
#include <iostream>
using std::cout;
using std::endl;
using std::ostream;
#include <iomanip>
using std::setw;
using std::setprecision;
using std::scientific;

namespace stable_distribution {

using boost::math::tools::toms748_solve;
using boost::math::students_t_distribution;
using boost::math::cdf;
using boost::math::quantile;
using std::pair;

//Functor which contains an approximation for stable p as a function of p for
//student t with df=alpha
class pt_solve {
private:
    double rplus;
    double rminus;
    double knot1;
    double knot2;
    double a0;
    double a1;
    double a2;
    double a3;
    double p;
    static const double pi;
    static const double pi2;
public:
    pt_solve(double pp, double alpha, double beta, int lower_tail, int log_p){
        double c_stable_plus = sin(pi2*alpha )*tgamma(alpha)/pi*alpha*(1+beta);
        double c_stable_minus = sin(pi2*alpha )*tgamma(alpha)/pi*alpha*(1-beta);
        double c_t=tgamma((alpha+1)/2)/(sqrt(alpha*pi)*tgamma(alpha/2))*pow(alpha,((alpha+1)/2));
        rplus=c_stable_plus/c_t;
        rminus=c_stable_minus/c_t;
        // construct a cubic spline for the mapping of pt to cdf
        students_t_distribution<double> st(alpha);
        if (alpha<1 && beta==1){
            knot1 = cdf(st,-tan(pi2*alpha));
        } else
            knot1 = .01;
        if (alpha<1 && beta==-1) {
            knot2 = cdf(st,tan(pi2*alpha));
        } else
            knot2 = .99;
        double dk=knot2-knot1;
        a0 = rminus*knot1;
        a1 = rminus;
        double b0 = 1-rplus*(1-knot2)-rminus*knot2;
        double b1 = rplus-rminus;
        a2 = -(dk*b1-3*b0)/(dk*dk);
        a3 = (dk*b1-2*b0)/(dk*dk*dk);
        /*    cout << "rminus = " << rminus << ",knot1 = " << knot1 << endl
         << "rplus = " << rplus << ", knot2 = " << knot2 << endl
         << "a0 = " << a0 << ", a1 = " << a1 << ", a2 = " << a2 << ", a3 = " << a3 << endl;
         */
        p = (log_p) ? exp(pp) : pp;
        p = (lower_tail) ? p : 1-p;
    }
    double operator () (double pt) {
        if (pt<=knot1)
            return pt*rminus-p;
        else if (pt>=knot2) {
            return 1-rplus*(1-pt)-p;
        } else {
            double ptmk1 = pt-knot1;
            double r = (((a3*ptmk1)+a2)*ptmk1+a1)*ptmk1+a0;
            return r-p;
        }
    }
	friend
	ostream& operator<< (ostream& os, pt_solve& solver) {
	os << "rplus  = " << setprecision(16) << scientific << solver.rplus << endl
		<< "rminus = " << setprecision(16) << scientific << solver.rminus << endl
		<< "knot1  = " << setprecision(16) << scientific << solver.knot1 << endl
		<< "knot2  = " << setprecision(16) << scientific << solver.knot2 << endl
		<< "a0     = " << setprecision(16) << scientific << solver.a0 << endl
		<< "a1     = " << setprecision(16) << scientific << solver.a1 << endl
		<< "a2     = " << setprecision(16) << scientific << solver.a2 << endl
		<< "a3     = " << setprecision(16) << scientific << solver.a3 << endl
		<< "p      = " << setprecision(16) << scientific << solver.p << endl
		<< "pi     = " << setprecision(16) << scientific << pt_solve::pi << endl
		<< "pi2    = " << setprecision(16) << scientific << pt_solve::pi2 << endl;

	return os;
	}
    
};

const double pt_solve::pi = boost::math::constants::pi<double>();
const double pt_solve::pi2 = boost::math::constants::half_pi<double>();

class rel_eps_tolerance
{


public:
    rel_eps_tolerance(double eps) : eps(eps) {};
    inline bool operator()(const double& a, const double& b)
    {
        return fabs(a - b) <= (eps * (std::min)(fabs(a), fabs(b)));
    }
private:
    double eps;
};


using boost::math::policies::policy;
using boost::math::policies::overflow_error;
using boost::math::policies::ignore_error;

typedef policy<overflow_error<ignore_error> > my_policy;

double q_guess(double p,double alpha,double beta,int lower_tail,int log_p){
    pt_solve pt_s(p,alpha,beta,lower_tail,log_p);
	students_t_distribution<double, my_policy> st(alpha);
    double lower = 0;
    double upper = 1;
    rel_eps_tolerance tol(1e-6);
    uintmax_t maxiter = 200;
    pair<double,double> r;
    r = toms748_solve(pt_s,lower,upper,tol,maxiter);
    double pt_=(r.first+r.second)/2;
    return quantile(st, fmax(1e-15,fmin(pt_,1-1e-15)));
}
} // namespace stable_distribution
