/// \file random_test.cpp
/// Low level unit test of random stable
/// \author Joseph Dunn
/// \copyright 2017 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <iomanip>
using std::setw;
using std::setprecision;
using std::scientific;
using std::fixed;
using std::right;
#include <sstream>
using std::stringstream;

#include <chrono>
using std::chrono::high_resolution_clock;
using std::chrono::duration;

struct auto_timer {
  high_resolution_clock::time_point start;
  std::ostream& os;
  auto_timer(std::ostream& os) : start(high_resolution_clock::now()),
  os(os){}
  auto_timer() : start(high_resolution_clock::now()),
  os(std::cout) {}
  ~auto_timer() {
    duration<double> elapsed = high_resolution_clock::now() - start;
    os << "Elapsed time = " << setprecision(3) << fixed << elapsed.count() << " seconds" << endl;
  }
};

#include <boost/filesystem.hpp>


#include "stable_distribution.h"
using namespace stable_distribution;

static void show_usage (string name){
    boost::filesystem::path p(name);
    cerr << "Usage: " << p.filename().string() << "alpha beta pm u1 u2" << endl;
}

int main(int argc, const char * argv[]) {
    // Check the number of parameters
    if (argc != 6) {
        // Tell the user how to run the program
        show_usage(string(argv[0]));
        return 1;
    }
    
    stringstream ss1((string(argv[1]))), ss2((string(argv[2]))), ss3((string(argv[3]))),
                 ss4((string(argv[4]))), ss5((string(argv[5])));
    double alpha; ss1 >> alpha;
    double beta; ss2 >> beta;
    int pm_in; ss3 >> pm_in;
    Parameterization pm = (pm_in==0) ? S0 : S1;
    double u1; ss4 >> u1;
    double u2; ss5 >> u2;
    double result = random_stable<double>(alpha, beta, u1, u2, pm);
    cout.precision(16);
    cout << "alpha = " << alpha << endl
         << "beta = " << beta << endl
         << "pm = " << ((pm==S0)?"S0":"S1")<< endl
         << "u1 = " << u1 << endl
         << "u2 = " << u2 << endl << endl
         << "result = " << result << endl;
    return 0;
}
