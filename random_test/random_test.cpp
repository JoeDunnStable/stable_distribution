//
//  main.cpp
//  random_test
//
//  Created by Joseph Dunn on 11/18/17.
//  Copyright © 2017 Joseph Dunn. All rights reserved.
//

#include <iostream>
using std::cout;
using std::endl;
#include <iomanip>
using std::setw;
using std::setprecision;
using std::scientific;
using std::fixed;
using std::right;
#include <sstream>
using std::stringstream;

#include <boost/timer/timer.hpp>
using boost::timer::auto_cpu_timer;

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
