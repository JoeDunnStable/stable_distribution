Skew Stable Distributions
=========================
Introduction
------------
This package contains C++ code to calculate the skew stable distribution, a
generalization of the normal distribution appropriate for distributions with
power law tails.  The package started out as a port of the stabledist package
from R, but the algorithms have been extensively modified to improve
numerical stability.
Motivation
----------
To date the only open source C or C++ software available to calculate the 
stable distribution is the FMStable package from Geoff Robinson which is limited
to the extremal, i.e., maximally skewed, cases of the stable distribution.
This package can calculate the general case, i.e., -1 <= beta <= 1.
Installation
------------
This package was developed on a Mac running OS X v10.12.6 and Xcode v8.3.3.  It 
should compile using any C++ compiler that supports C++11. You'll need the
following to install the package.

1. The meson build system available at http://mesonbuild.com/Getting-meson.html
2. While the meson system has support for several backends, the ninja one 
   works best.  It's available at https://github.com/ninja-build/ninja/releases
3. The Boost headers and libraries available at http://www.boost.org
4. The Eigen 3 headers available at http://eigen.tuxfamily.org
5. For multiprecision support, which is used in some of the test programs, 
   you'll need eithhe GNU mpfr library available at http://www.mpfr.org and
   the GNU gmp library available at https://gmplib.org

Most of these prerequisites can be downloaded using a package management system.
I used MacPorts to get them on my Mac.

Meson has built in dependency support for the Boost libraries, but I had
to define the BOOST_ROOT to get it to work.  For the other prequisites meson
works best using the pkconfig files.  For some reason mpfr and gmp came without
such files, but it's relatively easy to add appropriate entries.  Once that's
done Meson should do the rest.  I wanted to use the clang++ compiler so
I used the following commands from the root directory:

    $ mkdir build_double
    $ CXX=clang++ BOOST_ROOT=/opt/local meson build_double
    $ cd build_double
    $ ninja
    $ ninja test

There are also three options available for building a multiprecision version:
'mpreal', 'mpfr_float', or 'cpp_bin_float'.  I generally use the mpreal version, 
which involves the following commands starting from the parent directory:

    $ mkdir build_mpreal
    $ CXX=clang++ BOOST_ROOT=/opt/local meson build_mpreal
    $ cd build_mpreal
    $ mesonconf -D floating_type='mpreal'
    $ ninja
    $ ninja test

The ninja test suite includes only the more rapidly running tests.  Other 
testing programs can be found in the test subdirectory of each build director.
They should be run from the parent build directory.  For instance, the 
following command can be entered from the build_double directory:

    $ ./test/test_pstable

The output will be in a file in the output directory of the parent directory.

If you have a lot of time on your hands, you can run the following from the
build_mpreal directory:

    $ ./test/dump_stable

This will create a 224,000 line file 'stable_mpreal.out' in the output directory,
which can be used to check the double precision calculations by issuing the 
following command from the build_double directory:

    $ ./test/xcheck_to_file
   
The src files beginning with the word "trace" can be used as examples of how to
call the programs.

Documentation
-------------
The header files of the package are documented using the Doxygen system.  The 
documentation is generated automatically by ninja and can be found in the html
and latex subdirectories of the build directories.

The html version of the documentation can be accessed at html/index.html. 

Acknowledgements
----------------
This package builds on the work of a number of others.

1. John Nolan developed the integral representation of the skew stable
distribution that's used in this package
[density.pdf](http://fs2.american.edu/jpnolan/www/stable/density.pdf).
He has released a compiled version 3.12.02 of his Fortran 
program STABLE,
[stable.exe](http://academic2.american.edu/~jpnolan/stable/stable.exe).
It's not open source and he's not explicit about its license.  The two of the
test programs compare the output of the present package to the output of Nolan's
package.
2. Many of the routines in this package started out life as C++ translations of
the stabledist package for R, the source code for which is available at 
[stabledist_0.7-0](https://cran.r-project.org/src/contrib/stabledist_0.7-0.tar.gz)
3. The routines in the adaptive_integration routine started out life
as machine C++ translations of Fortran routines in QUADPACK, which is part of 
SLATEC and therefore in the public domain (http://en.wikipedia.org/wiki/QUADPACK).
The routines were then heavily modified to take advantage of the C++ language.
4. One of the modifications made to QUADPACK is the addition of the ability to calculate
the nodes and weights for the Gauss Kronrod integration on the fly.  For this purpose
Dirk Laurie's method is used [kronrod.ps](http://dip.sun.ac.za/~laurie/papers/kronrod/kronrod.ps).
The routines used here are C++ translations of Dirk Laurie's MATLAB code, which
is included in Walter Gautschi's OPQ suite [OPQ](https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html).
5. Some of the Boost header files are used.  They're available at http://www.boost.org.
6. The Eigen 3 package of headers is used.  They're available at http://eigen.tuxfamily.org.
7. The Pavel Holoborodko's Mpfr C++ header is used to wrap the GNU mpfr library. 
A copy of the header, mpreal.h, with modifications to work with Boost is included with 
this package.  The original is available at http://www.holoborodko.com/pavel/mpfr/.
8. The GNU mpfr package is used for the multiprecision version.  It's available at http://www.mpfr.org/mpfr-current/#download.
9. Three components, meta.h, Problem.h and NelderMeadSolver.h, of Patrick Wieschollek's CppNumericalSolver package are used 
by the stable_fit routine and are included in the present package.  The original
is available at https://github.com/PatWie/CppNumericalSolvers.

Details
-------
The function uses the approach of J.P. Nolan for general stable
distributions. Nolan (1997) derived expressions in form of integrals
based on the characteristic function for standardized stable random
variables. For `StandardStableDistribution::pdf` and
`StandardStableDistribution::cdf`, these integrals
are numerically evaluated.

There are several different parameterizations in use for the stable
distribution.  Nolan derived his integral representation using the parameterization
he labeled "S0" "[pm=0], which is based on the (M) representation
of Zolotarev.  Unlike the Zolotarev (M) parameterization, gamma and
delta are straightforward scale and shift parameters. This
representation is continuous in all 4 parameters.

Definition
----------
The random variable \f$ Y \f$ is said to be distriubuted according to the
\f$ S(\alpha, \beta, \gamma, \delta, 0) \f$ if its characteristic function is
given by the following:
\f{equation}{
{E \exp(itY)} =
\begin{cases}
\exp \{ -\gamma^\alpha |t|^\alpha \left [1+i\beta(\tan{\frac{\pi \alpha}{2}})(\text{sign } t)(|\gamma t|^{1-\alpha}-1) \right ]+i \delta t \} &\alpha \ne 1 \\
\exp \{-\gamma |t| \left[ 1+i \beta \frac {2}{\pi}(\text{sign }t)\log{(\gamma |t|)} \right ]+i \delta t \} & \alpha=1.
\end{cases}
\f}
where
\par
\f$ \alpha \f$, the shape parameter, satisfies \f$ 0 \lt \alpha \le 2 \f$,
\par
\f$ \beta \f$, the skewness parameter, satisfies \f$ -1 \le \beta \le 1 \f$,
\par
\f$ \gamma \f$, the scale parameter, satisfies \f$ 0 \lt \gamma \f$, and
\par
\f$ \delta \f$, the location parameter, is unconstrained.

Setting \f$\gamma = 1\f$ and \f$\delta =0\f$ gives the standard skew stable distribution whose characteristic function is therefore:
\f{equation}{
{E \exp(itY)} =
\begin{cases}
\exp \{ -|t|^\alpha \left [1+i\beta(\tan{\frac{\pi \alpha}{2}})(\text{sign } t)(|t|^{1-\alpha}-1) \right ] \} &\alpha \ne 1 \\
\exp \{-|t| \left[ 1+i \beta \frac {2}{\pi}(\text{sign }t)\log{(|t|)} \right ]\} & \alpha=1.
\end{cases}
\f}



