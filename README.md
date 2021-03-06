Skew Stable Distributions
=========================

Introduction
------------
This package contains C++ code to calculate the skew stable distribution, a
generalization of the normal distribution appropriate for distributions with
power law tails.  The package started out as a port of the stabledist package
from R, but the algorithms have been extensively modified to improve
numerical stability.  There is a accompanying RcppStable package which
uses this C++ library to implement an R package with the same interface
as the original stabledist package.

Motivation
----------
To date the only open source C or C++ software available to calculate the 
stable distribution is the FMStable package from Geoff Robinson which is limited
to the extremal, i.e., maximally skewed, cases of the stable distribution.
This package can calculate the general case, i.e., -1 <= beta <= 1.

Installation
------------
This package was developed on a Mac running OS X v10.13.1 and Xcode v9.1.  It 
should compile using any C++ compiler that supports C++11. You'll need the
following to install the package.

1. The meson build system available at http://mesonbuild.com/Getting-meson.html
2. While the meson system has support for several backends, the ninja one 
   works best.  It's available at https://github.com/ninja-build/ninja/releases
3. The Boost headers and libraries available at http://www.boost.org
4. The Eigen 3 headers available at http://eigen.tuxfamily.org
5. For multiprecision support, which is used in some of the test programs, 
   you'll need the GNU mpfr library available at http://www.mpfr.org and
   the GNU gmp library available at https://gmplib.org.
6. The Pavel Holoborodko's Mpfr C++ header is used to wrap the GNU mpfr library. 
   The header needs to be modified to allow it to work with Boost by commenting out the 
   line #define MPREAL_HAVE_DYNAMIC_STD_NUMERIC_LIMITS. 
   The original is available at http://www.holoborodko.com/pavel/mpfr/..

Most of these prerequisites can be downloaded using a package management system.
I used MacPorts to get them on my Mac.

Meson has built in dependency support for the Boost libraries, but I had
to define the BOOST_ROOT to get it to work.  
For the other prequisites meson
works best using the pkconfig files.  For some reason mpfr and gmp came without
such files, but it's relatively easy to add appropriate entries.  Once that's
done Meson should do the rest.  I wanted to use the clang++ compiler so
I used the following commands from the root directory:

    $ mkdir build
    $ CXX=clang++ BOOST_ROOT=/opt/local meson build
    $ cd build
    $ ninja
    $ ninja test

There are also three options available for adding a multiprecision support, 
which are implemented by defining MPREAL, CPP_BIN_FLOAT, and/or MPFR_FLOAT. As 
distributed they are all defined.

The output of the tests will be in files in the 
output-$PACKAGE_OS-$PACKAGE_COMPILER-$PACKAGE_VERSION- directory of the parent directory.

If you have a lot of time on your hands, you can run the following from the
build directory:

    $ ./stable_dump/stable_dump float_type

Where float_type is mpreal, mpfr_float or cpp_bin_float.  
This will create a 224,000 line file 'stable_float_type.out' in the output directory,
which can be used to check the double precision calculations by issuing the 
following command from the build directory:

    $ ./xcheck_to_file/xcheck_to_file float_type.
   
The src files beginning with the word "trace" can be used as examples of how to
call the programs.

Documentation
-------------
The header files of the package are documented using the Doxygen system.  The 
documentation is generated automatically by ninja and can be found in the doc/html
and doc/latex subdirectories of the build directories.

The html version of the documentation can be accessed at doc/html/index.html. 

Acknowledgements
----------------
This package builds on the work of a number of others.

1. John Nolan developed the integral representation of the skew stable
distribution that's used in this package
[density.pdf](http://fs2.american.edu/jpnolan/www/stable/density.pdf).
He has released a compiled version 3.12.02 of his Fortran 
program STABLE,
[stable.exe](http://academic2.american.edu/~jpnolan/stable/stable.exe).
It's not open source and he's not explicit about its license. 
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
6. The Eigen 3 package of headers is used.  It's available at http://eigen.tuxfamily.org.
7. Pavel Holoborodko's Mpfr C++ header is used to wrap the GNU mpfr library for the
mpreal float type.
The header needs to be modified to allow it to work with Boost by commenting out the 
line #define MPREAL_HAVE_DYNAMIC_STD_NUMERIC_LIMITS. 
The original is available at http://www.holoborodko.com/pavel/mpfr/.
8. The GNU mpfr package is used for the multiprecision version.  It's available at http://www.mpfr.org/mpfr-current/#download.
9. Three components, meta.h, Problem.h and NelderMeadSolver.h, of Patrick Wieschollek's CppNumericalSolver package are used 
by the stable_fit routine and are included in the present package.  The original
is available at https://github.com/PatWie/CppNumericalSolvers.
10. The documentation is prepared using Doxygen, which is available at http://doxygen.org, and 
mathjax



