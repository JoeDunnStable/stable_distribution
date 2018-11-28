///
///  bracket_and_solve.h
///  Modification of boost bracket_and_solve_root
/// \author Joseph Dunn
/// \copyright 2016, 2017, 2018 Joseph Dunn
/// \copyright Distributed under the terms of the GNU General Public License version 3

#ifndef bracket_and_solve_h
#define bracket_and_solve_h

#include <boost/math/tools/toms748_solve.hpp>

namespace boost { namespace math { namespace tools {
  
  /** modification of boost braket and solve root using -infinity as lower bound not 0 */
  template <class F, class T, class Tol, class Policy>
  std::pair<T, T> bracket_and_solve_root2(F& f, const T& guess, T factor, bool rising, Tol tol, boost::uintmax_t& max_iter, const Policy& pol)
  {
    BOOST_MATH_STD_USING
    static const char* function = "bracket_and_solve_root2<%1%>";
    //
    // Set up inital brackets:
    //
    T a = guess;
    T b = a;
    T fa = f(a);
    T fb = fa;
    //
    // Set up invocation count:
    //
    boost::uintmax_t count = max_iter - 1;
    
    int step = 32;
    
    if((fa < 0) == rising)
    {
      //
      // Zero is to the right of b, so walk upwards
      // until we find it:
      //
      while((boost::math::sign)(fb) == (boost::math::sign)(fa))
      {
        if(count == 0)
          return boost::math::detail::pair_from_single(policies::raise_evaluation_error(function, "Unable to bracket root, last nearest value was %1%", b, pol));
        //
        // Heuristic: normally it's best not to increase the step sizes as we'll just end up
        // with a really wide range to search for the root.  However, if the initial guess was *really*
        // bad then we need to speed up the search otherwise we'll take forever if we're orders of
        // magnitude out.  This happens most often if the guess is a small value (say 1) and the result
        // we're looking for is close to std::numeric_limits<T>::min().
        //
        if((max_iter - count) % step == 0)
        {
          factor *= 2;
          if(step > 1) step /= 2;
        }
        //
        // now go ahead and move our guess by "factor":
        //
        a = b;
        fa = fb;
        b += factor;
        fb = f(b);
        --count;
        BOOST_MATH_INSTRUMENT_CODE("a = " << a << " b = " << b << " fa = " << fa << " fb = " << fb << " count = " << count);
      }
    }
    else
    {
      //
      // Zero is to the left of a, so walk downwards
      // until we find it:
      //
      while((boost::math::sign)(fb) == (boost::math::sign)(fa))
      {
        if(count == 0)
          return boost::math::detail::pair_from_single(policies::raise_evaluation_error(function, "Unable to bracket root, last nearest value was %1%", a, pol));
        //
        // Heuristic: normally it's best not to increase the step sizes as we'll just end up
        // with a really wide range to search for the root.  However, if the initial guess was *really*
        // bad then we need to speed up the search otherwise we'll take forever if we're orders of
        // magnitude out.  This happens most often if the guess is a small value (say 1) and the result
        // we're looking for is close to numeric_limits<T>::min().
        //
        if((max_iter - count) % step == 0)
        {
          factor *= 2;
          if(step > 1) step /= 2;
        }
        //
        // n_gaussow go ahead and move are guess by "factor":
        //
        b = a;
        fb = fa;
        a -= factor;
        fa = f(a);
        --count;
        BOOST_MATH_INSTRUMENT_CODE("a = " << a << " b = " << b << " fa = " << fa << " fb = " << fb << " count = " << count);
      }
    }
    max_iter -= count;
    max_iter += 1;
    std::pair<T, T> r = toms748_solve(
                                      f, a, b, fa, fb,
                                      tol, count, pol);
    max_iter += count;
    BOOST_MATH_INSTRUMENT_CODE("max_iter = " << max_iter << " count = " << count);
    BOOST_MATH_LOG_COUNT(max_iter)
    return r;
  }
  
  template <class F, class T, class Tol>
  inline std::pair<T, T> bracket_and_solve_root2(F& f, const T& guess, const T& factor, bool rising, Tol tol, boost::uintmax_t& max_iter)
  {
    return bracket_and_solve_root2(f, guess, factor, rising, tol, max_iter, policies::policy<>());
  }
  
} // namespace tools
} // namespace math
} // namespace boost


#endif /* bracket_and_solve_h */
