// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.  
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info". 
// "http://www.gnu.org/licenses".
// ===================================================================================
#ifndef FCHEBROOTS_HPP
#define FCHEBROOTS_HPP

#include <cmath>
#include <limits>
#include <cassert>

#include "../../Utils/FNoCopyable.hpp"


/**
 * @author Matthias Messner (matthias.matthias@inria.fr)
 * Please read the license
 */

/**
 * @class FChebRoots
 *
 * The class @p FChebRoots provides the Chebyshev roots of order \f$\ell\f$
 * and the Chebyshev polynomials of first kind \f$T_n(x)\f$ and second kind
 * \f$U_{n-1}(x)\f$.
 *
 * @tparam ORDER interpolation order \f$\ell\f$
 */
template <int ORDER>
struct FChebRoots : FNoCopyable
{
	enum {order = ORDER}; //!< interpolation order

	/**
	 * Chebyshev roots in [-1,1] computed as \f$\bar x_n =
	 * \cos\left(\frac{\pi}{2}\frac{2n-1}{\ell}\right)\f$ for
	 * \f$n=1,\dots,\ell\f$
	 */
	const static double roots[];

  /**
   * Chebyshev polynomials of first kind \f$ T_n(x) = \cos(n \arccos(x)) \f$
	 *
	 * @param[in] n index
	 * @param[in] x coordinate in [-1,1]
	 * @return function value
   */
	static FReal T(const unsigned int n, FReal x)
	{
		//std::cout << x << std::endl;
    assert(std::fabs(x)-1.<10.*std::numeric_limits<FReal>::epsilon());
		if (std::fabs(x)>1.) {
			x = (x > FReal( 1.) ? FReal( 1.) : x);
      x = (x < FReal(-1.) ? FReal(-1.) : x);
    }

    return FReal(cos(n * acos(x)));
	}


  /**
	 * For the derivation of the Chebyshev polynomials of first kind \f$
   * \frac{\mathrm{d} T_n(x)}{\mathrm{d}x} = n U_{n-1}(x) \f$ the Chebyshev
   * polynomials of second kind \f$ U_{n-1}(x) = \frac{\sin(n
   * \arccos(x))}{\sqrt{1-x^2}} \f$ are needed.
	 *
	 * @param[in] n index
	 * @param[in] x coordinate in [-1,1]
	 * @return function value
   */
  static FReal U(const unsigned int n, FReal x)
  {
    assert(std::fabs(x)-1.<10.*std::numeric_limits<FReal>::epsilon());
		if (std::fabs(x)>1.) {
			x = (x > FReal( 1.) ? FReal( 1.) : x);
			x = (x < FReal(-1.) ? FReal(-1.) : x);
    }
  
    return FReal(n * (sin(n * acos(x))) / sqrt(1 - x*x));
  }
};

// We declare the roots here only once Please look to .cpp for definitions

// order 2
template<> const double FChebRoots<2>::roots[];

// order 3
template<> const double FChebRoots<3>::roots[];

// order 4
template<> const double FChebRoots<4>::roots[];

// order 5
template<> const double FChebRoots<5>::roots[];

// order 6
template<> const double FChebRoots<6>::roots[];

// order 7
template<> const double FChebRoots<7>::roots[];

// order 8
template<> const double FChebRoots<8>::roots[];

// order 9
template<> const double FChebRoots<9>::roots[];

// order 10
template<> const double FChebRoots<10>::roots[];

// order 11
template<> const double FChebRoots<11>::roots[];

// order 12
template<> const double FChebRoots<12>::roots[];

// order 13
template<> const double FChebRoots<13>::roots[];

// order 14
template<> const double FChebRoots<14>::roots[];


#endif
