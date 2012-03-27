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





// order 2
template<> const double FChebRoots<2>::roots[] = {-0.707106781186548,
																									0.707106781186547};

// order 3
template<> const double FChebRoots<3>::roots[] = {-8.66025403784439e-01,
																									0.0,
																									8.66025403784438e-01};

// order 4
template<> const double FChebRoots<4>::roots[] = {-0.923879532511287,
																									-0.382683432365090,
																									0.382683432365090,
																									0.923879532511287};

// order 5
template<> const double FChebRoots<5>::roots[] = {-9.51056516295154e-01,
																								 -5.87785252292473e-01,
																									0.0,
																									5.87785252292473e-01,
																									9.51056516295154e-01};

// order 6
template<> const double FChebRoots<6>::roots[] = {-0.965925826289068,
																									-0.707106781186548,
																									-0.258819045102521,
																									0.258819045102521,
																									0.707106781186547,
																									0.965925826289068};

// order 7
template<> const double FChebRoots<7>::roots[] = {-9.74927912181824e-01,
																									-7.81831482468030e-01,
																									-4.33883739117558e-01,
																									0.0,
																									4.33883739117558e-01,
																									7.81831482468030e-01,
																									9.74927912181824e-01};

// order 8
template<> const double FChebRoots<8>::roots[] = {-0.980785280403230,
																									-0.831469612302545,
																									-0.555570233019602,
																									-0.195090322016128,
																									0.195090322016128,
																									0.555570233019602,
																									0.831469612302545,
																									0.980785280403230};

// order 9
template<> const double FChebRoots<9>::roots[] = {-9.84807753012208e-01,
																									-8.66025403784439e-01,
																									-6.42787609686539e-01,
																									-3.42020143325669e-01,
																									0.0,
																									3.42020143325669e-01,
																									6.42787609686539e-01,
																									8.66025403784438e-01,
																									9.84807753012208e-01,};

// order 10
template<> const double FChebRoots<10>::roots[] = {-0.987688340595138,
																									 -0.891006524188368,
																									 -0.707106781186548,
																									 -0.453990499739547,
																									 -0.156434465040231,
																									 0.156434465040231,
																									 0.453990499739547,
																									 0.707106781186547,
																									 0.891006524188368,
																									 0.987688340595138};

// order 11
template<> const double FChebRoots<11>::roots[] = {-9.89821441880933e-01,
																									 -9.09631995354518e-01,
																									 -7.55749574354258e-01,
																									 -5.40640817455598e-01,
																									 -2.81732556841430e-01,
																									 0.0,
																									 2.81732556841430e-01,
																									 5.40640817455597e-01,
																									 7.55749574354258e-01,
																									 9.09631995354518e-01,
																									 9.89821441880933e-01};

// order 12
template<> const double FChebRoots<12>::roots[] = {-0.991444861373810,
																									 -0.923879532511287,
																									 -0.793353340291235,
																									 -0.608761429008721,
																									 -0.382683432365090,
																									 -0.130526192220052,
																									 0.130526192220052,
																									 0.382683432365090,
																									 0.608761429008721,
																									 0.793353340291235,
																									 0.923879532511287,
																									 0.991444861373810};


// order 13
template<> const double FChebRoots<13>::roots[] = {-9.92708874098054e-01,
																									 -9.35016242685415e-01,
																									 -8.22983865893656e-01,
																									 -6.63122658240795e-01,
																									 -4.64723172043769e-01,
																									 -2.39315664287558e-01,
																									 0.0,
																									 2.39315664287557e-01,
																									 4.64723172043769e-01,
																									 6.63122658240795e-01,
																									 8.22983865893656e-01,
																									 9.35016242685415e-01,
																									 9.92708874098054e-01};






#endif
