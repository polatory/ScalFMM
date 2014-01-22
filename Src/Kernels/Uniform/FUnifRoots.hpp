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
#ifndef FUNIFROOTS_HPP
#define FUNIFROOTS_HPP

#include <cmath>
#include <limits>
#include <cassert>

#include "../../Utils/FNoCopyable.hpp"


/**
 * @author Pierre Blanchard (pierre.blanchard@inria.fr)
 * Please read the license
 */

/**
 * @class FUnifRoots
 *
 * The class @p FUnifRoots provides the equispaced roots of order \f$\ell\f$
 * and the Lagrange polynomials \f$L_n(x)\f$.
 *
 * @tparam ORDER interpolation order \f$\ell\f$
 */
template <int ORDER>
struct FUnifRoots : FNoCopyable
{
  enum {order = ORDER}; //!< interpolation order

  /**
   * Lagrange roots in [-1,1] computed as \f$\bar x_n =
   * -1 + 2\frac{n-1}{\ell}\f$ for \f$n=1,\dots,\ell\f$
   */
  const static double roots[];

  /**
   * Lagrange polynomials \f$ L_n(x) = \Pi_{m=1 \atop m\neq n}^\ell \frac{x-\bar x_m}{\bar x_n-\bar x_m} \f$
   *
   * @param[in] n index
   * @param[in] x coordinate in [-1,1]
   * @return function value
   */
  static FReal L(const unsigned int n, FReal x)
  {
    //std::cout << x << std::endl;
    assert(std::fabs(x)-1.<10.*std::numeric_limits<FReal>::epsilon());
    if (std::fabs(x)>1.) {
      x = (x > FReal( 1.) ? FReal( 1.) : x);
      x = (x < FReal(-1.) ? FReal(-1.) : x);
    }

    FReal L=FReal(1.);
    for(unsigned int m=0;m<order;++m){
      if(m!=n)
        L *= (x-FUnifRoots<order>::roots[m])/(FUnifRoots<order>::roots[n]-FUnifRoots<order>::roots[m]);
      
    }

    return FReal(L);
  }


  /**
   * For the derivation of the Lagrange polynomials
   * \f$ \frac{\mathrm{d} L_n(x)}{\mathrm{d}x} = ... \f$
   *
   * @param[in] n index
   * @param[in] x coordinate in [-1,1]
   * @return function value
   */
  static FReal dL(const unsigned int n, FReal x)
  {
    assert(std::fabs(x)-1.<10.*std::numeric_limits<FReal>::epsilon());
    if (std::fabs(x)>1.) {
      x = (x > FReal( 1.) ? FReal( 1.) : x);
      x = (x < FReal(-1.) ? FReal(-1.) : x);
    }

    // optimized variant
    FReal NdL=FReal(0.);// init numerator
    FReal DdL=FReal(1.);// init denominator
    FReal tmpNdL;
    for(unsigned int p=0;p<order;++p){
      if(p!=n){
        tmpNdL=FReal(1.);
        for(unsigned int m=0;m<order;++m)
          if(m!=n && m!=p)
            tmpNdL*=x-FUnifRoots<order>::roots[m];
        NdL+=tmpNdL;
        DdL*=FUnifRoots<order>::roots[n]-FUnifRoots<order>::roots[p];
      }//endif
    }// p

    return FReal(NdL/DdL);

  }
};

// We declare the roots here only once Please look to .cpp for definitions

// order 2
template<> const double FUnifRoots<2>::roots[];

// order 3
template<> const double FUnifRoots<3>::roots[];

// order 4
template<> const double FUnifRoots<4>::roots[];

// order 5
template<> const double FUnifRoots<5>::roots[];

// order 6
template<> const double FUnifRoots<6>::roots[];

// order 7
template<> const double FUnifRoots<7>::roots[];

// order 8
template<> const double FUnifRoots<8>::roots[];

// order 9
template<> const double FUnifRoots<9>::roots[];

// order 10
template<> const double FUnifRoots<10>::roots[];

// order 11
template<> const double FUnifRoots<11>::roots[];

// order 12
template<> const double FUnifRoots<12>::roots[];

// order 13
template<> const double FUnifRoots<13>::roots[];

// order 14
template<> const double FUnifRoots<14>::roots[];


#endif
