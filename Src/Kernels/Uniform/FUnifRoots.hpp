// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas, Matthias Messner
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
// Keep in private GIT
// @SCALFMM_PRIVATE

#ifndef FUNIFROOTS_HPP
#define FUNIFROOTS_HPP

#include <cmath>
#include <limits>
#include <cassert>

#include "../../Utils/FNoCopyable.hpp"
#include "../../Utils/FMath.hpp"

#include <array>

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
template < class FReal, int ORDER>
struct FUnifRoots : FNoCopyable
{
    enum {order = ORDER}; //!< interpolation order

    /**
   * Lagrange roots in [-1,1] computed as \f$\bar x_n =
   * -1 + 2\frac{n-1}{\ell}\f$ for \f$n=1,\dots,\ell\f$
   */
    const static std::array<FReal,ORDER> roots;

    /**
   * Lagrange polynomials \f$ L_n(x) = \Pi_{m=0 \atop m\neq n}^{\ell-1} \frac{x-\bar x_m}{\bar x_n-\bar x_m} \f$
   * Expression with reduced roundoff errors:
   * \f$ L_n(x) = \frac{(-1)^(\ell-n-1)(\ell-1)^(\ell-1)}{(2h)^(\ell-1)n!(\ell-n-1)!} \Pi_{m=0 \atop m\neq n}^{\ell-1} (x-\bar x_m) \f$
   *
   * @param[in] n index
   * @param[in] x coordinate in [-1,1]
   * @return function value
   */
    static FReal L(const unsigned int n, FReal x)
    {
        assert(std::fabs(x)-1.<10.*std::numeric_limits<FReal>::epsilon());
        if (std::fabs(x)>1.) {
            //std::cout << "x=" << x << " out of bounds!" << std::endl;
            x = (x > FReal( 1.) ? FReal( 1.) : x);
            x = (x < FReal(-1.) ? FReal(-1.) : x);
        }

        // Specific precomputation of scale factor
        // in order to minimize round-off errors
        // NB: scale factor could be hardcoded (just as the roots)
        FReal scale;
        int omn = order-n-1;
        if(omn%2) scale=FReal(-1.); // (-1)^(n-1-(k+1)+1)=(-1)^(omn-1)
        else scale=FReal(1.);
        scale/=FReal(FMath::pow(FReal(2.),order-1)*FMath::factorial<FReal>(n)*FMath::factorial<FReal>(omn));

        // compute L
        FReal L=FReal(1.);
        for(unsigned int m=0;m<order;++m){
            if(m!=n){
                // previous version with risks of round-off error
                //L *= (x-FUnifRoots<order>::roots[m])/(FUnifRoots<order>::roots[n]-FUnifRoots<order>::roots[m]);

                // new version (reducing round-off)
                // regular grid on [-1,1] (h simplifies, only the size of the domain and a remains i.e. 2. and -1.)
                L *= (FReal(order-1)*(x+FReal(1.))-FReal(2.)*FReal(m));
            }
        }

        L*=scale;

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
                        tmpNdL*=x-FUnifRoots<FReal, order>::roots[m];
                NdL+=tmpNdL;
                DdL*=FUnifRoots<FReal, order>::roots[n]-FUnifRoots<FReal, order>::roots[p];
            }//endif
        }// p

        return FReal(NdL/DdL);

    }
};

template<int ORDER>
struct FUnifRootsCore{};

template<class FReal, int ORDER>
const std::array<FReal,ORDER> FUnifRoots<FReal,ORDER>::roots = FUnifRootsCore<ORDER>::template Build<FReal>();


// order 2
template<>
struct FUnifRootsCore<2>{
    template <class FReal>
    static std::array<FReal,2> Build(){
        return { { -1., 1.} };
    }
};

// order 3
template<>
struct FUnifRootsCore<3>{
    template <class FReal>
    static std::array<FReal,3> Build(){
        return { {-1.,           0.0,          1.} };
    }
};

// order 4
template<>
struct FUnifRootsCore<4>{
    template <class FReal>
    static std::array<FReal,4> Build(){
        return { {-1.,            -0.333333333333333,            0.333333333333333,            1.} };
    }
};

// order 5
template<>
struct FUnifRootsCore<5>{
    template <class FReal>
    static std::array<FReal,5> Build(){
        return { {-1.,       -0.5,         0.,          0.5,           1.} };
    }
};

// order 6
template<>
struct FUnifRootsCore<6>{
    template <class FReal>
    static std::array<FReal,6> Build(){
        return { {-1.,           -0.6,          -0.2,            0.2,           0.6,            1.} };
    }
};

// order 7
template<>
struct FUnifRootsCore<7>{
    template <class FReal>
    static std::array<FReal,7> Build(){
        return { {-1.,         -0.666666666666666,         -0.333333333333333,          0.,
            0.333333333333333,          0.666666666666666,           1.} };
    }
};

// order 8
template<>
struct FUnifRootsCore<8>{
    template <class FReal>
    static std::array<FReal,8> Build(){
        return { {-1.,
            -0.714285714285714,
            -0.428571428571429,
            -0.142857142857143,
            0.142857142857143,
            0.428571428571429,
            0.714285714285714,
            1.} };
    }
};

// order 9
template<>
struct FUnifRootsCore<9>{
    template <class FReal>
    static std::array<FReal,9> Build(){
        return { {-1.,
            -0.75,
            -0.5,
            -0.25,
            0.0,
            0.25,
            0.5,
            0.75,
            1.} };
    }
};

// order 10
template<>
struct FUnifRootsCore<10>{
    template <class FReal>
    static std::array<FReal,10> Build(){
        return { {-1.,
            -0.777777777777777,
            -0.555555555555555,
            -0.333333333333333,
            -0.111111111111111,
            0.111111111111111,
            0.333333333333333,
            0.555555555555555,
            0.777777777777777,
            1.} };
    }
};

// order 11
template<>
struct FUnifRootsCore<11>{
    template <class FReal>
    static std::array<FReal,11> Build(){
        return { {-1.,
            -0.8,
            -0.6,
            -0.4,
            -0.2,
            0.0,
            0.2,
            0.4,
            0.6,
            0.8,
            1.} };
    }
};

// order 12
template<>
struct FUnifRootsCore<12>{
    template <class FReal>
    static std::array<FReal,12> Build(){
        return { {-1.,
            -0.818181818181818,
            -0.636363636363636,
            -0.454545454545455,
            -0.272727272727273,
            -0.090909090909091,
            0.090909090909091,
            0.272727272727273,
            0.454545454545455,
            0.636363636363636,
            0.818181818181818,
            1.} };
    }
};


// order 13
template<>
struct FUnifRootsCore<13>{
    template <class FReal>
    static std::array<FReal,13> Build(){
        return { {-1.,
            -0.833333333333333,
            -0.666666666666666,
            -0.5,
            -0.333333333333333,
            -0.166666666666666,
            0.0,
            0.166666666666666,
            0.333333333333333,
            0.5,
            0.666666666666666,
            0.833333333333333,
            1.} };
    }
};

// order 14
template<>
struct FUnifRootsCore<14>{
    template <class FReal>
    static std::array<FReal,14> Build(){
        return { {-1.,
            -0.846153846153846,
            -0.692307692307692,
            -0.538461538461538,
            -0.384615384615385,
            -0.230769230769231,
            -0.076923076923077,
            0.076923076923077,
            0.230769230769231,
            0.384615384615385,
            0.538461538461538,
            0.692307692307692,
            0.846153846153846,
            1.} };
    }
};


// order 15
template<>
struct FUnifRootsCore<15>{
    template <class FReal>
    static std::array<FReal,15> Build(){
        return { {-1.0,
            -0.857142857142857,
            -0.714285714285714,
            -0.571428571428571,
            -0.428571428571429,
            -0.285714285714286,
            -0.142857142857143,
            0.0,
            0.142857142857143,
            0.285714285714286,
            0.428571428571429,
            0.571428571428571,
            0.714285714285714,
            0.857142857142857,
            1.0} };
    }
};


// order 20
template<>
struct FUnifRootsCore<16>{
    template <class FReal>
    static std::array<FReal,16> Build(){
        return { {-1.0,
            -0.8947368421052632,
            -0.7894736842105263,
            -0.6842105263157895,
            -0.5789473684210527,
            -0.4736842105263158,
            -0.3684210526315789,
            -0.2631578947368421,
            -0.1578947368421053,
            -0.0526315789473684,
            0.0526315789473684,
            0.1578947368421053,
            0.2631578947368421,
            0.3684210526315789,
            0.4736842105263158,
            0.5789473684210527,
            0.6842105263157895,
            0.7894736842105263,
            0.8947368421052632,
            1.0} };
    }
};

#endif
