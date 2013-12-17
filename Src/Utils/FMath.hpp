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
#ifndef FMATH_HPP
#define FMATH_HPP


#include <cmath>
#include <limits>

#include "FGlobal.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class
* Please read the license
*
* Propose basic math functions or indirections to std math.
*/
struct FMath{
    static const FReal FPi;        //< Pi constant
    static const FReal FTwoPi;     //< 2 Pi constant
    static const FReal FPiDiv2;    //< Pi/2 constant
    static const FReal Epsilon;    //< Epsilon

    /** To get absolute value */
    template <class NumType>
    static NumType Abs(const NumType inV){
        return (inV < 0 ? -inV : inV);
    }

    /** To get max between 2 values */
    template <class NumType>
    static NumType Max(const NumType inV1, const NumType inV2){
        return (inV1 > inV2 ? inV1 : inV2);
    }

    /** To get min between 2 values */
    template <class NumType>
    static NumType Min(const NumType inV1, const NumType inV2){
        return (inV1 < inV2 ? inV1 : inV2);
    }

    /** To know if 2 values seems to be equal */
    template <class NumType>
    static bool LookEqual(const NumType inV1, const NumType inV2){
        return (Abs(inV1-inV2) < std::numeric_limits<NumType>::epsilon());
        //const FReal relTol = FReal(0.00001);
        //const FReal absTol = FReal(0.00001);
        //return (Abs(inV1 - inV2) <= Max(absTol, relTol * Max(Abs(inV1), Abs(inV2))));
    }

    /** To know if 2 values seems to be equal */
    template <class NumType>
    static FReal RelatifDiff(const NumType inV1, const NumType inV2){
        return Abs(inV1 - inV2)*Abs(inV1 - inV2)/Max(Abs(inV1*inV1), Abs(inV2*inV2));
    }

    /** To get floor of a FReal */
    static float dfloor(const float inValue){
        return floorf(inValue);
    }
    static double dfloor(const double inValue){
        return floor(inValue);
    }

    /** To get ceil of a FReal */
    static float Ceil(const float inValue){
        return ceilf(inValue);
    }
    static double Ceil(const double inValue){
        return ceil(inValue);
    }

    /** To get pow */
    static double pow(double x, double y){
        return ::pow(x,y);
    }
    static double pow(float x, float y){
        return ::powf(x,y);
    }
    template <class NumType>
    static NumType pow(const NumType inValue, int power){
        if(power<0) power = -power;
        NumType result = 1;
        while(power-- > 0) result *= inValue;
        return result;
    }

    /** To get pow of 2 */
    static int pow2(const int power){
        return (1 << power);
    }

    /** To know if a value is between two others */
    template <class NumType>
    static bool Between(const NumType inValue, const NumType inMin, const NumType inMax){
        return ( inMin <= inValue && inValue < inMax );
    }

    /** To get sqrt of a FReal */
    static float Sqrt(const float inValue){
        return sqrtf(inValue);
    }
    static double Sqrt(const double inValue){
        return sqrt(inValue);
    }

    /** To get Log of a FReal */
    static float Log(const float inValue){
        return logf(inValue);
    }
    static double Log(const double inValue){
        return log(inValue);
    }

    /** To get Log2 of a FReal */
    static float Log2(const float inValue){
        return log2f(inValue);
    }
    static double Log2(const double inValue){
        return log2(inValue);
    }

    /** To get atan2 of a 2 FReal,  The return value is given in radians and is in the
 range -pi to pi, inclusive.  */
    static float Atan2(const float inValue1,const float inValue2){
        return atan2f(inValue1,inValue2);
    }
    static double Atan2(const double inValue1,const double inValue2){
        return atan2(inValue1,inValue2);
    }

    /** To get sin of a FReal */
    static float Sin(const float inValue){
        return sinf(inValue);
    }
    static double Sin(const double inValue){
        return sin(inValue);
    }

    /** To get cos of a FReal */
    static float Cos(const float inValue){
        return cosf(inValue);
    }
    static double Cos(const double inValue){
        return cos(inValue);
    }

    /** To get arccos of a float. The result is in the range [0, pi]*/
    static float ACos(const float inValue){
        return acosf(inValue);
    }
    /** To get arccos of a double. The result is in the range [0, pi]*/
    static double ACos(const double inValue){
        return acos(inValue);
    }

    /** To get atan2 of a 2 FReal */
    static float Fmod(const float inValue1,const float inValue2){
        return fmodf(inValue1,inValue2);
    }
    /** return the floating-point remainder of inValue1  / inValue2 */
    static double Fmod(const double inValue1,const double inValue2){
        return fmod(inValue1,inValue2);
    }

    /** To know if a variable is nan, based on the C++0x */
    template <class TestClass>
    static bool IsNan(const TestClass& value){
        //volatile const TestClass* const pvalue = &value;
        //return (*pvalue) != value;
        return std::isnan(value);
    }

    /** To know if a variable is not inf, based on the C++0x */
    template <class TestClass>
    static bool IsFinite(const TestClass& value){
        // return !(value <= std::numeric_limits<T>::min()) && !(std::numeric_limits<T>::max() <= value);
        return std::isfinite(value);
    }


    /** A class to compute accuracy */
    class FAccurater {
        FReal l2Dot;
        FReal l2Diff;
        FReal max;
        FReal maxDiff;
    public:
        FAccurater() : l2Dot(0), l2Diff(0), max(0), maxDiff(0) {
        }
        /** with inital values */
        FAccurater(const FReal inGood[], const FReal inBad[], const int nbValues)
            : l2Dot(0), l2Diff(0), max(0), maxDiff(0) {
            add(inGood, inBad, nbValues);
        }
        /** Add value to the current list */
        void add(const FReal inGood, const FReal inBad){
            l2Diff += (inBad - inGood) * (inBad - inGood);
            l2Dot  += inGood * inGood;

            max = Max(max , Abs(inGood));
            maxDiff = Max(maxDiff, Abs(inGood-inBad));
        }
        /** Add array of values */
        void add(const FReal inGood[], const FReal inBad[], const int nbValues){
            for(int idx = 0 ; idx < nbValues ; ++idx){
                add(inGood[idx],inBad[idx]);
            }
        }
        /** Get the L2 norm */
        FReal getL2Norm() const{
            return Sqrt(l2Diff / l2Dot);
        }
        /** Get the inf norm */
        FReal getInfNorm() const{
            return maxDiff / max;
        }
        /** Print */
        template <class StreamClass>
        friend StreamClass& operator<<(StreamClass& output, const FAccurater& inAccurater){
            output << "[Error] L2Norm = " << inAccurater.getL2Norm() << " \t Inf = " << inAccurater.getInfNorm();
            return output;
        }
    };
};


#endif //FMATH_HPP
