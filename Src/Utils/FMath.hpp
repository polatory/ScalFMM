// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
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
        //const FReal relTol = FReal(0.00001);
        //const FReal absTol = FReal(0.00001);
        //return (Abs(inV1 - inV2) <= Max(absTol, relTol * Max(Abs(inV1), Abs(inV2))));
        return Abs(inV1 - inV2) < std::numeric_limits<NumType>::epsilon();
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
    template <class NumType>
    static NumType pow(const NumType inValue, long power){
        NumType result = 1;
        while(power-- > 0) result *= inValue;
        return result;
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

    /** To get atan2 of a 2 FReal */
    static float Atan2(const float inValue1,const float inValue2){
        return atan2f(inValue1,inValue2);
    }
    static double Atan2(const double inValue1,const double inValue2){
        return atan2(inValue1,inValue2);
    }

    /** To get sqrt of a FReal */
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

    /** To get acos of a FReal */
    static float ACos(const float inValue){
        return acosf(inValue);
    }
    static double ACos(const double inValue){
        return acos(inValue);
    }

    /** To get atan2 of a 2 FReal */
    static float Fmod(const float inValue1,const float inValue2){
        return fmodf(inValue1,inValue2);
    }
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

    /** Compute a relative difference between two values */
    template <class ValueClass>
    static ValueClass RelativeDiff(const ValueClass& value1, const ValueClass& value2){
        if(Abs(value1) > Abs(value2)){
            return Abs((value2 - value1) / value1);
        }
        else{
            return Abs((value2 - value1) / value2);
        }
    }
};


#endif //FMATH_HPP
