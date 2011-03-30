#ifndef FMATH_HPP
#define FMATH_HPP
// /!\ Please, you must read the license at the bottom of this page

#include <math.h>

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
    static const FReal Epsilon;     //< Epsilon

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
        /*const FReal relTol = 0.00001;
		const FReal absTol = 0.00001;
                return (Abs(inV1 - inV2) <= Max(absTol, relTol * Max(Abs(inV1), Abs(inV2))));*/
        return Abs(inV1 - inV2) <= (Abs(inV1 < Abs(inV2) ? Abs(inV2) : Abs(inV1)) * 0.00001);
    }

    /** To get floor of a FReal */
    static FReal dfloor(const FReal inValue){
        return floor(inValue);
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
    static FReal Sqrt(const FReal inValue){
        return sqrt(inValue);
    }

    /** To get atan2 of a 2 FReal */
    static FReal Atan2(const FReal inValue1,const FReal inValue2){
        return atan2(inValue1,inValue2);
    }

    /** To get sqrt of a FReal */
    static FReal Sin(const FReal inValue){
        return sin(inValue);
    }

    /** To get cos of a FReal */
    static FReal Cos(const FReal inValue){
        return cos(inValue);
    }

    /** To get acos of a FReal */
    static FReal ACos(const FReal inValue){
        return acos(inValue);
    }

    /** To get atan2 of a 2 FReal */
    static FReal Fmod(const FReal inValue1,const FReal inValue2){
        return fmod(inValue1,inValue2);
    }
};

const FReal FMath::FPi = M_PI;
const FReal FMath::FPiDiv2 = M_PI_2;
const FReal FMath::Epsilon = 0.00000000000000000001;

#endif //FMATH_HPP

// [--LICENSE--]
