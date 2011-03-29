#ifndef FMATH_HPP
#define FMATH_HPP
// /!\ Please, you must read the license at the bottom of this page

#include <math.h>

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class 
* Please read the license
*
* Propose basic math functions or indirections
*/
struct FMath{
    static const double FPi;        //< Pi constant
    static const double FPiDiv2;    //< Pi/2 constant

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
        /*const double relTol = 0.00001;
		const double absTol = 0.00001;
                return (Abs(inV1 - inV2) <= Max(absTol, relTol * Max(Abs(inV1), Abs(inV2))));*/
        return Abs(inV1 - inV2) <= (Abs(inV1 < Abs(inV2) ? Abs(inV2) : Abs(inV1)) * 0.00001);
    }

    /** To get floor of a double */
    static double dfloor(const double inValue){
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

    /** To get sqrt of a double */
    static double Sqrt(const double inValue){
        return sqrt(inValue);
    }

    /** To get atan2 of a double */
    static double Atan2(const double inValue1,const double inValue2){
        return atan2(inValue1,inValue2);
    }

    /** To get sqrt of a double */
    static double Sin(const double inValue){
        return sin(inValue);
    }

    /** To get cos of a double */
    static double Cos(const double inValue){
        return cos(inValue);
    }

    /** To get acos of a double */
    static double ACos(const double inValue){
        return acos(inValue);
    }
};

const double FMath::FPi = M_PI;
const double FMath::FPiDiv2 = M_PI_2;

#endif //FMATH_HPP

// [--LICENSE--]
