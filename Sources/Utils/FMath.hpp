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
};


#endif //FMATH_HPP

// [--LICENSE--]
