#ifndef FCOMPLEXE_HPP
#define FCOMPLEXE_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "FMath.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class 
* Please read the license
*
* Propose basic complexe class
*/
class FComplexe {
    double imag;
    double real;
public:
    FComplexe()
        : imag(0), real(0){
    }
    FComplexe(const double inImag, const double inReal)
        : imag(inImag), real(inReal){
    }
    FComplexe(const FComplexe& other)
        : imag(other.imag), real(other.real){
    }
    FComplexe& operator=(const FComplexe& other){
        this->imag = other.imag;
        this->real = other.real;
        return *this;
    }
    bool operator==(const FComplexe& other){
        return FMath::LookEqual(this->imag,other.imag)
                       && FMath::LookEqual(this->real,other.real);
    }
    bool operator!=(const FComplexe& other){
        return !(*this == other);
    }

    double getImag() const{
        return this->imag;
    }

    double getReal() const{
        return this->real;
    }

    void setImag(const double inImag) {
        this->imag = inImag;
    }

    void setReal(const double inReal) {
        this->real = inReal;
    }
};


#endif //FCOMPLEXE_HPP

// [--LICENSE--]
