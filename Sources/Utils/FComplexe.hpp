#ifndef FCOMPLEXE_HPP
#define FCOMPLEXE_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "FMath.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class 
* Please read the license
*
* Propose basic complexe class.
*/
class FComplexe {
    FReal real;    //< Real
    FReal imag;    //< Imaginary

public:
    /** Default Constructor (set real&imaginary to 0) */
    FComplexe() : real(0),imag(0){
    }

    /** Constructor with values
      * @param inImag the imaginary
      * @param inReal the real
      */
    FComplexe(const FReal inImag, const FReal inReal)
        : real(inReal),imag(inImag){
    }

    /** Copy constructor */
    FComplexe(const FComplexe& other)
        : real(other.real), imag(other.imag){
    }

    /** Copy operator */
    FComplexe& operator=(const FComplexe& other){
        this->imag = other.imag;
        this->real = other.real;
        return *this;
    }

    /** Equality operator */
    bool operator==(const FComplexe& other){
        return FMath::LookEqual(this->imag,other.imag)
                       && FMath::LookEqual(this->real,other.real);
    }

    /** Different equal */
    bool operator!=(const FComplexe& other){
        return !(*this == other);
    }

    /** Get imaginary */
    FReal getImag() const{
        return this->imag;
    }

    /** Get real */
    FReal getReal() const{
        return this->real;
    }

    /** Set Imaginary */
    void setImag(const FReal inImag) {
        this->imag = inImag;
    }

    /** Set Real */
    void setReal(const FReal inReal) {
        this->real = inReal;
    }

    /**
     * Operator +=
     * in real with other real, same for imag
     * @param other the complexe to use data
     */
    FComplexe& operator+=(const FComplexe& other){
        this->real += other.real;
        this->imag += other.imag;
        return *this;
    }   

    /** Inc real and imaginary by values
      * @param inIncReal to inc the real
      * @param inIncImag to inc the imag
      */
    void inc(const FReal inIncReal, const FReal inIncImag){
        this->real += inIncReal;
        this->imag += inIncImag;
    }

    /** Inc real by FReal
      * @param inIncReal to inc the real
      */
    void incReal(const FReal inIncReal){
        this->real += inIncReal;
    }

    /** Inc imaginary by FReal
      * @param inIncImag to inc the imag
      */
    void incImag(const FReal inIncImag){
        this->imag += inIncImag;
    }

    /** Dec real by FReal
      * @param inDecReal to dec the real
      */
    void decReal(const FReal inIncReal){
        this->real -= inIncReal;
    }

    /** Dec imaginary by FReal
      * @param inDecImag to dec the imag
      */
    void decImag(const FReal inIncImag){
        this->imag -= inIncImag;
    }

    /** Mul real and imaginary by a FReal
      * @param inValue the coef to mul data
      */
    void mulRealAndImag(const FReal inValue){
        this->imag *= inValue;
        this->real *= inValue;
    }

    /** Mul a complexe by another "c*=c2" */
    FComplexe& operator*=(const FComplexe& other){
        const FReal tempReal = this->real;
        this->real = (tempReal * other.real) - (this->imag * other.imag);
        this->imag = (tempReal * other.imag) + (this->imag * other.real);
        return *this;
    }
};

/** Global operator Mul a complexe by another "c=c1*c2" */
FComplexe operator*=(const FComplexe& first, const FComplexe& second){
    const FComplexe result(
            (first.getReal() * second.getImag()) + (first.getImag() * second.getReal()),
            (first.getReal() * second.getReal()) - (first.getImag() * second.getImag())
            );
    return result;
}

#endif //FCOMPLEXE_HPP

// [--LICENSE--]
