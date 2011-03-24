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
    double imag;    //< Imaginary
    double real;    //< Real

public:
    /** Default Constructor */
    FComplexe() : imag(0), real(0){
    }

    /** Constructor with default value
      * @param inImag the imaginary
      * @param inReal the real
      */
    FComplexe(const double inImag, const double inReal)
        : imag(inImag), real(inReal){
    }

    /** Copy constructor */
    FComplexe(const FComplexe& other)
        : imag(other.imag), real(other.real){
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
    double getImag() const{
        return this->imag;
    }

    /** Get real */
    double getReal() const{
        return this->real;
    }

    /** Set Imaginary */
    void setImag(const double inImag) {
        this->imag = inImag;
    }

    /** Set Real */
    void setReal(const double inReal) {
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

    /** Mul real and imaginary by a double
      * @param inValue the coef to mul data
      */
    void mulRealAndImag(const double inValue){
        this->imag *= inValue;
        this->real *= inValue;
    }

    /** Inc real and imaginary by doubles
      * @param inIncReal to inc the real
      * @param inIncImag to inc the imag
      */
    void inc(const double inIncReal, const double inIncImag){
        this->real += inIncReal;
        this->imag += inIncImag;
    }

    /** Inc real by double
      * @param inIncReal to inc the real
      */
    void incReal(const double inIncReal){
        this->real += inIncReal;
    }

    /** Inc imaginary by double
      * @param inIncImag to inc the imag
      */
    void incImag(const double inIncImag){
        this->imag += inIncImag;
    }

    /** Dec real by double
      * @param inDecReal to dec the real
      */
    void decReal(const double inIncReal){
        this->real -= inIncReal;
    }

    /** Dec imaginary by double
      * @param inDecImag to dec the imag
      */
    void decImag(const double inIncImag){
        this->imag -= inIncImag;
    }
};


#endif //FCOMPLEXE_HPP

// [--LICENSE--]
