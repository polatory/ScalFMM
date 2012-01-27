#ifndef FBOOLARRAY_HPP
#define FBOOLARRAY_HPP
// [--License--]

// To get memcpy
#include <cstring>

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FBoolArray
* Please read the license
*
* A bool array is a dynamique allocated array that used 1 bit per value.
*/
class FBoolArray{
    /** Size of a unsigned long */
    const static int SizeOfLong = sizeof(unsigned long);

    /** Size of the array => number of real elements */
    const int size;
    /** The array to store bits */
    unsigned long* const array;
    /** Size of the memory allocated */
    const int memSize;

    /** get size to number of long */
    int LongFromSize(const int inSize){
        const int nbLong = (inSize / (SizeOfLong * 8) );
        return nbLong + 1;
    }

    /** Alloc an array */
    unsigned long * AllocArray(const int inSize){
        return new unsigned long[LongFromSize(inSize)];
    }

public :
    /** Constructor with size */
    FBoolArray(const int inSize) : size(inSize), array(AllocArray(inSize)), memSize(LongFromSize(inSize)*SizeOfLong) {
        setToZeros();
    }

    /** Constructor form another array */
    FBoolArray(const FBoolArray& other): size(other.size), array(AllocArray(other.size)), memSize(other.memSize){
        *this = other;
    }

    /** Destructor */
    ~FBoolArray(){
        delete [] array;
    }

    /**
    * Operator =
    * Array must have the same size
    */
    FBoolArray& operator=(const FBoolArray& other){
        memcpy(array, other.array, memSize);
        return *this;
    }

    /**
    * Operator ==
    * Array must have the same size
    */
    bool operator==(const FBoolArray& other){
        return memcmp(array, other.array, memSize) == 0;
    }

    /**
    * Operator !=
    * Array must have the same size
    */
    bool operator!=(const FBoolArray& other){
        return !(*this == other);
    }

    /** To get a value */
    bool get(const int inPos) const {
        const int posInArray = inPos / (SizeOfLong*8);
        const int bytePosition = inPos - (posInArray * 8);
        return (array[posInArray] >> bytePosition) & 1;
    }

    /** To set a value */
    void set(const int inPos, const bool inVal){
        const int posInArray = inPos / (SizeOfLong*8);
        const int bytePosition = inPos - (posInArray * 8);
        if(inVal) array[posInArray] |= (1UL << bytePosition);
        else array[posInArray] &= ~(1UL << bytePosition);
    }

    /** To get the size of the array */
    int getSize() const {
        return size;
    }

    /** Set all the memory to 0 */
    void setToZeros() const {
        memset( array, 0, memSize);
    }
};


#endif //FBOOLARRAY_HPP

// [--END--]
