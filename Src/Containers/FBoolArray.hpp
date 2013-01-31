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
#ifndef FBOOLARRAY_HPP
#define FBOOLARRAY_HPP


// To get memcpy
#include <cstring>

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FBoolArray
* Please read the license
*
* A bool array is a dynamique allocated array that used 1 bit per value.
* Under the wood, it use bit operations to acess/set value in an array of
* native type.
*/
class FBoolArray{
    /** Size of a unsigned long */
    const static int BytesInBlock = sizeof(unsigned long);
    const static int SizeOfBlock = BytesInBlock * 8;

    /** The array to store bits */
    unsigned long* const array;
    /** Size of the memory allocated */
    const int memSize;
    /** Size of the array => number of real elements */
    const int size;

    /** get size to number of long */
    int LongFromSize(const int inSize){
        return ((inSize + SizeOfBlock - 1) / SizeOfBlock);
    }

    /** Alloc an array */
    unsigned long * AllocArray(const int inSize){
        return new unsigned long[LongFromSize(inSize)];
    }

public :
    /** Constructor with size */
    FBoolArray(const int inSize) : array(AllocArray(inSize)), memSize(LongFromSize(inSize)*BytesInBlock), size(inSize) {
        setToZeros();
    }

    /** Constructor form another array */
    FBoolArray(const FBoolArray& other): array(AllocArray(other.size)), memSize(other.memSize), size(other.size){
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
        const int posInArray = inPos / SizeOfBlock;
        const int bytePosition = inPos - (posInArray * 8);
        return (array[posInArray] >> bytePosition) & 1;
    }

    /** To set a value */
    void set(const int inPos, const bool inVal){
        const int posInArray = inPos / SizeOfBlock;
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


