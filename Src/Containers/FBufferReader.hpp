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
#ifndef FBUFFERREADER_HPP
#define FBUFFERREADER_HPP

#include "FVector.hpp"


/** @author Berenger Bramas
  * This class provide a fast way to manage a memory and convert
  * the content to basic type.
  *
  * Specifie the needed space with reserve,
  * then fill it with data
  * finaly read and convert.
  */
class FBufferReader {
    FVector<char> buffer;   //< The memory buffer
    int index;              //< The current index reading position

public:
    /** Construct with a memory size init to 0 */
    explicit FBufferReader(const int inCapacity = 0) : buffer(inCapacity), index(0) {
        if(inCapacity){
            reserve(inCapacity);
        }
    }

    /** Destructor */
    virtual ~FBufferReader(){
    }

    /** Get the memory area */
    char* data(){
        return buffer.data();
    }

    /** Get the memory area */
    const char* data() const {
        return buffer.data();
    }

    /** Size of the meomry initialzed */
    int getSize() const{
        return buffer.getSize();
    }

    /** Move the read index to a position */
    void seek(const int inIndex){
        index = inIndex;
    }

    /** Get the read position */
    int tell() const {
        return index;
    }

    /** Reset and allocate nbBytes memory filled with 0 */
    void reserve(const int nbBytes){
        buffer.clear();
        buffer.set( 0, nbBytes);
    }

    /** Move the read index to 0 */
    void reset(){
        buffer.clear();
        index = 0;
    }

    /** Get a value with memory cast */
    template <class ClassType>
    ClassType getValue(){
        ClassType value = (*reinterpret_cast<ClassType*>(&buffer[index]));
        index += sizeof(ClassType);
        return value;
    }

    /** Fill a value with memory cast */
    template <class ClassType>
    void fillValue(ClassType* const inValue){
        (*inValue) = (*reinterpret_cast<ClassType*>(&buffer[index]));
        index += sizeof(ClassType);
    }

    /** Fill one/many value(s) with memcpy */
    template <class ClassType>
    void fillArray(ClassType* const inArray, const int inSize){
        memcpy( inArray, &buffer[index], sizeof(ClassType) * inSize);
        index += sizeof(ClassType) * inSize;
    }

    /** Same as fillValue */
    template <class ClassType>
    FBufferReader& operator>>(ClassType& object){
        fillValue(&object);
        return *this;
    }
};


#endif // FBUFFERREADER_HPP
