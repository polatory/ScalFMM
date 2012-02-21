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
#ifndef FBUFFERWRITER_HPP
#define FBUFFERWRITER_HPP

#include "FVector.hpp"

/** @author Berenger Bramas
  * This class provide a fast way to manage a memory and fill it
  *
  * Put some data
  * then insert back if needed
  * finaly use data pointer as you like
  */
class FBufferWriter {
private:
    FVector<char> buffer; //< The buffer

public:
    /** Constructor with a default capacity of 512 bytes */
    explicit FBufferWriter(const int inCapacity = 512) : buffer(inCapacity) {
    }

    /** Destructor */
    virtual ~FBufferWriter(){
    }

    /** Get allocated memory pointer */
    char* data(){
        return buffer.data();
    }

    /** Get allocated memory pointer */
    const char* data() const {
        return buffer.data();
    }

    /** Get the filled space */
    int getSize() const {
        return buffer.getSize();
    }

    /** Write data by mem cpy */
    template <class ClassType>
    void write(const ClassType& object){
        buffer.memocopy(reinterpret_cast<const char*>(&object), sizeof(ClassType));
    }

    /** Write back, position + sizeof(object) has to be < size */
    template <class ClassType>
    void writeAt(const int position, const ClassType& object){
        (*reinterpret_cast<ClassType*>(&buffer[position])) = object;
    }

    /** Write an array */
    template <class ClassType>
    void write(const ClassType* const objects, const int inSize){
        buffer.memocopy(reinterpret_cast<const char*>(objects), int(sizeof(ClassType)) * inSize);
    }

    /** Equivalent to write */
    template <class ClassType>
    FBufferWriter& operator<<(const ClassType& object){
        write(object);
        return *this;
    }

    /** Reset the writing index, but do not change the capacity */
    void reset(){
        buffer.clear();
    }
};


#endif // FBUFFERWRITER_HPP