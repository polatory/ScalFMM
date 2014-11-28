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
#ifndef FMPIBUFFERWRITER_HPP
#define FMPIBUFFERWRITER_HPP

#include <memory>
#include "../Utils/FMpi.hpp"
#include "FAbstractBuffer.hpp"
#include "../Utils/FAssert.hpp"

/** @author Cyrille Piacibello
 * This class provide the same features as FBufferWriter using MPI_Pack system
 *
 * Put some data
 * then insert back if needed
 * finally use data pointer as you like
 */
class FMpiBufferWriter : public FAbstractBufferWriter {
    MPI_Comm mpiComm;         //< Communicator needed by MPI_Pack functions
    int arrayCapacity;              //< Allocated Space
    std::unique_ptr<char[]> array;  //< Allocated Array
    int currentIndex;               //< Currently filled space

    /** Test and exit if not enought space */
    void expandIfNeeded(const size_t requestedSpace) {
        if( arrayCapacity < int(currentIndex + requestedSpace) ){
            arrayCapacity = int((currentIndex + requestedSpace + 1) * 1.5);
            char* arrayTmp = new char[arrayCapacity];
            memcpy(arrayTmp, array.get(), sizeof(char)*currentIndex);
            array.reset(arrayTmp);
        }
    }

public:
    /** Constructor with a default arrayCapacity of 512 bytes */
    explicit FMpiBufferWriter(const MPI_Comm inComm, const int inDefaultCapacity = 1024):
        mpiComm(inComm),
        arrayCapacity(inDefaultCapacity),
        array(new char[inDefaultCapacity]),
        currentIndex(0)
    {}


    /** Change the comm (or to set it later) */
    void setComm(const MPI_Comm inComm){
        mpiComm = inComm;
    }

    /** To change the capacity (but reset the head to 0 if size if lower) */
    void resize(const int newCapacity){
        if(newCapacity != arrayCapacity){
            arrayCapacity = newCapacity;
            char* arrayTmp = new char[arrayCapacity];
            currentIndex = (currentIndex < arrayCapacity ? currentIndex : arrayCapacity-1);
            memcpy(arrayTmp, array.get(), sizeof(char)*currentIndex);
            array.reset(arrayTmp);
        }
    }

    /** Destructor */
    virtual ~FMpiBufferWriter(){
    }

    /** Get allocated memory pointer */
    char* data(){
        return array.get();
    }

    /** Get allocated memory pointer */
    const char* data() const {
        return array.get();
    }

    /** Get the filled space */
    int getSize() const {
        return currentIndex;
    }

    /** Get the allocated space */
    int getCapacity() const {
        return arrayCapacity;
    }

    /** Write data by packing cpy */
    template <class ClassType>
    void write(const ClassType& object){
        expandIfNeeded(sizeof(ClassType));
        MPI_Pack(const_cast<ClassType*>(&object), 1, FMpi::GetType(object), array.get(), arrayCapacity, &currentIndex, mpiComm);
    }

    /**
   * Allow to pass rvalue to write
   */
    template <class ClassType>
    void write(const ClassType&& object){
        expandIfNeeded(sizeof(ClassType));
        MPI_Pack(const_cast<ClassType*>(&object), 1, FMpi::GetType(object), array.get(), arrayCapacity, &currentIndex, mpiComm);
    }

    /** Write back, position + sizeof(object) has to be < size */
    template <class ClassType>
    void writeAt(const int position, const ClassType& object){
        FAssertLF(int(position + sizeof(ClassType)) <= currentIndex)
        int noConstPosition = position;
        MPI_Pack(const_cast<ClassType*>(&object), 1, FMpi::GetType(object), array.get(), arrayCapacity, &noConstPosition, mpiComm);
    }

    /** Write an array
   * Warning : inSize is a number of ClassType object to write, not a size in bytes
   */
    template <class ClassType>
    void write(const ClassType* const objects, const int inSize){
        expandIfNeeded(sizeof(ClassType) * inSize);
        MPI_Pack( const_cast<ClassType*>(objects), inSize, FMpi::GetType(*objects), array.get(), arrayCapacity, &currentIndex, mpiComm);
    }

    /** Equivalent to write */
    template <class ClassType>
    FMpiBufferWriter& operator<<(const ClassType& object){
        write(object);
        return *this;
    }

    /** Reset the writing index, but do not change the arrayCapacity */
    void reset(){
        currentIndex = 0;
    }
};


#endif // FBUFFERWRITER_HPP
