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

/** @author Cyrille Piacibello
 * This class provide the same features as FBufferWriter using MPI_Pack system
 *
 * Put some data
 * then insert back if needed
 * finally use data pointer as you like
 */
class FMpiBufferWriter : public FAbstractBufferWriter {
    const MPI_Comm mpiComm;         //< Communicator needed by MPI_Pack functions
    const int arrayCapacity;        //< Allocated Space
    std::unique_ptr<char[]> array;  //< Allocated Array
    int currentIndex;               //< Currently filled space

    /** Test and exit if not enought space */
    void assertRemainingSpace(const size_t requestedSpace) const {
        if(int(currentIndex + requestedSpace) > arrayCapacity){
            printf("Error FMpiBufferWriter has not enough space\n");
            exit(0);
        }
    }

public:
    /** Constructor with a default arrayCapacity of 512 bytes */
    FMpiBufferWriter(const MPI_Comm inComm, const int inCapacity = 1024):
        mpiComm(inComm),
        arrayCapacity(inCapacity),
        array(new char[inCapacity]),
        currentIndex(0)
    {}

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
        assertRemainingSpace(sizeof(ClassType));
        MPI_Pack(const_cast<ClassType*>(&object), 1, FMpi::GetType(object), array.get(), arrayCapacity, &currentIndex, mpiComm);
    }

    /**
   * Allow to pass rvalue to write
   */
    template <class ClassType>
    void write(const ClassType&& object){
        assertRemainingSpace(sizeof(ClassType));
        MPI_Pack(const_cast<ClassType*>(&object), 1, FMpi::GetType(object), array.get(), arrayCapacity, &currentIndex, mpiComm);
    }

    /** Write back, position + sizeof(object) has to be < size */
    template <class ClassType>
    void writeAt(const int position, const ClassType& object){
        if(position + (int) sizeof(ClassType) > currentIndex){
            printf("Not enought space\n");
            exit(0);
        }
        int noConstPosition = position;
        MPI_Pack(const_cast<ClassType*>(&object), 1, FMpi::GetType(object), array.get(), arrayCapacity, &noConstPosition, mpiComm);
    }

    /** Write an array
   * Warning : inSize is a number of ClassType object to write, not a size in bytes
   */
    template <class ClassType>
    void write(const ClassType* const objects, const int inSize){
        assertRemainingSpace(sizeof(ClassType) * inSize);
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
