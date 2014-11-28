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
#ifndef FMPIBUFFERREADER_HPP
#define FMPIBUFFERREADER_HPP

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
class FMpiBufferReader : public FAbstractBufferReader {
    MPI_Comm comm;            //< Communicator needed by MPI_Pack functions
    int arrayCapacity;        //< Allocated space
    std::unique_ptr<char[]> array;  //< Allocated Array
    int currentIndex;

public :
    /*Constructor with a default arrayCapacity of 512 bytes */
    explicit FMpiBufferReader(const MPI_Comm inComm = MPI_COMM_WORLD, const int inDefaultCapacity = 512):
        comm(inComm),
        arrayCapacity(inDefaultCapacity),
        array(new char[inDefaultCapacity]),
        currentIndex(0){
        FAssertLF(array, "Cannot allocate array");
    }

    /** Change the comm (or to set it later) */
    void setComm(const MPI_Comm inComm){
        comm = inComm;
    }

    /** To change the capacity (but reset the head to 0) */
    void cleanAndResize(const int newCapacity){
        if(newCapacity != arrayCapacity){
            arrayCapacity = newCapacity;
            array.reset(new char[newCapacity]);
        }
        currentIndex = 0;
    }

    /** Destructor
   */
    virtual ~FMpiBufferReader(){
    }

    /** Get allocated memory pointer */
    char* data(){
        return array.get();
    }

    /** Get allocated memory pointer */
    const char* data() const {
        return array.get();
    }

    /** get the filled space */
    int getSize() const{
        return currentIndex;
    }

    /** Size of the memory initialized */
    int getCapacity() const{
        return arrayCapacity;
    }

    /** Move the read index to a position */
    void seek(const int inIndex){
      
      FAssertLF(inIndex <= arrayCapacity, "FMpiBufferReader :: Aborting :: Can't move index because buffer isn't long enough ",inIndex," ",arrayCapacity);
        currentIndex = inIndex;
    }

    /** Get the read position */
    int tell() const {
        return currentIndex;
    }

    /** Get a value with memory cast */
    template <class ClassType>
    ClassType getValue(){
        ClassType value;
        int previousIndex = currentIndex;
        seek(int(sizeof(value) + previousIndex));
        MPI_Unpack(array.get(),arrayCapacity,&previousIndex,&value,1,FMpi::GetType(value),comm);
        return value;
    }

    /** Get a value with memory cast at a specified index */
    template <class ClassType>
    ClassType getValue(const int ind){
        ClassType value;
        int previousIndex = ind;
        seek(int(sizeof(value)+ind));
        MPI_Unpack(array.get(),arrayCapacity,&previousIndex,&value,1,FMpi::GetType(value),comm);
        return value;
    }

    /** Fill a value with memory cast */
    template <class ClassType>
    void fillValue(ClassType* const inValue){
        int previousIndex = currentIndex;
        seek(int(sizeof(ClassType) + previousIndex));
        MPI_Unpack(array.get(),arrayCapacity,&previousIndex,inValue,1,FMpi::GetType(*inValue),comm);
    }

    /** Fill one/many value(s) with memcpy */
    template <class ClassType>
    void fillArray(ClassType* const inArray, const int inSize){
        int previousIndex = currentIndex;
        seek(int(sizeof(ClassType) * inSize + previousIndex));
        MPI_Unpack(array.get(),arrayCapacity,&previousIndex,inArray,inSize,FMpi::GetType(*inArray),comm);
    }

    /** Same as fillValue */
    template <class ClassType>
    FMpiBufferReader& operator>>(ClassType& object){
        fillValue(&object);
        return *this;
    }

};
#endif

