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

#include "../Utils/FMpi.hpp"


/** @author Cyrille Piacibello
 * This class provide the same features as FBufferWriter using MPI_Pack system
 *
 * Put some data
 * then insert back if needed
 * finally use data pointer as you like
 */
class FMpiBufferReader {

private:
  //Classe FMpiVector, used only by FMpiBuffer
  template<class ObjectType>
  class FMpiVector {
      protected :
    static const int DefaultSize = 10;
    ObjectType* array;                    //memory area
    int capacity;                         //Capacity of the array
    long int index;                       //current position in Byte !!
  
  public:
    FMpiVector():
      array(0),
      capacity(DefaultSize),
      index(0){
      array = reinterpret_cast<ObjectType*>(new char[sizeof(ObjectType)* DefaultSize]);
    }

    FMpiVector(const int inCapa):
      array(0),
      capacity(inCapa),
      index(0){
      array = reinterpret_cast<ObjectType*>(new char[sizeof(ObjectType)* inCapa]);
    }

    virtual ~FMpiVector(){
      delete[] reinterpret_cast< char* >(array);
    }

    //To get the capacity
    const int getCapacity() const{
      return this->capacity;
    }
    
    //To get the array
    ObjectType * data(){
      return array;
    }
    
    const ObjectType * data() const{
      return array;
    }
    
    //To delete all the element stored of the array
    void clear(){
      while(0 < index){
	(&array[--index])->~ObjectType();
      }
    }
    
    //To get how much space is used
    long int getSize() const{
      return this->index;
    }

    //To get how many objects are stored
    int getObjectsSize(){
      return (this->index / sizeof(ObjectType));
    }
    
    //To inc the index
    //Usually used with array.incIndex(sizeof(my_object_stored));
    void incIndex(const int inInc){
      if(index + inInc > capacity){
	fprintf(stderr,"Aborting : index array out of range\n");
	exit(0);
      }
      else{
	this->index+=inInc;
      }
    }
    
    //To set the index
    void setIndex(const int inInd){
      if(inInd>capacity){
	fprintf(stderr,"Aborting : index array out of range\n");
	exit(0);
      }
      else{
	this->index = inInd;
      }
    }
  };

  MPI_Comm comm;
  FMpiVector<char> array;

public :
  FMpiBufferReader(MPI_Comm inComm, const int inCapacity = 0):
    comm(inComm),
    array(inCapacity)
  {}
  
  /** Destructor
   */
  virtual ~FMpiBufferReader(){
  }
  
  /** Get the memory area */
  char* data(){
    return array.data();
  }
  
  /** Get the memory area */
  const char* data() const {
    return array.data();
  }
  
  /** Size of the memory initialized */
  int getSize() const{
    return array.getSize();
  }

  /** Move the read index to a position */
  void seek(const int inIndex){
    array.setIndex(inIndex);
  }

  /** Get the read position */
  int tell() const {
    return array.getSize();
  }
  
  /** Get a value with memory cast */
  template <class ClassType>
  ClassType getValue(){
    ClassType value;
    int currentIndex = array.getSize();
    array.incIndex(sizeof(value));
    MPI_Unpack(array.data(),sizeof(ClassType),&currentIndex,&value,1,FMpi::GetType(value),comm);
    return value;
  }

  /** Fill a value with memory cast */
  template <class ClassType>
  void fillValue(ClassType* const inValue){
    int currentIndex = array.getSize();
    array.incIndex(sizeof(ClassType));
    MPI_Pack(inValue,1,FMpi::GetType(*inValue),array.data(),array.getCapacity(),&currentIndex,comm);
  }

  /** Fill one/many value(s) with memcpy */
  template <class ClassType>
  void fillArray(ClassType* const inArray, const int inSize){
    int currentIndex = array.getSize();
    array.incIndex(sizeof(ClassType) * inSize);
    MPI_Pack(inArray,inSize,FMpi::GetType(*inArray),array.data(),array.getCapacity(),&currentIndex,comm);
  }

  /** Same as fillValue */
  template <class ClassType>
  FMpiBufferReader& operator>>(ClassType& object){
    fillValue(&object);
    return *this;
  }

};
#endif

