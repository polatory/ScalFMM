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


#include "../Utils/FMpi.hpp"

/** @author Cyrille Piacibello
 * This class provide the same features as FBufferWriter using MPI_Pack system
 *
 * Put some data
 * then insert back if needed
 * finally use data pointer as you like
 */
class FMpiBufferWriter {
  
private :
  //Classe FMpiVector, used only by FMpiBuffer.
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

  const MPI_Comm comm;      // Communicator needed by MPI_Pack functions
  FMpiVector<char> array;

public:
  /** Constructor with a default capacity of 512 bytes */
  FMpiBufferWriter(const MPI_Comm inComm, const int inCapacity = 512): 
    comm(inComm),
    array(inCapacity)
  {}
  
  /** Destructor */
  virtual ~FMpiBufferWriter(){
    delete &array;
  }
  
  /** Get allocated memory pointer */
  char* data(){
    return array.data();
  }

  /** Get allocated memory pointer */
  const char* data() const {
    return array.data();
  }

  /** Get the filled space */
  int getSize() const {
    return array.getSize();
  }

  /** Write data by packing cpy */
  template <class ClassType>
  void write(const ClassType& object){
    //        buffer.memocopy(reinterpret_cast<const char*>(&object), int(sizeof(ClassType)));
    int currentIndex = array.getSize();
    array.incIndex(sizeof(ClassType));
    MPI_Pack(const_cast<ClassType*>(&object),1,FMpi::GetType(object),array.data(),array.getCapacity(),&currentIndex,comm);
  }

  /**
   * Allow to pass rvalue to write
   */
  template <class ClassType>
  void write(const ClassType&& object){
    //        buffer.memocopy(reinterpret_cast<const char*>(&object), int(sizeof(ClassType)));
    int currentIndex = array.getSize();
    array.incIndex(sizeof(ClassType));
    MPI_Pack(const_cast<ClassType*>(&object),1,FMpi::GetType(object),array.data(),array.getCapacity(),&currentIndex,comm);
  }
  
  /** Write back, position + sizeof(object) has to be < size */
  template <class ClassType>
  void writeAt(const int position, const ClassType& object){
    //(*reinterpret_cast<ClassType*>(&buffer[position])) = object;
    if(position < array.getSize()){
      fprintf(stderr,"Aborting : writeAt is overwritting data\n");
    }
    else{
      int temp = position;
      if(position + FMpi::GetType(object) < array.getCapacity()){
	MPI_Pack(&object,1,FMpi::GetType(object),array.data(),array.getCapacity(),&temp,comm);
      }
      array.setIndex(temp);
    }
  }
  
  /** Write an array 
   * Warning : inSize is a number of ClassType object to write, not a size in bytes
   */
  template <class ClassType>
  void write(const ClassType* const objects, const int inSize){
    int currentIndex = array.getSize();
    array.incIndex(sizeof(ClassType) * inSize);
    MPI_Pack(objects,inSize,FMpi::GetType(*objects),array.data(),array.getCapacity(),&currentIndex,comm);
  }
  
  /** Equivalent to write */
  template <class ClassType>
  FMpiBufferWriter& operator<<(const ClassType& object){
    write(object);
    return *this;
  }
  
  /** Reset the writing index, but do not change the capacity */
  void reset(){
    array.clear();
  }
};


#endif // FBUFFERWRITER_HPP
