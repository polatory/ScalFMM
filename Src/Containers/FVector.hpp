#ifndef FVECTOR_HPP
#define FVECTOR_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Utils/FGlobal.hpp"

// To get memcpy
#include <string.h>

/**
 * @author Berenger Bramas (berenger.bramas@inria.fr)
 * @class FVector
 * Please read the license
 *
 * This class is a vector container.
 * It is a very basic vector to enable strong performance.
 *
 * Please refere to unit test fvectorUTest.cpp to know more.
 */
template<class T>
class FVector {
protected:
    T* array;        /**< memory area*/

    int capacity;    /**< memory capacity*/
    int index;       /**< index in array  */

    static const int DefaultSize = 10;      /**< Default size */
    static const int SizeOfT = sizeof(T);   /**< size of the object*/

public:
    /**
    *@brief constructor
    */
    FVector() : array(reinterpret_cast< T* >( new char[SizeOfT*DefaultSize] )), capacity(DefaultSize), index(0) {
    }

    /**
    *@brief constructor
    *@param inSize the buffer size
    *@param inPointOfStart the point of start [0;1]
    */
    FVector(const int inSize): array(reinterpret_cast< T* >( new char[SizeOfT*inSize]))
        , capacity(inSize), index(0) {
    }

    /**
     * Copy constructor
     * @param other original vector
     */
    FVector(const FVector& other): array(reinterpret_cast< T* >( new char[SizeOfT*other.capacity]))
        , capacity(other.capacity), index(other.index) {
        for(int idx = 0 ; idx < other.index ; ++idx){
            new((void*)&this->array[idx]) T;
            this->array[idx] = other.array[idx];
        }
    }

    /** Copy operator
      * @param other the original vector
      * @return this after copying data
      * Objects of the current vector are deleted before
      * copying data.
      */
    FVector& operator=(const FVector& other){
        clear();
        if(other.getSize() > this->capacity){
            delete [] reinterpret_cast< char* >(this->array);

            this->capacity = int(other.getSize() * 1.5);
            array = reinterpret_cast< T* >( new char[SizeOfT*this->capacity]);
        }
        this->index = other.index;

        for(int idx = 0 ; idx < other.index ; ++idx){
            new((void*)&this->array[idx]) T;
            this->array[idx] = other.array[idx];
        }

        return *this;
    }

    /**
    *@brief destructor
    */
    virtual ~FVector(){
        clear();
        delete [] reinterpret_cast< char* >(this->array);
    }

    /**
    * @brief get the buffer capacity
    * @return the buffer capacity
    * The capacity is the current memory size allocated.
    */
    int getCapacity() const{
        return this->capacity;
    }

    /**
    *@brief set the buffer capacity
    *@param inthis->capacity to change the capacity
    * Warning : the memocy is duplicated using memcpy
    */
    void setCapacity(int inCapacity){
        if( inCapacity != this->capacity ){

            if( inCapacity < this->index){
                inCapacity = this->index;
            }

            T* const buffer = reinterpret_cast< T* >( new char[SizeOfT*inCapacity]);

            memcpy(buffer,this->array,SizeOfT*this->index);

            delete [] reinterpret_cast< char* >(this->array);
            this->array = buffer;

            this->capacity = inCapacity;
        }
    }


    /**
    * @return end->data
    * This function return the last insted data
    */
    const T& head() const {
        return this->array[this->index - 1];
    }

    /**
    *@return end->data
    * This function return the last insted data
    */
    T& head() {
        return this->array[this->index - 1];
    }

    /**
    *@brief delete all, then size = 0 but capacity is unchanged
    */
    void clear(){
        while(this->index > 0){
            (&this->array[--this->index])->~T();
        }
    }

    /**
    *@return The number of element added into the vector
    */
    int getSize() const{
        return this->index;
    }

    /**
    *@brief pop the first value
    */
    void pop(){
        (&this->array[--this->index])->~T();
    }

    /**
    *@param inValue the new value
    */
    void push( const T & inValue ){
        if( this->index == this->capacity ){
            setCapacity(int(this->capacity * 1.5));
        }
        new((void*)&this->array[this->index]) T;
        this->array[this->index] = inValue;
        ++this->index;
    }

    /**
    *@brief get a const reference of a given value
    *@param inPosition the query position
    *@return the value
    */
    T& operator[](const int inPosition ){
            return this->array[inPosition];
    }

    /**
    *@brief get a const reference of a given value
    *@param inPosition the query position
    *@return the value
    */
    const T& operator[](const int inPosition ) const{
            return this->array[inPosition];
    }

    /** To get the entire array
      * @return the array allocated by the vector
      */
    T* data(){
        return this->array;
    }

    /** To get the entire array
      * @return the array allocated by the vector
      */
    const T* data() const{
        return this->array;
    }


    /** This class is a basic iterator
      * <code>
      *  typename FVector<int>::ConstBasicIterator iter(myVector);<br>
      *  while( iter.hasNotFinished() ){<br>
      *      printf("%d\n",iter.data());<br>
      *      iter.gotoNext();<br>
      *  } <br>
      * </code>
      */
    class BasicIterator{
    protected:
        FVector* const vector;  /**< the vector to work on*/
        int index;              /**< the current node*/

    public:

        virtual ~BasicIterator(){}

        BasicIterator(FVector<T>& inVector) : vector(&inVector), index(0){}

        void gotoNext(){
            ++this->index;
        }

        bool hasNotFinished() const{
            return this->index < this->vector->index;
        }

        T& data(){
            return this->vector->array[this->index];
        }

        const T& data() const{
            return this->vector->array[this->index];
        }

    };
    friend class BasicIterator;

    /** This class is a basic const iterator
      * it uses a const vector to work on
      */
    class ConstBasicIterator{
    protected:
        const FVector* const vector;  /**< the vector to work on*/
        int index;              /**< the current node*/

    public:

        virtual ~ConstBasicIterator(){}

        ConstBasicIterator(const FVector<T>& inVector) : vector(&inVector), index(0){}

        void gotoNext(){
            ++this->index;
        }

        bool hasNotFinished() const{
            return this->index < this->vector->index;
        }

        const T& data() const{
            return this->vector->array[this->index];
        }

    };
    friend class ConstBasicIterator;

};



#endif // FVECTOR_HPP

