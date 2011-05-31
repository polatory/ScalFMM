#ifndef FVECTOR_HPP
#define FVECTOR_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Utils/FGlobal.hpp"

#include <string.h>

/**
 * @author Berenger Bramas (berenger.bramas@inria.fr)
 * @class FList
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

    FVector(const FVector& other): array(reinterpret_cast< T* >( new char[SizeOfT*other.capacity]))
        , capacity(other.capacity), index(other.index) {
        for(int idx = 0 ; idx < other.index ; ++idx){
            new((void*)&this->array[idx]) T;
            this->array[idx] = other.array[idx];
        }
    }

    FVector& operator=(const FVector& other){
        clear();
        if(other.getSize() > this->capacity){
            delete [] reinterpret_cast< char* >(this->array);

            this->capacity = other.getSize() * 1.5;
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
    *@brief get the buffer capacity
    *@return the buffer capacity
    */
    int getCapacity() const{
        return this->capacity;
    }

    /**
    *@brief set the buffer capacity
    *@param inthis->capacity to change the capacity
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
    *@brief last data
    *@return end->data
    */
    const T& head() const {
        return this->array[this->index - 1];
    }

    /**
    *@brief last data
    *@return end->data
    */
    T& head() {
        return this->array[this->index - 1];
    }

    /**
    *@brief delete all, then size = 0
    */
    void clear(){
        while(this->index > 0){
            (&this->array[--this->index])->~T();
        }
    }

    /**
    *@brief count the value
    *@param inValue the value to test
    *@return the value occured number
    */
    int getSize() const{
        return this->index;
    }

    /**
    *@brief pop the first node
    */
    void pop(){
        (&this->array[--this->index])->~T();
    }

    /**
    *@brief push a new node in head
    *@param inValue the new value
    */
    void push( const T & inValue ){
        if( this->index == this->capacity ){
            setCapacity(this->capacity * 1.5);
        }
        new((void*)&this->array[this->index]) T;
        this->array[this->index] = inValue;
        ++this->index;
    }


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

