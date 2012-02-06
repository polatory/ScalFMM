// ===================================================================================
// Ce LOGICIEL "ScalFmm" est couvert par le copyright Inria 20xx-2012.
// Inria détient tous les droits de propriété sur le LOGICIEL, et souhaite que
// la communauté scientifique l'utilise afin de le tester et de l'évaluer.
// Inria donne gracieusement le droit d'utiliser ce LOGICIEL. Toute utilisation
// dans un but lucratif ou à des fins commerciales est interdite sauf autorisation
// expresse et préalable d'Inria.
// Toute utilisation hors des limites précisées ci-dessus et réalisée sans l'accord
// expresse préalable d'Inria constituerait donc le délit de contrefaçon.
// Le LOGICIEL étant un produit en cours de développement, Inria ne saurait assurer
// aucune responsabilité et notamment en aucune manière et en aucun cas, être tenu
// de répondre d'éventuels dommages directs ou indirects subits par l'utilisateur.
// Tout utilisateur du LOGICIEL s'engage à communiquer à Inria ses remarques
// relatives à l'usage du LOGICIEL
// ===================================================================================
#ifndef FVECTOR_HPP
#define FVECTOR_HPP


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
    typedef T ValueType; /**< data type of data in FVector */

    /**
    *@brief constructor
    */
    FVector() : array(0), capacity(DefaultSize), index(0) {
        array = reinterpret_cast< T* >( new char[SizeOfT*DefaultSize] );
    }

    /**
    *@brief constructor
    *@param inSize the buffer size
    *@param inPointOfStart the point of start [0;1]
    */
    FVector(const int inSize): array(0), capacity(inSize), index(0) {
        array = reinterpret_cast< T* >( new char[SizeOfT*inSize]);
    }

    /**
     * Copy constructor
     * @param other original vector
     */
    FVector(const FVector& other): array(0), capacity(other.capacity), index(other.index) {
        array = reinterpret_cast< T* >( new char[SizeOfT*other.capacity]);
        // Copy each element
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
        if(&other != this){
            // clear current element
            clear();
            // alloc bigger if needed
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
    * if size == 10 and setCapactity(5), then it will
    * set the capacity to 10.
    */
    void setCapacity(int inCapacity){
        if( inCapacity != this->capacity ){
            // if the capacity is not enought
            if( inCapacity < this->index){
                inCapacity = this->index;
            }

            // Copy elements
            T* const buffer = reinterpret_cast< T* >( new char[SizeOfT*inCapacity]);
            memcpy( buffer, this->array, SizeOfT*this->index);
            delete [] reinterpret_cast< char* >(this->array);
            this->array = buffer;

            this->capacity = inCapacity;
        }
    }


    /**
    * @return end->data
    * This function return the last inserted data
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
    * this is not the capcity
    */
    int getSize() const{
        return this->index;
    }

    /**
    *@brief pop the first value
    * Warning, FVector do not check that there is an element before poping
    */
    void pop(){
        (&this->array[--this->index])->~T();
    }

    /**
    *@param inValue the new value
    */
    void push( const T & inValue ){
        // if needed, increase the vector
        if( this->index == this->capacity ){
            setCapacity(int(this->capacity * 1.5));
        }
        // add the new element
        new((void*)&this->array[this->index]) T;
        this->array[this->index] = inValue;
        ++this->index;
    }

    /**
    *@brief get a reference of a given value
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
        /** Empty destructor */
        virtual ~BasicIterator(){}

        /** Constructor need a vector */
        explicit BasicIterator(FVector<T>& inVector) : vector(&inVector), index(0){}

        /** Go to next vector element */
        void gotoNext(){
            ++this->index;
        }

        /** is it over
          * @return true if we are over the vector
          */
        bool hasNotFinished() const{
            return this->index < this->vector->index;
        }

        /** Get current data */
        T& data(){
            return this->vector->array[this->index];
        }

        /** Get current data */
        const T& data() const{
            return this->vector->array[this->index];
        }

        /** Set the data */
        void setData(const T& inData){
            this->vector->array[this->index] = inData;
        }

        /** Remove current data
          * It will move all the data after to their previous position
          */
        void remove(){
            if( hasNotFinished() ){
                for(int idxMove = this->index + 1; idxMove < this->vector->index ; ++idxMove){
                    vector->array[idxMove - 1] = vector->array[idxMove];
                }
                this->vector->index -= 1;
            }
        }

    };
    friend class BasicIterator;

    /** This class is a basic const iterator
      * it uses a const vector to work on
      */
    class ConstBasicIterator{
    protected:
        const FVector* const vector;  /**< the vector to work on*/
        int index;                    /**< the current node*/

    public:
        /** Empty destructor */
        virtual ~ConstBasicIterator(){}

        /** Constructor need a vector */
        explicit ConstBasicIterator(const FVector<T>& inVector) : vector(&inVector), index(0){}

        /** Go to next vector element */
        void gotoNext(){
            ++this->index;
        }

        /** is it over
          * @return true if we are over the vector
          */
        bool hasNotFinished() const{
            return this->index < this->vector->index;
        }

        /** Get current data */
        const T& data() const{
            return this->vector->array[this->index];
        }

    };
    friend class ConstBasicIterator;

};



#endif // FVECTOR_HPP

