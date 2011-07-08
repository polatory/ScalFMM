#ifndef FBUFFERVECTOR_HPP
#define FBUFFERVECTOR_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Utils/FAssertable.hpp"
// To get memcpy
#include <cstring>

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FBufferVector
* Please read the license
* This class is a buffer to send big data via MPI.
* It has a capacity and alloc an array only if needed.
* You have to be sure that the buffer is not full before adding any data.
* This can be done by checking :
* <code>
* if(!sendBuffer.addEnoughSpace(object->bytesToSendUp())){<br>
*    app.sendData(idxReceiver,sendBuffer.getSize(),sendBuffer.getData(),globalTag);<br>
*    sendBuffer.clear();<br>
* }<br>
* sendBuffer.addDataUp(tag,object);<br>
* </code>
*/
template <int Capacity>
        class FBufferVector : protected FAssertable {
    // the Buffer (can be null)
    char* buffer;
    // current ocupped size
    int occuped;

    // check if buffer is allocated if not does it
    void checkBuffer(){
        if(!buffer){
            buffer = new char[Capacity];
        }
    }

public:
    // Constructor, set buffer to null
    FBufferVector() : buffer(0), occuped(0) {
        FAssertable::assert(Capacity > 0, "Capacity has to be positive", __LINE__, __FILE__);
    }

    // Dealloc buffer if needed
    ~FBufferVector(){
        if(buffer) delete(buffer);
    }

    // To check if there is enough free space
    bool addEnoughSpace(const int neededSpace) const{
        return occuped + neededSpace + sizeof(int) <= Capacity;
    }

    // get the capacity of the buffer
    int getCapacity() const{
        return Capacity;
    }

    // get the current size
    int getSize(){
        return occuped;
    }

    // get a pointer to the buffer
    void* getData() {
        return buffer;
    }

    // ask to an object to write himself into the buffer
    // and update occuped space
    template < class T >
    void addDataUp(const int tag, const T& object){
        checkBuffer();
        memcpy(&buffer[occuped],&tag,sizeof(int));
        occuped += sizeof(int);
        occuped += object.writeUp(&buffer[occuped], Capacity - occuped);
    }

    // ask to an object to write himself into the buffer
    // and update occuped space
    template < class T >
    void addDataDown(const int tag, const T& object){
        checkBuffer();
        memcpy(&buffer[occuped],&tag,sizeof(int));
        occuped += sizeof(int);
        occuped += object.writeDown(&buffer[occuped], Capacity - occuped);
    }

    // reset occuped memory
    void clear(){
        occuped = 0;
    }
};


#endif //FBUFFERVECTOR_HPP

// [--LICENSE--]
