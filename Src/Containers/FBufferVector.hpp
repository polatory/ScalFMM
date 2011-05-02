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
*/
template <int Capacity>
        class FBufferVector : protected FAssertable {
    char* buffer;
    int occuped;

    void checkBuffer(){
        if(!buffer){
            buffer = new char[Capacity];
        }
    }

public:
    FBufferVector() : buffer(0), occuped(0) {
        assert(Capacity > 0, "Capacity has to be positive", __LINE__, __FILE__);
    }

    ~FBufferVector(){
        if(buffer) delete(buffer);
    }

    bool addEnoughSpace(const int neededSpace) const{
        return occuped + neededSpace + sizeof(int) <= Capacity;
    }

    int getCapacity() const{
        return Capacity;
    }

    int getSize(){
        return occuped;
    }

    void* getData() {
        return buffer;
    }

    template < class T >
    void addDataUp(const int tag, const T& object){
        checkBuffer();
        memcpy(&buffer[occuped],&tag,sizeof(int));
        occuped += sizeof(int);
        occuped += object.writeUp(&buffer[occuped], Capacity - occuped);
    }

    template < class T >
    void addDataDown(const int tag, const T& object){
        checkBuffer();
        memcpy(&buffer[occuped],&tag,sizeof(int));
        occuped += sizeof(int);
        occuped += object.writeDown(&buffer[occuped], Capacity - occuped);
    }

    void clear(){
        occuped = 0;
    }
};


#endif //FBUFFERVECTOR_HPP

// [--LICENSE--]
