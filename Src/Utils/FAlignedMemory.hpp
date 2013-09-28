#ifndef FALIGNEDMEMORY_HPP
#define FALIGNEDMEMORY_HPP

namespace FAlignedMemory {

inline void* Allocate16BAligned(const size_t inSize){
    unsigned char* memoryBlock = new unsigned char[inSize + 15 + 16];
    {
        unsigned char** storeRealAddress = reinterpret_cast<unsigned char**>( (reinterpret_cast<long long int>(memoryBlock) + 15) & ~0xFLL);
        (*storeRealAddress) = memoryBlock;
    }
    return reinterpret_cast<void*>( (reinterpret_cast<long long int>(memoryBlock) + 16 + 15) & ~0xFLL);
}

template <class ObjectType>
inline ObjectType* AllocateType16BAligned(const int inNbElements){
    return reinterpret_cast<ObjectType*>(Allocate16BAligned(sizeof(ObjectType) * inNbElements));
}

/** delete a block allocated with allocate16BAligned
*/
inline void Dealloc16BAligned(const void*const memoryBlock){
    if( memoryBlock ){
        const unsigned char*const* storeRealAddress = reinterpret_cast<const unsigned char*const *>(reinterpret_cast<const unsigned char*>(memoryBlock) - 16);
        delete[] reinterpret_cast<const unsigned char*>(*storeRealAddress);
    }
}

}


#endif // FALIGNEDMEMORY_HPP
