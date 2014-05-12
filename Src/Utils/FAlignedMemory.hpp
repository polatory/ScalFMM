// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas, Matthias Messner
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
#ifndef FALIGNEDMEMORY_HPP
#define FALIGNEDMEMORY_HPP

#include <cstdint>

/**
 * This should be used to allocate and deallocate aligned memory.
 */
namespace FAlignedMemory {

/**
 * @brief Allocate16BAligned
 * @param inSize in Bytes
 * @return the address of the allocated block of size inSize
 */
inline void* Allocate16BAligned(const size_t inSize){
    unsigned char* memoryBlock = new unsigned char[inSize + 15 + 16];
    {
        unsigned char** storeRealAddress = reinterpret_cast<unsigned char**>( (reinterpret_cast<long long int>(memoryBlock) + 15) & ~0xFLL);
        (*storeRealAddress) = memoryBlock;
    }
    return reinterpret_cast<void*>( (reinterpret_cast<long long int>(memoryBlock) + 16 + 15) & ~0xFLL);
}

  /**
   * @brief Allocate32BAligned
   * @param inSize in Bytes
   * @return the address of the allocated block of size inSize
   */
  inline void* Allocate32BAligned(const size_t inSize){
    unsigned char* memoryBlock;
    int resMemAl = posix_memalign((void**)&memoryBlock,32,inSize);
    if(resMemAl != 0){
      fprintf(stderr,"Allocation failed : Error Code : %d",resMemAl);
    }
    return memoryBlock;
  }

/**
 * Allocate an array of inNbElements elements of type  ObjectType.
 * The objects are not initialized!
 */
template <class ObjectType>
inline ObjectType* AllocateType16BAligned(const int inNbElements){
    return reinterpret_cast<ObjectType*>(Allocate16BAligned(sizeof(ObjectType) * inNbElements));
}

/**
 * Delete a block allocated with allocate16BAligned
 */
inline void Dealloc16BAligned(const void*const memoryBlock){
    if( memoryBlock ){
        const unsigned char*const* storeRealAddress = reinterpret_cast<const unsigned char*const *>(reinterpret_cast<const unsigned char*>(memoryBlock) - 16);
        delete[] reinterpret_cast<const unsigned char*>(*storeRealAddress);
    }
}

/**
 * Delete a block allocated with Allocate32BAligned
 */
  inline void Dealloc32BAligned(const void*const memoryBlock){
    const void * toBeFreed = memoryBlock;
    free(const_cast<void *>(toBeFreed));
  }

}


#endif // FALIGNEDMEMORY_HPP
