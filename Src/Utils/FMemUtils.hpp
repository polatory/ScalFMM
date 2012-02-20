// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
// ===================================================================================
#ifndef FMEMUTILS_HPP
#define FMEMUTILS_HPP

#include "FGlobal.hpp"

// To get memcpy
#include <cstring>
#include <climits>


/** The memory utils class proposes some methods
  * to copy/set memory with an size bigger than size_t
  */
namespace FMemUtils {
    static const FSize MaxSize_t = UINT_MAX; //std::numeric_limits<std::size_t>::max();

    /** memcpy */
    static void* memcpy(void* const dest, const void* const source, const FSize nbBytes){
        if( nbBytes < MaxSize_t){
            return ::memcpy(dest, source, size_t(nbBytes));
        }
        else{
            char* iterDest          = static_cast<char*>(dest);
            const char* iterSource  = static_cast<const char*>(source);

            for(FSize idx = 0 ; idx < nbBytes - MaxSize_t ; idx += MaxSize_t ){
                ::memcpy(iterDest, iterSource, size_t(MaxSize_t));
                iterDest += MaxSize_t;
                iterSource += MaxSize_t;
            }
            ::memcpy(iterDest, iterSource, size_t(nbBytes%MaxSize_t));

            return dest;
        }
    }

    /** memset */
    static void* memset(void* const dest, const int val, const FSize nbBytes){
        if( nbBytes < MaxSize_t){
            return ::memset(dest, val, size_t(nbBytes));
        }
        else{
            char* iterDest  = static_cast<char*>(dest);

            for(FSize idx = 0 ; idx < nbBytes - MaxSize_t ; idx += MaxSize_t ){
                ::memset(iterDest, val, size_t(MaxSize_t));
                iterDest += MaxSize_t;
            }
            ::memset(iterDest, val, size_t(nbBytes%MaxSize_t));

            return dest;
        }
    }

    /** copy all value from one vector to the other */
    template <class TypeClass>
    void copyall(TypeClass*const dest, const TypeClass*const source, const int nbElements){
        for(int idx = 0 ; idx < nbElements ; ++idx){
            dest[idx] = source[idx];
        }
    }

    /** copy all value from one vector to the other */
    template <class TypeClass>
    void addall(TypeClass*const dest, const TypeClass*const source, const int nbElements){
        for(int idx = 0 ; idx < nbElements ; ++idx){
            dest[idx] += source[idx];
        }
    }

    /** Delete all */
    template <class TypeClass>
    void DeleteAll(TypeClass*const array[], const int size){
        for(int idx = 0 ; idx < size ; ++idx){
            delete[] array[idx];
        }
    }
}

#endif // FMEMUTILS_HPP
