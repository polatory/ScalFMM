#ifndef FMEMUTILS_HPP
#define FMEMUTILS_HPP

#include "FGlobal.hpp"

// To get memcpy
#include <cstring>
#include <limits>

namespace FMemUtils {
    static const FSize MaxSize_t = std::numeric_limits<std::size_t>::max();

    static void* memcpy(void* const dest, const void* const source, const FSize nbBytes){
        if( nbBytes < MaxSize_t){
            return ::memcpy(dest, source, size_t(nbBytes));
        }
        else{
            char* iterDest          = static_cast<char*>(dest);
            const char* iterSource  = static_cast<const char*>(source);

            for(FSize idx = 0 ; idx < nbBytes ; idx += MaxSize_t ){
                ::memcpy(iterDest, iterSource, size_t(MaxSize_t));
                iterDest += MaxSize_t;
                iterSource += MaxSize_t;
            }
            ::memcpy(iterDest, iterSource, size_t(nbBytes%MaxSize_t));

            return dest;
        }
    }

    static void* memset(void* const dest, const int val, const FSize nbBytes){
        if( nbBytes < MaxSize_t){
            return ::memset(dest, val, size_t(nbBytes));
        }
        else{
            char* iterDest          = static_cast<char*>(dest);

            for(FSize idx = 0 ; idx < nbBytes ; idx += MaxSize_t ){
                ::memset(iterDest, val, size_t(MaxSize_t));
                iterDest += MaxSize_t;
            }
            ::memset(iterDest, val, size_t(nbBytes%MaxSize_t));

            return dest;
        }
    }

};

#endif // FMEMUTILS_HPP
