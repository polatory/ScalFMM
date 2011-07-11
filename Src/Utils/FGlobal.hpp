#ifndef FGLOBAL_HPP
#define FGLOBAL_HPP

#include "ScalFmmConfig.h"
#include <stdlib.h>

///////////////////////////////////////////////////////
// Operating System
///////////////////////////////////////////////////////

#if defined(_WIN32) || defined(ming)
    #define WINDOWS
#else
    #define POSIX
#endif

///////////////////////////////////////////////////////
// Debug
///////////////////////////////////////////////////////

// Uncomment the next line to use debug mode
#define SCALFMM_USE_DEBUG

///////////////////////////////////////////////////////
// Debug
///////////////////////////////////////////////////////

// Uncomment the next line to use trace mode
//#define SCALFMM_USE_TRACE

///////////////////////////////////////////////////////
// Types
///////////////////////////////////////////////////////

typedef float FReal;

///////////////////////////////////////////////////////
// Restrict
///////////////////////////////////////////////////////

static const int MaxTreeHeight = 20;


///////////////////////////////////////////////////////
// Morton index
///////////////////////////////////////////////////////

typedef long long MortonIndex;


///////////////////////////////////////////////////////
// Restrict
///////////////////////////////////////////////////////

#ifdef WINDOWS
    #define FRestrict __restrict
#else
    #define FRestrict __restrict__
#endif

///////////////////////////////////////////////////////
// Prefetch
///////////////////////////////////////////////////////

#ifdef __GNUC__
    #define Prefetch_Read(X) __builtin_prefetch(X)
    #define Prefetch_Write(X) __builtin_prefetch(X,1,1)
#else
    #ifdef __INTEL_COMPILER
        #define Prefetch_Read(X) _mm_prefetch(X,_MM_HINT_T0)
        #define Prefetch_Write(X) _mm_prefetch(X,_MM_HINT_T0)
    #else
        #warning compiler is not defined
        #define Prefetch_Read(X)
        #define Prefetch_Write(X)
    #endif
#endif

#endif //FGLOBAL_HPP

