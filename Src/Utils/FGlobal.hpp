#ifndef FGLOBAL_HPP
#define FGLOBAL_HPP

#include "ScalFmmConfig.h"

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
// MPI
///////////////////////////////////////////////////////

#define SCALFMM_USE_MPI

///////////////////////////////////////////////////////
// Threads
///////////////////////////////////////////////////////

static const int FThreadNumbers = 1;

///////////////////////////////////////////////////////
// Types
///////////////////////////////////////////////////////

typedef float FReal;

///////////////////////////////////////////////////////
// Restrict
///////////////////////////////////////////////////////

#ifdef WINDOWS
    #define FRestrict __restrict
#else
    #define FRestrict __restrict__
#endif

#endif //FGLOBAL_HPP

