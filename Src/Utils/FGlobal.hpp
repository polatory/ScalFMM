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
#ifndef FGLOBAL_HPP
#define FGLOBAL_HPP

#include "ScalFmmConfig.h"

///////////////////////////////////////////////////////
// Memory profiling
///////////////////////////////////////////////////////

#include "FMemStats.h"

///////////////////////////////////////////////////////
// Stdlib
///////////////////////////////////////////////////////

#include <cstdlib>

///////////////////////////////////////////////////////
// Operating System
///////////////////////////////////////////////////////

#if defined(_WIN32) || defined(ming)
    #define WINDOWS
#else
    #define POSIX
#endif


///////////////////////////////////////////////////////
// Types
///////////////////////////////////////////////////////

typedef long long int FSize;

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
    #define Prefetch_Read(X)  __builtin_prefetch(X)
    #define Prefetch_Write(X) __builtin_prefetch(X,1,1)
#else
    #ifdef __INTEL_COMPILER
        #define Prefetch_Read(X)  _mm_prefetch(X,_MM_HINT_T0)
        #define Prefetch_Write(X) _mm_prefetch(X,_MM_HINT_T0)
    #else
        #warning compiler is not defined
        #define Prefetch_Read(X)
        #define Prefetch_Write(X)
    #endif
#endif


///////////////////////////////////////////////////////
// Test OMP4
///////////////////////////////////////////////////////

#if _OPENMP >= 201307
#ifndef __INTEL_COMPILER
#define SCALFMM_USE_OMP4
#endif
#endif


///////////////////////////////////////////////////////
// Default P2P Alignement
///////////////////////////////////////////////////////

static const int FP2PDefaultAlignement = 64;


#endif //FGLOBAL_HPP

