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

typedef long long FSize;

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

