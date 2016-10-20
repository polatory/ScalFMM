// ===================================================================================
// Copyright ScalFmm 2016 INRIA
//
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by Mozilla Public License Version 2.0 (MPL 2.0) and
// abiding by the rules of distribution of free software.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// Mozilla Public License Version 2.0 (MPL 2.0) for more details.
// https://www.mozilla.org/en-US/MPL/2.0/
// ===================================================================================
#ifndef FAVX2_HPP
#define FAVX2_HPP

#include "FGlobal.hpp"
#ifndef SCALFMM_USE_AVX2
#error The AVX header is included while SCALFMM_USE_AVX is turned OFF
#endif

#include "immintrin.h"

#ifdef __MIC__

//Side effect operators DOUBLE
inline __m512d& operator+=(__m512d & a, const __m512d & b){
  return (a = _mm512_add_pd (a,b));
}

inline __m512d& operator-=(__m512d& a, const __m512d& b){
  return (a = _mm512_sub_pd (a,b));
}

inline __m512d& operator*=(__m512d& a, const __m512d& b){
  return (a = _mm512_mul_pd (a,b));
}

inline __m512d& operator/=(__m512d& a, const __m512d& b){
  return (a = _mm512_div_pd (a,b));
}

//No side effect operators DOUBLE
inline __m512d operator+(const __m512d& a,const  __m512d& b){
  return _mm512_add_pd (a,b);
}

inline __m512d operator-(const __m512d& a, const __m512d& b){
  return _mm512_sub_pd (a,b);
}

inline __m512d operator*(const __m512d& v1, const __m512d& v2){
    return _mm512_mul_pd(v1, v2);
}

inline __m512d operator/(const __m512d& v1, const __m512d& v2){
    return _mm512_div_pd(v1, v2);
}

//Side effect operators SINGLE
inline __m512& operator+=(__m512 & a, const __m512 & b){
  return (a = _mm512_add_ps (a,b));
}

inline __m512& operator-=(__m512& a, const __m512& b){
  return (a = _mm512_sub_ps (a,b));
}

inline __m512& operator*=(__m512& a, const __m512& b){
  return (a = _mm512_mul_ps (a,b));
}

inline __m512& operator/=(__m512& a, const __m512& b){
  return (a = _mm512_div_ps (a,b));
}

//No side effect operators SINGLE
inline __m512 operator+(const __m512& a,const  __m512& b){
  return _mm512_add_ps (a,b);
}

inline __m512 operator-(const __m512& a, const __m512& b){
  return _mm512_sub_ps (a,b);
}

inline __m512 operator*(const __m512& v1, const __m512& v2){
    return _mm512_mul_ps(v1, v2);
}

inline __m512 operator/(const __m512& v1, const __m512& v2){
    return _mm512_div_ps(v1, v2);
}

#endif

#endif

