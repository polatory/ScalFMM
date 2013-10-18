#ifndef FSSE_HPP
#define FSSE_HPP

#include <xmmintrin.h>  // SSE
#include <emmintrin.h>  //SSE2
#include <pmmintrin.h> //SSE3
#include <tmmintrin.h>  //SSSE3
//#include <smmintrin.h> // SSE4

#ifdef __INTEL_COMPILER

inline __m128d& operator+=(__m128d& v1, const __m128d& v2){
    return (v1 = _mm_add_pd(v1, v2));
}

inline __m128d& operator-=(__m128d& v1, const __m128d& v2){
    return (v1 = _mm_sub_pd(v1, v2));
}

inline __m128d& operator*=(__m128d& v1, const __m128d& v2){
    return (v1 = _mm_mul_pd(v1, v2));
}

inline __m128d& operator/=(__m128d& v1, const __m128d& v2){
    return (v1 = _mm_div_pd(v1, v2));
}

inline __m128d operator+(const __m128d& v1, const __m128d& v2){
    return _mm_add_pd(v1, v2);
}

inline __m128d operator-(const __m128d& v1, const __m128d& v2){
    return _mm_sub_pd(v1, v2);
}

inline __m128d operator*(const __m128d& v1, const __m128d& v2){
    return _mm_mul_pd(v1, v2);
}

inline __m128d operator/(const __m128d& v1, const __m128d& v2){
    return _mm_div_pd(v1, v2);
}

inline __m128& operator+=(__m128& v1, const __m128& v2){
    return (v1 = _mm_add_ps(v1, v2));
}

inline __m128& operator-=(__m128& v1, const __m128& v2){
    return (v1 = _mm_sub_ps(v1, v2));
}

inline __m128& operator*=(__m128& v1, const __m128& v2){
    return (v1 = _mm_mul_ps(v1, v2));
}

inline __m128& operator/=(__m128& v1, const __m128& v2){
    return (v1 = _mm_div_ps(v1, v2));
}

inline __m128 operator+(const __m128& v1, const __m128& v2){
    return _mm_add_ps(v1, v2);
}

inline __m128 operator-(const __m128& v1, const __m128& v2){
    return _mm_sub_ps(v1, v2);
}

inline __m128 operator*(const __m128& v1, const __m128& v2){
    return _mm_mul_ps(v1, v2);
}

inline __m128 operator/(const __m128& v1, const __m128& v2){
    return _mm_div_ps(v1, v2);
}

#endif

#endif // FSSE_HPP
