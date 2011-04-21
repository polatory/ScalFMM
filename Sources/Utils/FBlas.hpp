#ifndef FBLAS_HPP
#define FBLAS_HPP
#include "ScalFMM_config.h"
#ifdef  SCALFMM_USE_MKL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

///////////////////////////////////////////////////////
// GEMV
///////////////////////////////////////////////////////

template <typename T>
void cblas_gemv(const CBLAS_ORDER order ,
                   const CBLAS_TRANSPOSE TransA , const int M , const int N ,
                   const void *alpha , const void *A , const int lda ,
                   const void *X , const int incX , const void *beta ,
                   void *Y , const int incY){
    T t;
    t.you_cannot_use_this_function_with_this_type();
}

template <>
void cblas_gemv<double>(const CBLAS_ORDER order ,
        const CBLAS_TRANSPOSE TransA , const int M , const int N ,
        const void *alpha , const void *A , const int lda ,
        const void *X , const int incX , const void *beta ,
        void *Y , const int incY){
    cblas_zgemv(order,TransA,M,N,alpha,A,lda,X,incX,beta,Y,incY);
}

template <>
void cblas_gemv<float>(const CBLAS_ORDER order ,
        const CBLAS_TRANSPOSE TransA , const int M , const int N ,
        const void *alpha , const void *A , const int lda ,
        const void *X , const int incX , const void *beta ,
        void *Y , const int incY){
    cblas_cgemv(order,TransA,M,N,alpha,A,lda,X,incX,beta,Y,incY);
}

#endif //FBLAS_HPP

