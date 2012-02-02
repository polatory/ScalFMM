#ifndef FBLAS_HPP
#define FBLAS_HPP

// This file interfaces the blas functions
// to enable a generic use.
// If no blas has been enabled in the cmake,
// the function will be empty

///////////////////////////////////////////////////////
// Manage Blas Version
///////////////////////////////////////////////////////

#include "FGlobal.hpp"



#ifdef SCALFMM_USE_CBLAS

    #ifdef  SCALFMM_USE_MKL_AS_BLAS
    #include <mkl_cblas.h>
    #else
    #include <cblas.h>
    #endif
#else
    enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
    enum CBLAS_TRANSPOSE {CblasNoTrans, CblasTrans, CblasConjTrans};
#endif

template <typename T>
void cblas_gemv(const CBLAS_ORDER order ,
                   const CBLAS_TRANSPOSE TransA , const int M , const int N ,
                   const void *alpha , const void *A , const int lda ,
                   const void *X , const int incX , const void *beta ,
                   void *Y , const int incY){
}
template <typename T>
void cblas_dotu_sub( const int N , const void *X , const int incX ,
                     const void *Y , const int incY , void *dotu){
}

#ifdef SCALFMM_USE_CBLAS
    ///////////////////////////////////////////////////////
    // GEMV
    ///////////////////////////////////////////////////////

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


    ///////////////////////////////////////////////////////
    // Dotu
    ///////////////////////////////////////////////////////

    template <>
    void cblas_dotu_sub<double>(const int N , const void *X , const int incX ,
                            const void *Y , const int incY , void *dotu){
        cblas_zdotu_sub(N,X,incX,Y,incY,dotu);
    }

    template <>
    void cblas_dotu_sub<float>(const int N , const void *X , const int incX ,
                           const void *Y , const int incY , void *dotu){
        cblas_cdotu_sub(N,X,incX,Y,incY,dotu);
    }

#endif //SCALFMM_USE_CBLAS



#endif //FBLAS_HPP

