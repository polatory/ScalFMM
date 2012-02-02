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

    ///////////////////////////////////////////////////////
    // GEMV
    ///////////////////////////////////////////////////////

    void cblas_gemv(const CBLAS_ORDER order ,
            const CBLAS_TRANSPOSE TransA , const int M , const int N ,
            const double *alpha , const double *A , const int lda ,
            const double *X , const int incX , const double *beta ,
            double *Y , const int incY){
        cblas_zgemv(order,TransA,M,N,alpha,A,lda,X,incX,beta,Y,incY);
    }

    void cblas_gemv(const CBLAS_ORDER order ,
            const CBLAS_TRANSPOSE TransA , const int M , const int N ,
            const float *alpha , const float *A , const int lda ,
            const float *X , const int incX , const float *beta ,
            float *Y , const int incY){
        cblas_cgemv(order,TransA,M,N,alpha,A,lda,X,incX,beta,Y,incY);
    }


    ///////////////////////////////////////////////////////
    // Dotu
    ///////////////////////////////////////////////////////

    void cblas_dotu_sub(const int N , const double *X , const int incX ,
                            const double *Y , const int incY , double *dotu){
        cblas_zdotu_sub(N,X,incX,Y,incY,dotu);
    }

    void cblas_dotu_sub(const int N , const float *X , const int incX ,
                           const float *Y , const int incY , float *dotu){
        cblas_cdotu_sub(N,X,incX,Y,incY,dotu);
    }

#else
    enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
    enum CBLAS_TRANSPOSE {CblasNoTrans, CblasTrans, CblasConjTrans};

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

#endif //SCALFMM_USE_CBLAS

#endif //FBLAS_HPP

