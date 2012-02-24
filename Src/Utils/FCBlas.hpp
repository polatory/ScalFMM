// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
// ===================================================================================
#ifndef FCBLAS_HPP
#define FCBLAS_HPP

#include "FGlobal.hpp"

#ifdef  SCALFMM_USE_MKL_AS_CBLAS
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

namespace FCBlas {
    ///////////////////////////////////////////////////////
    // GEMV
    ///////////////////////////////////////////////////////

    void Gemv(const CBLAS_ORDER order ,
                    const CBLAS_TRANSPOSE TransA , const int M , const int N ,
                    const double *alpha , const double *A , const int lda ,
                    const double *X , const int incX , const double *beta ,
                    double *Y , const int incY){
            cblas_zgemv(order,TransA,M,N,alpha,A,lda,X,incX,beta,Y,incY);
    }

    void Gemv(const CBLAS_ORDER order ,
                    const CBLAS_TRANSPOSE TransA , const int M , const int N ,
                    const float *alpha , const float *A , const int lda ,
                    const float *X , const int incX , const float *beta ,
                    float *Y , const int incY){
            cblas_cgemv(order,TransA,M,N,alpha,A,lda,X,incX,beta,Y,incY);
    }

    ///////////////////////////////////////////////////////
    // GEMM
    ///////////////////////////////////////////////////////

    void Gemm(const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB,
              const int m, const int n, const int k,
              const double*const alpha, const double*const A, const int ldA,
              const double*const B, const int ldB, const double*const beta, double*const C, const int ldC){
        cblas_zgemm(  order, TransA, TransB,
                         m, n, k,
                         alpha, A,ldA,
                         B,ldB,beta,
                         C, ldC);
    }

    void Gemm(const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB,
              const int m, const int n, const int k,
              const float*const alpha, const float*const A, const int ldA,
              const float*const B, const int ldB, const float*const beta, float*const C, const int ldC){
        cblas_cgemm(  order, TransA, TransB,
                         m, n, k,
                         alpha, A,ldA,
                         B,ldB,beta,
                         C, ldC);
    }


    ///////////////////////////////////////////////////////
    // Dotu
    ///////////////////////////////////////////////////////

    void Dotu_sub(const int N , const double *X , const int incX , const double *Y , const int incY , double *dotu){
            cblas_zdotu_sub(N,X,incX,Y,incY,dotu);
    }

    void Dotu_sub(const int N , const float *X , const int incX , const float *Y , const int incY , float *dotu){
            cblas_cdotu_sub(N,X,incX,Y,incY,dotu);
    }
}

#endif // FCBLAS_HPP
