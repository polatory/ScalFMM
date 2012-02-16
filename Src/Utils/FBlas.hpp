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



//#ifdef SCALFMM_USE_CBLAS
//
//#ifdef  SCALFMM_USE_MKL_AS_BLAS
//#include <mkl_cblas.h>
//#else
//#include <cblas.h>
//#include <clapack.h>
//#endif
//#else
//    enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
//    enum CBLAS_TRANSPOSE {CblasNoTrans, CblasTrans, CblasConjTrans};
//#endif
//
/////////////////////////////////////////////////////////
//// GEMV
/////////////////////////////////////////////////////////
//
//void cblas_gemv(const CBLAS_ORDER order ,
//								const CBLAS_TRANSPOSE TransA , const int M , const int N ,
//								const double *alpha , const double *A , const int lda ,
//								const double *X , const int incX , const double *beta ,
//								double *Y , const int incY){
//	cblas_zgemv(order,TransA,M,N,alpha,A,lda,X,incX,beta,Y,incY);
//}
//
//void cblas_gemv(const CBLAS_ORDER order ,
//								const CBLAS_TRANSPOSE TransA , const int M , const int N ,
//								const float *alpha , const float *A , const int lda ,
//								const float *X , const int incX , const float *beta ,
//								float *Y , const int incY){
//	cblas_cgemv(order,TransA,M,N,alpha,A,lda,X,incX,beta,Y,incY);
//}
//
//
/////////////////////////////////////////////////////////
//// Dotu
/////////////////////////////////////////////////////////
//
//void cblas_dotu_sub(const int N , const double *X , const int incX ,
//										const double *Y , const int incY , double *dotu){
//	cblas_zdotu_sub(N,X,incX,Y,incY,dotu);
//}
//
//void cblas_dotu_sub(const int N , const float *X , const int incX ,
//										const float *Y , const int incY , float *dotu){
//	cblas_cdotu_sub(N,X,incX,Y,incY,dotu);
//}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

// for real
const double D_ZERO =  0.0;
const double D_ONE  =  1.0;
const double D_MONE = -1.0;
const float  S_ZERO =  0.0;
const float  S_ONE  =  1.0;
const float  S_MONE = -1.0;
// for complex
const double Z_ZERO =  0.0;
const double Z_ONE  =  1.0;
const double Z_MONE = -1.0;
const float  C_ZERO =  0.0;
const float  C_ONE  =  1.0;
const float  C_MONE = -1.0;

//const double D_PREC = 1e-16;

const unsigned N_ONE = 1;
const int N_MONE = -1;
const char JOB_STR[] = "NTOSVULCR";


extern "C"
{
	// double //////////////////////////////////////////////////////////
	// blas 1
	double ddot_(const unsigned*, const double*, const unsigned*, const double*, const unsigned*);
	void dscal_(const unsigned*, const double*, const double*, const unsigned*);
	// blas 2
	void dgemv_(const char*, const unsigned*, const unsigned*, const double*,
							const double*, const unsigned*, const double*, const unsigned*,
							const double*, double*, const unsigned*);
	void dcopy_(const unsigned*, const double*, const unsigned*, double*, const unsigned*);
	// lapack
	void dgesvd_(const char*, const char*, const unsigned*, const unsigned*,
							 double*, const unsigned*, double*, double*, const unsigned*,
							 double*, const unsigned*, double*, const unsigned*, int*);
	
	// single //////////////////////////////////////////////////////////
	// blas 1
	float sdot_(const unsigned*, const float*, const unsigned*,	const float*, const unsigned*);
	void sscal_(const unsigned*, const float*, const float*, const unsigned*);
	// blas 2
	void sgemv_(const char*, const unsigned*, const unsigned*, const float*,
							const float*, const unsigned*, const float*, const unsigned*,
							const float*, float*, const unsigned*);
	void scopy_(const unsigned*, const float*, const unsigned*,	float*, const unsigned*);

	// lapack
	void sgesvd_(const char*, const char*, const unsigned*, const unsigned*,
							 float*, const unsigned*, float*, float*, const unsigned*,
							 float*, const unsigned*, float*, const unsigned*, int*);

	// double complex //////////////////////////////////////////////////
	void zgemv_(const char*, const unsigned*, const unsigned*, const double*,
							const double*, const unsigned*, const double*, const unsigned*,
							const double*, double*, const unsigned*);
//	void zcopy_(const unsigned*, const double*, const unsigned*,
//							double*, const unsigned*);

	// single complex //////////////////////////////////////////////////
	void cgemv_(const char*, const unsigned*, const unsigned*, const float*,
							const float*, const unsigned*, const float*, const unsigned*,
							const float*, float*, const unsigned*);
//	void ccopy_(const unsigned*, const float*, const unsigned*,
//							float*, const unsigned*);
}


namespace FBlas {

	// copy
	inline void copy(const unsigned n, double* orig, double* dest)
	{	dcopy_(&n, orig, &N_ONE, dest, &N_ONE);	}
	inline void copy(const unsigned n, float* orig, float* dest)
	{	scopy_(&n, orig, &N_ONE, dest, &N_ONE);	}
	// copy (different increment)
	inline void copy(const unsigned n, double* orig, const unsigned inco, double* dest, const unsigned incd)
	{	dcopy_(&n, orig, &inco, dest, &incd);	}
	inline void copy(const unsigned n, float* orig, const unsigned inco, float* dest, const unsigned incd)
	{	scopy_(&n, orig, &inco, dest, &incd);	}

	// scale
	inline void scal(const unsigned n, const double d, double* const x)
	{	dscal_(&n, &d, x, &N_ONE); }
	inline void scal(const unsigned n, const float d, float* const x)
	{	sscal_(&n, &d, x, &N_ONE); }


//	// y = d Ax
//	inline void gemv(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
//	{	cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, d, A, m, x, N_ONE, D_ZERO, y, N_ONE); }
//	inline void gemv(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
//	{	cblas_sgemv(CblasColMajor, CblasNoTrans, m, n, d, A, m, x, N_ONE, S_ZERO, y, N_ONE); }
	// y = d Ax
	inline void gemv(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
	{	dgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &D_ZERO, y, &N_ONE); }
	inline void gemv(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
	{	sgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &S_ZERO, y, &N_ONE); }
	inline void c_gemv(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
	{	cgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &C_ZERO, y, &N_ONE); }
	inline void c_gemv(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
	{	zgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &Z_ZERO, y, &N_ONE); }

//	// y += d Ax
//	inline void gemva(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
//	{	cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, d, A, m, x, N_ONE, D_ONE, y, N_ONE); }
//	inline void gemva(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
//	{	cblas_sgemv(CblasColMajor, CblasNoTrans, m, n, d, A, m, x, N_ONE, S_ONE, y, N_ONE); }
	// y += d Ax
	inline void gemva(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
	{	dgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &D_ONE, y, &N_ONE);	}
	inline void gemva(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
	{	sgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &S_ONE, y, &N_ONE);	}
	inline void c_gemva(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
	{	cgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &C_ONE, y, &N_ONE);	}
	inline void c_gemva(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
	{	zgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &Z_ONE, y, &N_ONE);	}

//	// y = d A^T x
//	inline void gemtv(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
//	{ cblas_dgemv(CblasColMajor, CblasTrans, m, n, d, A, m, x, N_ONE, D_ZERO, y, N_ONE); }
//	inline void gemtv(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
//	{	cblas_sgemv(CblasColMajor, CblasTrans, m, n, d, A, m, x, N_ONE, S_ZERO, y, N_ONE); }
	// y = d A^T x
	inline void gemtv(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
	{	dgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &D_ZERO, y, &N_ONE); }
	inline void gemtv(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
	{	sgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &S_ZERO, y, &N_ONE); }
	inline void c_gemtv(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
	{	cgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &C_ZERO, y, &N_ONE); }
	inline void c_gemtv(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
	{	zgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &Z_ZERO, y, &N_ONE); }
	inline void c_gemhv(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
	{	cgemv_(JOB_STR+7, &m, &n, &d, A, &m, x, &N_ONE, &C_ZERO, y, &N_ONE); } // hermitian transposed
	inline void c_gemhv(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
	{	zgemv_(JOB_STR+7, &m, &n, &d, A, &m, x, &N_ONE, &Z_ZERO, y, &N_ONE); } // hermitian transposed

//	// y += d A^T x
//	inline void gemtva(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
//	{	cblas_dgemv(CblasColMajor, CblasTrans, m, n, d, A, m, x, N_ONE, D_ONE, y, N_ONE); }
//	inline void gemtva(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
//	{	cblas_sgemv(CblasColMajor, CblasTrans, m, n, d, A, m, x, N_ONE, S_ONE, y, N_ONE); }
	// y += d A^T x
	inline void gemtva(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
	{	dgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &D_ONE, y, &N_ONE);	}
	inline void gemtva(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
	{	sgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &S_ONE, y, &N_ONE);	}
	inline void c_gemtva(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
	{	cgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &C_ONE, y, &N_ONE);	}
	inline void c_gemtva(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
	{	zgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &Z_ONE, y, &N_ONE); }
	inline void c_gemhva(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
	{	cgemv_(JOB_STR+7, &m, &n, &d, A, &m, x, &N_ONE, &C_ONE, y, &N_ONE);	} // hermitian transposed
	inline void c_gemhva(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
	{	zgemv_(JOB_STR+7, &m, &n, &d, A, &m, x, &N_ONE, &Z_ONE, y, &N_ONE);	} // hermitian transposed

	// singular value decomposition
	inline int gesvd(unsigned m, unsigned n, double* A, double* S, double* VT, unsigned ldVT, unsigned nwk, double* wk)
	{
		int INF;
		dgesvd_(JOB_STR+2, JOB_STR+3, &m, &n, A, &m, S, A, &m, VT, &ldVT,	wk, &nwk, &INF);
		return INF;
	}
	inline int gesvd(unsigned m, unsigned n, float* A, float* S, float* VT, unsigned ldVT, unsigned nwk, float* wk)
	{
		int INF;
		sgesvd_(JOB_STR+2, JOB_STR+3, &m, &n, A, &m, S, A, &m, VT, &ldVT,	wk, &nwk, &INF);
		return INF;
	}

  // singular value decomposition (SO)
  inline int gesvdSO(unsigned m, unsigned n, double* A, double* S, double* U, unsigned ldU, unsigned nwk, double* wk)
  {
    int INF;
    dgesvd_(JOB_STR+3, JOB_STR+2, &m, &n, A, &m, S, U, &m, A, &ldU, wk, &nwk, &INF);
    return INF;
  }
  inline int gesvdSO(unsigned m, unsigned n, float* A, float* S, float* U, unsigned ldU, unsigned nwk, float* wk)
  {
    int INF;
    sgesvd_(JOB_STR+3, JOB_STR+2, &m, &n, A, &m, S, U, &m, A, &ldU, wk, &nwk, &INF);
    return INF;
  }

	// Scalar product v1'*v2
	inline double scpr(const unsigned n, const double* const v1, const double* const v2)
	{	return ddot_(&n, v1, &N_ONE, v2, &N_ONE); }
	inline float scpr(const unsigned n, const float* const v1, const float* const v2)
	{	return sdot_(&n, v1, &N_ONE, v2, &N_ONE);	}

} // end namespace FCBlas

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


//#else
//enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
//enum CBLAS_TRANSPOSE {CblasNoTrans, CblasTrans, CblasConjTrans};
//
//template <typename T>
//void cblas_gemv(const CBLAS_ORDER order ,
//								const CBLAS_TRANSPOSE TransA , const int M , const int N ,
//								const void *alpha , const void *A , const int lda ,
//								const void *X , const int incX , const void *beta ,
//								void *Y , const int incY){
//}
//template <typename T>
//void cblas_dotu_sub( const int N , const void *X , const int incX ,
//										 const void *Y , const int incY , void *dotu){
//}


#endif //FBLAS_HPP

