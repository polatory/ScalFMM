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

#include "FGlobal.hpp"

// This file interfaces the blas functions
// to enable a generic use.
// If no blas has been enabled in the cmake,
// the function will be empty


// for real
const double D_ZERO =  0.0;
const double D_ONE  =  1.0;
const double D_MONE = -1.0;
const float  S_ZERO =  0.0;
const float  S_ONE  =  1.0;
const float  S_MONE = -1.0;
// for complex
const double Z_ZERO[2] =  {0.0,0.0};
const double Z_ONE[2]  =  {1.0,0.0};
const double Z_MONE[2] =  {-1.0,0.0};
const float  C_ZERO[2] =  {0.0,0.0};
const float  C_ONE[2]  =  {1.0,0.0};
const float  C_MONE[2] =  {-1.0,0.0};

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
	void dcopy_(const unsigned*, const double*, const unsigned*, double*, const unsigned*);
	void daxpy_(const unsigned*, const double*, const double*, const unsigned*, double*, const unsigned*);
	// blas 2
	void dgemv_(const char*, const unsigned*, const unsigned*, const double*,
							const double*, const unsigned*, const double*, const unsigned*,
							const double*, double*, const unsigned*);
	// blas 3
	void dgemm_(const char*, const char*, const unsigned*, const unsigned*,
							const unsigned*, const double*, double*, const unsigned*,
							double*, const unsigned*, const double*, double*,	const unsigned*);
	// lapack
	void dgesvd_(const char*, const char*, const unsigned*, const unsigned*,
							 double*, const unsigned*, double*, double*, const unsigned*,
							 double*, const unsigned*, double*, const unsigned*, int*);
	
	// single //////////////////////////////////////////////////////////
	// blas 1
	float sdot_(const unsigned*, const float*, const unsigned*,	const float*, const unsigned*);
	void sscal_(const unsigned*, const float*, const float*, const unsigned*);
	void scopy_(const unsigned*, const float*, const unsigned*,	float*, const unsigned*);
	void saxpy_(const unsigned*, const float*, const float*, const unsigned*, float*, const unsigned*);
	// blas 2
	void sgemv_(const char*, const unsigned*, const unsigned*, const float*,
							const float*, const unsigned*, const float*, const unsigned*,
							const float*, float*, const unsigned*);
	// blas 3
	void sgemm_(const char*, const char*, const unsigned*, const unsigned*,
							const unsigned*, const float*, float*, const unsigned*,
							float*, const unsigned*, const float*, float*, const unsigned*);
	// lapack
	void sgesvd_(const char*, const char*, const unsigned*, const unsigned*,
							 float*, const unsigned*, float*, float*, const unsigned*,
							 float*, const unsigned*, float*, const unsigned*, int*);

	// double complex //////////////////////////////////////////////////
	// blas 1
	void zscal_(const unsigned*, const double*, const double*, const unsigned*);
	void zcopy_(const unsigned*, const double*, const unsigned*, double*, const unsigned*);
	void zaxpy_(const unsigned*, const double*, const double*, const unsigned*, double*, const unsigned*);
	// blas 2
	void zgemv_(const char*, const unsigned*, const unsigned*, const double*,
							const double*, const unsigned*, const double*, const unsigned*,
							const double*, double*, const unsigned*);
	// blas 3
	void zgemm_(const char*, const char*, const unsigned*, const unsigned*,
							const unsigned*, const double*, double*, const unsigned*,
							double*, const unsigned*, const double*, double*, const unsigned*);


	// single complex //////////////////////////////////////////////////
	// blas 1
	void cscal_(const unsigned*, const float*, const float*, const unsigned*);
	void ccopy_(const unsigned*, const float*, const unsigned*,	float*, const unsigned*);
	void caxpy_(const unsigned*, const float*, const float*, const unsigned*, float*, const unsigned*);
	// blas 2
	void cgemv_(const char*, const unsigned*, const unsigned*, const float*,
							const float*, const unsigned*, const float*, const unsigned*,
							const float*, float*, const unsigned*);
	// blas 3
	void cgemm_(const char*, const char*, const unsigned*, const unsigned*,
							const unsigned*, const float*, float*, const unsigned*,
							float*, const unsigned*, const float*, float*, const unsigned*);

}


namespace FBlas {

	// copy
	inline void copy(const unsigned n, double* orig, double* dest)
	{	dcopy_(&n, orig, &N_ONE, dest, &N_ONE);	}
	inline void copy(const unsigned n, float* orig, float* dest)
	{	scopy_(&n, orig, &N_ONE, dest, &N_ONE);	}
	inline void c_copy(const unsigned n, double* orig, double* dest)
	{	zcopy_(&n, orig, &N_ONE, dest, &N_ONE);	}
	inline void c_copy(const unsigned n, float* orig, float* dest)
	{	ccopy_(&n, orig, &N_ONE, dest, &N_ONE);	}

	// copy (variable increment)
	inline void copy(const unsigned n, double* orig, const unsigned inco, double* dest, const unsigned incd)
	{	dcopy_(&n, orig, &inco, dest, &incd);	}
	inline void copy(const unsigned n, float* orig, const unsigned inco, float* dest, const unsigned incd)
	{	scopy_(&n, orig, &inco, dest, &incd);	}
	inline void c_copy(const unsigned n, double* orig, const unsigned inco, double* dest, const unsigned incd)
	{	zcopy_(&n, orig, &inco, dest, &incd);	}
	inline void c_copy(const unsigned n, float* orig, const unsigned inco, float* dest, const unsigned incd)
	{	ccopy_(&n, orig, &inco, dest, &incd);	}

	// scale
	inline void scal(const unsigned n, const double d, double* const x)
	{	dscal_(&n, &d, x, &N_ONE); }
	inline void scal(const unsigned n, const float d, float* const x)
	{	sscal_(&n, &d, x, &N_ONE); }
	inline void c_scal(const unsigned n, const double d, double* const x)
	{	zscal_(&n, &d, x, &N_ONE); }
	inline void c_scal(const unsigned n, const float d, float* const x)
	{	cscal_(&n, &d, x, &N_ONE); }

	// scale (variable increment)
	inline void scal(const unsigned n, const double d, double* const x, const unsigned incd)
	{	dscal_(&n, &d, x, &incd); }
	inline void scal(const unsigned n, const float d, float* const x, const unsigned incd)
	{	sscal_(&n, &d, x, &incd); }
	inline void c_scal(const unsigned n, const double d, double* const x, const unsigned incd)
	{	zscal_(&n, &d, x, &incd); }
	inline void c_scal(const unsigned n, const float d, float* const x, const unsigned incd)
	{	cscal_(&n, &d, x, &incd); }

	// set zero
	inline void setzero(const unsigned n, double* const x)
	{	dscal_(&n, &D_ZERO, x, &N_ONE); }
	inline void setzero(const unsigned n, float* const x)
	{	sscal_(&n, &S_ZERO, x, &N_ONE); }
	inline void c_setzero(const unsigned n, double* const x)
        {	zscal_(&n, Z_ZERO, x, &N_ONE); }
	inline void c_setzero(const unsigned n, float* const x)
        {	cscal_(&n, C_ZERO, x, &N_ONE); }

	// y += x
	inline void add(const unsigned n, double* const x, double* const y)
	{	daxpy_(&n, &D_ONE, x, &N_ONE, y, &N_ONE);	}
	inline void add(const unsigned n, float* const x, float* const y)
	{	saxpy_(&n, &S_ONE, x, &N_ONE, y, &N_ONE);	}
	inline void c_add(const unsigned n, float* const x, float* const y)
        {	caxpy_(&n, C_ONE, x, &N_ONE, y, &N_ONE);	}
	inline void c_add(const unsigned n, double* const x,double* const y)
        {	zaxpy_(&n, Z_ONE, x, &N_ONE, y, &N_ONE);	}

	// y += d x
	inline void axpy(const unsigned n, const double d, const double* const x, double* const y)
	{	daxpy_(&n, &d, x, &N_ONE, y, &N_ONE);	}
	inline void axpy(const unsigned n, const float d, const float* const x, float* const y)
	{	saxpy_(&n, &d, x, &N_ONE, y, &N_ONE);	}
        inline void c_axpy(const unsigned n, const float* d, const float* const x, float* const y)
        {	caxpy_(&n, d, x, &N_ONE, y, &N_ONE);	}
        inline void c_axpy(const unsigned n, const double* d, const double* const x, double* const y)
        {	zaxpy_(&n, d, x, &N_ONE, y, &N_ONE);	}



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
        inline void c_gemv(const unsigned m, const unsigned n, float* d, float* A, float *x, float *y)
        {	cgemv_(JOB_STR, &m, &n, d, A, &m, x, &N_ONE, C_ZERO, y, &N_ONE); }
        inline void c_gemv(const unsigned m, const unsigned n, double* d, double* A, double *x, double *y)
        {	zgemv_(JOB_STR, &m, &n, d, A, &m, x, &N_ONE, Z_ZERO, y, &N_ONE); }

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
        inline void c_gemva(const unsigned m, const unsigned n, const float* d, const float* A, const float *x, float *y)
        {	cgemv_(JOB_STR, &m, &n, d, A, &m, x, &N_ONE, C_ONE, y, &N_ONE);	}
        inline void c_gemva(const unsigned m, const unsigned n, const double* d, const double* A, const double *x, double *y)
        {	zgemv_(JOB_STR, &m, &n, d, A, &m, x, &N_ONE, Z_ONE, y, &N_ONE);	}

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
        inline void c_gemtv(const unsigned m, const unsigned n, float* d, float* A, float *x, float *y)
        {	cgemv_(JOB_STR+1, &m, &n, d, A, &m, x, &N_ONE, C_ZERO, y, &N_ONE); }
        inline void c_gemtv(const unsigned m, const unsigned n, double* d, double* A, double *x, double *y)
        {	zgemv_(JOB_STR+1, &m, &n, d, A, &m, x, &N_ONE, Z_ZERO, y, &N_ONE); }
        inline void c_gemhv(const unsigned m, const unsigned n, float* d, float* A, float *x, float *y)
        {	cgemv_(JOB_STR+7, &m, &n, d, A, &m, x, &N_ONE, C_ZERO, y, &N_ONE); } // hermitian transposed
        inline void c_gemhv(const unsigned m, const unsigned n, double* d, double* A, double *x, double *y)
        {	zgemv_(JOB_STR+7, &m, &n, d, A, &m, x, &N_ONE, Z_ZERO, y, &N_ONE); } // hermitian transposed

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
        inline void c_gemtva(const unsigned m, const unsigned n, float* d, float* A, float *x, float *y)
        {	cgemv_(JOB_STR+1, &m, &n, d, A, &m, x, &N_ONE, C_ONE, y, &N_ONE);	}
        inline void c_gemtva(const unsigned m, const unsigned n, double* d, double* A, double *x, double *y)
        {	zgemv_(JOB_STR+1, &m, &n, d, A, &m, x, &N_ONE, Z_ONE, y, &N_ONE); }
        inline void c_gemhva(const unsigned m, const unsigned n, float* d, float* A, float *x, float *y)
        {	cgemv_(JOB_STR+7, &m, &n, d, A, &m, x, &N_ONE, C_ONE, y, &N_ONE);	} // hermitian transposed
        inline void c_gemhva(const unsigned m, const unsigned n, double* d, double* A, double *x, double *y)
        {	zgemv_(JOB_STR+7, &m, &n, d, A, &m, x, &N_ONE, Z_ONE, y, &N_ONE);	} // hermitian transposed




	// C = d A B, A is m x p, B is p x n
	inline void gemm(unsigned m, unsigned p, unsigned n, double d,
									 double* A, unsigned ldA, double* B, unsigned ldB, double* C, unsigned ldC)
	{	dgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB, &D_ZERO, C, &ldC);	}
	inline void gemm(unsigned m, unsigned p, unsigned n, float d,
									 float* A, unsigned ldA, float* B, unsigned ldB, float* C, unsigned ldC)
	{	sgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB, &S_ZERO, C, &ldC);	}
        inline void c_gemm(const unsigned m, const unsigned p, const unsigned n, const float* d,
                                                                                 float* A, const unsigned ldA, float* B, const unsigned ldB, float* C, const unsigned ldC)
        {
            cgemm_(JOB_STR, JOB_STR, &m, &n, &p, d, A, &ldA, B, &ldB, C_ZERO, C, &ldC);	}
        inline void c_gemm(const unsigned m, const unsigned p, const unsigned n, const double* d,
                                                                                 double* A, const unsigned ldA, double* B, const unsigned ldB, double* C, const unsigned ldC)
        {
            zgemm_(JOB_STR, JOB_STR, &m, &n, &p, d, A, &ldA, B, &ldB, Z_ZERO, C, &ldC);	}

	// C += d A B, A is m x p, B is p x n
	inline void gemma(unsigned m, unsigned p, unsigned n, double d,
										double* A, unsigned ldA, double* B, unsigned ldB,	double* C, unsigned ldC)
	{	dgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB, &D_ONE, C, &ldC); }
	inline void gemma(unsigned m, unsigned p, unsigned n, float d,
										float* A, unsigned ldA, float* B, unsigned ldB,	float* C, unsigned ldC)
	{	sgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB, &S_ONE, C, &ldC); }
        inline void c_gemma(unsigned m, unsigned p, unsigned n, float* d,
											float* A, unsigned ldA, float* B, unsigned ldB,	float* C, unsigned ldC)
        {	cgemm_(JOB_STR, JOB_STR, &m, &n, &p, d, A, &ldA, B, &ldB, C_ONE, C, &ldC); }
        inline void c_gemma(unsigned m, unsigned p, unsigned n, double* d,
											double* A, unsigned ldA, double* B, unsigned ldB,	double* C, unsigned ldC)
        {	zgemm_(JOB_STR, JOB_STR, &m, &n, &p, d, A, &ldA, B, &ldB, Z_ONE, C, &ldC); }

	// C = d A^T B, A is m x p, B is m x n
	inline void gemtm(unsigned m, unsigned p, unsigned n, double d,
										double* A, unsigned ldA, double *B, unsigned ldB,	double* C, unsigned ldC)
	{	dgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB, &D_ZERO, C, &ldC);	}
	inline void gemtm(unsigned m, unsigned p, unsigned n, float d,
										float* A, unsigned ldA, float *B, unsigned ldB,	float* C, unsigned ldC)
	{	sgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB, &S_ZERO, C, &ldC);	}
        inline void c_gemtm(unsigned m, unsigned p, unsigned n, float* d,
											float* A, unsigned ldA, float *B, unsigned ldB,	float* C, unsigned ldC)
        {	cgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, d, A, &ldA, B, &ldB, C_ZERO, C, &ldC);	}
        inline void c_gemtm(unsigned m, unsigned p, unsigned n, double* d,
											double* A, unsigned ldA, double *B, unsigned ldB,	double* C, unsigned ldC)
        {	zgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, d, A, &ldA, B, &ldB, Z_ZERO, C, &ldC);	}
        inline void c_gemhm(unsigned m, unsigned p, unsigned n, float* d, // hermitialn transposed
											float* A, unsigned ldA, float *B, unsigned ldB,	float* C, unsigned ldC)
        {	cgemm_(JOB_STR+7, JOB_STR, &p, &n, &m, d, A, &ldA, B, &ldB, C_ZERO, C, &ldC);	}
        inline void c_gemhm(unsigned m, unsigned p, unsigned n, double* d, // hermitian transposed
											double* A, unsigned ldA, double *B, unsigned ldB,	double* C, unsigned ldC)
        {	zgemm_(JOB_STR+7, JOB_STR, &p, &n, &m, d, A, &ldA, B, &ldB, Z_ZERO, C, &ldC);	}

	// C += d A^T B, A is m x p, B is m x n
	inline void gemtma(unsigned m, unsigned p, unsigned n, double d,
										 double* A, unsigned ldA, double *B, unsigned ldB, double* C, unsigned ldC)
	{	dgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB, &D_ONE, C, &ldC); }
	inline void gemtma(unsigned m, unsigned p, unsigned n, float d,
										 float* A, unsigned ldA, float *B, unsigned ldB, float* C, unsigned ldC)
	{	sgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB, &S_ONE, C, &ldC); }
        inline void c_gemtma(unsigned m, unsigned p, unsigned n, float* d,
											 float* A, unsigned ldA, float *B, unsigned ldB, float* C, unsigned ldC)
        {	cgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, d, A, &ldA, B, &ldB, C_ONE, C, &ldC); }
        inline void c_gemtma(unsigned m, unsigned p, unsigned n, double* d,
											 double* A, unsigned ldA, double *B, unsigned ldB, double* C, unsigned ldC)
        {	zgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, d, A, &ldA, B, &ldB, Z_ONE, C, &ldC); }
        inline void c_gemhma(unsigned m, unsigned p, unsigned n, float* d, // hermitian transposed
											 float* A, unsigned ldA, float *B, unsigned ldB, float* C, unsigned ldC)
        {	cgemm_(JOB_STR+7, JOB_STR, &p, &n, &m, d, A, &ldA, B, &ldB, C_ONE, C, &ldC); }
        inline void c_gemhma(unsigned m, unsigned p, unsigned n, double* d, // hermitian transposed
											 double* A, unsigned ldA, double *B, unsigned ldB, double* C, unsigned ldC)
        {	zgemm_(JOB_STR+7, JOB_STR, &p, &n, &m, d, A, &ldA, B, &ldB, Z_ONE, C, &ldC); }





	// singular value decomposition
	inline int gesvd(unsigned m, unsigned n, double* A, double* S, double* VT, unsigned ldVT,
									 unsigned nwk, double* wk)
	{
		int INF;
		dgesvd_(JOB_STR+2, JOB_STR+3, &m, &n, A, &m, S, A, &m, VT, &ldVT,	wk, &nwk, &INF);
		return INF;
	}
	inline int gesvd(unsigned m, unsigned n, float* A, float* S, float* VT, unsigned ldVT,
									 unsigned nwk, float* wk)
	{
		int INF;
		sgesvd_(JOB_STR+2, JOB_STR+3, &m, &n, A, &m, S, A, &m, VT, &ldVT,	wk, &nwk, &INF);
		return INF;
	}

  // singular value decomposition (SO)
  inline int gesvdSO(unsigned m, unsigned n, double* A, double* S, double* U, unsigned ldU,
										 unsigned nwk, double* wk)
  {
    int INF;
    dgesvd_(JOB_STR+3, JOB_STR+2, &m, &n, A, &m, S, U, &m, A, &ldU, wk, &nwk, &INF);
    return INF;
  }
  inline int gesvdSO(unsigned m, unsigned n, float* A, float* S, float* U, unsigned ldU,
										 unsigned nwk, float* wk)
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

