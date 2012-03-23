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

// ==== CMAKE =====
// @FUSE_BLAS
// ================

#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FBlas.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Utils/FAssertable.hpp"
#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Kernels/Chebyshev/FChebTensor.hpp"
#include "../../Src/Kernels/Chebyshev/FChebInterpolator.hpp"

void applyM2M(FReal *const S,	FReal *const w, const unsigned int n,	FReal *const W, const unsigned int N)
{ FBlas::gemtva(n, N, FReal(1.), S,	w, W); }

void applym2m(FReal *const S,	FReal *const w, const unsigned int n,	FReal *const W, const unsigned int N)
{ FBlas::gemtm(n, n, n*n, FReal(1.), S, n, w, n, W, n); }

void applyL2L(FReal *const S,	FReal *const F, const unsigned int n,	FReal *const f, const unsigned int N)
{ FBlas::gemva(n, n, FReal(1.), S, F, f);	}

void applyl2l(FReal *const S,	FReal *const F, const unsigned int n,	FReal *const f, const unsigned int N)
{ FBlas::gemm(n, n, n*n, FReal(1.), S, n, F, n, f, n); }


/**
* In this file we show how to use octree
*/

int main(int argc, char* argv[])
{
	const unsigned int ORDER = 2;
	const unsigned int nnodes = TensorTraits<ORDER>::nnodes;
	FPoint X[nnodes];
	FChebTensor<ORDER>::setRoots(FPoint(0.,0.,0.), FReal(2.), X);
	FPoint x[nnodes];
//	const FPoint cx(-.5, -.5, -.5);
//	const FPoint cx(-.5, -.5,  .5);
//	const FPoint cx(-.5,  .5, -.5);
//	const FPoint cx(-.5,  .5,  .5);
//	const FPoint cx( .5, -.5, -.5);
//	const FPoint cx( .5, -.5,  .5);
//	const FPoint cx( .5,  .5, -.5);
	const FPoint cx( .5,  .5,  .5);
	const FReal  wx(1.);
	FChebTensor<ORDER>::setRoots(cx, wx, x);

	FReal w[nnodes], f[nnodes];
	for (unsigned int n=0; n<nnodes; ++n) {
		w[n] = f[n] = FReal(n);
//		std::cout << w[n] << "\t" << X[n] << "\t" << x[n] << std::endl;
	}



	FReal coords[3][ORDER];
	FChebTensor<ORDER>::setChebyshevRoots(cx, wx, coords);
	
//	for (unsigned int n=0; n<ORDER; ++n) {
//		std::cout << coords[0][n] << "\t"
//							<< coords[1][n] << "\t"
//							<< coords[2][n] << std::endl;
//	}

	FReal S[3][ORDER*ORDER];
	FChebInterpolator<ORDER> Interpolator;
	Interpolator.assembleInterpolator(ORDER, coords[0], S[0]);
	Interpolator.assembleInterpolator(ORDER, coords[1], S[1]);
	Interpolator.assembleInterpolator(ORDER, coords[2], S[2]);


//	std::cout << std::endl;
//	for (unsigned int i=0; i<ORDER; ++i) {
//		for (unsigned int j=0; j<ORDER; ++j)
//			std::cout << S[0][j*ORDER + i] << " ";
//		std::cout << std::endl;
//	}
//
//	std::cout << std::endl;
//	for (unsigned int i=0; i<ORDER; ++i) {
//		for (unsigned int j=0; j<ORDER; ++j)
//			std::cout << S[1][j*ORDER + i] << " ";
//		std::cout << std::endl;
//	}
//
//	std::cout << std::endl;
//	for (unsigned int i=0; i<ORDER; ++i) {
//		for (unsigned int j=0; j<ORDER; ++j)
//			std::cout << S[2][j*ORDER + i] << " ";
//		std::cout << std::endl;
//	}


	FReal Skron[nnodes * nnodes];
	Interpolator.assembleInterpolator(nnodes, x, Skron);
//	std::cout << std::endl;
//	for (unsigned int i=0; i<nnodes; ++i) {
//		for (unsigned int j=0; j<nnodes; ++j)
//			std::cout << Skron[j*nnodes + i] << " ";
//		std::cout << std::endl;
//	}


	FReal W0[nnodes]; FBlas::setzero(nnodes, W0);
	applyM2M(Skron, w, nnodes, W0, nnodes);
	for (unsigned int n=0; n<nnodes; ++n)
		std::cout << w[n] << "\t"	<< W0[n] << std::endl;

	FReal F0[nnodes]; FBlas::setzero(nnodes, F0);
	applyL2L(Skron, f, nnodes, F0, nnodes);
	for (unsigned int n=0; n<nnodes; ++n)
		std::cout << f[n] << "\t"	<< F0[n] << std::endl;

//	std::cout << std::endl;
//	for (unsigned int i=0; i<ORDER; ++i) {
//		for (unsigned int j=0; j<ORDER*ORDER; ++j)
//			std::cout << w[j*ORDER + i] << " ";
//		std::cout << std::endl;
//	}




	unsigned int perm[3][nnodes];
	for (unsigned int i=0; i<ORDER; ++i) {
		for (unsigned int j=0; j<ORDER; ++j) {
			for (unsigned int k=0; k<ORDER; ++k) {
				const unsigned int index = k*ORDER*ORDER + j*ORDER + i;
				perm[0][index] = k*ORDER*ORDER + j*ORDER + i;
				perm[1][index] = i*ORDER*ORDER + k*ORDER + j;
				perm[2][index] = j*ORDER*ORDER + i*ORDER + k;
			}
		}
	}

//	for (unsigned int n=0; n<nnodes; ++n)
//		std::cout << perm[0][n] << "\t" << perm[1][n] << "\t" << perm[2][n] << std::endl;


	FReal W[nnodes]; FBlas::setzero(nnodes, W);
	applym2m(S[0], w, ORDER, W, ORDER);
//	std::cout << std::endl;
//	for (unsigned int i=0; i<ORDER; ++i) {
//		for (unsigned int j=0; j<ORDER*ORDER; ++j)
//			std::cout << W[j*ORDER + i] << " ";
//		std::cout << std::endl;
//	}
//	std::cout << std::endl;

	for (unsigned int n=0; n<nnodes; ++n)
		w[n] = W[perm[1][n]];
	applym2m(S[2], w, ORDER, W, ORDER);
//	for (unsigned int i=0; i<ORDER; ++i) {
//		for (unsigned int j=0; j<ORDER*ORDER; ++j)
//			std::cout << W[j*ORDER + i] << " ";
//		std::cout << std::endl;
//	}
//	std::cout << std::endl;

	for (unsigned int n=0; n<nnodes; ++n)
		w[perm[1][n]] = W[perm[2][n]];
	applym2m(S[1], w, ORDER, W, ORDER);
//	for (unsigned int i=0; i<ORDER; ++i) {
//		for (unsigned int j=0; j<ORDER*ORDER; ++j)
//			std::cout << W[j*ORDER + i] << " ";
//		std::cout << std::endl;
//	}
//	std::cout << std::endl;

	for (unsigned int n=0; n<nnodes; ++n)
		w[perm[2][n]] = W[n];
//	for (unsigned int i=0; i<ORDER; ++i) {
//		for (unsigned int j=0; j<ORDER*ORDER; ++j)
//			std::cout << w[j*ORDER + i] << " ";
//		std::cout << std::endl;
//	}
//	std::cout << std::endl;

	std::cout << std::endl;
	FReal error(0.);
	for (unsigned int n=0; n<nnodes; ++n) {
		std::cout << n << "\t" << w[n] << " - " << W0[n] << " = " << w[n]-W0[n] << std::endl;
		error += w[n] - W0[n];
	}

	std::cout << "\nERROR = " << error << std::endl << std::endl;



	FReal F[nnodes]; FBlas::setzero(nnodes, F);
	applyl2l(S[0], f, ORDER, F, ORDER);

	for (unsigned int n=0; n<nnodes; ++n)
		f[n] = F[perm[1][n]];
	applyl2l(S[2], f, ORDER, F, ORDER);

	for (unsigned int n=0; n<nnodes; ++n)
		f[perm[1][n]] = F[perm[2][n]];
	applyl2l(S[1], f, ORDER, F, ORDER);

	for (unsigned int n=0; n<nnodes; ++n)
		f[perm[2][n]] = F[n];

	std::cout << std::endl;
	FReal ferror(0.);
	for (unsigned int n=0; n<nnodes; ++n) {
		std::cout << n << "\t" << f[n] << " - " << F0[n] << " = " << f[n]-F0[n] << std::endl;
		ferror += f[n] - F0[n];
	}

	std::cout << "\nERROR = " << ferror << std::endl << std::endl;

	return 0;
}


// [--END--]
