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

#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Kernels/Chebyshev/FChebMatrixKernel.hpp"
#include "../../Src/Kernels/Chebyshev/FChebRoots.hpp"
#include "../../Src/Kernels/Chebyshev/FChebTensor.hpp"
#include "../../Src/Kernels/Chebyshev/FChebSymM2LHandler.hpp"




FReal computeL2norm(unsigned int N, FReal *const u, FReal *const v)
{
	FReal      dot = FReal(0.);
	FReal diff_dot = FReal(0.);
	for (unsigned int i=0; i<N; ++i) {
		FReal w = v[i] - u[i];
		diff_dot += w    * w;
		dot      += u[i] * u[i];
	}
	return FMath::Sqrt(diff_dot / dot);
}



FReal computeINFnorm(unsigned int N, FReal *const u, FReal *const v)
{
	FReal      max = FReal(0.);
	FReal diff_max = FReal(0.);
	for (unsigned int n=0; n<N; ++n) {
		if (     max<std::abs(u[n]))           max = std::abs(u[n]);
		if (diff_max<std::abs(u[n]-v[n])) diff_max = std::abs(u[n]-v[n]);
	}
	return diff_max / max;
}




/**
* In this file we show how to use octree
*/

int main(int argc, char* argv[])
{
	typedef FChebMatrixKernelR MatrixKernelClass;
	MatrixKernelClass MatrixKernel;

	const unsigned int ORDER = 10;
	const FReal epsilon = 1e-10;
	const unsigned int nnodes = TensorTraits<ORDER>::nnodes;
	
	// cell size
	FReal width = FReal(2.);

	// centers of cells X and Y
	FPoint cx(0., 0., 0.);
	FPoint cy(FReal(2.)*width, 0., 0.);

	// compute Cheb points in cells X and Y
	FPoint rootsX[nnodes], rootsY[nnodes];
	FChebTensor<ORDER>::setRoots(cx, width, rootsX);
	FChebTensor<ORDER>::setRoots(cy, width, rootsY);

	
	// initialize timer
	FTic time;


	// fully pivoted ACA ///////////////////////////
	std::cout << "Fully pivoted ACA of Acc(" << ORDER << ", " << epsilon << ")" << std::endl;
	std::cout << "|- Assembling K" << std::flush;
	FReal *const K = new FReal [nnodes * nnodes];
	time.tic();
	EntryComputer<MatrixKernelClass> Computer(nnodes, rootsX, nnodes, rootsY);
	Computer(0, nnodes, 0, nnodes, K);
	std::cout << ", finished in " << time.tacAndElapsed() << "s." << std::endl;
	
	// generate right hand side vector
	FReal *const w = new FReal [nnodes];
	const FReal FRandMax = FReal(RAND_MAX);
	for (unsigned int j=0; j<nnodes; ++j)
		w[j] = FReal(rand())/FRandMax;
	
	// compute f0 = Kw
	FReal *const f0 = new FReal [nnodes];
	FBlas::gemv(nnodes, nnodes, FReal(1.), K, w, f0);
	
	// call fACA /////////////////////////////////
	std::cout << "|- Computing fACA" << std::flush;
	time.tic();
	FReal *U, *V;
	unsigned int k;
	fACA(K, nnodes, nnodes, epsilon, U, V, k);
	std::cout << " (k = " << k << "), finished in " << time.tacAndElapsed() << "s." << std::endl;
	delete [] K;
	
	// compute f1 = UV'w
	FReal *const f1 = new FReal [nnodes];
	FReal *const c1 = new FReal [k];
	FBlas::gemtv(nnodes, k, FReal(1.), V, w,  c1);
	FBlas::gemv( nnodes, k, FReal(1.), U, c1, f1);
	delete [] c1;
	std::cout << "   |- L2 error = " << computeL2norm(nnodes, f0, f1) << std::endl;
	delete [] U;
	delete [] V;

	// call pACA ///////////////////////////////////
	std::cout << "|- Computing pACA" << std::flush;
	time.tic();
	pACA(Computer, nnodes, nnodes, epsilon, U, V, k);
	std::cout << " (k = " << k << "), finished in " << time.tacAndElapsed() << "s." << std::endl;
	// compute f1 = UV'w
	FReal *const f2 = new FReal [nnodes];
	FReal *const c2 = new FReal [k];
	FBlas::gemtv(nnodes, k, FReal(1.), V, w,  c2);
	FBlas::gemv( nnodes, k, FReal(1.), U, c2, f2);
	delete [] c2;
	std::cout << "   |- L2 error = " << computeL2norm(nnodes, f0, f2) << std::endl;
	delete [] U;
	delete [] V;
	

	delete [] w;
	delete [] f0;
	delete [] f1;
	delete [] f2;

	return 0;
}


// [--END--]
