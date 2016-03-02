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
// Keep in private GIT
// @SCALFMM_PRIVATE


#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FMath.hpp"

#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "../../Src/Kernels/Chebyshev/FChebRoots.hpp"
#include "../../Src/Kernels/Chebyshev/FChebTensor.hpp"
#include "../../Src/Kernels/Chebyshev/FChebSymM2LHandler.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

#include "../../Src/Utils/FAca.hpp"

/**
* In this file we show how to use octree
*/

int main(int argc, char* argv[])
{
    FHelpDescribeAndExit(argc, argv, "Look to the code.");

    typedef double FReal;

    typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
	MatrixKernelClass MatrixKernel;

	const unsigned int ORDER = 9;
	const FReal epsilon = 1e-9;
	const unsigned int nnodes = TensorTraits<ORDER>::nnodes;

	// initialize timer
	FTic time;

	/*	
	//// width of cell X and Y
	//FReal wx = FReal(2.);
	//FReal wy = FReal(2.);
	//// centers of cells X and Y
    //FPoint<FReal> cx( 0., 0., 0.);
    //FPoint<FReal> cy( 4., 0., 0.);

	// width of cell X and Y
	FReal wx = FReal(2.);
	FReal wy = FReal(4.);
	// centers of cells X and Y
    FPoint<FReal> cx( 1.,-1.,-1.);
    //FPoint<FReal> cy(-4., 0., 0.);
    FPoint<FReal> cy(-4., 4., 4.);


	std::cout << "[cx = " << cx << ", wx = " << wx
						<< "] [cy = " << cy << ", wy = " << wy
                        << "] -> dist = " << FPoint<FReal>(cx-cy).norm() << std::endl;

	// compute Cheb points in cells X and Y
    FPoint<FReal> rootsX[nnodes], rootsY[nnodes];
    FChebTensor<FReal,ORDER>::setRoots(cx, wx, rootsX);
    FChebTensor<FReal,ORDER>::setRoots(cy, wy, rootsY);

	


	// fully pivoted ACA ///////////////////////////
	std::cout << "Fully pivoted ACA of Acc(" << ORDER << ", " << epsilon << ")" << std::endl;
	std::cout << "|- Assembling K" << std::flush;
	FReal *const K = new FReal [nnodes * nnodes];
	time.tic();
    EntryComputer<FReal,MatrixKernelClass> Computer(nnodes, rootsX, nnodes, rootsY);
	Computer(0, nnodes, 0, nnodes, K);
	std::cout << ", finished in " << time.tacAndElapsed() << "s." << std::endl;
	
	// generate right hand side vector
    FReal *const w = new FReal [nnodes];
	for (unsigned int j=0; j<nnodes; ++j)
        w[j] = FReal(drand48());
	
	// compute f0 = Kw
	FReal *const f0 = new FReal [nnodes];
	FBlas::gemv(nnodes, nnodes, FReal(1.), K, w, f0);

	// allocate UV and k
	FReal *U, *V;
	unsigned int k;


	// call fACA /////////////////////////////////
	std::cout << "|- Computing fACA" << std::flush;
	time.tic();
	fACA(K, nnodes, nnodes, epsilon, U, V, k);
	std::cout << " (k = " << k << "), finished in " << time.tacAndElapsed() << "s." << std::endl;
	delete [] K;
	
	// compute f1 = UV'w
	FReal *const f1 = new FReal [nnodes];
	FReal *const c1 = new FReal [k];
	FBlas::gemtv(nnodes, k, FReal(1.), V, w,  c1);
	FBlas::gemv( nnodes, k, FReal(1.), U, c1, f1);
	delete [] c1;
	std::cout << "|  |- L2 error = " << computeL2norm(nnodes, f0, f1) << std::endl;
	delete [] U;
	delete [] V;
	delete [] f1;


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
	std::cout << "|  |- L2 error = " << computeL2norm(nnodes, f0, f2) << std::endl;
	delete [] U;
	delete [] V;
	

	delete [] w;
	delete [] f0;
	delete [] f2;
*/

	// test different levels of multipole expansions
	const FReal pw = FReal(4.);
	const FReal cw = pw / FReal(2.);

    const FPoint<FReal> cpx( 0., 0., 0.);
    const FPoint<FReal> ccx(-1.,-1.,-1.);


	// compute Cheb points in cells X and Y
    FPoint<FReal> rootsX[nnodes], rootsY[nnodes];
    FChebTensor<FReal,ORDER>::setRoots(ccx, cw, rootsX);

	unsigned int all_counter = 0;
	unsigned int ccounter = 0;
	unsigned int pcounter = 0;

	unsigned int overall_rank = 0;
	unsigned int all_overall_rank = 0;

	// allocate UV and k
	FReal *U, *V;
	unsigned int rank;


	for (int i=-1; i<=1; ++i)
		for (int j=-1; j<=1; ++j)
			for (int k=-1; k<=1; ++k) {
                const FPoint<FReal> cpy(FReal(i)*pw, FReal(j)*pw, FReal(k)*pw);
				
				// exclude target cell itself
				if ((FMath::Abs(i) + FMath::Abs(j) + FMath::Abs(k)) != 0) {
					
					//////////////////////////////////////////////////////////////////////
					// use all multipole expansions form child level
					for (int ci=-1; ci<=1; ci+=2)
						for (int cj=-1; cj<=1; cj+=2)
							for (int ck=-1; ck<=1; ck+=2) {
                                const FPoint<FReal> ccy(cpy.getX() + FReal(ci)*cw/2., cpy.getY() + FReal(cj)*cw/2., cpy.getZ() + FReal(ck)*cw/2.);
								
								// exclude near-field 
                                if (FPoint<FReal>(ccx-ccy).norm() > FReal(3.5)) {
                                    FChebTensor<FReal,ORDER>::setRoots(ccy, cw, rootsY);
                                    EntryComputer<FReal,MatrixKernelClass> Computer(&MatrixKernel,nnodes, rootsX, nnodes, rootsY);
                                    FAca::pACA(Computer, nnodes, nnodes, epsilon, U, V, rank);
									//std::cout << "- Compress " << ccy << "\tof width " << cw << " to rank " << rank << std::endl;
									
									all_counter++;
									all_overall_rank += rank;
									
									delete [] U;
									delete [] V;
								}
							}
					
					//////////////////////////////////////////////////////////////////////
					// exclude source cells whose multipole expansion are taken from parent level
					if (cpy.getX()*ccx.getX() >= 0 &&	cpy.getY()*ccx.getY() >= 0 &&	cpy.getZ()*ccx.getZ() >= 0) {
						
						for (int ci=-1; ci<=1; ci+=2)
							for (int cj=-1; cj<=1; cj+=2)
								for (int ck=-1; ck<=1; ck+=2) {
                                    const FPoint<FReal> ccy(cpy.getX() + FReal(ci)*cw/2., cpy.getY() + FReal(cj)*cw/2., cpy.getZ() + FReal(ck)*cw/2.);
									
									// exclude near-field 
                                    if (FPoint<FReal>(ccx-ccy).norm() > FReal(3.5)) {
                                        FChebTensor<FReal,ORDER>::setRoots(ccy, cw, rootsY);
                                        EntryComputer<FReal,MatrixKernelClass> Computer(&MatrixKernel,nnodes, rootsX, nnodes, rootsY);
                                        FAca::pACA(Computer, nnodes, nnodes, epsilon, U, V, rank);
										//std::cout << "- Compress " << ccy << "\tof width " << cw << " to rank " << rank << std::endl;
										
										ccounter++;
										overall_rank += rank;
										
										delete [] U;
										delete [] V;
									}
								}
						
					}
					// remaining far-field: source cells whose multipole expansion is taken from parent level
					else {
                        FChebTensor<FReal,ORDER>::setRoots(cpy, pw, rootsY);
                        EntryComputer<FReal,MatrixKernelClass> Computer(&MatrixKernel,nnodes, rootsX, nnodes, rootsY);
                        FAca::pACA(Computer, nnodes, nnodes, epsilon, U, V, rank);
						//std::cout << "- Compress " << cpy << "\tof width " << pw << " to rank " << rank << std::endl;

						pcounter++;
						overall_rank += rank;

						delete [] U;
						delete [] V;
					}
					
				}
			}

	std::cout << "-- Overall low rank " << overall_rank
						<< ", child cells = " << ccounter
						<< ", parent cells = " << pcounter << std::endl;

	std::cout << "-- Overall low rank " << all_overall_rank
						<< ", child cells = " << all_counter << std::endl;
	
	std::cout << "--- Ratio = " << FReal(overall_rank) / FReal(all_overall_rank) << std::endl;

	return 0;
}


// [--END--]
