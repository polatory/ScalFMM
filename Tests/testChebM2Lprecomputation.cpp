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


// system includes
#include <iostream>
#include <stdexcept>


#include "../Src/Utils/FGlobal.hpp"
#include "../Src/Utils/F3DPosition.hpp"
#include "../Src/Utils/FMath.hpp"
#include "../Src/Utils/FTic.hpp"


#include "../Src/Chebyshev/FChebTensor.hpp"
#include "../Src/Chebyshev/FChebM2LHandler.hpp"
#include "../Src/Chebyshev/FChebMatrixKernel.hpp"






int main(int argc, char* argv[])
{
	// start timer /////////////////////////////////
	FTic time;
	
	// define set matrix kernel
	typedef FChebMatrixKernelR MatrixKernelClass;
	MatrixKernelClass MatrixKernel;

	// constants
  const FReal epsilon     = FReal(atof(argv[1]));
  const unsigned int order = 4;
	
	// number of interpolation points per cell
	const unsigned int nnodes = TensorTraits<order>::nnodes;

	// interpolation points of source (Y) and target (X) cell
	F3DPosition X[nnodes], Y[nnodes];
	// set roots of target cell (X)
	FChebTensor<order>::setRoots(F3DPosition(0.,0.,0.), FReal(2.), X);
	

	/*
	// allocate memory
	FReal *Qu, *C, *Qb;
	Qu = Qb = C = NULL;
	unsigned int ninteractions = 0;

	////////////////////////////////////////////////
	std::cout << "\nAssembly of 316 times "
						<< nnodes << "x" << nnodes << " M2L operators";
	time.tic();
	// compute 316 m2l operators
	ninteractions = 316;
	C = new FReal [nnodes*nnodes * ninteractions];
	unsigned int counter = 0;
	for (int i=-3; i<=3; ++i) {
		for (int j=-3; j<=3; ++j) {
			for (int k=-3; k<=3; ++k) {
				if (abs(i)>1 || abs(j)>1 || abs(k)>1) {
					// set roots of source cell (Y)
					const F3DPosition cy(FReal(2.*i), FReal(2.*j), FReal(2.*k));
					FChebTensor<order>::setRoots(cy, FReal(2.), Y);
					// evaluate m2l operator
					for (unsigned int n=0; n<nnodes; ++n)
						for (unsigned int m=0; m<nnodes; ++m)
							C[counter*nnodes*nnodes + n*nnodes + m] = MatrixKernel.evaluate(X[m], Y[n]);
					// increment interaction counter
					counter++;
				}
			}
		}
	}
	if (counter != 316)
		std::runtime_error("Number of interactions must correspond to 316");
	std::cout << " took " << time.tacAndElapsed() << " sec." << std::endl;
	////////////////////////////////////////////////

	////////////////////////////////////////////////
	std::cout << "\nSVD compression ";
	time.tic();
	const unsigned int rank = Compress<order>(epsilon, ninteractions, Qu, C, Qb);
	std::cout << "to low rank = " << rank << " (eps = " << epsilon
						<< ") took " << time.tacAndElapsed() << " sec." << std::endl;
	////////////////////////////////////////////////

	// free memory
	if (C  != NULL) delete [] C;
	if (Qu != NULL) delete [] Qu;
	if (Qb != NULL) delete [] Qb;
	*/

	////////////////////////////////////////////////
	// allocate memory
	FReal *Qu1, *C1, *Qb1;
	Qu1 = Qb1 = C1 = NULL;
	////////////////////////////////////////////////
	std::cout << "\nAssembly of an " << nnodes << "x" << nnodes << " M2L operator";
	time.tic();
	// compute 316 m2l operators
	C1 = new FReal [nnodes*nnodes];
	const unsigned int i = 2;
	const unsigned int j = 0;
	const unsigned int k = 0;
	const F3DPosition cy(FReal(2.*i), FReal(2.*j), FReal(2.*k));
	FChebTensor<order>::setRoots(cy, FReal(2.), Y);
	// evaluate m2l operator
	for (unsigned int n=0; n<nnodes; ++n)
		for (unsigned int m=0; m<nnodes; ++m)
			C1[n*nnodes + m] = MatrixKernel.evaluate(X[m], Y[n]);
	std::cout << " took " << time.tacAndElapsed() << " sec." << std::endl;
	////////////////////////////////////////////////

	////////////////////////////////////////////////
	// get a copy C2 of the M2L operator C1
	FReal *Qu2, *C2, *Qb2;
	Qu2 = Qb2 = C2 = NULL;
	C2 = new FReal [nnodes * nnodes];
	FBlas::copy(nnodes*nnodes, C1, C2);
	// Omega_x^{1/2} C2 Omega_y^{1/2}
	FReal weights[nnodes];
	FChebTensor<order>::setRootOfWeights(weights);
	for (unsigned int n=0; n<nnodes; ++n) {
		FBlas::scal(nnodes, weights[n], C2+n, nnodes); // scale rows
		FBlas::scal(nnodes, weights[n], C2+n*nnodes);  // scale cols
	}
	////////////////////////////////////////////////

	////////////////////////////////////////////////
	std::cout << "\nSVD compression of K ";
	time.tic();
	const unsigned int rank1 = Compress<order>(epsilon, 1, Qu1, C1, Qb1);
	std::cout << "to low rank = " << rank1 << " (eps = " << epsilon
						<< ") took " << time.tacAndElapsed() << " sec." << std::endl;
	////////////////////////////////////////////////

	////////////////////////////////////////////////
	std::cout << "SVD compression of Omega_x^{1/2} K Omega_y^{1/2} ";
	time.tic();
	const unsigned int rank2 = Compress<order>(epsilon, 1, Qu2, C2, Qb2);
	std::cout << "to low rank = " << rank2 << " (eps = " << epsilon
						<< ") took " << time.tacAndElapsed() << " sec." << std::endl;
	////////////////////////////////////////////////

	// free memory
	if (C1  != NULL) delete [] C1;
	if (Qu1 != NULL) delete [] Qu1;
	if (Qb1 != NULL) delete [] Qb1;

	if (C2  != NULL) delete [] C2;
	if (Qu2 != NULL) delete [] Qu2;
	if (Qb2 != NULL) delete [] Qb2;


	return 0;
} 	
