// [--License--]

#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../Src/Utils/FTic.hpp"
#include "../Src/Utils/FMath.hpp"

#include "../Src/Containers/FVector.hpp"

#include "../Src/Utils/FAssertable.hpp"
#include "../Src/Utils/F3DPosition.hpp"

#include "../Src/Chebyshev/FChebParticle.hpp"
#include "../Src/Chebyshev/FChebLeaf.hpp"
#include "../Src/Chebyshev/FChebInterpolator.hpp"
#include "../Src/Chebyshev/FChebM2LHandler.hpp"
#include "../Src/Chebyshev/FChebMatrixKernel.hpp"




const FReal computeL2norm(unsigned int N, FReal *const u, FReal *const v)
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



const FReal computeINFnorm(unsigned int N, FReal *const u, FReal *const v)
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
	typedef FChebParticle ParticleClass;
	typedef FVector<FChebParticle> ContainerClass;
	typedef FChebLeaf<ParticleClass,ContainerClass> LeafClass;
	typedef FChebMatrixKernelR MatrixKernelClass;

	///////////////////////What we do/////////////////////////////
	std::cout << "\nTask: Compute interactions between source particles in leaf Y and target\n";
	std::cout << " particles in leaf X. Compare the fast summation K ~ Sx K Sy' with the\n";
	std::cout << " direct computation.\n" << std::endl;
	//////////////////////////////////////////////////////////////

	MatrixKernelClass MatrixKernel;
	const FReal FRandMax = FReal(RAND_MAX);
	FTic time;

	
	// Leaf size
	FReal width(FReal(rand()) / FRandMax * FReal(10.));

	////////////////////////////////////////////////////////////////////
	LeafClass X;
	F3DPosition cx(0., 0., 0.);
	const long M = 10000;
	std::cout << "Fill the leaf X of width " << width
						<< " centered at cx=[" << cx.getX() << "," << cx.getY() << "," << cx.getZ()
						<< "] with M=" << M << " target particles" << std::endl;
	{
		FChebParticle particle;
		for(long i=0; i<M; ++i){
			FReal x = (FReal(rand())/FRandMax - FReal(.5)) * width + cx.getX();
			FReal y = (FReal(rand())/FRandMax - FReal(.5)) * width + cx.getY();
			FReal z = (FReal(rand())/FRandMax - FReal(.5)) * width + cx.getZ();
			particle.setPosition(x, y, z);
			X.push(particle);
		}
	}


	////////////////////////////////////////////////////////////////////
	LeafClass Y;
	F3DPosition cy(FReal(2.)*width, 0., 0.);
	const long N = 10000;
	std::cout << "Fill the leaf Y of width " << width
						<< " centered at cy=[" << cy.getX() << "," << cy.getY() << "," << cy.getZ()
						<< "] with N=" << N << " target particles" << std::endl;
	{
		FChebParticle particle;
		for(long i=0; i<N; ++i){
			FReal x = (FReal(rand())/FRandMax - FReal(.5)) * width + cy.getX();
			FReal y = (FReal(rand())/FRandMax - FReal(.5)) * width + cy.getY();
			FReal z = (FReal(rand())/FRandMax - FReal(.5)) * width + cy.getZ();
			particle.setPosition(x, y, z);
			particle.setPhysicalValue(FReal(rand())/FRandMax);
			Y.push(particle);
		}
	}



	////////////////////////////////////////////////////////////////////
	// approximative computation
	const unsigned int ORDER = 5;
	const FReal epsilon = FReal(atof(argv[1]));
	const unsigned int nnodes = TensorTraits<ORDER>::nnodes;
	typedef FChebInterpolator<ORDER> InterpolatorClass;
	InterpolatorClass S;

	// set compressed M2L operators
	FChebM2LHandler<ORDER,MatrixKernelClass> M2L(epsilon); 
	//M2L.ComputeAndCompressAndSet();  // <- precompute M2L operators
	M2L.ReadFromBinaryFileAndSet();    // <- read precomputed M2L operators from binary file

	std::cout << "\nCompute interactions approximatively, ACC = (" << ORDER << ", " << epsilon << ") ..." << std::endl;
	time.tic();
	FReal W[nnodes * 2]; // multipole expansion
	FReal F[nnodes * 2]; // local expansion
	for (unsigned int n=0; n<nnodes*2; ++n)
		W[n] = F[n] = FReal(0.);

	// Anterpolate: W_n = \sum_j^N S(y_j,\bar y_n) * w_j
	S.applyP2M(cy, width, W, Y.getSrc()); // the multipole expansions are set to 0 in S.applyP2M

	
	// M2L (compressed)
	const int diffx = int((cy.getX()-cx.getX()) / width);
	const int diffy = int((cy.getY()-cx.getY()) / width);
	const int diffz = int((cy.getZ()-cx.getZ()) / width);
	const int transfer[3] = {diffx, diffy, diffz};
	M2L.applyB(W, W+nnodes);
	M2L.applyC(transfer, width, W+nnodes, F+nnodes);
	M2L.applyU(F+nnodes, F);
	
	//for (unsigned int n=0; n<nnodes; ++n) F[n] *= MatrixKernel.getScaleFactor(width);

	/*
	// M2L (direct)
	F3DPosition rootsX[nnodes], rootsY[nnodes];
	FChebTensor<ORDER>::setRoots(cx, width, rootsX);
	FChebTensor<ORDER>::setRoots(cy, width, rootsY);
	for (unsigned int i=0; i<nnodes; ++i) {
		F[i] = FReal(0.);
		for (unsigned int j=0; j<nnodes; ++j)
			F[i] += MatrixKernel.evaluate(rootsX[i], rootsY[j]) * W[j];
	}
	*/

	// Interpolate f_i = \sum_m^L S(x_i,\bar x_m) * F_m
	S.applyL2P(cx, width, F, X.getTargets());

	time.tac();
	std::cout << "Done in " << time.elapsed() << "sec." << std::endl;
	// -----------------------------------------------------

	////////////////////////////////////////////////////////////////////
	// direct computation
	std::cout << "Compute interactions directly ..." << std::endl;
	time.tic();

	FReal* approx_f = new FReal[M];
	FReal*        f = new FReal[M];
	for (unsigned int i=0; i<M; ++i) f[i] = FReal(0.);
	ContainerClass::ConstBasicIterator iterY(*(Y.getSrc()));
	while(iterY.hasNotFinished()){
		const F3DPosition& y = iterY.data().getPosition();
		const FReal        w = iterY.data().getPhysicalValue();
		unsigned int counter = 0;
		ContainerClass::ConstBasicIterator iterX(*(X.getSrc()));
		while(iterX.hasNotFinished()){
			const F3DPosition& x = iterX.data().getPosition();
			f[counter++] += MatrixKernel.evaluate(x,y) * w;
			iterX.gotoNext();
		}
		iterY.gotoNext();
	}
	time.tac();
	std::cout << "Done in " << time.elapsed() << "sec." << std::endl;


	////////////////////////////////////////////////////////////////////
	ContainerClass::ConstBasicIterator iterX(*(X.getSrc()));
	unsigned int counter = 0;
	while(iterX.hasNotFinished()) {
		approx_f[counter++] = iterX.data().getPotential();
		iterX.gotoNext();
	}

	std::cout << "\nRelative L2 error  = " << computeL2norm( M, f, approx_f) << std::endl;
	std::cout << "Relative Lmax error = "  << computeINFnorm(M, f, approx_f) << "\n" << std::endl;

	// free memory
	delete [] approx_f;
	delete [] f;


	return 0;
}


// [--END--]
