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

#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "../Src/Chebyshev/FChebParticle.hpp"
#include "../Src/Chebyshev/FChebLeaf.hpp"
#include "../Src/Chebyshev/FChebCell.hpp"
#include "../Src/Chebyshev/FChebMatrixKernel.hpp"
#include "../Src/Chebyshev/FChebKernels.hpp"

//#include "../Src/Utils/FTic.hpp"
#include "../Src/Utils/FParameters.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Core/FFmmAlgorithm.hpp"
#include "../Src/Core/FFmmAlgorithmThread.hpp"


/** This program show an example of use of
  * the fmm basic algo
  * it also check that each particles is little or longer
  * related that each other
  */



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




// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
	const unsigned int ORDER = 3;
	const FReal epsilon              = FParameters::getValue(argc, argv, "-eps", FReal(1e-3));
	const long NbPart                = FParameters::getValue(argc, argv, "-num", 100000);
	const unsigned int TreeHeight    = FParameters::getValue(argc, argv, "-h", 5);
	const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, "-sh", 2);

	const FReal Width = 10.;
	// init random fun
	const FReal FRandMax = FReal(RAND_MAX) / Width;
	srand( static_cast<unsigned int>(time(NULL)) );
	// init timer
	FTic time;

	// typedefs
	typedef FChebParticle ParticleClass;
	typedef FVector<FChebParticle> ContainerClass;
	typedef FChebLeaf<ParticleClass,ContainerClass> LeafClass;
	typedef FChebMatrixKernelRR MatrixKernelClass;
	typedef FChebCell<ORDER> CellClass;
	typedef FOctree<ParticleClass,CellClass,ContainerClass,LeafClass> OctreeClass;
	typedef FChebKernels<ParticleClass,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
	//typedef FFmmAlgorithm<OctreeClass,ParticleClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
	typedef FFmmAlgorithmThread<OctreeClass,ParticleClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;

	// What we do //////////////////////////////////////////////////////
	std::cout << ">> Testing the Chebyshev interpolation base FMM algorithm.\n";

	
	const F3DPosition BoxCenter(.5*Width, .5*Width, .5*Width);
	OctreeClass tree(TreeHeight, SubTreeHeight, Width, BoxCenter);
	
	// -----------------------------------------------------
	std::cout << "Creating and inserting " << NbPart << " particles ..." << std::endl;
	time.tic();
	ParticleClass particle;
	for(long i=0; i<NbPart; ++i) {
		particle.setPosition(FReal(rand())/FRandMax, FReal(rand())/FRandMax, FReal(rand())/FRandMax);
		particle.setPhysicalValue(FReal(rand())/FReal(RAND_MAX));
		tree.insert(particle);
	}
	std::cout << "Done  " << "(" << time.tacAndElapsed() << ")." << std::endl;
	// -----------------------------------------------------

	
	// -----------------------------------------------------
	std::cout << "\nChebyshev FMM ... " << std::endl;
	time.tic();
	KernelClass kernels(TreeHeight, BoxCenter, Width, epsilon);
	FmmClass algorithm(&tree,&kernels);
	algorithm.execute();
	std::cout << "completed in " << time.tacAndElapsed() << "sec." << std::endl;
	// -----------------------------------------------------
	

	// -----------------------------------------------------
	// find first non empty leaf cell 
	OctreeClass::Iterator iLeafs(&tree);
	iLeafs.gotoBottomLeft();
	
	const ContainerClass *const Targets = iLeafs.getCurrentListTargets();
	const unsigned int NumTargets = Targets->getSize();
	FReal* Potential = new FReal [NumTargets];
	FBlas::scal(NumTargets, FReal(0.), Potential);
	
	std::cout << "\nDirect computation of " << NumTargets << " target particles ..." << std::endl;
	const MatrixKernelClass MatrixKernel;
	do {
			const ContainerClass *const Sources = iLeafs.getCurrentListSrc();
			unsigned int counter = 0;
			ContainerClass::ConstBasicIterator iTarget(*Targets);
			while(iTarget.hasNotFinished()) {
				ContainerClass::ConstBasicIterator iSource(*Sources);
				while(iSource.hasNotFinished()) {
					if (&iTarget.data() != &iSource.data())
						Potential[counter] += MatrixKernel.evaluate(iTarget.data().getPosition(),
																												iSource.data().getPosition())
							* iSource.data().getPhysicalValue();
					iSource.gotoNext();
				}
				counter++;
				iTarget.gotoNext();
			}
	} while(iLeafs.moveRight());

	
	FReal* ApproxPotential = new FReal [NumTargets];
	unsigned int counter = 0;
	ContainerClass::ConstBasicIterator iTarget(*Targets);
	while(iTarget.hasNotFinished()) {
		ApproxPotential[counter] = iTarget.data().getPotential();
		//std::cout << Potential[counter] << " - " << ApproxPotential[counter] << " = "
		//					<< Potential[counter]-ApproxPotential[counter] << "\t rel error = "
		//					<< (Potential[counter]-ApproxPotential[counter]) / Potential[counter]
		//					<< std::endl;
		counter++;
		iTarget.gotoNext();
	}
	
	std::cout << "\nRelative L2 error  = " << computeL2norm( NumTargets, Potential, ApproxPotential)
						<< std::endl;
	std::cout << "Relative Lmax error = "  << computeINFnorm(NumTargets, Potential, ApproxPotential)
						<< "\n" << std::endl;

	// free memory
	delete [] Potential;
	delete [] ApproxPotential;
	

	/*
	// Check if particles are strictly within its containing cells
	const FReal BoxWidthLeaf = BoxWidth / FReal(FMath::pow(2, TreeHeight-1));
	OctreeClass::Iterator octreeIterator(&tree);
	octreeIterator.gotoBottomLeft();
	do{
		const CellClass *const LeafCell = octreeIterator.getCurrentCell();
		const F3DPosition& LeafCellCenter = LeafCell -> getPosition();
		const ContainerClass *const Particles = octreeIterator.getCurrentListSrc();
		ContainerClass::ConstBasicIterator particleIterator(*Particles);
		while(particleIterator.hasNotFinished()) {
			const F3DPosition distance(LeafCellCenter-particleIterator.data().getPosition());
			std::cout << "center - particle = " << distance << " < " << BoxWidthLeaf/FReal(2.) << std::endl;
			if (std::abs(distance.getX())>BoxWidthLeaf/FReal(2.) ||
					std::abs(distance.getY())>BoxWidthLeaf/FReal(2.) ||
					std::abs(distance.getZ())>BoxWidthLeaf/FReal(2.)) {
				std::cout << "stop" << std::endl;
				exit(-1);
			}
			particleIterator.gotoNext();
		}
	} while(octreeIterator.moveRight());
	*/


    return 0;
}



