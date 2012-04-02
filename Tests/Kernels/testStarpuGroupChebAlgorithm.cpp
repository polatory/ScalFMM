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
// @FUSE_STARPU
// ================

#include <iostream>

#include <cstdio>
#include <cstdlib>

#include <starpu.h>

#include "../../Src/Files/FFmaScanfLoader.hpp"

#include "../../Src/Kernels/Chebyshev/FChebParticle.hpp"
#include "../../Src/Kernels/Chebyshev/FChebLeaf.hpp"
#include "../../Src/Kernels/Chebyshev/FChebCell.hpp"
#include "../../Src/Kernels/Chebyshev/FChebMatrixKernel.hpp"
#include "../../Src/Kernels/Chebyshev/FChebKernel.hpp"
#include "../../Src/Kernels/Chebyshev/FChebSymKernel.hpp"

#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"
#include "../../Src/Core/FFmmAlgorithmStarpuGroup.hpp"

template <int ORDER>
class GChebCell : public FChebCell<ORDER> {
public:
    void intialCopy(const GChebCell*const other){
        FMemUtils::copyall(FChebCell<ORDER>::getLocal(), other->FChebCell<ORDER>::getLocal(), TensorTraits<ORDER>::nnodes * 2);
        FMemUtils::copyall(FChebCell<ORDER>::getMultipole(), other->FChebCell<ORDER>::getMultipole(), TensorTraits<ORDER>::nnodes * 2);
    }
    void copyUp(const GChebCell*const other){
        FMemUtils::copyall(FChebCell<ORDER>::getMultipole(), other->FChebCell<ORDER>::getMultipole(), TensorTraits<ORDER>::nnodes * 2);
    }
    void restoreCopy(const GChebCell*const other){
        FMemUtils::copyall(FChebCell<ORDER>::getLocal(), other->FChebCell<ORDER>::getLocal(), TensorTraits<ORDER>::nnodes * 2);
        FMemUtils::copyall(FChebCell<ORDER>::getMultipole(), other->FChebCell<ORDER>::getMultipole(), TensorTraits<ORDER>::nnodes * 2);
    }
};


template <class ParticleClass>
class Container : public FVector<ParticleClass> {
public:
    void reduce(const Container*const other){
        for( int idx = 0 ; idx < FVector<ParticleClass>::getSize() ; ++idx){
            FVector<ParticleClass>::data()[idx].incForces(other->FVector<ParticleClass>::data()[idx].getForces());
            FVector<ParticleClass>::data()[idx].incPotential(other->FVector<ParticleClass>::data()[idx].getPotential());
        }
    }
};


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




// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
	const char* const filename       = FParameters::getStr(argc,argv,"-f", "../Data/test20k.fma");
	const unsigned int TreeHeight    = FParameters::getValue(argc, argv, "-h", 5);
	const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, "-sh", 2);
	const unsigned int NbThreads     = FParameters::getValue(argc, argv, "-t", 1);
        const unsigned int BlockSize     = FParameters::getValue(argc, argv, "-bs", 250);

        const unsigned int ORDER = 7;
        const FReal epsilon = FReal(1e-7);

	// set threads
	omp_set_num_threads(NbThreads); 
	std::cout << "Using " << omp_get_max_threads() << " threads." << std::endl;

	// init timer
	FTic time;

	// typedefs for STARPU
	typedef FChebParticle ParticleClass;
        typedef Container<ParticleClass> ContainerClass;
	typedef FChebLeaf<ParticleClass,ContainerClass> LeafClass;
	typedef FChebMatrixKernelR MatrixKernelClass;
        typedef GChebCell<ORDER> CellClass;
	typedef FOctree<ParticleClass,CellClass,ContainerClass,LeafClass> OctreeClass;
        //typedef FChebKernel<ParticleClass,CellClass,RealContainerClass,MatrixKernelClass,ORDER> KernelClass;
        typedef FChebSymKernel<ParticleClass,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
        typedef FFmmAlgorithmStarpuGroup<OctreeClass,ParticleClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;

	// What we do //////////////////////////////////////////////////////
	std::cout << ">> Testing the Chebyshev interpolation base FMM algorithm.\n";
	
	// open particle file
	FFmaScanfLoader<ParticleClass> loader(filename);
	if(!loader.isOpen()) throw std::runtime_error("Particle file couldn't be opened!");
	
	// init oct-tree
	OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());

	// -----------------------------------------------------
	std::cout << "Creating and inserting " << loader.getNumberOfParticles() << " particles in a octree of height " << TreeHeight
						<< " ..." << std::endl;
	time.tic();
	loader.fillTree(tree);
	std::cout << "Done  " << "(" << time.tacAndElapsed() << ")." << std::endl;
	// -----------------------------------------------------

	
	// -----------------------------------------------------
	std::cout << "\nChebyshev FMM ... " << std::endl;
	time.tic();
	KernelClass kernels(TreeHeight, loader.getCenterOfBox(), loader.getBoxWidth(), epsilon);
        FmmClass algorithm(&tree,&kernels, BlockSize);
        algorithm.buildGroups(NbThreads);
        std::cout << "init fmm & kernel " << time.tacAndElapsed() << "sec." << std::endl;

        time.tic();
	algorithm.execute();
	std::cout << "completed in " << time.tacAndElapsed() << "sec." << std::endl;

        algorithm.releaseGroups();
	// -----------------------------------------------------
	

	// -----------------------------------------------------
	// find first non empty leaf cell 
	OctreeClass::Iterator iLeafs(&tree);
	iLeafs.gotoBottomLeft();
	
	const ContainerClass *const Targets = iLeafs.getCurrentListTargets();
	const unsigned int NumTargets = Targets->getSize();

	FReal* Potential = new FReal [NumTargets];
	FBlas::setzero(NumTargets, Potential);

	FReal* Force = new FReal [NumTargets * 3];
	FBlas::setzero(NumTargets * 3, Force);
	
	std::cout << "\nDirect computation of " << NumTargets << " target particles ..." << std::endl;
	const MatrixKernelClass MatrixKernel;
	do {
			const ContainerClass *const Sources = iLeafs.getCurrentListSrc();
			unsigned int counter = 0;
			ContainerClass::ConstBasicIterator iTarget(*Targets);
			while(iTarget.hasNotFinished()) {
				const FReal wt = iTarget.data().getPhysicalValue();
				ContainerClass::ConstBasicIterator iSource(*Sources);
				while(iSource.hasNotFinished()) {
					if (&iTarget.data() != &iSource.data()) {
					  const FReal one_over_r = MatrixKernel.evaluate(iTarget.data().getPosition(),
																													 iSource.data().getPosition());
						const FReal ws = iSource.data().getPhysicalValue();
						// potential
						Potential[counter] += one_over_r * ws;
						// force
						FPoint force(iSource.data().getPosition() - iTarget.data().getPosition());
						force *= ((ws*wt) * (one_over_r*one_over_r*one_over_r));
						Force[counter*3 + 0] += force.getX();
						Force[counter*3 + 1] += force.getY();
						Force[counter*3 + 2] += force.getZ();
					}
					iSource.gotoNext();
				}
				counter++;
				iTarget.gotoNext();
			}
	} while(iLeafs.moveRight());

	
	FReal* ApproxPotential = new FReal [NumTargets];
	FReal* ApproxForce     = new FReal [NumTargets * 3];

	unsigned int counter = 0;
	ContainerClass::ConstBasicIterator iTarget(*Targets);
	while(iTarget.hasNotFinished()) {
		ApproxPotential[counter]   = iTarget.data().getPotential();
		ApproxForce[counter*3 + 0] = iTarget.data().getForces().getX();
		ApproxForce[counter*3 + 1] = iTarget.data().getForces().getY();
		ApproxForce[counter*3 + 2] = iTarget.data().getForces().getZ();
		counter++;
		iTarget.gotoNext();
	}

	std::cout << "\nPotential error:" << std::endl;
	std::cout << "Relative L2 error   = " << computeL2norm( NumTargets, Potential, ApproxPotential)
						<< std::endl;
	std::cout << "Relative Lmax error = " << computeINFnorm(NumTargets, Potential, ApproxPotential)
						<< std::endl;

	std::cout << "\nForce error:" << std::endl;
	std::cout << "Relative L2 error   = " << computeL2norm( NumTargets*3, Force, ApproxForce)
						<< std::endl;
	std::cout << "Relative Lmax error = " << computeINFnorm(NumTargets*3, Force, ApproxForce)
						<< std::endl;
	std::cout << std::endl;

	// free memory
	delete [] Potential;
	delete [] ApproxPotential;
	delete [] Force;
	delete [] ApproxForce;
	

	/*
	// Check if particles are strictly within its containing cells
	const FReal BoxWidthLeaf = BoxWidth / FReal(FMath::pow(2, TreeHeight-1));
	OctreeClass::Iterator octreeIterator(&tree);
	octreeIterator.gotoBottomLeft();
	do{
		const CellClass *const LeafCell = octreeIterator.getCurrentCell();
		const FPoint& LeafCellCenter = LeafCell -> getPosition();
		const ContainerClass *const Particles = octreeIterator.getCurrentListSrc();
		ContainerClass::ConstBasicIterator particleIterator(*Particles);
		while(particleIterator.hasNotFinished()) {
			const FPoint distance(LeafCellCenter-particleIterator.data().getPosition());
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



