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

#include <cstdio>
#include <cstdlib>

#include "../../Src/Files/FFmaScanfLoader.hpp"
#include "../../Src/Files/FFmaBinLoader.hpp"

#include "../../Src/Kernels/Chebyshev/FChebCell.hpp"
#include "../../Src/Kernels/Interpolation/FInterpMatrixKernel.hpp"

#include "../../Src/Kernels/Chebyshev/FChebFlopsSymKernel.hpp"

#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"




int main(int argc, char* argv[])
{
	const char* const filename       = FParameters::getStr(argc,argv,"-f", "../Data/test20k.fma");
	const unsigned int TreeHeight    = FParameters::getValue(argc, argv, "-h", 5);
	const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, "-sh", 2);

	const unsigned int ORDER = 8;
	const FReal epsilon = FReal(1e-8);

	// init timer
	FTic time;

    // typedefs
    typedef FP2PParticleContainer ContainerClass;
    typedef FSimpleLeaf<ContainerClass> LeafClass;
	typedef FInterpMatrixKernelR MatrixKernelClass;
	typedef FChebCell<ORDER> CellClass;
    typedef FOctree<CellClass,ContainerClass,LeafClass> OctreeClass;
    typedef FChebFlopsSymKernel<CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
    typedef FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;


	// What we do //////////////////////////////////////////////////////
	std::cout << ">> Testing the Chebyshev interpolation base FMM algorithm.\n";
	
	// open particle file
    FFmaScanfLoader loader(filename);
    //FFmaBinLoader loader(filename);
	if(!loader.isOpen()) throw std::runtime_error("Particle file couldn't be opened!");
	
	// init oct-tree
	OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());

	// -----------------------------------------------------
	std::cout << "Creating and inserting " << loader.getNumberOfParticles()
						<< " particles in a octree of height " << TreeHeight << " ..." << std::endl;
	time.tic();

    {
        FPoint particlePosition;
        FReal physicalValue = 0.0;
        for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(&particlePosition,&physicalValue);
            tree.insert(particlePosition, physicalValue);
        }
    }

	std::cout << "Done  " << "(" << time.tacAndElapsed() << ")." << std::endl;
	// -----------------------------------------------------

	
	// -----------------------------------------------------
	std::cout << "\nChebyshev FMM ... " << std::endl;
	KernelClass kernels(TreeHeight, loader.getCenterOfBox(), loader.getBoxWidth(), epsilon);
	FmmClass algorithm(&tree,&kernels);
	time.tic();
	algorithm.execute();
	std::cout << "completed in " << time.tacAndElapsed() << "sec." << std::endl;
	// -----------------------------------------------------
	
    return 0;
}



