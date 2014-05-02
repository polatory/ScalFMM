// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info".
// "http://www.gnu.org/licenses".
// ===================================================================================


/**
 *@author Matthias Messner
 *
 * **/
// ==== CMAKE =====
// @FUSE_BLAS
// ================

#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "../../Src/Files/FFmaScanfLoader.hpp"

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
	const unsigned int TreeHeight    = FParameters::getValue(argc, argv, "-depth", 5);
	const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, "-subdepth", 2);

	const unsigned int ORDER = 8;
	const FReal epsilon = FReal(1e-8);

	// init timer
	FTic time;

	// typedefs
	typedef FP2PParticleContainer<> ContainerClass;
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
	//
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
	KernelClass kernels(TreeHeight, loader.getBoxWidth(),loader.getCenterOfBox(), epsilon);
	FmmClass algorithm(&tree,&kernels);
	time.tic();
	algorithm.execute();
	std::cout << "completed in " << time.tacAndElapsed() << "sec." << std::endl;
	// -----------------------------------------------------

	return 0;
}



