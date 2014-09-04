// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas
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

// ==== CMAKE =====
// @FUSE_BLAS
// ================

#include <iostream>
#include <fstream>

#include <cstdio>
#include <cstdlib>
#include <string>

#include "ScalFmmConfig.h"

#include "Files/FFmaGenericLoader.hpp"
#include "Files/FTreeIO.hpp"

#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "Utils/FParameters.hpp"

#include "Containers/FOctree.hpp"


#include "Core/FFmmAlgorithmNoP2P.hpp"


/**
 * This program runs the FMM Algorithm with the Chebyshev kernel and compares the results with a direct computation.
 */
/// \file  ChebyshevInterpolationFMM.cpp
//!
//! \brief This program runs the FMM Algorithm with the interpolation kernel based on Chebyshev interpolation (1/r kernel)
//!  \authors B. Bramas, O. Coulaud
//!
//!  This code is a short example to use the Chebyshev Interpolation approach for the 1/r kernel
//!
//!@Algorithm
//!  <b> General arguments:</b>
//!     \param   -help(-h)      to see the parameters available in this driver
//!     \param   -depth          The depth of the octree
//!     \param   -subdepth     Specifies the size of the sub octree
//!     \param   -t                   The number of threads
//!
//!     \param   -f name          Name of the particles file. The file have to be in our FMA format
//!     \param   -bin                 if the file is in binary mode

//!
//

void usage() {
	std::cout << "Driver for Chebyshev interpolation kernel  (1/r kernel)" << std::endl;
	std::cout <<	 "Options  "<< std::endl
			<<     "      -help         to see the parameters    " << std::endl
			<<	  "      -depth       the depth of the octree   "<< std::endl
			<<	  "      -subdepth  specifies the size of the sub octree   " << std::endl
			<<     "      -fin   name    name specifies the name of the particle distribution" << std::endl
			<<     "      -t  n  specifies the number of threads used in the computations" << std::endl
			<<     "      -savetree  file  specify the file name to save the tree" << std::endl
			<<     "      -fout   name    name specifies the name for the computational values" << std::endl ;

}

// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
	const std::string defaultFile(/*SCALFMMDataPath+*/"../Data/test20k.fma" );
	std::string filename                = FParameters::getStr(argc,argv,"-fin", defaultFile.c_str());
	const unsigned int TreeHeight        = FParameters::getValue(argc, argv, "-depth", 5);
	const unsigned int SubTreeHeight  = FParameters::getValue(argc, argv, "-subdepth", 2);
	const unsigned int NbThreads        = FParameters::getValue(argc, argv, "-t", 1);
	if(FParameters::existParameter(argc, argv, "-h")||FParameters::existParameter(argc, argv, "-help")){
		usage() ;
		exit(EXIT_SUCCESS);
	}
#ifdef _OPENMP
	omp_set_num_threads(NbThreads);
	std::cout << "\n>> Using " << omp_get_max_threads() << " threads.\n" << std::endl;
#else
	std::cout << "\n>> Sequential version.\n" << std::endl;
#endif
	//
	std::cout <<	 "Parameters  "<< std::endl
			<<     "      Octree Depth      "<< TreeHeight <<std::endl
			<<	  "      SubOctree depth " << SubTreeHeight <<std::endl
			<<     "      Input file  name: " <<filename <<std::endl
			<<     "      Thread number:  " << NbThreads <<std::endl
			<<std::endl;
	//
	// init timer
	FTic time;

	// open particle file
	////////////////////////////////////////////////////////////////////
	//
	FFmaGenericLoader loader(filename);
	//
	////////////////////////////////////////////////////////////////////
	// begin Chebyshev kernel

	// accuracy
	const unsigned int ORDER = 3;
	// typedefs
	typedef FP2PParticleContainerIndexed<>                     ContainerClass;
	typedef FSimpleLeaf< ContainerClass >                        LeafClass;
	typedef FChebCell<ORDER>                                         CellClass;
	typedef FOctree<CellClass,ContainerClass,LeafClass>  OctreeClass;
	//
	typedef FInterpMatrixKernelR                                                                              MatrixKernelClass;
	typedef FChebSymKernel<CellClass,ContainerClass,MatrixKernelClass,ORDER>  KernelClass;
	//
//#ifdef _OPENMP
//	typedef FFmmAlgorithmThread<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
//#else
//	typedef FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
//#endif
	typedef FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;

	// init oct-tree
	OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());


	{ // -----------------------------------------------------
		std::cout << "Creating & Inserting " << loader.getNumberOfParticles()
																															<< " particles ..." << std::endl;
		std::cout << "\tHeight : " << TreeHeight << " \t sub-height : " << SubTreeHeight << std::endl;
		time.tic();
		//
		FPoint position;
		FReal physicalValue = 0.0;
		//
		for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
			//
			// Read particle per particle from file
			loader.fillParticle(&position,&physicalValue);
			//
			// put particle in octree
			tree.insert(position, idxPart, physicalValue);
		}

		time.tac();
		std::cout << "Done  " << "(@Creating and Inserting Particles = "
				<< time.elapsed() << " s) ." << std::endl;
	} // -----------------------------------------------------

	{ // -----------------------------------------------------
		std::cout << "\nChebyshev FMM (ORDER="<< ORDER << ") ... " << std::endl;

		time.tic();
		//
		  const MatrixKernelClass MatrixKernel;
		KernelClass kernels(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox(),&MatrixKernel);
		//
		FmmClass algorithm(&tree, &kernels);
		//
		algorithm.execute();   // Here the call of the FMM algorithm
		//
		time.tac();
		std::cout << "Done  " << "(@Algorithm = " << time.elapsed() << " s) ." << std::endl;
	}
	// -----------------------------------------------------
	//
	// Some output
	//
	//
	{ // -----------------------------------------------------
		long int N1=0, N2= loader.getNumberOfParticles()/2, N3= loader.getNumberOfParticles() -1; ;
		FReal energy =0.0 ;
		//
		//   Loop over all leaves
		//
		std::cout <<std::endl<<" &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& "<<std::endl;
		std::cout << std::scientific;
		std::cout.precision(10) ;

		tree.forEachLeaf([&](LeafClass* leaf){
			const FReal*const posX = leaf->getTargets()->getPositions()[0];
			const FReal*const posY = leaf->getTargets()->getPositions()[1];
			const FReal*const posZ = leaf->getTargets()->getPositions()[2];

			const FReal*const potentials = leaf->getTargets()->getPotentials();
			const FReal*const forcesX = leaf->getTargets()->getForcesX();
			const FReal*const forcesY = leaf->getTargets()->getForcesY();
			const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
			const int nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
			const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();

			const FVector<int>& indexes = leaf->getTargets()->getIndexes();

			for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
				const int indexPartOrig = indexes[idxPart];
				if ((indexPartOrig == N1) || (indexPartOrig == N2) || (indexPartOrig == N3)  ) {
					std::cout << "Index "<< indexPartOrig <<"  potential  " << potentials[idxPart]
					                                                                      << " Pos "<<posX[idxPart]<<" "<<posY[idxPart]<<" "<<posZ[idxPart]
					                                                                                                                               << "   Forces: " << forcesX[idxPart] << " " << forcesY[idxPart] << " "<< forcesZ[idxPart] <<std::endl;
				}
				energy += potentials[idxPart]*physicalValues[idxPart] ;
			}
		});
		std::cout <<std::endl<<"Energy: "<< energy<<std::endl;
		std::cout <<std::endl<<" &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& "<<std::endl<<std::endl;

	}
	// -----------------------------------------------------
	if(FParameters::existParameter(argc, argv, "-savetree")){

		const std::string filename  = FParameters::getStr(argc,argv,"-savetree", "tree.bin");
		FTreeIO::Save<OctreeClass, CellClass, LeafClass, ContainerClass >(filename.c_str(), tree);

		//////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////
		int idCell=0 ;
		OctreeClass::Iterator octreeIterator(&tree);
		std::ofstream file("aa-REF.tree", std::ofstream::out );
		file << "Number of particles per leaf"<< std::endl;
		tree.forEachCellLeaf([&](CellClass* cell, LeafClass* leaf){
			file << "Cell Id " << ++idCell<< " Nb particles "<<  leaf->getSrc()->getNbParticles()<<std::endl;
		});
		file << std::endl << "//////////////////////////////////////////////////////////////////////////////////"
				<< std::endl;

		octreeIterator.gotoTop() ;  // here we are at level 1 (first child)
		//	octreeIterator.moveDown() ;
		octreeIterator.gotoLeft();
		//	octreeIterator.moveDown() ; // We are at the levell 2
		std::cout << " Number of Cells: " << idCell <<std::endl;
		//
		idCell=0 ;
		std::cout << "Top of the octree " << octreeIterator.level() << std::endl ;
		for(int idxLevel = 1; idxLevel < TreeHeight ;  ++idxLevel){
			file << std::endl << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<< std::endl;
			file << "  Level " << idxLevel <<"  Level  "<<  octreeIterator.level()<<  "  -- leave level " <<   std::boolalpha <<  octreeIterator.isAtLeafLevel() << std::endl;
			do{
				file <<"Cell id  "<< ++idCell << "   " <<*(octreeIterator.getCurrentCell())<< std::endl ;

			} while(octreeIterator.moveRight());
			octreeIterator.moveDown() ;
			octreeIterator.gotoLeft();
		}
		std::cout << "   END    " << std::endl;

	}

	if(FParameters::existParameter(argc, argv, "-fout")){

		filename  = FParameters::getStr(argc,argv,"-fout", defaultFile.c_str());
		FFmaGenericWriter writer(filename) ;
		//
		FSize NBparticles = loader.getNumberOfParticles() ;
		FmaRWParticle<8,8>	*particles = new  FmaRWParticle<8,8>[NBparticles];
		//
		tree.forEachLeaf([&](LeafClass* leaf){
			const FReal*const posX = leaf->getTargets()->getPositions()[0];
			const FReal*const posY = leaf->getTargets()->getPositions()[1];
			const FReal*const posZ = leaf->getTargets()->getPositions()[2];

			const FReal*const potentials = leaf->getTargets()->getPotentials();
			const FReal*const forcesX = leaf->getTargets()->getForcesX();
			const FReal*const forcesY = leaf->getTargets()->getForcesY();
			const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
			const int nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
			const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();

			const FVector<int>& indexes = leaf->getTargets()->getIndexes();

			for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
				const int indexPartOrig = indexes[idxPart];
				particles[indexPartOrig].setPosition(posX[idxPart],posY[idxPart],posZ[idxPart]);
				particles[indexPartOrig].setPhysicalValue(physicalValues[idxPart]) ;
				particles[indexPartOrig].setPotential( potentials[idxPart]);
				particles[indexPartOrig].setForces(forcesX[idxPart] ,forcesY[idxPart], forcesZ[idxPart]);
			}
		});

		writer.writeHeader(loader.getCenterOfBox(), loader.getBoxWidth(), NBparticles,*particles) ;

		writer.writeArrayOfParticles(particles,  NBparticles);
	}
	return 0;
}
