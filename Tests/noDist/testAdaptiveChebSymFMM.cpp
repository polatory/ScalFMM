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

// ==== CMAKE =====
// @FUSE_BLAS
// ================

#include <iostream>
#include <cstdio>


#include "Utils/FParameters.hpp"
#include "Utils/FTic.hpp"

#include "Containers/FOctree.hpp"
//#include "Containers/FVector.hpp"

//#include "Components/FSimpleLeaf.hpp"

#include "Utils/FPoint.hpp"

#include "Files/FFmaGenericLoader.hpp"
#include "Files/FRandomLoader.hpp"

#include "Components/FBasicKernels.hpp"
#include "Components/FSimpleIndexedLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebCell.hpp"

#include "Adaptive/FAdaptiveCell.hpp"
#include "Adaptive/FAdaptiveKernelWrapper.hpp"
#include "Adaptive/FAbstractAdaptiveKernel.hpp"
//
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Adaptive/FAdaptChebSymKernel.hpp"
#include "Adaptive/FAdaptTools.hpp"
//#include "AdaptiveTree/FAdaptChebSymKernel.hpp"

//
//
#include "Core/FFmmAlgorithm.hpp"
//#include "Core/FFmmAlgorithmThread.hpp"
//#include "Core/FFmmAlgorithmTask.hpp"
#include "Utils/FParameterNames.hpp"

/** This program show an example of use of the fmm basic algo
 * it also check that each particles is impacted each other particles
 */

void usage() {
	std::cout << "Driver to obtain statistics on the octree" << std::endl;
	std::cout <<	 "Options  "<< std::endl
			<<     "      -help       to see the parameters    " << std::endl
			<<	     "      -depth        the depth of the octree   "<< std::endl
			<<	     "      -subdepth   specifies the size of the sub octree   " << std::endl
			<<     "      -fin name specifies the name of the particle distribution" << std::endl
			<<     "      -sM    s_min^M threshold for Multipole (l+1)^2 for Spherical harmonics"<<std::endl
			<<     "      -sL    s_min^L threshold for Local  (l+1)^2 for Spherical harmonics"<<std::endl;
}
// Simply create particles and try the kernels
int main(int argc, char ** argv){
	//
	const FParameterNames LocalOptionMinMultipoleThreshod {
		{"-sM"},
		" s_min^M threshold for Multipole (l+1)^2 for Spherical harmonic."
	};
	const FParameterNames LocalOptionMinLocalThreshod {
		{"-SL"},
		" s_min^L threshold for Local  (l+1)^2 for Spherical harmonics."
	};

	FHelpDescribeAndExit(argc, argv,
			"Test Adaptive kernel and compare it with the direct computation.",
			FParameterDefinitions::OctreeHeight,FParameterDefinitions::NbThreads,
			FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::InputFile,
			LocalOptionMinMultipoleThreshod,LocalOptionMinLocalThreshod);

	for (int i = 0 ; i< argc ; ++i){
		std::cout << argv[i] << "  " ;
	}
	std::cout << std::endl<< std::endl;
	// ---------------------------------------------
	const std::string fileName(FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options,   "../Data/noDistprolate50.out.fma"));
	const unsigned int TreeHeight      = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 3);
	const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeSubHeight.options, 2);
	//
	// accuracy
	const unsigned int P = 5 ;
	//    typedef FTestCell                   CellClass;
	//    typedef FAdaptiveTestKernel< CellClass, ContainerClass >         KernelClass;
	typedef FChebCell<P>                                        CellClass;
	typedef FP2PParticleContainerIndexed<>            ContainerClass;
	typedef FSimpleIndexedLeaf<ContainerClass>    LeafClass;
	typedef FInterpMatrixKernelR                               MatrixKernelClass;
	//
	typedef FAdaptiveChebSymKernel<CellClass,ContainerClass,MatrixKernelClass,P> KernelClass;
	//
	//
	typedef FAdaptiveCell< CellClass, ContainerClass >                                        CellWrapperClass;
	typedef FAdaptiveKernelWrapper< KernelClass, CellClass, ContainerClass >   KernelWrapperClass;
	typedef FOctree< CellWrapperClass, ContainerClass , LeafClass >                  OctreeClass;

	// FFmmAlgorithmTask FFmmAlgorithmThread
	typedef FFmmAlgorithm<OctreeClass, CellWrapperClass, ContainerClass, KernelWrapperClass, LeafClass >     FmmClass;

	///////////////////////What we do/////////////////////////////
	std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
	//////////////////////////////////////////////////////////////
	//
	const int sminM    = FParameters::getValue(argc,argv,LocalOptionMinMultipoleThreshod.options, P*P*P);
	const int sminL     = FParameters::getValue(argc,argv,LocalOptionMinLocalThreshod.options, P*P*P);
	//


	FTic counter;

	//////////////////////////////////////////////////////////////////////////////////
	// Not Random Loader
	//////////////////////////////////////////////////////////////////////////////////
	FFmaGenericLoader loader(fileName);
	const long int NbPart  = loader.getNumberOfParticles() ;
	// Random Loader
	//const int NbPart       = FParameters::getValue(argc,argv,"-nb", 2000000);
	//	FRandomLoader loader(NbPart, 1, FPoint(0.5,0.5,0.5), 1);
	//////////////////////////////////////////////////////////////////////////////////

	OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());

	//////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////

	std::cout << "Creating & Inserting " << NbPart << " particles ..." << std::endl;
	std::cout << "\tHeight : " << TreeHeight << " \t sub-height : " << SubTreeHeight << std::endl;
	std::cout 		<< "         criteria SM:  "<< sminM     <<std::endl
			<< "         criteria SL:  "<< sminL     <<std::endl <<std::endl;
	//

	counter.tic();
	FReal L= loader.getBoxWidth();
	//
	FmaRWParticle<8,8>* const particles = new FmaRWParticle<8,8>[NbPart];

	FPoint minPos(L,L,L), maxPos(-L,-L,-L);
	//
	loader.fillParticle(particles,NbPart);

	for(int idxPart = 0 ; idxPart < NbPart; ++idxPart){
		const FPoint PP(particles[idxPart].getPosition() ) ;
		//
		minPos.setX(FMath::Min(minPos.getX(),PP.getX())) ;
		minPos.setY(FMath::Min(minPos.getY(),PP.getY())) ;
		minPos.setZ(FMath::Min(minPos.getZ(),PP.getZ())) ;
		maxPos.setX(FMath::Max(maxPos.getX(),PP.getX())) ;
		maxPos.setY(FMath::Max(maxPos.getY(),PP.getY())) ;
		maxPos.setZ(FMath::Max(maxPos.getZ(),PP.getZ())) ;
		//
		tree.insert(PP, idxPart, particles[idxPart].getPhysicalValue());

	}

	counter.tac();
	std::cout << "Data are inside the box delimited by "<<std::endl
			<< "         Min corner:  "<< minPos<<std::endl
			<< "         Max corner:  "<< maxPos<<std::endl <<std::endl;
	std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << " s)." << std::endl;

	//////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////

	std::cout << "Working on particles ..." << std::endl;
	// For debug purpose
	//  Set Global id
	//
	long int idCell  = setGlobalID(tree);
	counter.tic();
	const MatrixKernelClass MatrixKernel;
	KernelWrapperClass kernels(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox(),&MatrixKernel,sminM,sminL);            // FTestKernels FBasicKernels
	FmmClass algo(&tree,&kernels);  //FFmmAlgorithm FFmmAlgorithmThread
	algo.execute();

	counter.tac();
	std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << " s)." << std::endl;
	//
	FReal energy= 0.0 , energyD = 0.0 ;
	/////////////////////////////////////////////////////////////////////////////////////////////////
	// Compute direct energy
	/////////////////////////////////////////////////////////////////////////////////////////////////

	for(int idx = 0 ; idx <  loader.getNumberOfParticles()  ; ++idx){
		energyD +=  particles[idx].getPotential()*particles[idx].getPhysicalValue() ;
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////
	// Compare
	/////////////////////////////////////////////////////////////////////////////////////////////////
	FMath::FAccurater potentialDiff;
	FMath::FAccurater fx, fy, fz;
	{ // Check that each particle has been summed with all other

		//    std::cout << "indexPartOrig || DIRECT V fx || FMM V fx" << std::endl;

		tree.forEachLeaf([&](LeafClass* leaf){
			const FReal*const potentials        = leaf->getTargets()->getPotentials();
			const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();
			const FReal*const forcesX            = leaf->getTargets()->getForcesX();
			const FReal*const forcesY            = leaf->getTargets()->getForcesY();
			const FReal*const forcesZ            = leaf->getTargets()->getForcesZ();
			const int nbParticlesInLeaf           = leaf->getTargets()->getNbParticles();
			const FVector<int>& indexes      = leaf->getTargets()->getIndexes();

			for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
				const int indexPartOrig = indexes[idxPart];
				potentialDiff.add(particles[indexPartOrig].getPotential(),potentials[idxPart]);
				fx.add(particles[indexPartOrig].getForces()[0],forcesX[idxPart]);
				fy.add(particles[indexPartOrig].getForces()[1],forcesY[idxPart]);
				fz.add(particles[indexPartOrig].getForces()[2],forcesZ[idxPart]);
				energy   += potentials[idxPart]*physicalValues[idxPart];

				//          std::cout << indexPartOrig
				//                    << " " << particles[indexPartOrig].getPotential() << " " << particles[indexPartOrig].getForces()[0]
				//                    << " " << potentials[idxPart] << " " << forcesX[idxPart]
				//                    << std::endl;

			}
		});
	}

	delete[] particles;

	// Print for information
	std::cout << "Energy [relative L2 error] "  << FMath::Abs(energy-energyD) /energyD << std::endl;
	std::cout << "Potential " << potentialDiff << std::endl;
	std::cout << "Fx " << fx << std::endl;
	std::cout << "Fy " << fy << std::endl;
	std::cout << "Fz " << fz << std::endl;


	OctreeClass::Iterator octreeIterator(&tree);
	std::ofstream file("aa.tree", std::ofstream::out );

	//
	////////////////////////////////////////////////////////////////////
	//              Export adaptive tree in our format
	////////////////////////////////////////////////////////////////////
	//
	// -----------------------------------------------------
	//
	//
	//  Set Global id
	//
//	long int idCell  = setGlobalID(tree);
	//////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////

	tree.forEachCellLeaf([&](CellWrapperClass* cell, LeafClass* leaf){
		file << "Cell Id " << cell->getGlobalId( ) << " Nb particles "<<  leaf->getSrc()->getNbParticles()<<std::endl;
	});

	octreeIterator.gotoTop() ;  // here we are at level 1 (first child)
	//	octreeIterator.moveDown() ;
	octreeIterator.gotoLeft();
	//	octreeIterator.moveDown() ; // We are at the levell 2
	std::cout << " Number of Cells: " << idCell <<std::endl;
	//
	std::cout << "Top of the octree " << octreeIterator.level() << std::endl ;
	for(int idxLevel = 1; idxLevel < TreeHeight ;  ++idxLevel){
		file << std::endl << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<< std::endl;
		file << "  Level " << idxLevel <<"  Level  "<<  octreeIterator.level()<<  "  -- leave level " <<   std::boolalpha <<  octreeIterator.isAtLeafLevel() << std::endl;
		do{
			if(octreeIterator.getCurrentCell()->hasDevelopment()){
				file <<"Cell id  "<< octreeIterator.getCurrentCell()->getGlobalId( ) << "   "<<*(octreeIterator.getCurrentCell())<< std::endl ;
			}
		} while(octreeIterator.moveRight());
		octreeIterator.moveDown() ;
		octreeIterator.gotoLeft();
	}
	std::cout << "   END    " << std::endl;

	// Check
	octreeIterator.gotoBottomLeft();
	do {
		std::cout << " Level " <<octreeIterator.level() <<std::endl;
	}while(octreeIterator.moveUp() );
	std::cout << "   RETURN 0  " << std::endl;

	return 0;
}



