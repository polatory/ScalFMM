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
// @FUSE_FFT
//  ==== Git =====
// @SCALFMM_PRIVATE
// ================

#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <string>

#include "ScalFmmConfig.h"

#include "Files/FFmaGenericLoader.hpp"

#include "Kernels/Uniform/FUnifCell.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Uniform/FUnifKernel.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "Utils/FParameters.hpp"

#include "Containers/FOctree.hpp"

#ifdef _OPENMP
#include "Core/FFmmAlgorithmThread.hpp"
#else
#include "Core/FFmmAlgorithm.hpp"
#endif

#include "Utils/FParameterNames.hpp"

/**
 * This program runs the FMM Algorithm with the uniform interpolation kernel and compares the results with a direct computation.
 */
/// \file  LagrangeInterpolationFMM.cpp
//!
//! \brief This program runs the FMM Algorithm with the interpolation kernel based on uniform (grid points) interpolation (1/r kernel)
//!  \authors B. Bramas, O. Coulaud
//!
//!  This code is a short example to use the  Interpolation approach for the 1/r kernel


// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
	FHelpDescribeAndExit(argc, argv,
			"Driver for Lagrange interpolation kernel  (1/r kernel).",
			FParameterDefinitions::InputFile, FParameterDefinitions::OctreeHeight,
			FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::InputFile,FParameterDefinitions::OutputFile,
			FParameterDefinitions::NbThreads);


	const std::string defaultFile(/*SCALFMMDataPath+*/"Data/unitCubeXYZQ100.bfma" );
	const std::string filename                = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, defaultFile.c_str());
	const unsigned int TreeHeight        = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 5);
	const unsigned int SubTreeHeight  = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeSubHeight.options, 2);
	const unsigned int NbThreads        = FParameters::getValue(argc, argv, FParameterDefinitions::NbThreads.options, 1);


#ifdef _OPENMP
	omp_set_num_threads(NbThreads);
	std::cout << "\n>> Using " << omp_get_max_threads() << " threads.\n" << std::endl;
#else
	std::cout << "\n>> Sequential version.\n" << std::endl;
#endif
	//
	std::cout <<     "Parameters  "<< std::endl
			<<     "      Octree Depth      "<< TreeHeight <<std::endl
			<<        "      SubOctree depth " << SubTreeHeight <<std::endl
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
	// begin Lagrange kernel

	// accuracy
	const unsigned int ORDER = 7;
	// typedefs
	typedef FP2PParticleContainerIndexed<>                     ContainerClass;
	typedef FSimpleLeaf< ContainerClass >                        LeafClass;
	typedef FUnifCell<ORDER>                                            CellClass;
	typedef FOctree<CellClass,ContainerClass,LeafClass>  OctreeClass;
	//
	typedef FInterpMatrixKernelR                                                                              MatrixKernelClass;
	const MatrixKernelClass MatrixKernel;
	typedef FUnifKernel<CellClass,ContainerClass,MatrixKernelClass,ORDER>  KernelClass;
	//
#ifdef _OPENMP
	typedef FFmmAlgorithmThread<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
#else
	typedef FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
#endif
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
		std::cout << "\nLagrange FMM (ORDER="<< ORDER << ") ... " << std::endl;

		time.tic();
		//
		KernelClass *kernels = new KernelClass(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox(),&MatrixKernel);
		//
		FmmClass algo(&tree, kernels);
		//
		algo.execute();   // Here the call of the FMM algorithm
		//
		time.tac();
		std::cout << "Timers Far Field \n"
				<< "P2M " << algo.getTime(FAlgorithmTimers::P2MTimer) << " seconds\n"
				<< "M2M " << algo.getTime(FAlgorithmTimers::M2MTimer) << " seconds\n"
				<< "M2L " << algo.getTime(FAlgorithmTimers::M2LTimer) << " seconds\n"
				<< "L2L " << algo.getTime(FAlgorithmTimers::L2LTimer) << " seconds\n"
				<< "P2P and L2P " << algo.getTime(FAlgorithmTimers::NearTimer) << " seconds\n"
				<< std::endl;


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
	if(FParameters::existParameter(argc, argv, FParameterDefinitions::OutputFile.options)){
		std::string name(FParameters::getStr(argc,argv,FParameterDefinitions::OutputFile.options,   "output.fma"));
		FFmaGenericWriter writer(name) ;
		//
		int NbPoints = loader.getNumberOfParticles();
		FReal * particles ;
		particles = new FReal[8*NbPoints] ;
		memset(particles,0,8*NbPoints*sizeof(FReal));
		int j = 0 ;
		tree.forEachLeaf([&](LeafClass* leaf){
			//
			// Input
			const FReal*const posX = leaf->getTargets()->getPositions()[0];
			const FReal*const posY = leaf->getTargets()->getPositions()[1];
			const FReal*const posZ = leaf->getTargets()->getPositions()[2];
			const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();
			const FVector<int>& indexes = leaf->getTargets()->getIndexes();
			//
			// Computed data
			const FReal*const potentials = leaf->getTargets()->getPotentials();
			const FReal*const forcesX = leaf->getTargets()->getForcesX();
			const FReal*const forcesY = leaf->getTargets()->getForcesY();
			const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
			//
			const int nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
			for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
				j = 8*indexes[idxPart];
				particles[j]      = posX[idxPart] ;
				particles[j+1]  = posY[idxPart] ;
				particles[j+2]  = posZ[idxPart] ;
				particles[j+3]  = physicalValues[idxPart] ;
				particles[j+4]  = potentials[idxPart] ;
				particles[j+5]  =  forcesX[idxPart] ;
				particles[j+6]  =  forcesY[idxPart] ;
				particles[j+7]  =  forcesZ[idxPart] ;
			}
		});

		writer.writeHeader( loader.getCenterOfBox(), loader.getBoxWidth() ,  NbPoints, sizeof(FReal), 8) ;
		writer.writeArrayOfReal(particles,  8 , NbPoints);
	}

	return 0;
}