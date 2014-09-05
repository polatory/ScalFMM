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

#include <cstdio>
#include <cstdlib>
#include <string>

#include "ScalFmmConfig.h"
#include "Utils/FGlobal.hpp"

#include "Files/FFmaGenericLoader.hpp"

#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "Utils/FParameters.hpp"

#include "Containers/FOctree.hpp"

#include "Core/FFmmAlgorithmPeriodic.hpp"

template <class Output>
void Print(const Output& value){
	std::cout<< "--- Output from program : " << value << "\n";
}

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
//!     \param   -f name          Name of the particles file with extension (.fma or .bfma). The data in  file have to be in our FMA format
//!
//

void usage() {
	std::cout << "Driver for Chebyshev interpolation kernel  (1/r kernel)" << std::endl;
	std::cout <<	 "Options  "<< std::endl
			<<     "      -help         to see the parameters    " << std::endl
			<<	  "      -depth       the depth of the octree   "<< std::endl
			<<	  "      -subdepth  specifies the size of the sub octree   " << std::endl
			<<     "      -fin   name    name specifies the name of the particle distribution" << std::endl
			<<     "      -fout   name    to write the computation in a file" << std::endl
			<<     "      -t  n  specifies the number of threads used in the computations" << std::endl;
}

// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
	const std::string defaultFile(/*SCALFMMDataPath+*/"../Data/test20k.fma" );
	const std::string filename                = FParameters::getStr(argc,argv,"-fin", defaultFile.c_str());
	const std::string filenameOut          = FParameters::getStr(argc,argv,"-fout", "resultPer.fma");
	const unsigned int TreeHeight        = FParameters::getValue(argc, argv, "-depth", 5);
	const unsigned int SubTreeHeight  = FParameters::getValue(argc, argv, "-subdepth", 2);
	const unsigned int NbThreads        = FParameters::getValue(argc, argv, "-t", 1);
	const int PeriodicDeep                     = FParameters::getValue(argc,argv,"-per", 3);

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
	FSize nbParticles = loader.getNumberOfParticles() ;
	FmaRWParticle<8,8>* const particles = new FmaRWParticle<8,8>[nbParticles];

	//
	////////////////////////////////////////////////////////////////////
	// begin Chebyshev kernel

	// accuracy
	const unsigned int ORDER = 7;
	// typedefs
	typedef FP2PParticleContainerIndexed<>                     ContainerClass;
	typedef FSimpleLeaf< ContainerClass >                        LeafClass;
	typedef FChebCell<ORDER>                                         CellClass;
	typedef FOctree<CellClass,ContainerClass,LeafClass>  OctreeClass;
	//
	typedef FInterpMatrixKernelR                                                                              MatrixKernelClass;
	const MatrixKernelClass MatrixKernel;
	typedef FChebSymKernel<CellClass,ContainerClass,MatrixKernelClass,ORDER>  KernelClass;
	//
	typedef FFmmAlgorithmPeriodic<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;

	// init oct-tree
	OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());


	{ // -----------------------------------------------------
		std::cout << "Creating & Inserting " << loader.getNumberOfParticles()
																																	<< " particles ..." << std::endl;
		std::cout << "\tHeight : " << TreeHeight << " \t sub-height : " << SubTreeHeight << std::endl;
		time.tic();
		//
		FPoint position;
		//
		loader.fillParticle(particles,nbParticles);

		for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
			tree.insert(particles[idxPart].getPosition(), idxPart, particles[idxPart].getPhysicalValue() );
		}

		time.tac();
		std::cout << "Done  " << "(@Creating and Inserting Particles = "
				<< time.elapsed() << " s) ." << std::endl;
	} // -----------------------------------------------------

	/////////////////////////////////////////////////////////////////////////////////////////////////
	{ // -----------------------------------------------------
		std::cout << "\nChebyshev FMM (ORDER="<< ORDER << ") ... " << std::endl;

		time.tic();

		FmmClass algo(&tree,PeriodicDeep );
		const MatrixKernelClass MatrixKernel;
		KernelClass kernels(algo.extendedTreeHeight(), algo.extendedBoxWidth(), algo.extendedBoxCenter(),&MatrixKernel);
		algo.setKernel(&kernels);
		algo.execute();
		//
		time.tac();
		std::cout << "Done  " << "(@Algorithm = " << time.elapsed() << " s) ." << std::endl;

	}
	// -----------------------------------------------------
	//
	// Some output
	//
	//
	FmaRWParticle<8,8>* const particlesOut = new FmaRWParticle<8,8>[nbParticles];

	{ // -----------------------------------------------------
		//
		FReal energy= 0.0 , energyD = 0.0 ;
		/////////////////////////////////////////////////////////////////////////////////////////////////
		// Compute direct energy
		/////////////////////////////////////////////////////////////////////////////////////////////////

		for(int idx = 0 ; idx <  loader.getNumberOfParticles()  ; ++idx){
			energyD +=  particles[idx].getPotential()*particles[idx].getPhysicalValue() ;
		}
		//
		//   Loop over all leaves
		//
		std::cout <<std::endl<<" &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& "<<std::endl;
		std::cout << std::scientific;
		std::cout.precision(10) ;
		/////////////////////////////////////////////////////////////////////////////////////////////////
		// Compare
		/////////////////////////////////////////////////////////////////////////////////////////////////
		FMath::FAccurater potentialDiff;
		FMath::FAccurater fx, fy, fz;

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
				//
				particlesOut[indexPartOrig].setPosition(posX[idxPart],posY[idxPart],posZ[idxPart]) ;
				particlesOut[indexPartOrig].setPhysicalValue(physicalValues[idxPart]) ;
				particlesOut[indexPartOrig].setPotential (potentials[idxPart]) ;
				particlesOut[indexPartOrig].setForces(forcesX[idxPart],forcesY[idxPart],forcesZ[idxPart]) ;

				potentialDiff.add(particles[indexPartOrig].getPotential(),potentials[idxPart]);
				fx.add(particles[indexPartOrig].getForces()[0],forcesX[idxPart]);
				fy.add(particles[indexPartOrig].getForces()[1],forcesY[idxPart]);
				fz.add(particles[indexPartOrig].getForces()[2],forcesZ[idxPart]);
				energy+=potentials[idxPart]*physicalValues[idxPart];
			}
		});
		energy *= 0.5;
		std::cout <<std::endl<<"Energy: "<< energy<<std::endl;
		std::cout <<std::endl<<" &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& "<<std::endl<<std::endl;


		// remove index
		Print("Test1 - Error Relative L2 norm Potential ");
		printf("         Pot L2Norm     %e\n",potentialDiff.getL2Norm());
		printf("         Pot RL2Norm   %e\n",potentialDiff.getRelativeL2Norm());
		printf("         Pot RMSError   %e\n",potentialDiff.getRMSError());
		Print("Fx diff is = ");
		printf("         Fx L2Norm     %e\n",fx.getL2Norm());
		printf("         Fx RL2Norm   %e\n",fx.getRelativeL2Norm());
		printf("         Fx RMSError   %e\n",fx.getRMSError());
		Print("Fy diff is = ");
		printf("        Fy L2Norm     %e\n",fy.getL2Norm());
		printf("        Fy RL2Norm   %e\n",fy.getRelativeL2Norm());
		printf("        Fy RMSError   %e\n",fy.getRMSError());
		Print("Fz diff is = ");
		printf("        Fz L2Norm     %e\n",fz.getL2Norm());
		printf("        Fz RL2Norm   %e\n",fz.getRelativeL2Norm());
		printf("        Fz RMSError   %e\n",fz.getRMSError());
		FReal L2error = (fx.getRelativeL2Norm()*fx.getRelativeL2Norm() + fy.getRelativeL2Norm()*fy.getRelativeL2Norm()  + fz.getRelativeL2Norm() *fz.getRelativeL2Norm()  );
		printf(" Total L2 Force Error= %e\n",FMath::Sqrt(L2error)) ;
		printf("  Energy Error  =   %.12e\n",FMath::Abs(energy-energyD));
		printf("  Energy FMM    =   %.12e\n",FMath::Abs(energy));
		printf("  Energy DIRECT =   %.12e\n",FMath::Abs(energyD));

	}
	// -----------------------------------------------------
	if(FParameters::existParameter(argc, argv, "-fout") ){

		std::cout << "Generate " << filenameOut <<"  for output file" << std::endl;
		//
		std::cout << " numberofParticles: " << nbParticles <<"  " << sizeof(nbParticles) <<std::endl;
		std::cout << " Box size: " << loader.getBoxWidth() << "  " << sizeof(loader.getBoxWidth())<<std::endl;
		//
		FFmaGenericWriter writer(filenameOut) ;
		writer.writeHeader(loader.getCenterOfBox(), loader.getBoxWidth() , nbParticles,*particlesOut) ;
		writer.writeArrayOfParticles(particlesOut, nbParticles);

	}
	delete [] particlesOut;
	return 0;
}
