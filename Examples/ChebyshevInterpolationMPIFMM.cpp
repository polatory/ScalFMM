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
// @FUSE_MPI
// @FUSE_BLAS
// ================

#include <iostream>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>


#include "ScalFmmConfig.h"
#include "Containers/FOctree.hpp"
#include "Utils/FMpi.hpp"
#include "Core/FFmmAlgorithmThreadProc.hpp"

#include "Files/FFmaGenericLoader.hpp"
#include "Files/FMpiFmaGenericLoader.hpp"
#include "Files/FMpiTreeBuilder.hpp"

#include "BalanceTree/FLeafBalance.hpp"

#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"
#include "Kernels/Chebyshev/FChebCell.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"


/// \file ChebyshevInterpolationMPIFMM
//!
//! \brief This program runs the MPI FMM with Chebyshev interpolation of 1/r kernel
//!  \authors B. Bramas, O. Coulaud
//!
//!  This code is a short example to use the FMM Algorithm Proc with Chebyshev Interpolation for the 1/r kernel


// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
    FHelpDescribeAndExit(argc, argv,
			 "Driver for Chebyshev Interpolation kernel using MPI  (1/r kernel). "
			 "Usully run using : mpirun -np nb_proc_needed ./ChebyshevInterpolationAlgorithm [params].",
			 FParameterDefinitions::InputFile, FParameterDefinitions::OctreeHeight,
			 FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::InputFile,
			 FParameterDefinitions::NbThreads);

    const std::string defaultFile(/*SCALFMMDataPath+*/"../Data/test20k.fma");
    const std::string filename                = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, defaultFile.c_str());
    const unsigned int TreeHeight       = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 5);
    const unsigned int SubTreeHeight  = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeSubHeight.options, 2);
    const unsigned int NbThreads        = FParameters::getValue(argc, argv, FParameterDefinitions::NbThreads.options, 1);


#ifdef _OPENMP
	omp_set_num_threads(NbThreads);
	std::cout << "\n>> Using " << omp_get_max_threads() << " threads.\n" << std::endl;
#else
	std::cout << "\n>> Sequential version.\n" << std::endl;
#endif
	//
	std::cout <<	 "Parameters  "<< std::endl
		  <<     "      Octree Depth      "<< TreeHeight <<std::endl
		  <<	 "      SubOctree depth "<< SubTreeHeight <<std::endl
		  <<     "      Input file  name: " <<filename <<std::endl
		  <<     "      Thread number:  " << NbThreads <<std::endl
		  <<std::endl;
	//init values for MPI
	FMpi app(argc,argv);
	//
	// init timer
	FTic time;

    typedef double FReal;
    FMpiFmaGenericLoader<FReal> loader(filename,app.global());

	if(!loader.isOpen()) throw std::runtime_error("Particle file couldn't be opened!") ;
	////////////////////////////////////////////////////////////////////



    // begin spherical kernel
	// accuracy
	const unsigned int ORDER = 7;
	// typedefs
    typedef FP2PParticleContainerIndexed<FReal>                      ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass >                       LeafClass;
    typedef FChebCell<FReal,ORDER>                                    CellClass;
    typedef FOctree<FReal,CellClass,ContainerClass,LeafClass>         OctreeClass;
    typedef FInterpMatrixKernelR<FReal>                                MatrixKernelClass;
	const MatrixKernelClass MatrixKernel;
    typedef FChebSymKernel<FReal, CellClass,ContainerClass,MatrixKernelClass,ORDER>  KernelClass;

	//
	typedef FFmmAlgorithmThreadProc<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClassProc;

	// init oct-tree
	OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());


	{ // -----------------------------------------------------
	  if(app.global().processId() == 0){
	    std::cout << "Creating & Inserting " << loader.getNumberOfParticles()
		      << " particles ..." << std::endl;
	    std::cout << "\tHeight : " << TreeHeight << " \t sub-height : " << SubTreeHeight << std::endl;
	  }
	  time.tic();
	  //

	  struct TestParticle{
	    FSize index;
        FPoint<FReal> position;
	    FReal physicalValue;
        const FPoint<FReal>& getPosition(){
	      return position;
	    }
	  };

	  TestParticle* particles = new TestParticle[loader.getMyNumberOfParticles()];
	  memset(particles, 0, (unsigned int) (sizeof(TestParticle) * loader.getMyNumberOfParticles()));

	  //idx (in file) of the first part that will be used by this proc.
	  FSize idxStart = loader.getStart();
	  printf("Proc %d idxStart %lld \n",app.global().processId(),idxStart);

	  for(FSize idxPart = 0 ; idxPart < loader.getMyNumberOfParticles() ; ++idxPart){
	    //Storage of the index (in the original file) of each part.
	    particles[idxPart].index = idxPart + idxStart;
	    // Read particles from file
	    loader.fillParticle(&particles[idxPart].position,&particles[idxPart].physicalValue);
	  }

	  FVector<TestParticle> finalParticles;
      FLeafBalance balancer;
	  // FMpiTreeBuilder< FReal,TestParticle >::ArrayToTree(app.global(), particles, loader.getMyNumberOfParticles(),
	  //					 tree.getBoxCenter(),
	  //					 tree.getBoxWidth(),
	  //					 tree.getHeight(), &finalParticles,&balancer);
      FMpiTreeBuilder< FReal, TestParticle >::DistributeArrayToContainer(app.global(),particles,
								      loader.getMyNumberOfParticles(),
								      tree.getBoxCenter(),
								      tree.getBoxWidth(),tree.getHeight(),
								      &finalParticles, &balancer);

	  for(int idx = 0 ; idx < finalParticles.getSize(); ++idx){
	    tree.insert(finalParticles[idx].position,finalParticles[idx].index,finalParticles[idx].physicalValue);
	  }
	  printf("%d parts have been inserted in Tree \n",finalParticles.getSize());
	  delete[] particles;

	  time.tac();
	  double timeUsed = time.elapsed();
	  double minTime,maxTime;
	  std::cout << "Done  " << "(@Reading and Inserting Particles = "  << time.elapsed() << "s)." << std::endl;
	  MPI_Reduce(&timeUsed,&minTime,1,MPI_DOUBLE,MPI_MIN,0,app.global().getComm());
	  MPI_Reduce(&timeUsed,&maxTime,1,MPI_DOUBLE,MPI_MAX,0,app.global().getComm());
	  if(app.global().processId() == 0){
	      printf("Reading and Inserting Particles Time used : \t MIN : %f \t MAX %f in s\n",minTime,maxTime);
	  }
	} // -----------------------------------------------------

	{ // -----------------------------------------------------
	  std::cout << "\nChebyshev Interpolation  FMM Proc (P="<< ORDER << ") ... " << std::endl;

	  time.tic();
	  //
	  // Here we use a pointer due to the limited size of the stack
	  //
	  KernelClass *kernels = new KernelClass(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox(),&MatrixKernel);
	  //
	  FmmClassProc algorithm(app.global(),&tree, kernels);
	  //
	  algorithm.execute();   // Here the call of the FMM algorithm
	  //
	  time.tac();
	  double timeUsed = time.elapsed();
	  double minTime,maxTime;
	  std::cout << "Done  " << "(@Algorithm = " << time.elapsed() << " s)." << std::endl;
	  MPI_Reduce(&timeUsed,&minTime,1,MPI_DOUBLE,MPI_MIN,0,app.global().getComm());
	  MPI_Reduce(&timeUsed,&maxTime,1,MPI_DOUBLE,MPI_MAX,0,app.global().getComm());
	  if(app.global().processId() == 0){
	      printf("Time used : \t MIN : %f \t MAX %f \n",minTime,maxTime);
	  }
	}
	// -----------------------------------------------------
	//
	// Some output
	//
	//
	{ // -----------------------------------------------------
	  long int N1=0, N2= loader.getNumberOfParticles()/2, N3= (loader.getNumberOfParticles()-1); ;
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
		  std::cout << "Proc "<< app.global().processId() << " Index "<< indexPartOrig <<"  potential  " << potentials[idxPart]
			    << " Pos "<<posX[idxPart]<<" "<<posY[idxPart]<<" "<<posZ[idxPart]
			    << "   Forces: " << forcesX[idxPart] << " " << forcesY[idxPart] << " "<< forcesZ[idxPart] <<std::endl;
		}
		energy += potentials[idxPart]*physicalValues[idxPart] ;
	      }
	    });
	  FReal gloEnergy = app.global().reduceSum(energy);
	  if(0 == app.global().processId()){
	    std::cout <<std::endl << "Proc "<< app.global().processId() << " Energy: "<< gloEnergy <<std::endl;
	    std::cout <<std::endl <<" &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& "<<std::endl<<std::endl;
	  }
	}
	// -----------------------------------------------------


	return 0;
}
