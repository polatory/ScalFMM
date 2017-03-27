// See LICENCE file at project root

// ==== CMAKE =====
// @FUSE_BLAS
//  ==== Git =====
//
// ================

#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <string>

#include "ScalFmmConfig.h"

#include "Files/FFmaGenericLoader.hpp"

#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "Utils/FParameters.hpp"

#include "Containers/FOctree.hpp"

#ifdef _OPENMP
#include "Core/FFmmAlgorithmThread.hpp"
#else
#include "Core/FFmmAlgorithm.hpp"
#endif

#include "Utils/FTemplate.hpp"
#include "Utils/FParameterNames.hpp"


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
	    <<     "      -f   name    name specifies the name of the particle distribution" << std::endl
	    <<     "      -t  n  specifies the number of threads used in the computations" << std::endl;
}

// Simply create particles and try the kernels
struct TempMainStruct{
  template <const unsigned int ORDER>
  static void Run(int argc, char* argv[])
  {
    //
    const std::string defaultFile(/*SCALFMMDataPath+*/"../Data/UTest/unitCubeRef20kDouble.bfma" );
    const std::string filename                = FParameters::getStr(argc,argv,"-f", defaultFile.c_str());
    const unsigned int TreeHeight        = FParameters::getValue(argc, argv, "-depth", 5);
    const unsigned int SubTreeHeight  = FParameters::getValue(argc, argv, "-subdepth", 2);
    const unsigned int NbThreads        = FParameters::getValue(argc, argv, "-t", 1);

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
    typedef double FReal;
    FFmaGenericLoader<FReal> loader(filename);
    //
    FSize nbParticles = loader.getNumberOfParticles() ;
    FmaRWParticle<FReal,8,8>* const particles = new FmaRWParticle<FReal, 8,8>[nbParticles];

    loader.fillParticle(particles,nbParticles);
    FReal  energyD = 0.0 ;
    /////////////////////////////////////////////////////////////////////////////////////////////////
    // Compute direct energy
    /////////////////////////////////////////////////////////////////////////////////////////////////

    for(int idx = 0 ; idx <  nbParticles  ; ++idx){
      energyD +=  particles[idx].getPotential()*particles[idx].getPhysicalValue() ;
    }

    //
    ////////////////////////////////////////////////////////////////////
    // begin Chebyshev kernel

    // accuracy
    // typedefs
    typedef FP2PParticleContainerIndexed<FReal>                                 ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass >                                    LeafClass;

    typedef FChebCell<FReal,ORDER>                                                    CellClass;
    typedef FOctree<FReal,CellClass,ContainerClass,LeafClass>                       OctreeClass;
    //
    typedef FInterpMatrixKernelR<FReal>                                        MatrixKernelClass;
    typedef FChebSymKernel<FReal, CellClass,ContainerClass,MatrixKernelClass,ORDER>  KernelClass;
    //
#ifdef _OPENMP
    typedef FFmmAlgorithmThread<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
#else
    typedef FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
#endif
    // init oct-tree
    OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());


    { // -----------------------------------------------------
      std::cout << "Creating & Inserting " <<nbParticles
		<< " particles ..." << std::endl;
      std::cout << "\tHeight : " << TreeHeight << " \t sub-height : " << SubTreeHeight << std::endl;
      time.tic();
      //
      for(FSize idxPart = 0 ; idxPart < nbParticles; ++idxPart){
	//
	// Read particle per particle from file
	//
	// put particle in octree
	tree.insert(particles[idxPart].getPosition() , idxPart, particles[idxPart].getPhysicalValue() );
      }

      time.tac();
      std::cout << "Done  " << "(@Creating and Inserting Particles = "
		<< time.elapsed() << " s) ." << std::endl;
    } // -----------------------------------------------------

    { // -----------------------------------------------------
      std::cout << "\nChebyshev FMM (ORDER="<< ORDER << ") ... " << std::endl;
      const MatrixKernelClass matrixClass;
      time.tic();
      //
      KernelClass kernels(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox(),&matrixClass);
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
      FReal energy =0.0;

      //
      //   Loop over all leaves
      //
      std::cout <<std::endl<<" &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& "<<std::endl;
      std::cout << std::scientific;
      std::cout.precision(10) ;

      /////////////////////////////////////////////////////////////////////////////////////////////////
      // Compare
      /////////////////////////////////////////////////////////////////////////////////////////////////
      printf("Compute Diff...");
      FMath::FAccurater<FReal> potentialDiff;
      FMath::FAccurater<FReal> fx, fy, fz, f;
      { // Check that each particle has been summed with all other

	tree.forEachLeaf([&](LeafClass* leaf){
	    const FReal*const potentials        = leaf->getTargets()->getPotentials();
	    const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();
	    const FReal*const forcesX            = leaf->getTargets()->getForcesX();
	    const FReal*const forcesY            = leaf->getTargets()->getForcesY();
	    const FReal*const forcesZ            = leaf->getTargets()->getForcesZ();
	    const FSize nbParticlesInLeaf           = leaf->getTargets()->getNbParticles();
        const FVector<FSize>& indexes      = leaf->getTargets()->getIndexes();

	    for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
	      const FSize indexPartOrig = indexes[idxPart];
	      potentialDiff.add(particles[indexPartOrig].getPotential(),potentials[idxPart]);
	      fx.add(particles[indexPartOrig].getForces()[0],forcesX[idxPart]);
	      //	const std::string outputFile("accuracyChebyschev.txt") ;
	      fy.add(particles[indexPartOrig].getForces()[1],forcesY[idxPart]);
	      fz.add(particles[indexPartOrig].getForces()[2],forcesZ[idxPart]);
	      f.add(particles[indexPartOrig].getForces()[0],forcesX[idxPart]);
	      f.add(particles[indexPartOrig].getForces()[1],forcesY[idxPart]);
	      f.add(particles[indexPartOrig].getForces()[2],forcesZ[idxPart]);
	      energy   += potentials[idxPart]*physicalValues[idxPart];
	    }
	  });
      }
      std::cout << energy << " " << energyD << std::endl;
      delete[] particles;

      f.setNbElements(nbParticles);
      std::cout << "FChebSymKernel Energy "  << FMath::Abs(energy-energyD) <<  "  Relative     "<< FMath::Abs(energy-energyD) / FMath::Abs(energyD) <<std::endl;
      std::cout << "FChebSymKernel Potential " << potentialDiff << std::endl;
      std::cout << "FChebSymKernel Fx " << fx << std::endl;
      std::cout << "FChebSymKernel Fy " << fy << std::endl;
      std::cout << "FChebSymKernel Fz " << fz << std::endl;
      std::cout << "FChebSymKernel F  " << f << std::endl;

      const std::string outputFile("accuracyChebyschev.txt") ;
      std::ofstream output(outputFile,std::fstream::out | std::fstream::app);

      output << ORDER << "   " << FMath::Abs(energy-energyD) / FMath::Abs(energyD)
	     << "  " << potentialDiff.getRelativeL2Norm() << "  "<< potentialDiff.getRelativeInfNorm()<< "  "<< potentialDiff.getRMSError()
	     << "  " << f.getRelativeL2Norm() << "  "<< f.getRelativeInfNorm() <<"  "<< f.getRMSError()
	     << std::endl;
    }
    // -----------------------------------------------------

  }
};


int main(int argc, char** argv){
  //	const std::string outputFile("accuracyChebyschev.txt") ;
  //	std::ofstream output(outputFile,std::fstream::out | std::fstream::app);
  FHelpDescribeAndExit(argc, argv,
		       "Test for different accuracy).",
		       FParameterDefinitions::InputFile, FParameterDefinitions::OctreeHeight,
		       FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::InputFile,
		       FParameterDefinitions::NbThreads);
  //
  const unsigned int order = FParameters::getValue(argc, argv, "-order", 5);
  std::cout << "Order given by user is : " << order << "\n";

  // il faudrait un For ici
  FRunIf::Run<unsigned int, 2, 13, 1, TempMainStruct>(order,argc, argv);
  return 0;
}
