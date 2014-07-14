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
 *@author Pierre Blanchard
 *
 * **/
// ==== CMAKE =====
// @FUSE_FFT
// ================

#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "Files/FFmaGenericLoader.hpp"


#include "Kernels/Uniform/FUnifCell.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Uniform/FUnifKernel.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "Utils/FParameters.hpp"
#include "Utils/FMemUtils.hpp"

#include "Containers/FOctree.hpp"
#include "Containers/FVector.hpp"

#include "Core/FFmmAlgorithm.hpp"
#include "Core/FFmmAlgorithmThread.hpp"

/**
 * This program runs the FMM Algorithm with the Uniform kernel and compares the results with a direct computation.
 */

// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
  const char* const filename       = FParameters::getStr(argc,argv,"-f", "../Data/test20k.fma");
  const unsigned int TreeHeight    = FParameters::getValue(argc, argv, "-depth", 3);
  const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, "-subdepth", 2);
  const unsigned int NbThreads     = FParameters::getValue(argc, argv, "-t", 1);

#ifdef _OPENMP
  omp_set_num_threads(NbThreads);
  std::cout << "\n>> Using " << omp_get_max_threads() << " threads.\n" << std::endl;
#else
  std::cout << "\n>> Sequential version.\n" << std::endl;
#endif

    // init timer
    FTic time;

  // interaction kernel evaluator
//  typedef FInterpMatrixKernelLJ MatrixKernelClass;
  typedef FInterpMatrixKernelR MatrixKernelClass;
//  typedef FInterpMatrixKernelRR MatrixKernelClass;
  const MatrixKernelClass MatrixKernel;

  // init particles position and physical value
  struct TestParticle{
    FPoint position;
    FReal forces[3];
    FReal physicalValue;
    FReal potential;
  };

  // open particle file
  FFmaGenericLoader loader(filename);
  if(!loader.isOpen()) throw std::runtime_error("Particle file couldn't be opened!");

  TestParticle* const particles = new TestParticle[loader.getNumberOfParticles()];
  for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
    FPoint position;
    FReal physicalValue = 0.0;
    loader.fillParticle(&position,&physicalValue);
    // get copy
    particles[idxPart].position       = position;
    particles[idxPart].physicalValue  = physicalValue;
    particles[idxPart].potential      = 0.0;
    particles[idxPart].forces[0]      = 0.0;
    particles[idxPart].forces[1]      = 0.0;
    particles[idxPart].forces[2]      = 0.0;
  }

  ////////////////////////////////////////////////////////////////////

  { // begin direct computation
    std::cout << "\nDirect Computation ... " << std::endl;
    time.tic();
    {
      for(int idxTarget = 0 ; idxTarget < loader.getNumberOfParticles() ; ++idxTarget){
        for(int idxOther = idxTarget + 1 ; idxOther < loader.getNumberOfParticles() ; ++idxOther){
          FP2P::MutualParticles(particles[idxTarget].position.getX(), particles[idxTarget].position.getY(),
                                particles[idxTarget].position.getZ(), particles[idxTarget].physicalValue,
                                &particles[idxTarget].forces[0], &particles[idxTarget].forces[1],
                                &particles[idxTarget].forces[2], &particles[idxTarget].potential,
                                particles[idxOther].position.getX(), particles[idxOther].position.getY(),
                                particles[idxOther].position.getZ(), particles[idxOther].physicalValue,
                                &particles[idxOther].forces[0], &particles[idxOther].forces[1],
                                &particles[idxOther].forces[2], &particles[idxOther].potential,
                                &MatrixKernel);
        }
      }
    }
    time.tac();
    std::cout << "Done  " << "(@Direct Computation = "
              << time.elapsed() << "s)." << std::endl;

  } // end direct computation

  ////////////////////////////////////////////////////////////////////

  {	// begin Lagrange kernel

    // accuracy
    const unsigned int ORDER = 5;

    // typedefs
    typedef FP2PParticleContainerIndexed<> ContainerClass;
    typedef FSimpleLeaf< ContainerClass >  LeafClass;
    typedef FUnifCell<ORDER> CellClass;
    typedef FOctree<CellClass,ContainerClass,LeafClass> OctreeClass;
    typedef FUnifKernel<CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
    typedef FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
    //  typedef FFmmAlgorithmThread<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;

    // init oct-tree
    OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());


    { // -----------------------------------------------------
      std::cout << "Creating & Inserting " << loader.getNumberOfParticles()
                << " particles ..." << std::endl;
      std::cout << "\tHeight : " << TreeHeight << " \t sub-height : " << SubTreeHeight << std::endl;
      time.tic();

      for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        // put in tree
        tree.insert(particles[idxPart].position, idxPart, particles[idxPart].physicalValue);
      }

      time.tac();
      std::cout << "Done  " << "(@Creating and Inserting Particles = "
                << time.elapsed() << "s)." << std::endl;
    } // -----------------------------------------------------

    { // -----------------------------------------------------
      std::cout << "\nLagrange/Uniform grid FMM (ORDER="<< ORDER << ") ... " << std::endl;
      time.tic();
      KernelClass kernels(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());
      FmmClass algorithm(&tree, &kernels);
      algorithm.execute();
      time.tac();
      std::cout << "Done  " << "(@Algorithm = " << time.elapsed() << "s)." << std::endl;
    } // -----------------------------------------------------


    { // -----------------------------------------------------
      std::cout << "\nError computation ... " << std::endl;
      FMath::FAccurater potentialDiff;
      FMath::FAccurater fx, fy, fz;

      FReal checkPotential[20000];

      { // Check that each particle has been summed with all other

        tree.forEachLeaf([&](LeafClass* leaf){
            const FReal*const potentials = leaf->getTargets()->getPotentials();
            const FReal*const forcesX = leaf->getTargets()->getForcesX();
            const FReal*const forcesY = leaf->getTargets()->getForcesY();
            const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
            const int nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
            const FVector<int>& indexes = leaf->getTargets()->getIndexes();

            for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
              const int indexPartOrig = indexes[idxPart];
              //PB: store potential in nbParticles array
              checkPotential[indexPartOrig]=potentials[idxPart];              

              potentialDiff.add(particles[indexPartOrig].potential,potentials[idxPart]);
              fx.add(particles[indexPartOrig].forces[0],forcesX[idxPart]);
              fy.add(particles[indexPartOrig].forces[1],forcesY[idxPart]);
              fz.add(particles[indexPartOrig].forces[2],forcesZ[idxPart]);
            }
          });
      }

      // Print for information
      std::cout << "Potential " << potentialDiff << std::endl;
      std::cout << "Fx " << fx << std::endl;
      std::cout << "Fy " << fy << std::endl;
      std::cout << "Fz " << fz << std::endl;
    } // -----------------------------------------------------

  } // end Lagrange kernel

  return 0;
}
