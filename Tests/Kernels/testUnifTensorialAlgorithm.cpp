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
// @FUSE_FFT
// ================

#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "../../Src/Files/FFmaScanfLoader.hpp"
#include "../../Src/Files/FFmaBinLoader.hpp"



#include "../../Src/Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "../../Src/Kernels/Uniform/FUnifCell.hpp"
#include "../../Src/Kernels/Uniform/FUnifTensorialKernel.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FMemUtils.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"

/**
 * This program runs the FMM Algorithm with the Uniform kernel and compares the results with a direct computation.
 */

// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
  const char* const filename       = FParameters::getStr(argc,argv,"-f", "../Data/test20k.fma");
  const unsigned int TreeHeight    = FParameters::getValue(argc, argv, "-h", 3);
  const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, "-sh", 2);
  const unsigned int NbThreads     = FParameters::getValue(argc, argv, "-t", 1);

#ifdef _OPENMP
  omp_set_num_threads(NbThreads);
  std::cout << "\n>> Using " << omp_get_max_threads() << " threads.\n" << std::endl;
#else
  std::cout << "\n>> Sequential version.\n" << std::
#endif

    // init timer
    FTic time;


  // typedefs
  typedef FInterpMatrixKernel_R_IJ MatrixKernelClass;
//  typedef FInterpMatrixKernel_IOR MatrixKernelClass;

  const KERNEL_FUNCTION_IDENTIFIER MK_ID = MatrixKernelClass::Identifier;
  const unsigned int NRHS = MatrixKernelClass::NRHS;
  const unsigned int NLHS = MatrixKernelClass::NLHS;

  // init particles position and physical value
  struct TestParticle{
    FPoint position;
    FReal forces[3][NLHS];
    FReal physicalValue[NRHS];
    FReal potential[NLHS];
  };

  // open particle file
  FFmaScanfLoader loader(filename);
  if(!loader.isOpen()) throw std::runtime_error("Particle file couldn't be opened!");

  TestParticle* const particles = new TestParticle[loader.getNumberOfParticles()];
  for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
    FPoint position;
    FReal physicalValue = 0.0;
    loader.fillParticle(&position,&physicalValue);
    // get copy
    particles[idxPart].position       = position;
    for(unsigned idxRhs = 0; idxRhs<NRHS;++idxRhs)
      particles[idxPart].physicalValue[idxRhs]  = physicalValue; // copy same physical value in each component
    for(unsigned idxLhs = 0; idxLhs<NLHS;++idxLhs){
      particles[idxPart].potential[idxLhs]      = 0.0;
      particles[idxPart].forces[0][idxLhs]      = 0.0;
      particles[idxPart].forces[1][idxLhs]      = 0.0;
      particles[idxPart].forces[2][idxLhs]      = 0.0;
    }
  }

  ////////////////////////////////////////////////////////////////////

  { // begin direct computation
    std::cout << "\nDirect Computation ... " << std::endl;
    time.tic();
    {
      for(int idxTarget = 0 ; idxTarget < loader.getNumberOfParticles() ; ++idxTarget){
        for(int idxOther = idxTarget + 1 ; idxOther < loader.getNumberOfParticles() ; ++idxOther){
          if(MK_ID == R_IJ)
            FP2P::MutualParticlesRIJ(particles[idxTarget].position.getX(), particles[idxTarget].position.getY(),
                                     particles[idxTarget].position.getZ(), particles[idxTarget].physicalValue,
                                     particles[idxTarget].forces[0], particles[idxTarget].forces[1],
                                     particles[idxTarget].forces[2], particles[idxTarget].potential,
                                     particles[idxOther].position.getX(), particles[idxOther].position.getY(),
                                     particles[idxOther].position.getZ(), particles[idxOther].physicalValue,
                                     particles[idxOther].forces[0], particles[idxOther].forces[1],
                                     particles[idxOther].forces[2], particles[idxOther].potential);
          else if(MK_ID == ID_OVER_R)
            FP2P::MutualParticlesIOR(particles[idxTarget].position.getX(), particles[idxTarget].position.getY(),
                                     particles[idxTarget].position.getZ(), particles[idxTarget].physicalValue,
                                     particles[idxTarget].forces[0], particles[idxTarget].forces[1],
                                     particles[idxTarget].forces[2], particles[idxTarget].potential,
                                     particles[idxOther].position.getX(), particles[idxOther].position.getY(),
                                     particles[idxOther].position.getZ(), particles[idxOther].physicalValue,
                                     particles[idxOther].forces[0], particles[idxOther].forces[1],
                                     particles[idxOther].forces[2], particles[idxOther].potential);
          else 
            std::runtime_error("Provide a valid matrix kernel!");
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
    const unsigned int ORDER = 7;

    typedef FP2PParticleContainerIndexed<NRHS,NLHS> ContainerClass;

    typedef FSimpleLeaf< ContainerClass >  LeafClass;
    typedef FUnifCell<ORDER,NRHS,NLHS> CellClass;
    typedef FOctree<CellClass,ContainerClass,LeafClass> OctreeClass;
    typedef FUnifTensorialKernel<CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
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
        // PB: here we have to know the number of NRHS...
        tree.insert(particles[idxPart].position, idxPart, particles[idxPart].physicalValue[0], particles[idxPart].physicalValue[1], particles[idxPart].physicalValue[2]);
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
      FMath::FAccurater potentialDiff[NLHS];
      FMath::FAccurater fx[NLHS], fy[NLHS], fz[NLHS];

      FReal checkPhysVal[20000][NRHS];
      FReal checkPotential[20000][NLHS];
      FReal checkfx[20000][NLHS];

      { // Check that each particle has been summed with all other

        tree.forEachLeaf([&](LeafClass* leaf){
            for(unsigned idxLhs = 0; idxLhs<NLHS;++idxLhs){

              const FReal*const physVals = leaf->getTargets()->getPhysicalValues(idxLhs);
              const FReal*const potentials = leaf->getTargets()->getPotentials(idxLhs);
              const FReal*const forcesX = leaf->getTargets()->getForcesX(idxLhs);
              const FReal*const forcesY = leaf->getTargets()->getForcesY(idxLhs);
              const FReal*const forcesZ = leaf->getTargets()->getForcesZ(idxLhs);
              const int nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
              const FVector<int>& indexes = leaf->getTargets()->getIndexes();

              for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                const int indexPartOrig = indexes[idxPart];

                //PB: store potential in array[nbParticles]
                checkPhysVal[indexPartOrig][idxLhs]=physVals[idxPart];              
                checkPotential[indexPartOrig][idxLhs]=potentials[idxPart];              
                checkfx[indexPartOrig][idxLhs]=forcesX[idxPart];              

                potentialDiff[idxLhs].add(particles[indexPartOrig].potential[idxLhs],potentials[idxPart]);
                fx[idxLhs].add(particles[indexPartOrig].forces[0][idxLhs],forcesX[idxPart]);
                fy[idxLhs].add(particles[indexPartOrig].forces[1][idxLhs],forcesY[idxPart]);
                fz[idxLhs].add(particles[indexPartOrig].forces[2][idxLhs],forcesZ[idxPart]);
              }
            }// NLHS
          });
      }

//      std::cout << "Check Potential, forceX " << std::endl;
//      for(int idxPart = 0 ; idxPart < 20 ; ++idxPart)
//        for(unsigned idxLhs = 0; idxLhs<NLHS;++idxLhs){
//          std::cout << checkPhysVal[idxPart][idxLhs] << ", "<< particles[idxPart].physicalValue[idxLhs]<< "|| ";
//          std::cout << checkPotential[idxPart][idxLhs] << ", "<< particles[idxPart].potential[idxLhs]<< "|| ";
//          std::cout << checkfx[idxPart][idxLhs] << ", "<< particles[idxPart].forces[0][idxLhs] << std::endl;
//        }
//      std::cout << std::endl;

      // Print for information
      std::cout << "\nAbsolute errors: " << std::endl;
      std::cout << "Potential: ";
      for(unsigned idxLhs = 0; idxLhs<NLHS;++idxLhs) std::cout << potentialDiff[idxLhs] << ", " ;
      std::cout << std::endl;
      std::cout << "Fx: "; 
      for(unsigned idxLhs = 0; idxLhs<NLHS;++idxLhs) std::cout << fx[idxLhs] << ", " ;
      std::cout  << std::endl;
      std::cout << "Fy: "; 
      for(unsigned idxLhs = 0; idxLhs<NLHS;++idxLhs) std::cout << fy[idxLhs] << ", " ;
      std::cout  << std::endl;
      std::cout << "Fz: "; 
      for(unsigned idxLhs = 0; idxLhs<NLHS;++idxLhs) std::cout << fz[idxLhs] << ", " ;
      std::cout << std::endl;

    } // -----------------------------------------------------

  } // end Lagrange kernel

  return 0;
}
