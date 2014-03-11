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



#include "../../Src/Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "../../Src/Kernels/Chebyshev/FChebCell.hpp"
#include "../../Src/Kernels/Chebyshev/FChebTensorialKernel.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FMemUtils.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"

/**
 * This program runs the FMM Algorithm with the Chebyshev kernel and compares the results with a direct computation.
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
//  typedef FInterpMatrixKernel_R_IJ MatrixKernelClass;
//  typedef FInterpMatrixKernel_IOR MatrixKernelClass;
//  typedef FInterpMatrixKernelR MatrixKernelClass;
  typedef FInterpMatrixKernel_R_IJK MatrixKernelClass;

  const KERNEL_FUNCTION_IDENTIFIER MK_ID = MatrixKernelClass::Identifier;
  const unsigned int NPV  = MatrixKernelClass::NPV;
  const unsigned int NPOT = MatrixKernelClass::NPOT;
  const unsigned int NRHS = MatrixKernelClass::NRHS;
  const unsigned int NLHS = MatrixKernelClass::NLHS;

  const double CoreWidth = 0.1;
  const MatrixKernelClass DirectMatrixKernel(CoreWidth);
  std::cout<< "CoreWidth2 = "<< DirectMatrixKernel.getCoreWidth2()<<std::endl;

  // init particles position and physical value
  struct TestParticle{
    FPoint position;
    FReal forces[3][NPOT];
    FReal physicalValue[NPV];
    FReal potential[NPOT];
  };

  const FReal FRandMax = FReal(RAND_MAX);

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
    // Set physical values
    for(unsigned idxPV = 0; idxPV<NPV;++idxPV){
    //   Either copy same physical value in each component
      particles[idxPart].physicalValue[idxPV]  = physicalValue; 
    // ... or set random value
//      particles[idxPart].physicalValue[idxPV]  = physicalValue*FReal(rand())/FRandMax;
    }

    for(unsigned idxPot = 0; idxPot<NPOT;++idxPot){
      particles[idxPart].potential[idxPot]      = 0.0;
      particles[idxPart].forces[0][idxPot]      = 0.0;
      particles[idxPart].forces[1][idxPot]      = 0.0;
      particles[idxPart].forces[2][idxPot]      = 0.0;
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
                                     particles[idxOther].forces[2], particles[idxOther].potential,
                                     &DirectMatrixKernel);
          else if(MK_ID == ID_OVER_R)
            FP2P::MutualParticlesIOR(particles[idxTarget].position.getX(), particles[idxTarget].position.getY(),
                                     particles[idxTarget].position.getZ(), particles[idxTarget].physicalValue,
                                     particles[idxTarget].forces[0], particles[idxTarget].forces[1],
                                     particles[idxTarget].forces[2], particles[idxTarget].potential,
                                     particles[idxOther].position.getX(), particles[idxOther].position.getY(),
                                     particles[idxOther].position.getZ(), particles[idxOther].physicalValue,
                                     particles[idxOther].forces[0], particles[idxOther].forces[1],
                                     particles[idxOther].forces[2], particles[idxOther].potential);
          else if(MK_ID == ONE_OVER_R)
            FP2P::MutualParticles(particles[idxTarget].position.getX(), particles[idxTarget].position.getY(),
                                  particles[idxTarget].position.getZ(), particles[idxTarget].physicalValue[0],
                                  particles[idxTarget].forces[0], particles[idxTarget].forces[1],
                                  particles[idxTarget].forces[2], particles[idxTarget].potential,
                                  particles[idxOther].position.getX(), particles[idxOther].position.getY(),
                                  particles[idxOther].position.getZ(), particles[idxOther].physicalValue[0],
                                  particles[idxOther].forces[0], particles[idxOther].forces[1],
                                  particles[idxOther].forces[2], particles[idxOther].potential);
          else if(MK_ID == R_IJK)
            FP2P::MutualParticlesRIJK(particles[idxTarget].position.getX(), particles[idxTarget].position.getY(),
                                      particles[idxTarget].position.getZ(), particles[idxTarget].physicalValue,
                                      particles[idxTarget].forces[0], particles[idxTarget].forces[1],
                                      particles[idxTarget].forces[2], particles[idxTarget].potential,
                                      particles[idxOther].position.getX(), particles[idxOther].position.getY(),
                                      particles[idxOther].position.getZ(), particles[idxOther].physicalValue,
                                      particles[idxOther].forces[0], particles[idxOther].forces[1],
                                      particles[idxOther].forces[2], particles[idxOther].potential,
                                      &DirectMatrixKernel);
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

  {	// begin Chebyshev kernel

    // accuracy
const unsigned int ORDER = 7 ;
    const FReal epsilon = FReal(1e-7);

    typedef FP2PParticleContainerIndexed<NRHS,NLHS> ContainerClass;

    typedef FSimpleLeaf< ContainerClass >  LeafClass;
    typedef FChebCell<ORDER,NRHS,NLHS> CellClass;
    typedef FOctree<CellClass,ContainerClass,LeafClass> OctreeClass;
    typedef FChebTensorialKernel<CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
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
        // PB: here we have to know NPV...
        if(NPV==1)
          tree.insert(particles[idxPart].position, idxPart, particles[idxPart].physicalValue[0]);
        else if(NPV==3)
          tree.insert(particles[idxPart].position, idxPart, particles[idxPart].physicalValue[0], particles[idxPart].physicalValue[1], particles[idxPart].physicalValue[2]);
        else if(NPV==9) // R_IJK
          tree.insert(particles[idxPart].position, idxPart, 
                      particles[idxPart].physicalValue[0], particles[idxPart].physicalValue[1], particles[idxPart].physicalValue[2],
                      particles[idxPart].physicalValue[3], particles[idxPart].physicalValue[4], particles[idxPart].physicalValue[5],
                      particles[idxPart].physicalValue[6], particles[idxPart].physicalValue[7], particles[idxPart].physicalValue[8]);
        else 
          std::runtime_error("Insert correct number of physical values!");

      }

      time.tac();
      std::cout << "Done  " << "(@Creating and Inserting Particles = "
                << time.elapsed() << "s)." << std::endl;
    } // -----------------------------------------------------

    { // -----------------------------------------------------
      std::cout << "\nChebyshev FMM (ORDER="<< ORDER << ") ... " << std::endl;
      time.tic();
      KernelClass kernels(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox(), epsilon,CoreWidth);
      FmmClass algorithm(&tree, &kernels);
      algorithm.execute();
      time.tac();
      std::cout << "Done  " << "(@Algorithm = " << time.elapsed() << "s)." << std::endl;
    } // -----------------------------------------------------


    { // -----------------------------------------------------
      std::cout << "\nError computation ... " << std::endl;
      FMath::FAccurater potentialDiff[NPOT];
      FMath::FAccurater fx[NPOT], fy[NPOT], fz[NPOT];

      FReal checkPotential[20000][NPOT];
      FReal checkfx[20000][NPOT];

      { // Check that each particle has been summed with all other

        tree.forEachLeaf([&](LeafClass* leaf){
            for(unsigned idxPot = 0; idxPot<NPOT;++idxPot){

              const FReal*const potentials = leaf->getTargets()->getPotentials(idxPot);
              const FReal*const forcesX = leaf->getTargets()->getForcesX(idxPot);
              const FReal*const forcesY = leaf->getTargets()->getForcesY(idxPot);
              const FReal*const forcesZ = leaf->getTargets()->getForcesZ(idxPot);
              const int nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
              const FVector<int>& indexes = leaf->getTargets()->getIndexes();

              for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                const int indexPartOrig = indexes[idxPart];

                //PB: store potential in array[nbParticles]
                checkPotential[indexPartOrig][idxPot]=potentials[idxPart];
                checkfx[indexPartOrig][idxPot]=forcesX[idxPart];

                // update accuracy
                potentialDiff[idxPot].add(particles[indexPartOrig].potential[idxPot],potentials[idxPart]);
                fx[idxPot].add(particles[indexPartOrig].forces[0][idxPot],forcesX[idxPart]);
                fy[idxPot].add(particles[indexPartOrig].forces[1][idxPot],forcesY[idxPart]);
                fz[idxPot].add(particles[indexPartOrig].forces[2][idxPot],forcesZ[idxPart]);
              }
            }// NPOT
          });
      }

//      std::cout << "Check Potential, forceX " << std::endl;
//      for(int idxPart = 0 ; idxPart < 20 ; ++idxPart)
//        for(unsigned idxPot = 0; idxPot<NPOT;++idxPot){
//          std::cout << checkPotential[idxPart][idxPot] << ", "<< particles[idxPart].potential[idxPot]<< "|| ";
//          std::cout << checkfx[idxPart][idxPot] << ", "<< particles[idxPart].forces[0][idxPot] << std::endl;
//        }
//      std::cout << std::endl;

      // Print for information
      std::cout << "\nRelative Inf/L2 errors: " << std::endl;
      std::cout << "  Potential: " << std::endl;
      for(unsigned idxPot = 0; idxPot<NPOT;++idxPot) {
        std::cout << "    " << idxPot << ": "
                  << potentialDiff[idxPot].getRelativeInfNorm() << ", " 
                  << potentialDiff[idxPot].getRelativeL2Norm() 
                  << std::endl;
      }
      std::cout << std::endl;
      std::cout << "  Fx: " << std::endl; 
      for(unsigned idxPot = 0; idxPot<NPOT;++idxPot) {
        std::cout << "    " << idxPot << ": "
                  << fx[idxPot].getRelativeInfNorm() << ", " 
                  << fx[idxPot].getRelativeL2Norm()
                  << std::endl;
      }
      std::cout  << std::endl;
      std::cout << "  Fy: " << std::endl; 
      for(unsigned idxPot = 0; idxPot<NPOT;++idxPot) {
        std::cout << "    " << idxPot << ": "
                  << fy[idxPot].getRelativeInfNorm() << ", " 
                  << fy[idxPot].getRelativeL2Norm()
                  << std::endl;
      }
      std::cout  << std::endl;
      std::cout << "  Fz: " << std::endl; 
      for(unsigned idxPot = 0; idxPot<NPOT;++idxPot) {
        std::cout << "    " << idxPot << ": "
                  << fz[idxPot].getRelativeInfNorm() << ", " 
                  << fz[idxPot].getRelativeL2Norm()
                  << std::endl;
      }
      std::cout << std::endl;

      { // -----------------------------------------------------


        std::cout << "\nStore results in file ... "<<std::endl;
        std::ostringstream sstream;
        if(MK_ID == R_IJ)
          sstream << "testChebRij_h"<< TreeHeight << "_a" << 1000*CoreWidth;
        else if(MK_ID == R_IJK)
          sstream << "testChebRijk_h"<< TreeHeight << "_a" << 1000*CoreWidth;

        if(TreeHeight>2) sstream <<"_o"<<ORDER;
        const std::string para_ext = sstream.str();
        std::string outname = "../tmp/"+para_ext+".dat";
        std::ofstream fout(outname.c_str());
        fout.precision(15);
        
        for(unsigned idxLhs = 0; idxLhs<NLHS;++idxLhs) 
          fout << potentialDiff[idxLhs].getL2Norm() << " " 
               << potentialDiff[idxLhs].getInfNorm() << " "
               << potentialDiff[idxLhs].getRelativeL2Norm() << " " 
               << potentialDiff[idxLhs].getRelativeInfNorm() << std::endl;
        for(unsigned idxLhs = 0; idxLhs<NLHS;++idxLhs) 
          fout << fx[idxLhs].getL2Norm() << " " 
               << fx[idxLhs].getInfNorm() << " "
               << fx[idxLhs].getRelativeL2Norm() << " " 
               << fx[idxLhs].getRelativeInfNorm() << std::endl;
        for(unsigned idxLhs = 0; idxLhs<NLHS;++idxLhs) 
          fout << fy[idxLhs].getL2Norm() << " " 
               << fy[idxLhs].getInfNorm() << " "
               << fy[idxLhs].getRelativeL2Norm() << " " 
               << fy[idxLhs].getRelativeInfNorm() << std::endl;
        for(unsigned idxLhs = 0; idxLhs<NLHS;++idxLhs) 
          fout << fz[idxLhs].getL2Norm() << " " 
               << fz[idxLhs].getInfNorm() << " "
               << fz[idxLhs].getRelativeL2Norm() << " " 
               << fz[idxLhs].getRelativeInfNorm() << std::endl;
        fout << std::flush;
        fout.close();

      } // -----------------------------------------------------

    } // -----------------------------------------------------

  } // end Chebyshev kernel

  return 0;
}
