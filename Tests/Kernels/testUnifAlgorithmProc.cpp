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
// @FUSE_MPI
// ================

#include <iostream>

#include <cstdio>
#include <cstdlib>


#include "../../Src/Kernels/Uniform/FUnifCell.hpp"
#include "../../Src/Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "../../Src/Kernels/Uniform/FUnifKernel.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FMemUtils.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Files/FFmaScanfLoader.hpp"
#include "../../Src/Files/FFmaBinLoader.hpp"
#include "../../Src/Files/FMpiFmaLoader.hpp"
#include "../../Src/Files/FMpiTreeBuilder.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"
#include "../../Src/Core/FFmmAlgorithmThreadProc.hpp"

/**
 * This program runs the FMM Algorithm Distributed with the Uniform kernel
 */

// Simply create particles and try the kernel
int main(int argc, char* argv[])
{
  
  const unsigned int ORDER = 7;
  const FReal epsilon = FReal(1e-7);

  typedef FP2PParticleContainerIndexed<> ContainerClass;
  typedef FSimpleLeaf< ContainerClass >  LeafClass;

  typedef FInterpMatrixKernelR MatrixKernelClass;
  typedef FUnifCell<ORDER> CellClass;
  typedef FOctree<CellClass,ContainerClass,LeafClass> OctreeClass;

  typedef FUnifKernel<CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
  typedef FFmmAlgorithmThreadProc<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
  
  FMpi app(argc,argv);
  
  const char* const filename       = FParameters::getStr(argc,argv,"-f", "../Data/test20k.fma");
  const unsigned int TreeHeight    = FParameters::getValue(argc, argv, "-h", 5);
  const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, "-sh", 2);
  const unsigned int NbThreads     = FParameters::getValue(argc, argv, "-t", 1);

  std::cout << ">> This executable has to be used to test Proc Uniform Algorithm. \n";
  

#ifdef _OPENMP
  omp_set_num_threads(NbThreads);
  std::cout << "\n>> Using " << omp_get_max_threads() << " threads.\n" << std::endl;
#else
  std::cout << "\n>> Sequential version.\n" << std::endl;
#endif
    
  std::cout << "Opening : " <<filename << "\n" << std::endl;
    // init timer
    FTic time;
  
  // init particles position and physical value
  struct TestParticle{
    FPoint position;
    FReal physicalValue;
    const FPoint& getPosition(){
      return position;
    }
  };

  // open particle file
  FMpiFmaLoader loader(filename,app.global());
  if(!loader.isOpen()) throw std::runtime_error("Particle file couldn't be opened!");
  
  OctreeClass tree(TreeHeight, SubTreeHeight,loader.getBoxWidth(),loader.getCenterOfBox());
  time.tic();
  TestParticle* particles = new TestParticle[loader.getNumberOfParticles()];
  memset(particles,0,(unsigned int) (sizeof(TestParticle)* loader.getNumberOfParticles()));
  for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
    loader.fillParticle(&particles[idxPart].position,&particles[idxPart].physicalValue);
  }
  FVector<TestParticle> finalParticles;
  FMpiTreeBuilder< TestParticle >::ArrayToTree(app.global(), particles, loader.getNumberOfParticles(),
					       tree.getBoxCenter(),
					       tree.getBoxWidth(),
					       tree.getHeight(), &finalParticles);		 

  { // -----------------------------------------------------
    std::cout << "Creating & Inserting " << loader.getNumberOfParticles()
	      << " particles ..." << std::endl;
    std::cout << "\tHeight : " << TreeHeight << " \t sub-height : " << SubTreeHeight << std::endl;
    time.tic();

    for(int idxPart = 0 ; idxPart < finalParticles.getSize() ; ++idxPart){
      // put in tree
      tree.insert(finalParticles[idxPart].position, idxPart, finalParticles[idxPart].physicalValue);
    }

    time.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = "
	      << time.elapsed() << "s)." << std::endl;
  } // -----------------------------------------------------

  { // -----------------------------------------------------
    std::cout << "\nChebyshev FMM (ORDER="<< ORDER << ",EPS="<< epsilon <<") ... " << std::endl;
    time.tic();
    KernelClass kernels(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());
    FmmClass algorithm(app.global(),&tree, &kernels);
    algorithm.execute();
    time.tac();
    std::cout << "Done  " << "(@Algorithm = " << time.elapsed() << "s)." << std::endl;
  } // -----------------------------------------------------


  { // -----------------------------------------------------
    std::cout << "\nError computation ... " << std::endl;
    FReal potential;
    FReal fx, fy, fz;
    { // Check that each particle has been summed with all other
      tree.forEachLeaf([&](LeafClass* leaf){
	  const FReal*const potentials = leaf->getTargets()->getPotentials();
	  const FReal*const forcesX = leaf->getTargets()->getForcesX();
	  const FReal*const forcesY = leaf->getTargets()->getForcesY();
	  const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
	  const int nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
	  for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
	    potential += potentials[idxPart];
	    fx += forcesX[idxPart];
	    fy += forcesY[idxPart];
	    fz += forcesZ[idxPart];
	  }
	});
    }

    // Print for information
    std::cout << "Potential " << potential << std::endl;
    std::cout << "Fx " << fx << std::endl;
    std::cout << "Fy " << fy << std::endl;
    std::cout << "Fz " << fz << std::endl;
    
  } // end Chebyshev kernel


  return 0;
}
