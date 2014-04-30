// ===================================================================================
// Ce LOGICIEL "ScalFmm" est couvert par le copyright Inria 20xx-2012.
// Inria d��tient tous les droits de propri��t�� sur le LOGICIEL, et souhaite que
// la communaut�� scientifique l'utilise afin de le tester et de l'��valuer.
// Inria donne gracieusement le droit d'utiliser ce LOGICIEL. Toute utilisation
// dans un but lucratif ou �� des fins commerciales est interdite sauf autorisation
// expresse et pr��alable d'Inria.
// Toute utilisation hors des limites pr��cis��es ci-dessus et r��alis��e sans l'accord
// expresse pr��alable d'Inria constituerait donc le d��lit de contrefa��on.
// Le LOGICIEL ��tant un produit en cours de d��veloppement, Inria ne saurait assurer
// aucune responsabilit�� et notamment en aucune mani��re et en aucun cas, ��tre tenu
// de r��pondre d'��ventuels dommages directs ou indirects subits par l'utilisateur.
// Tout utilisateur du LOGICIEL s'engage �� communiquer �� Inria ses remarques
// relatives �� l'usage du LOGICIEL
// ===================================================================================

// ==== CMAKE =====
// @FUSE_BLAS
// ================

#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "../../Src/Files/FFmaScanfLoader.hpp"
#include "../../Src/Files/FFmaBinLoader.hpp"


#include "../../Src/Kernels/Chebyshev/FChebCell.hpp"
#include "../../Src/Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "../../Src/Kernels/Chebyshev/FChebKernel.hpp"
#include "../../Src/Kernels/Chebyshev/FChebSymKernel.hpp"

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
  const char* const filename             = FParameters::getStr(argc,argv,"-f", "../Data/test20k.fma");
  const unsigned int TreeHeight       = FParameters::getValue(argc, argv, "-h", 5);
  const unsigned int SubTreeHeight  = FParameters::getValue(argc, argv, "-sh", 2);
  const unsigned int NbThreads        = FParameters::getValue(argc, argv, "-t", 1);

#ifdef _OPENMP
  omp_set_num_threads(NbThreads);
  std::cout << "\n>> Using " << omp_get_max_threads() << " threads.\n" << std::endl;
#else
  std::cout << "\n>> Sequential version.\n" << std::
#endif

  // init timer
  FTic time;

  // init particles position and physical value
  struct TestParticle{
    FPoint position;
    FReal forces[3];
    FReal physicalValue;
    FReal potential;
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
    particles[idxPart].physicalValue  = physicalValue;
    particles[idxPart].potential      = 0.0;
    particles[idxPart].forces[0]      = 0.0;
    particles[idxPart].forces[1]      = 0.0;
    particles[idxPart].forces[2]      = 0.0;
  }

  ////////////////////////////////////////////////////////////////////

  { // begin direct computation

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
                                &particles[idxOther].forces[2], &particles[idxOther].potential);
        }
      }
    }
    time.tac();
    printf("Elapsed Time for direct computation: %f\n",time.elapsed());

  } // end direct computation

  ////////////////////////////////////////////////////////////////////

  {	// begin Chebyshev kernel

    // accuracy
    const unsigned int ORDER = 7;
    const FReal epsilon = FReal(1e-7);
    // typedefs
    typedef FP2PParticleContainerIndexed<> ContainerClass;
    typedef FSimpleLeaf< ContainerClass >  LeafClass;
    //typedef FInterpMatrixKernelLJ MatrixKernelClass;
    typedef FInterpMatrixKernelR MatrixKernelClass;
    typedef FChebCell<ORDER> CellClass;
    typedef FOctree<CellClass,ContainerClass,LeafClass> OctreeClass;
    //typedef FChebKernel<CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
    typedef FChebSymKernel<CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
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
      std::cout << "\nChebyshev FMM (ORDER="<< ORDER << ",EPS="<< epsilon <<") ... " << std::endl;
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

  } // end Chebyshev kernel




  //// -----------------------------------------------------
  //{ // cost of symmetric m2l opertors, weighted rank, etc.
  //	const unsigned int nnodes = ORDER*ORDER*ORDER;
  //	const SymmetryHandler<ORDER> *const SymHandler = kernels.getPtrToSymHandler();
  //	unsigned int expansionCounter[343];
  //	for (unsigned i=0; i<343; ++i) expansionCounter[i] = 0;
  //	for (unsigned i=0; i<343; ++i) if (SymHandler->pindices[i]) expansionCounter[SymHandler->pindices[i]]++;
  //
  //	unsigned int overallCost = 0;
  //	unsigned int overallWeightedRank = 0;
  //	unsigned int nbExpansions = 0;
  //	for (unsigned i=0; i<343; ++i)
  //		if (expansionCounter[i]) {
  //			const unsigned int cost = (2*nnodes*SymHandler->LowRank[i]) * expansionCounter[i];
  //			const unsigned int weightedRank = SymHandler->LowRank[i] * expansionCounter[i];
  //			overallCost += cost;
  //			overallWeightedRank += weightedRank;
  //			nbExpansions += expansionCounter[i];
  //			std::cout << "expansionCounter[" << i << "] = " << expansionCounter[i]
  //								<< "\tlow rank = " << SymHandler->LowRank[i]
  //								<< "\t(2*nnodes*rank) * nb_exp = " << cost
  //								<< std::endl;
  //		}
  //	std::cout << "=== Overall cost = " << overallCost << "\t Weighted rank = " << (double)overallWeightedRank / (double)nbExpansions << std::endl;
  //	if (nbExpansions!=316) std::cout << "Something went wrong, number of counted expansions = " << nbExpansions << std::endl;
  //}
  //// -----------------------------------------------------


  // -----------------------------------------------------
  // find first non empty leaf cell
  /*if (FParameters::findParameter(argc,argv,"-dont_check_accuracy") == FParameters::NotFound) {
    OctreeClass::Iterator iLeafs(&tree);
    iLeafs.gotoBottomLeft();

    const ContainerClass *const Targets = iLeafs.getCurrentListTargets();
    const unsigned int NumTargets = Targets->getSize();

    FReal* Potential = new FReal [NumTargets];
    FBlas::setzero(NumTargets, Potential);

    FReal* Force = new FReal [NumTargets * 3];
    FBlas::setzero(NumTargets * 3, Force);

    std::cout << "\nDirect computation of " << NumTargets << " target particles ..." << std::endl;
    const MatrixKernelClass MatrixKernel;
    do {
    const ContainerClass *const Sources = iLeafs.getCurrentListSrc();
    unsigned int counter = 0;
    ContainerClass::ConstBasicIterator iTarget(*Targets);
    while(iTarget.hasNotFinished()) {
    const FReal wt = iTarget.data().getPhysicalValue();
    ContainerClass::ConstBasicIterator iSource(*Sources);
    while(iSource.hasNotFinished()) {
    if (&iTarget.data() != &iSource.data()) {
    const FReal potential_value = MatrixKernel.evaluate(iTarget.data().getPosition(),
    iSource.data().getPosition());
    const FReal ws = iSource.data().getPhysicalValue();
    // potential
    Potential[counter] += potential_value * ws;
    // force
    if (MatrixKernelClass::Identifier == ONE_OVER_R) { // laplace force
    FPoint force(iSource.data().getPosition() - iTarget.data().getPosition());
    force *= ((ws*wt) * (potential_value*potential_value*potential_value));
    Force[counter*3 + 0] += force.getX();
    Force[counter*3 + 1] += force.getY();
    Force[counter*3 + 2] += force.getZ();
    } else if (MatrixKernelClass::Identifier == LEONARD_JONES_POTENTIAL) { // lenard-jones force
    FPoint force(iSource.data().getPosition() - iTarget.data().getPosition());
    const FReal one_over_r = FReal(1.) / FMath::Sqrt(force.getX()*force.getX() +
    force.getY()*force.getY() +
    force.getZ()*force.getZ()); // 1 + 15 + 5 = 21 flops
    const FReal one_over_r3 = one_over_r * one_over_r * one_over_r;
    const FReal one_over_r6 = one_over_r3 * one_over_r3;
    const FReal one_over_r4 = one_over_r3 * one_over_r;
    force *= ((ws*wt) * (FReal(12.)*one_over_r6*one_over_r4*one_over_r4 - FReal(6.)*one_over_r4*one_over_r4));
    Force[counter*3 + 0] += force.getX();
    Force[counter*3 + 1] += force.getY();
    Force[counter*3 + 2] += force.getZ();
    }
    }
    iSource.gotoNext();
    }
    counter++;
    iTarget.gotoNext();
    }
    } while(iLeafs.moveRight());


    FReal* ApproxPotential = new FReal [NumTargets];
    FReal* ApproxForce     = new FReal [NumTargets * 3];

    unsigned int counter = 0;
    ContainerClass::ConstBasicIterator iTarget(*Targets);
    while(iTarget.hasNotFinished()) {
    ApproxPotential[counter]   = iTarget.data().getPotential();
    ApproxForce[counter*3 + 0] = iTarget.data().getForces().getX();
    ApproxForce[counter*3 + 1] = iTarget.data().getForces().getY();
    ApproxForce[counter*3 + 2] = iTarget.data().getForces().getZ();
    counter++;
    iTarget.gotoNext();
    }

    std::cout << "\nPotential error:" << std::endl;
    std::cout << "Relative L2 error   = " << computeL2norm( NumTargets, Potential, ApproxPotential)
    << std::endl;
    std::cout << "Relative Lmax error = " << computeINFnorm(NumTargets, Potential, ApproxPotential)
    << std::endl;

    std::cout << "\nForce error:" << std::endl;
    std::cout << "Relative L2 error   = " << computeL2norm( NumTargets*3, Force, ApproxForce)
    << std::endl;
    std::cout << "Relative Lmax error = " << computeINFnorm(NumTargets*3, Force, ApproxForce)
    << std::endl;
    std::cout << std::endl;

    // free memory
    delete [] Potential;
    delete [] ApproxPotential;
    delete [] Force;
    delete [] ApproxForce;
    }*/

  /*
  // Check if particles are strictly within its containing cells
  const FReal BoxWidthLeaf = loader.getBoxWidth() / FReal(FMath::pow(2, TreeHeight-1));
  OctreeClass::Iterator octreeIterator(&tree);
  octreeIterator.gotoBottomLeft();
  do{
  const CellClass *const LeafCell = octreeIterator.getCurrentCell();
  const FPoint LeafCellCenter(LeafCell->getCoordinate().getX() * BoxWidthLeaf + BoxWidthLeaf/2 + loader.getCenterOfBox().getX(),
  LeafCell->getCoordinate().getY() * BoxWidthLeaf + BoxWidthLeaf/2 + loader.getCenterOfBox().getY(),
  LeafCell->getCoordinate().getZ() * BoxWidthLeaf + BoxWidthLeaf/2 + loader.getCenterOfBox().getZ());
  const ContainerClass *const Particles = octreeIterator.getCurrentListSrc();
  ContainerClass::ConstBasicIterator particleIterator(*Particles);
  while(particleIterator.hasNotFinished()) {
  const FPoint distance(LeafCellCenter-particleIterator.data().getPosition());
  std::cout << "center - particle = " << distance << " < " << BoxWidthLeaf/FReal(2.) << std::endl;
  if (std::abs(distance.getX())>BoxWidthLeaf/FReal(2.) ||
  std::abs(distance.getY())>BoxWidthLeaf/FReal(2.) ||
  std::abs(distance.getZ())>BoxWidthLeaf/FReal(2.)) {
  std::cout << "stop" << std::endl;
  exit(-1);
  }
  particleIterator.gotoNext();
  }
  } while(octreeIterator.moveRight());
  */


  return 0;
}
