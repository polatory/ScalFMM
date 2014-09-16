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

#include <cstdio>
#include <cstdlib>


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

#include "../../Src/Files/FFmaGenericLoader.hpp"
#include "../../Src/Files/FMpiFmaGenericLoader.hpp"
#include "../../Src/Files/FMpiTreeBuilder.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"
#include "../../Src/Core/FFmmAlgorithmThreadProc.hpp"

#include "../../Src/BalanceTree/FLeafBalance.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

/**
 * This program runs the FMM Algorithm Distributed with the Chebyshev kernel
 */

// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
    FHelpDescribeAndExit(argc, argv,
                         "Test with MPI the chebyshev FMM and compare it to the direct computation for debugging purpose.",
                         FParameterDefinitions::InputFile, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::NbThreads);
  
  const unsigned int ORDER = 7;
  const FReal epsilon = FReal(1e-7);

  typedef FP2PParticleContainerIndexed<> ContainerClass;
  typedef FSimpleLeaf< ContainerClass >  LeafClass;

  typedef FInterpMatrixKernelR MatrixKernelClass;
  typedef FChebCell<ORDER> CellClass;
  typedef FOctree<CellClass,ContainerClass,LeafClass> OctreeClass;

  typedef FChebSymKernel<CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
  typedef FFmmAlgorithmThreadProc<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
  
  FMpi app(argc,argv);
  
  const char* const filename       = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");
  const unsigned int TreeHeight    = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 5);
  const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeSubHeight.options, 2);
  const unsigned int NbThreads     = FParameters::getValue(argc, argv, FParameterDefinitions::NbThreads.options, 1);

  std::cout << ">> This executable has to be used to test Proc Chebyshev Algorithm. \n";
  

#ifdef _OPENMP
  omp_set_num_threads(NbThreads);
  std::cout << "\n>> Using " << omp_get_max_threads() << " threads.\n" << std::endl;
#else
  std::cout << "\n>> Sequential version.\n" << std::endl;
#endif
    
  std::cout << "Opening : " <<filename << "\n" << std::endl;
  // init timer
  FTic time;
  
  // Create Matrix Kernel
  const MatrixKernelClass MatrixKernel;

  // init particles position and physical value
  struct TestParticle{
    FPoint position;
    FReal physicalValue;
    const FPoint& getPosition(){
      return position;
    }
  };

  // open particle file
	FMpiFmaGenericLoader loader(filename,app.global());

  if(!loader.isOpen()) throw std::runtime_error("Particle file couldn't be opened!");
  
  OctreeClass tree(TreeHeight, SubTreeHeight,loader.getBoxWidth(),loader.getCenterOfBox());
  time.tic();
  TestParticle* particles = new TestParticle[loader.getNumberOfParticles()];
  memset(particles,0,(unsigned int) (sizeof(TestParticle)* loader.getNumberOfParticles()));
  for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
    loader.fillParticle(&particles[idxPart].position,&particles[idxPart].physicalValue);
  }
  FVector<TestParticle> finalParticles;
  FLeafBalance balancer;
  // FMpiTreeBuilder< TestParticle >::ArrayToTree(app.global(), particles, loader.getNumberOfParticles(),
  // 					       tree.getBoxCenter(),
  // 					       tree.getBoxWidth(),
  // 					       tree.getHeight(), &finalParticles,&balancer);
  FMpiTreeBuilder< TestParticle >::DistributeArrayToContainer(app.global(),particles, 
							      loader.getMyNumberOfParticles(),
							      tree.getBoxCenter(),
							      tree.getBoxWidth(),tree.getHeight(),
							      &finalParticles, &balancer);
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
    KernelClass kernels(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox(), &MatrixKernel, epsilon);
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
	  const FVector<int>& indexes = leaf->getTargets()->getIndexes();
    
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