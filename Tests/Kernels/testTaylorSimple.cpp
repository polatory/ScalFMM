#include <limits>
#include <iostream>
#include <string>

#include "../../Src/Utils/FPoint.hpp"
#include "../../Src/Utils/FDebug.hpp"
#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Containers/FVector.hpp"
#include "../../Src/Containers/FOctree.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"
#include "../../Src/Core/FFmmAlgorithmTask.hpp"

#include "../../Src/Kernels/Taylor/FTaylorCell.hpp"
#include "../../Src/Kernels/Taylor/FTaylorKernel.hpp"

#include "../../Src/Files/FFmaLoader.hpp"

int main(int argc,char* argv[]){
  static const int P = 3;
  static const int order = 1;
  FPoint rootCenter(FReal(0),FReal(0),FReal(0));
  FReal boxWidth = FReal(8);

  typedef FTaylorCell<P,order> CellClass;
  typedef FP2PParticleContainer ContainerClass;
  
  typedef FSimpleLeaf< ContainerClass > LeafClass;
  typedef FOctree< CellClass, ContainerClass , LeafClass > OctreeClass;
  typedef FTaylorKernel<CellClass,ContainerClass,P,order> KernelClass;
  
  typedef FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
  typedef FFmmAlgorithmThread<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClassThread;
  typedef FFmmAlgorithmTask<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClassTask;


  const int NbLevels = 3;
  const int SizeSubLevels = 1;
  FTic counter;
  
  OctreeClass tree(NbLevels, SizeSubLevels, boxWidth, rootCenter);
  
  FPoint part1Pos = FPoint(FReal(3.75),FReal(0.25),FReal(0.25));
  FReal physVal1 = 1;

  FPoint part2Pos = FPoint(FReal(-3.75),FReal(0.25),FReal(0.25));
  FReal physVal2 = 100;
  
  tree.insert(part1Pos,physVal1);
  tree.insert(part2Pos,physVal2);

  KernelClass kernels(NbLevels, boxWidth, rootCenter);

  FmmClassThread algo(&tree,&kernels);
  algo.execute();
  
  { // get sum forces&potential
    FReal potential = 0;
    
    
    tree.forEachLeaf([&](LeafClass* leaf){
	FReal fx = 0.0, fy = 0.0, fz = 0.0;
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
	  printf("part : %f, fx = %f, fy = %f, fz = %f\n",leaf->getTargets()->getPositions()[0][0],fx,fy,fz);
	}
	
      }); 
    
  }

  return 0;
}