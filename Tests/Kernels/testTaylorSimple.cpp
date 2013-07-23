#include <limits>
#include <iostream>
#include <string>

#include "../../Src/Utils/FPoint.hpp"
#include "../../Src/Utils/FDebug.hpp"
#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Components/FTypedLeaf.hpp"

#include "../../Src/Containers/FVector.hpp"
#include "../../Src/Containers/FOctree.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"
#include "../../Src/Core/FFmmAlgorithmTask.hpp"

#include "../../Src/Kernels/Taylor/FTaylorCell.hpp"
#include "../../Src/Kernels/Taylor/FTaylorKernel_Opt.hpp"

#include "../../Src/Files/FFmaLoader.hpp"

int main(int argc,char* argv[]){
  const int P             = 10;
  const int NbLevels      = FParameters::getValue(argc, argv, "-h",  3);
  const int SizeSubLevels = FParameters::getValue(argc, argv, "-sh", 2);
  static const int order = 1;
  FPoint rootCenter(FReal(0.0),FReal(0.0),FReal(0.0));
  FReal boxWidth = FReal(8);

  typedef FTaylorCell<P,order> CellClass;
  typedef FP2PParticleContainer ContainerClass;
  
  typedef FSimpleLeaf< ContainerClass > LeafClass;
  typedef FOctree< CellClass, ContainerClass , LeafClass > OctreeClass;
  typedef FTaylorKernel<CellClass,ContainerClass,P,order> KernelClass;
  
  typedef FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
  typedef FFmmAlgorithmThread<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClassThread;
  typedef FFmmAlgorithmTask<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClassTask;


  FTic counter;
  
  OctreeClass tree(NbLevels, SizeSubLevels, boxWidth, rootCenter);

  FPoint part1Pos = FPoint(FReal(3.25),FReal(0.75),FReal(0.75));
  FReal physVal1 = -1;
  FPoint part3Pos = FPoint(FReal(2.25),FReal(1.75),FReal(0.75));
  FReal physVal3 = -1;

  FPoint part2Pos = FPoint(FReal(-3.25),FReal(1.75),FReal(0.5));
  FReal physVal2 = 2;
  
  tree.insert(part1Pos,physVal1);
  tree.insert(part3Pos,physVal3);
  tree.insert(part2Pos,physVal2);


  KernelClass kernels(NbLevels, boxWidth, rootCenter);


  FmmClass algo(&tree,&kernels);
  algo.execute();
  

 
  { // get sum forces&potential
    FReal Energy = 0.0;
    
    
    tree.forEachLeaf([&](LeafClass* leaf){
	FReal fx = 0.0, fy = 0.0, fz = 0.0;
	const FReal * FRestrict charges = leaf->getTargets()->getPhysicalValues();
	const FReal*const potentials = leaf->getTargets()->getPotentials();
	const FReal*const forcesX = leaf->getTargets()->getForcesX();
	const FReal*const forcesY = leaf->getTargets()->getForcesY();
	const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
	const int nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
	
	for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
	  Energy   += charges[idxPart]*potentials[idxPart];
	  fx        += forcesX[idxPart];
	  fy        += forcesY[idxPart];
	  fz        += forcesZ[idxPart];
	  printf(" part : %f, p= %f fx = %f, fy = %f, fz = %f\n",leaf->getTargets()->getPositions()[0][0],Energy,fx,fy,fz);
	}
	
      }); 
    //
    FReal potTheo ;
    //
    FReal dx = part1Pos.getX() - part2Pos.getX();
    FReal dy = part1Pos.getY() - part2Pos.getY();
    FReal dz = part1Pos.getZ() - part2Pos.getZ();
    FReal dx1 = part3Pos.getX() - part2Pos.getX();
    FReal dy1 = part3Pos.getY() - part2Pos.getY();
    FReal dz1 = part3Pos.getZ() - part2Pos.getZ();
    FReal dx2 = part1Pos.getX() - part3Pos.getX();
    FReal dy2 = part1Pos.getY() - part3Pos.getY();
    FReal dz2 = part1Pos.getZ() - part3Pos.getZ();
    std::cout << dx  <<"  " <<dy  <<"  "<< dz <<"  " <<std::endl;
    std::cout << dx1 <<"  " <<dy1 <<"  "<< dz1 <<"  " <<std::endl;
    std::cout << dx2 <<"  " <<dy2 <<"  "<< dz2 <<"  " <<std::endl;
    //
    potTheo = physVal2*physVal1 / FMath::Sqrt(dx*dx + dy*dy + dz*dz)+ physVal2*physVal3 / FMath::Sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1) + physVal1*physVal3/ FMath::Sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
    FReal coeffa   = physVal2*physVal1*(FMath::Sqrt(FReal(1.0) / (dx*dx + dy*dy + dz*dz)) ) / (dx*dx + dy*dy + dz*dz);
    Energy *=FReal(0.5) ;
    printf("Exact potential : %f     Computed  potential : %f  Error: %e \n",potTheo, Energy,std::abs(potTheo- Energy));
    printf("Exact Force : %f  %f : %f  \n",dx*coeffa,dy*coeffa,dz*coeffa);
    std::cout << Energy/potTheo << std::endl;
  }
  return 0;
}
