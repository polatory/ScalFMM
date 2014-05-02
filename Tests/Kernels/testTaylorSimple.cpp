#include <limits>
#include <iostream>
#include <string>

#include "../../Src/Utils/FPoint.hpp"
#include "../../Src/Utils/FLog.hpp"
#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Components/FTypedLeaf.hpp"

#include "../../Src/Containers/FVector.hpp"
#include "../../Src/Containers/FOctree.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
//#include "../../Src/Core/FFmmAlgorithmThread.hpp"
//#include "../../Src/Core/FFmmAlgorithmTask.hpp"

#include "../../Src/Kernels/Taylor/FTaylorCell.hpp"
#include "../../Src/Kernels/Taylor/FTaylorKernel.hpp"


//#include "../../Src/Files/FFmaGenericLoader.hpp"

 int main(int argc,char* argv[]){
   const int P             = 3;
   const int NbLevels      = FParameters::getValue(argc, argv, "-depth",  4);
   const int SizeSubLevels = FParameters::getValue(argc, argv, "-subdepth", 2);
   static const int order = 1;
   FPoint rootCenter(FReal(0.0),FReal(0.0),FReal(0.0));
   FReal boxWidth = FReal(4);

   typedef FTaylorCell<P,order> CellClass;
   typedef FP2PParticleContainer<> ContainerClass;

   typedef FSimpleLeaf< ContainerClass > LeafClass;
   typedef FOctree< CellClass, ContainerClass , LeafClass > OctreeClass;
   typedef FTaylorKernel<CellClass,ContainerClass,P,order> KernelClass;

   typedef FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
//   typedef FFmmAlgorithmThread<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClassThread;
 //  typedef FFmmAlgorithmTask<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClassTask;


   FTic counter;

   OctreeClass tree(NbLevels, SizeSubLevels, boxWidth, rootCenter);

   int nbPart = 3;

   FPoint tabPart[nbPart];
   tabPart[0] = FPoint(FReal(-1.9),FReal(-1.9),FReal(1.6));
   tabPart[1] = FPoint(FReal(0.9),FReal(1.1),FReal(1.6));
   tabPart[2] = FPoint(FReal(1.9),FReal(0.1),FReal(1.6));
   //tabPart[2] = FPoint(FReal(0.1),FReal(1.9),FReal(1.6));
   //tabPart[4] = FPoint(FReal(0.1),FReal(0.1),FReal(1.6));

   FReal tabPhyValue[nbPart];
   tabPhyValue[0] = FReal(1);
   tabPhyValue[1] = FReal(1);
   tabPhyValue[2] = FReal(1);
   // tabPhyValue[3] = FReal(1);
   //tabPhyValue[4] = FReal(1);

   for(int l=0 ; l<nbPart ; ++l){
     tree.insert(tabPart[l],tabPhyValue[l]);
   }

   KernelClass kernels(NbLevels, boxWidth, rootCenter);
   FReal tabResDirect[nbPart*3];

   for(int r=0 ; r<nbPart ; ++r){
     tabResDirect[r*3]   = FReal(0);
     tabResDirect[r*3+1] = FReal(0);
     tabResDirect[r*3+2] = FReal(0);
   }

   FReal potTheoric = FReal(0.0);

   //iteration over particles src
   for(int t=0 ; t<nbPart ; ++t)
     {   
       for(int u=t+1 ; u<nbPart ; ++u)
	 {
	   //Calcul of interaction between particules
	   FReal dx = tabPart[t].getX()-tabPart[u].getX();
	   FReal dy = tabPart[t].getY()-tabPart[u].getY();
	   FReal dz = tabPart[t].getZ()-tabPart[u].getZ();

	   FReal dist =  FMath::Sqrt(dx*dx+dy*dy+dz*dz);
	   FReal dist2 = dx*dx+dy*dy+dz*dz;

	   FReal fx = tabPhyValue[t]*tabPhyValue[u]*dx/(dist2*dist);
	   FReal fy = tabPhyValue[t]*tabPhyValue[u]*dy/(dist2*dist);
	   FReal fz = tabPhyValue[t]*tabPhyValue[u]*dz/(dist2*dist);


	   //Computation of force on part[t]
	   tabResDirect[t*3]   += fx;
	   tabResDirect[t*3+1] += fy;
	   tabResDirect[t*3+2] += fz;

	   //Computation of force on part[u]
	   tabResDirect[u*3]   += -fx;
	   tabResDirect[u*3+1] += -fy;
	   tabResDirect[u*3+2] += -fz;

	   //Computation of potential
	   potTheoric += tabPhyValue[t]*tabPhyValue[u]/dist;
	 }
     }


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
	   printf(" part : (%f,%F,%F), p= %f fx = %f, fy = %f, fz = %f\n",leaf->getTargets()->getPositions()[0][0],leaf->getTargets()->getPositions()[1][0],leaf->getTargets()->getPositions()[2][0],Energy,fx,fy,fz);
	 }

       }); 

     FReal dx =tabPart[0].getX() -tabPart[1].getX();
     FReal dy =tabPart[0].getY() -tabPart[1].getY();
     FReal dz =tabPart[0].getZ() -tabPart[1].getZ();
     FReal dx1 =tabPart[2].getX() -tabPart[1].getX();
     FReal dy1 =tabPart[2].getY() -tabPart[1].getY();
     FReal dz1 =tabPart[2].getZ() -tabPart[1].getZ();
     FReal dx2 =tabPart[0].getX() -tabPart[2].getX();
     FReal dy2 =tabPart[0].getY() -tabPart[2].getY();
     FReal dz2 =tabPart[0].getZ() -tabPart[2].getZ();
     std::cout << dx  <<"  " <<dy  <<"  "<< dz <<"  " <<std::endl;
     std::cout << dx1 <<"  " <<dy1 <<"  "<< dz1 <<"  " <<std::endl;
     std::cout << dx2 <<"  " <<dy2 <<"  "<< dz2 <<"  " <<std::endl;
     
     Energy *=FReal(0.5) ;
     printf("Exact potential : %f     Computed  potential : %f  Error: %e \n",potTheoric, Energy,std::abs(potTheoric- Energy));
     //printf("Exact Force : %f  %f : %f  \n",dx*coeffa,dy*coeffa,dz*coeffa);
     //std::cout << Energy/potTheo << std::endl;
     for (int j =0; j<nbPart ; ++j)
       {
	 printf("particule : (%f,%f,%f) fx = %f, fy = %f, fz = %f\n",tabPart[j].getX(),tabPart[j].getY(),tabPart[j].getZ(),
		tabResDirect[j*3],tabResDirect[j*3+1],tabResDirect[j*3+2]);
       }
     
   }
  return 0;
 }
