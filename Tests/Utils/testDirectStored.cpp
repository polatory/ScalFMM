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
 *@author Cyrille Piacibello
 * 
 * This exec will be used to store the result of direct computation in
 * Binary Files. Then, a kernel (Chebyshev) is called on the same file, and
 * results are compared.
 * 
 * Format of result file :
 * Each data is a FReal : posX,posY,posZ,physicalValue,forceX,forceY,forceZ,potential
 */

// ==== CMAKE =====
// // @FUSE_BLAS
// // ================
//
//
#include <iostream>

#include <cstdio>
#include <cstdlib>
#include "ScalFmmConfig.h"

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Containers/FOctree.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Kernels/P2P/FP2P.hpp"
#include "../../Src/Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "../../Src/Files/FFmaGenericLoader.hpp"
#include "../../Src/Files/FFmaBinLoaderResult.hpp"

#include "../../Src/Utils/FTic.hpp"

#include "../../Src/Kernels/Chebyshev/FChebCell.hpp"
#include "../../Src/Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "../../Src/Kernels/Chebyshev/FChebKernel.hpp"
#include "../../Src/Kernels/Chebyshev/FChebSymKernel.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"

int main(int argc, char* argv[])
{
  // get info from commandline
  const char* const defaultFilename = (sizeof(FReal) == sizeof(float))?
    "../Data/test20k.bin.fma.single":
    "../Data/test20k.bin.fma.double";
  const char* const filename        = FParameters::getStr(argc,argv,"-f", defaultFilename);
  const char* const fileresult      = FParameters::getStr(argc,argv,"-fr","../Data/test20k.bin.fma.double.result");;
  
  //For Fmm Computation
  const unsigned int TreeHeight    = FParameters::getValue(argc, argv, "-h", 5);
  const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, "-sh", 2);
  
  FTic time;
  struct TestParticle{
    FPoint position;
    FReal forces[3];
    FReal physicalValue;
    FReal potential;
  };
  
  printf("Input is %s, \n Results will be stored in %s \n",filename,fileresult);
  
  {
    //First Part 
    //Direct Computation and Storage of result

    //open particle file
    FFmaGenericLoader loader(filename,true);
    if(!loader.isOpen()){
      std::cout << "Loader Error, " << filename << " is missing\n";
      return 1;
    }

  
    time.tic();
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

    time.tac();
    printf("Elapsed Time for Loading %lld particles from File: \t %f\n",loader.getNumberOfParticles(),time.elapsed());
  
    //Direct Computation
    time.tic();
    for(int idxTarget = 0 ; idxTarget < loader.getNumberOfParticles() ; ++idxTarget){
      for(int idxOther = idxTarget + 1 ; idxOther < loader.getNumberOfParticles() ; ++idxOther){
	FP2P::MutualParticles(particles[idxTarget].position.getX(), particles[idxTarget].position.getY(),
			      particles[idxTarget].position.getZ(),particles[idxTarget].physicalValue,
			      &particles[idxTarget].forces[0],&particles[idxTarget].forces[1],
			      &particles[idxTarget].forces[2],&particles[idxTarget].potential,
			      particles[idxOther].position.getX(), particles[idxOther].position.getY(),
			      particles[idxOther].position.getZ(),particles[idxOther].physicalValue,
			      &particles[idxOther].forces[0],&particles[idxOther].forces[1],
			      &particles[idxOther].forces[2],&particles[idxOther].potential);
      }
    }
    time.tac();
    printf("Elapsed Time for Direct Computation: \t %f\n",time.elapsed());
  
    //Write results in output file
    FILE * fd = fopen(fileresult,"w");
  
    //Size of elements to be written
    
  
    int realSize = sizeof(FReal);
    FSize nbPart   = loader.getNumberOfParticles();
    FReal boxWidth = loader.getBoxWidth() / 2;
    FReal centerX  = loader.getCenterOfBox().getX();
    FReal centerY  = loader.getCenterOfBox().getY();
    FReal centerZ  = loader.getCenterOfBox().getZ();

    fwrite(&realSize,   sizeof(int),   1, fd);
    fwrite(&nbPart,     sizeof(FSize), 1, fd);
    fwrite(&boxWidth,   sizeof(FReal), 1, fd);
    fwrite(&centerX,    sizeof(FReal), 1, fd);
    fwrite(&centerY,    sizeof(FReal), 1, fd);
    fwrite(&centerZ,    sizeof(FReal), 1, fd);
  
  
    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart)
      {
	fwrite( particles[idxPart].position.getDataValue(),realSize,3,fd);
	fwrite(&particles[idxPart].physicalValue,          realSize,1,fd);
	fwrite( particles[idxPart].forces,                 realSize,3,fd);
	fwrite(&particles[idxPart].potential,              realSize,1,fd);
      }
    
    delete [] particles; 
    fclose(fd);
  }
  
  //Second Part
  {
    //Fmm Computation and comparison to Stored results

    // begin Chebyshev kernel

    // accuracy
    const unsigned int ORDER = 7;

    // typedefs
    typedef FP2PParticleContainerIndexed<> ContainerClass;
    typedef FSimpleLeaf< ContainerClass >  LeafClass;

    typedef FInterpMatrixKernelR MatrixKernelClass;
    typedef FChebCell<ORDER> CellClass;
    typedef FOctree<CellClass,ContainerClass,LeafClass> OctreeClass;

    typedef FChebSymKernel<CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
    typedef FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;


    time.tic();
    //Open result file to fill octree for Chebyshev kernel AND compare results
    FFmaBinLoaderResult loader2(fileresult);
    
    // init oct-tree
    OctreeClass tree(TreeHeight, SubTreeHeight, loader2.getBoxWidth(), loader2.getCenterOfBox());
  
    TestParticle* const particles2 = new TestParticle[loader2.getNumberOfParticles()];  
    for(int idxPart = 0 ; idxPart < loader2.getNumberOfParticles() ; ++idxPart)
      {
	loader2.fillParticle(&particles2[idxPart].position,
			     &particles2[idxPart].physicalValue,
			     &particles2[idxPart].forces[0],
			     &particles2[idxPart].forces[1],
			     &particles2[idxPart].forces[2],
			     &particles2[idxPart].potential);
	tree.insert(FPoint(particles2[idxPart].position),idxPart,particles2[idxPart].physicalValue);
      }
    time.tac();
    printf("Elapsed Time for Reading File: \t %f\n",time.elapsed());  
    
    time.tic();
    KernelClass kernels(TreeHeight, loader2.getBoxWidth(), loader2.getCenterOfBox());
    FmmClass algorithm(&tree, &kernels);
    algorithm.execute();
    time.tac();
    printf("Elapsed Time for Fmm Computation: \t %f\n",time.elapsed());  
 

    FMath::FAccurater potentialDiff;
    FMath::FAccurater fx, fy, fz;

    //Compare the kernel to the stored elements
    {
      tree.forEachLeaf([&](LeafClass* leaf){
	  const FReal*const potentials = leaf->getTargets()->getPotentials();
	  const FReal*const forcesX = leaf->getTargets()->getForcesX();
	  const FReal*const forcesY = leaf->getTargets()->getForcesY();
	  const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
	  const int nbParticlesInLeaf = leaf->getTargets()->getNbParticles();

	  const FVector<int>& indexes = leaf->getTargets()->getIndexes();
	
	  for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
	    const int indexPartOrig = indexes[idxPart];
	    potentialDiff.add(particles2[indexPartOrig].potential,potentials[idxPart]);
	    fx.add(particles2[indexPartOrig].forces[0],forcesX[idxPart]);
	    fy.add(particles2[indexPartOrig].forces[1],forcesY[idxPart]);
	    fz.add(particles2[indexPartOrig].forces[2],forcesZ[idxPart]);
	  }
	});
      // Print for information
      std::cout << "Potential " << potentialDiff << std::endl;
      std::cout << "Fx " << fx << std::endl;
      std::cout << "Fy " << fy << std::endl;
      std::cout << "Fz " << fz << std::endl;
    
    }

   }
  return 0;
}
