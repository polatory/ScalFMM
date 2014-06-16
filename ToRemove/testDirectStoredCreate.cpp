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
 * @brief This exec can be used to create a result file from direct
 * computation of a data file. The result file will contain forces,
 * positions, physical values, and potential for each particles.
 * 
 * This file can be loaded by FFmaBinLoaderResult.
 * 
 * Format of result file :
 * Each data is a FReal : posX,posY,posZ,physicalValue,forceX,forceY,forceZ,potential
 */


#include <iostream>

#include <cstdio>
#include <cstdlib>
#include "ScalFmmConfig.h"

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Kernels/P2P/FP2P.hpp"
#include "../../Src/Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "../../Src/Files/FFmaGenericLoader.hpp"

#include "../../Src/Utils/FTic.hpp"

int main(int argc, char* argv[])
{
  // get info from commandline
  const char* const defaultFilename = (sizeof(FReal) == sizeof(float))?
    "../Data/test20k.bin.fma.single":
    "../Data/test20k.bin.fma.double";
  const char* const filename        = FParameters::getStr(argc,argv,"-f", defaultFilename);
  const char* const fileresult      = FParameters::getStr(argc,argv,"-fr","../Data/test20k.bin.fma.double.result");;
  
  printf("Input is %s \n, Results will be stored in %s \n",filename,fileresult);
  

  //Direct Computation and Storage of result
  //open particle file
  FFmaGenericLoader loader(filename);
  if(!loader.isOpen()){
    std::cout << "Loader Error, " << filename << " is missing\n";
    return 1;
  }

  struct TestParticle{
    FPoint position;
    FReal forces[3];
    FReal physicalValue;
    FReal potential;
  };
  
  FTic time;
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
  printf("Elapsed Time for Loading File: \t %f\n",time.elapsed());
  
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

  time.tic();
  //Write results in output file
  FILE * fd = fopen(fileresult,"w");
  
  //Size of elements to be written
  size_t realSize = sizeof(FReal);
  
  fwrite(&realSize, sizeof(int),   1, fd);
  FSize nbPart = loader.getNumberOfParticles();
  FReal boxWidth = loader.getBoxWidth();
  FReal centerX = loader.getCenterOfBox().getX();
  FReal centerY = loader.getCenterOfBox().getY();
  FReal centerZ = loader.getCenterOfBox().getZ();

  fwrite(&nbPart,     sizeof(FSize), 1, fd);
  fwrite(&boxWidth,   sizeof(FReal), 1, fd);
  fwrite(&centerX,    sizeof(FReal), 1, fd);
  fwrite(&centerY,    sizeof(FReal), 1, fd);
  fwrite(&centerZ,    sizeof(FReal), 1, fd);
    
  
  for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart)
    {
      fwrite( particles[idxPart].position.getDataValue(),realSize,3,fd);
      fwrite( particles[idxPart].forces,realSize,3,fd);
      fwrite(&particles[idxPart].physicalValue,realSize,1,fd);
      fwrite(&particles[idxPart].potential,realSize,1,fd);
    }
  
  time.tac();
  printf("Elapsed Time for Writing results into file: \t %f\n",time.elapsed());

  delete [] particles; 
  fclose(fd);
  

  return 0;
}
