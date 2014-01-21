// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
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

#ifndef FFMABINLOADERRESULT_HPP
#define FFMABINLOADERRESULT_HPP

#include "./FFmaBinLoader.hpp"

/**
 * @author Cyrille Piacibello (cyrille.piacibello@inria.fr)
 * @class FFmaBinLoaderResult
 *
 * Please read the license
 *
 * @brief The aim of this class is to load particles and forces from a
 * file, to compare with results from any Fmm Kernel.
 */

class FFmaBinLoaderResult : public  FFmaBinLoader {

protected:

  typedef FFmaBinLoader Parent;
  
  size_t removeWarning;
  
public:

  FFmaBinLoaderResult(const char * filename):
    Parent::FFmaBinLoader(filename)
  {
  }

  void fillParticle(FPoint*const inParticlePosition, FReal*const physicalValue, FReal* forceX, FReal* forceY, FReal* forceZ, FReal * const potential){
    FReal x,y,z,data,fx,fy,fz,pot;
    
    //Same datas as particle files
    removeWarning += fread(&x,    sizeof(FReal), 1, file);
    removeWarning += fread(&y,    sizeof(FReal), 1, file);
    removeWarning += fread(&z,    sizeof(FReal), 1, file);
    removeWarning += fread(&data, sizeof(FReal), 1, file);

    //results datas
    removeWarning += fread(&fx,  sizeof(FReal), 1, file);
    removeWarning += fread(&fy,  sizeof(FReal), 1, file);
    removeWarning += fread(&fz,  sizeof(FReal), 1, file);
    removeWarning += fread(&pot, sizeof(FReal), 1, file);
    
    //return
    inParticlePosition->setPosition(x,y,z);
    (*physicalValue) = data;
    *forceX = fx;
    *forceY = fy;
    *forceZ = fz;
    *potential = pot;
  }
};


#endif //FFMABINLOADERRESULT_HPP
