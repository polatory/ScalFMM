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

#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <time.h>

#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FTic.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Utils/FAssert.hpp"
#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Components/FBasicCell.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Files/FFmaBinLoaderResult.hpp"

#include "../../Src/Components/FBasicParticleContainer.hpp"


/**
  * In this file we show an example of FBasicLoader use
* Inserting 2000000 particles ...
* Done  (5.77996).
* Deleting particles ...
* Done  (0.171918).
  */

int main(int argc, char ** argv ){
    typedef FBasicParticleContainer<1>      ContainerClass;
    typedef FSimpleLeaf< ContainerClass >                     LeafClass;
    typedef FOctree< FBasicCell, ContainerClass , LeafClass >  OctreeClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable is useless to execute.\n";
    std::cout << ">> It is only interesting to wath the code to understand\n";
    std::cout << ">> how to use the FMA loader\n";
    //////////////////////////////////////////////////////////////

    // Use testLoaderFMABinCreate.exe to create this file
    FTic counter;
    const char* const filename = FParameters::getStr(argc,argv,"-f", "../Data/test20k.bin.fma.double.res");
    std::cout << "Opening : " << filename << "\n";

    // open basic result particles loader
    FFmaBinLoaderResult loader(filename);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << "is missing\n";
        return 1;
    }

    {
        // octree
        OctreeClass tree(FParameters::getValue(argc,argv,"-h", 5), FParameters::getValue(argc,argv,"-sh", 3),
                         loader.getBoxWidth(), loader.getCenterOfBox());

        // -----------------------------------------------------
        std::cout << "Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
        counter.tic();

        FPoint particlePosition;
        FReal physicalValue = 0.0;
	FReal forces[3];
	FReal potential;
        for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
	  loader.fillParticle(&particlePosition,&physicalValue,forces,&forces[1],&forces[2],&potential);
	  tree.insert(particlePosition,physicalValue);
        }

        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << ")." << std::endl;

        // -----------------------------------------------------
        std::cout << "Deleting particles ..." << std::endl;
        counter.tic();
    }
    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << ")." << std::endl;
    // -----------------------------------------------------

    return 0;
}


