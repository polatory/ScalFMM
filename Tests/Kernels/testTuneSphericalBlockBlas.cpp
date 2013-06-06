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

// ==== CMAKE =====
// @FUSE_BLAS
// ================

#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Components/FBasicCell.hpp"

#include "../../Src/Kernels/Spherical/FSphericalBlockBlasKernel.hpp"
#include "../../Src/Kernels/Spherical/FSphericalCell.hpp"

#include "../../Src/Files/FFmaScanfLoader.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

/** This program find the best block blas size
  */



// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef FSphericalCell                 CellClass;
    typedef FP2PParticleContainer         ContainerClass;

    typedef FSimpleLeaf< ContainerClass >                     LeafClass;
    typedef FOctree< CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FSphericalBlockBlasKernel< CellClass, ContainerClass > KernelClass;

    typedef FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test Spherical algorithm.\n";
    //////////////////////////////////////////////////////////////
    const int DevP = FParameters::getValue(argc,argv,"-p", 8);
    const int NbLevels = FParameters::getValue(argc,argv,"-h", 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    const int MaxBlockSize = FParameters::getValue(argc,argv,"-mbz", 10000);
    FTic counter;

    const char* const filename = FParameters::getStr(argc,argv,"-f", "../Data/test20k.fma");
    std::cout << "Opening : " << filename << "\n";

    // -----------------------------------------------------
    double minTime = 999999999999999999.0;
    int bestSize = -1;

    CellClass::Init(DevP, true);
    for(int idxBlockSize = 1 ; idxBlockSize < MaxBlockSize ; idxBlockSize *= 2){
        FFmaScanfLoader loader(filename);
        if(!loader.isOpen()){
            std::cout << "Loader Error, " << filename << " is missing\n";
            return 1;
        }

        OctreeClass tree(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());

        // -----------------------------------------------------

        for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            FPoint particlePosition;
            FReal physicalValue = 0.0;
            loader.fillParticle(&particlePosition,&physicalValue);
            tree.insert(particlePosition, physicalValue );
        }

        // -----------------------------------------------------

        KernelClass kernels(DevP, NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), idxBlockSize);
        FmmClass algo(&tree,&kernels);

        // -----------------------------------------------------

        counter.tic();

        algo.execute();

        counter.tac();
        std::cout << "For block size = " << idxBlockSize << ", time is = " << counter.elapsed() << "s" << std::endl;

        // update best
        if(counter.elapsed() < minTime){
            minTime = counter.elapsed();
            bestSize = idxBlockSize;
        }
    }

    std::cout << "Best block size is = " << bestSize << std::endl;

    return 0;
}



