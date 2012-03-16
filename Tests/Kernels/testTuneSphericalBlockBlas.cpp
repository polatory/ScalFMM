
// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
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
#include "../../Src/Kernels/Spherical/FSphericalParticle.hpp"
#include "../../Src/Kernels/Spherical/FSphericalCell.hpp"

#include "../../Src/Files/FFmaScanfLoader.hpp"

/** This program find the best block blas size
  */



// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef FSphericalParticle             ParticleClass;
    typedef FSphericalCell                 CellClass;
    typedef FVector<ParticleClass>         ContainerClass;

    typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
    typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FSphericalBlockBlasKernel<ParticleClass, CellClass, ContainerClass > KernelClass;

    typedef FFmmAlgorithm<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
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
        FFmaScanfLoader<ParticleClass> loader(filename);
        if(!loader.isOpen()){
            std::cout << "Loader Error, " << filename << " is missing\n";
            return 1;
        }

        OctreeClass tree(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());

        // -----------------------------------------------------

        loader.fillTree(tree);

        // -----------------------------------------------------
        counter.tic();

        KernelClass kernels(DevP, NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), idxBlockSize);

        FmmClass algo(&tree,&kernels);
        algo.execute();

        counter.tac();
        std::cout << "For block size = " << idxBlockSize << ", time is = " << counter.elapsed() << "s" << std::endl;

        if(counter.elapsed() < minTime){
            minTime = counter.elapsed();
            bestSize = idxBlockSize;
        }
    }

    std::cout << "Best block size is = " << bestSize << std::endl;

    return 0;
}



