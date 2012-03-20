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
// @FUSE_STARPU
// ================

#include <starpu.h>


#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Components/FTestKernels.hpp"
#include "../../Src/Components/FTestParticle.hpp"
#include "../../Src/Components/FTestCell.hpp"

#include "../../Src/Core/FFmmAlgorithmStarpu.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Components/FFmaParticle.hpp"
#include "../../Src/Extensions/FExtendForces.hpp"
#include "../../Src/Extensions/FExtendPotential.hpp"

#include "../../Src/Components/FBasicCell.hpp"

#include "../../Src/Files/FFmaLoader.hpp"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


// export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/
// Compile With openmp : g++ testFmbAlgorithm.cpp ../../Src/Utils/FDebug.cpp ../../Src/Utils/FTrace.cpp -lgomp -fopenmp -lstarpu -O2 -o testFmbAlgorithm.exe
//
// g++ -L../starpu/lib/ -I../starpu/include testFmbAlgorithmNoProc.cpp ../../Src/Utils/FDebug.cpp ../../Src/Utils/FTrace.cpp ../../Src/Utils/FMath.cpp ../../Src/Utils/F3DPosition.cpp -lgomp -fopenmp -lstarpu -O2 -o testFmbAlgorithm.exe

////////////////////////////////////////////////////////////////
// Define classes
////////////////////////////////////////////////////////////////

// just to be able to load a fma file
class TestParticle : public FTestParticle, public FExtendPhysicalValue{
};

////////////////////////////////////////////////////////////////
// Typedefs
////////////////////////////////////////////////////////////////
typedef TestParticle             ParticleClass;
typedef StarVector<ParticleClass> ContainerClass;
typedef DataVector<ParticleClass> RealContainerClass;

typedef FTestCell                RealCellClass;
typedef FStarCell<RealCellClass> CellClass;


typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;

typedef FTestKernels<ParticleClass, RealCellClass, RealContainerClass >          KernelClass;

typedef FFmmAlgorithmStarpu<OctreeClass, ParticleClass, CellClass, RealCellClass, ContainerClass,KernelClass,LeafClass>  AlgorithmClass;

////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test fmb algorithm.\n";
    //////////////////////////////////////////////////////////////
    const int NbLevels = FParameters::getValue(argc,argv,"-h", 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    FTic counter;
    const char* const filename = FParameters::getStr(argc,argv,"-f", "../Data/test20k.fma");

    std::cout << "Opening : " << filename << "\n";

    FFmaLoader<ParticleClass> loader(filename);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // -----------------------------------------------------
    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    loader.fillTree(tree);

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    KernelClass kernel;
    AlgorithmClass algo( &tree, &kernel);
    std::cout << "There are " << starpu_worker_get_count() << " workers" << std::endl;
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    // Check result
    ValidateFMMAlgo<OctreeClass, ParticleClass, CellClass, ContainerClass, LeafClass>(&tree);

    return 0;
}
