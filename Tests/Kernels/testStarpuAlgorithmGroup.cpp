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

#include "../../Src/Core/FFmmAlgorithmStarpuGroup.hpp"
#include "../../Src/Core/FFmmAlgorithm.hpp"

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


////////////////////////////////////////////////////////////////
// Define classes
////////////////////////////////////////////////////////////////

class TestCell : public FTestCell {
public:
    void intialCopy(const TestCell*const other){
        setDataUp( other->getDataUp() );
        setDataDown( other->getDataDown() );
    }
    void copyUp(const TestCell*const other){
        setDataUp( other->getDataUp() );
    }
    void restoreCopy(const TestCell*const other){
        setDataUp( other->getDataUp() );
        setDataDown( other->getDataDown() );
    }

};

// just to be able to load a fma file
class TestParticle : public FTestParticle, public FExtendPhysicalValue{
};

template <class ParticleClass>
class Container : public FVector<ParticleClass> {
public:
    void reduce(const Container*const other){
        for( int idx = 0 ; idx < FVector<ParticleClass>::getSize() ; ++idx){
            FVector<ParticleClass>::data()[idx].setDataDown(FVector<ParticleClass>::data()[idx].getDataDown() +
                                                            other->FVector<ParticleClass>::data()[idx].getDataDown());
        }
    }
};

////////////////////////////////////////////////////////////////
// Typedefs
////////////////////////////////////////////////////////////////
typedef TestParticle             ParticleClass;
typedef Container<ParticleClass>   ContainerClass;

typedef TestCell                CellClass;

typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;

typedef FTestKernels<ParticleClass, CellClass, ContainerClass >          KernelClass;

typedef FFmmAlgorithmStarpuGroup<OctreeClass, ParticleClass, CellClass, ContainerClass,KernelClass,LeafClass>  AlgorithmClass;

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
    const int BlockSize = FParameters::getValue(argc,argv,"-bs", 250);
    const bool usePerfModel = (FParameters::findParameter(argc, argv, "-perf") != FParameters::NotFound);
    const bool useReductionPart = (FParameters::findParameter(argc, argv, "-reducepart") != FParameters::NotFound);
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
    AlgorithmClass algo( &tree, &kernel, BlockSize, usePerfModel, useReductionPart);

    // -----------------------------------------------------

    std::cout << "Build gouped tree..." << std::endl;
    counter.tic();
    algo.buildGroups();
    counter.tac();
    std::cout << "Done  in " << counter.elapsed() << "s." << std::endl;

    // -----------------------------------------------------

    std::cout << "Execute Fmm..." << std::endl;
    counter.tic();
    algo.execute();
    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Release gouped tree..." << std::endl;
    counter.tic();
    algo.releaseGroups();
    counter.tac();
    std::cout << "Done in " << counter.elapsed() << "s." << std::endl;

    // -----------------------------------------------------

    // Check result
    ValidateFMMAlgo<OctreeClass, ParticleClass, CellClass, ContainerClass, LeafClass>(&tree);

    return 0;
}
