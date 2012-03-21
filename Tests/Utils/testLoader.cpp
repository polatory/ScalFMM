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

#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <time.h>

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Utils/FAssertable.hpp"
#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Components/FBasicParticle.hpp"
#include "../../Src/Components/FBasicCell.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Files/FBasicLoader.hpp"

/**
  * In this file we show an example of FBasicLoader use
* Inserting 2000000 particles ...
* Done  (5.77996).
* Deleting particles ...
* Done  (0.171918).
  */

int main(int argc, char ** argv){
    typedef FVector<FBasicParticle>      ContainerClass;
    typedef FSimpleLeaf<FBasicParticle, ContainerClass >                     LeafClass;
    typedef FOctree<FBasicParticle, FBasicCell, ContainerClass , LeafClass >  OctreeClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable is useless to execute.\n";
    std::cout << ">> It is only interesting to wath the code to understand\n";
    std::cout << ">> how to use the basic loader.\n";
    //////////////////////////////////////////////////////////////

    // Use testLoaderCreate.exe to create this file
    FTic counter;

    const char* const filename = FParameters::getStr(argc,argv,"-f", "../Data/test20k.basic");
    std::cout << "Opening : " << filename << "\n";

    // open basic particles loader
    FBasicLoader<FBasicParticle> loader(filename);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << "is missing\n";
        return 1;
    }

    {
        // otree
        OctreeClass tree(FParameters::getValue(argc,argv,"-h", 5), FParameters::getValue(argc,argv,"-sh", 3),
                         loader.getBoxWidth(), loader.getCenterOfBox());

        // -----------------------------------------------------
        std::cout << "Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
        counter.tic();

        loader.fillTree(tree);

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


