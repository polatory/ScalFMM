// [--License--]

#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <time.h>

#include "../Src/Utils/FTic.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Utils/FAssertable.hpp"
#include "../Src/Utils/F3DPosition.hpp"

#include "../Src/Components/FFmaParticle.hpp"
#include "../Src/Components/FBasicCell.hpp"

#include "../Src/Extensions/FExtendParticleType.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"

#include "../Src/Files/FFmaTsmLoader.hpp"


#include "../Src/Utils/FParameters.hpp"

/**
  * In this file we show an example of FFmaLoader use
  */

class ParticleTsm : public FFmaParticle, public FExtendParticleType {
};



int main(int argc, char ** argv ){
    typedef FVector<ParticleTsm>      ContainerClass;
    typedef FSimpleLeaf<ParticleTsm, ContainerClass >                     LeafClass;
    typedef FOctree<ParticleTsm, FBasicCell, ContainerClass , LeafClass >  OctreeClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable is useless to execute.\n";
    std::cout << ">> It is only interesting to wath the code to understand\n";
    std::cout << ">> how to use the Tsm loader\n";
    //////////////////////////////////////////////////////////////

    // Use testLoaderCreate.exe to create this file
    FTic counter;
    const char* const defaultFilename = "../Data/test20k.tsm.fma";
    const char* filename;

    if(argc == 1){
        std::cout << "You have to give a .tos.fma file in argument.\n";
        std::cout << "The program will try a default file : " << defaultFilename << "\n";
        filename = defaultFilename;
    }
    else{
        filename = argv[1];
        std::cout << "Opening : " << filename << "\n";
    }

    // open basic particles loader
    FFmaTsmLoader<ParticleTsm> loader(filename);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << "is missing\n";
        return 1;
    }
    {
        // otree
        OctreeClass tree(FParameters::getValue(argc,argv,"-h", 5), FParameters::getValue(argc,argv,"-sh", 3),
                         loader.getBoxWidth(),loader.getCenterOfBox());

        // -----------------------------------------------------
        std::cout << "Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
        counter.tic();

        ParticleTsm part;
        for(int idx = 0 ; idx < loader.getNumberOfParticles() ; ++idx){
            loader.fillParticle(part);
            tree.insert(part);
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


// [--END--]
