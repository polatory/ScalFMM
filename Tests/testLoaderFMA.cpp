// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../Src/Utils/FParameters.hpp"
#include "../Src/Utils/FTic.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Utils/FAssertable.hpp"
#include "../Src/Utils/F3DPosition.hpp"

#include "../Src/Components/FFmaParticle.hpp"
#include "../Src/Components/FBasicCell.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"

#include "../Src/Files/FFmaLoader.hpp"

// Compile by : g++ testLoaderFMA.cpp ../Src/Utils/FDebug.cpp ../Src/Utils/FTrace.cpp -O2 -o testLoaderFMA.exe


/**
  * In this file we show an example of FBasicLoader use
* Démarrage de /home/berenger/Dropbox/Personnel/FMB++/FMB++-build-desktop/FMB++...
* Inserting 2000000 particles ...
* Done  (5.77996).
* Deleting particles ...
* Done  (0.171918).
  */

int main(int argc, char ** argv ){
    typedef FVector<FFmaParticle>      ContainerClass;
    typedef FSimpleLeaf<FFmaParticle, ContainerClass >                     LeafClass;
    typedef FOctree<FFmaParticle, FBasicCell, ContainerClass , LeafClass >  OctreeClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable is useless to execute.\n";
    std::cout << ">> It is only interesting to wath the code to understand\n";
    std::cout << ">> how to use the FMA loader\n";
    //////////////////////////////////////////////////////////////

    // Use testLoaderCreate.exe to create this file
    FTic counter;
    const char* const defaultFilename = "../../Data/testLoaderFMA.fma";
    const char* filename;

    if(argc == 1){
        std::cout << "You have to give a .fma file in argument.\n";
        std::cout << "The program will try a default file : " << defaultFilename << "\n";
        filename = defaultFilename;
    }
    else{
        filename = argv[1];
        std::cout << "Opening : " << filename << "\n";
    }

    // open basic particles loader
    FFmaLoader<FFmaParticle> loader(filename);
    if(!loader.hasNotFinished()){
        std::cout << "Loader Error, " << filename << "is missing\n";
        return 1;
    }

    {
        // otree
        OctreeClass tree(10, 3,loader.getBoxWidth(),loader.getCenterOfBox());

        // -----------------------------------------------------
        std::cout << "Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
        counter.tic();

        FFmaParticle part;
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


// [--LICENSE--]
