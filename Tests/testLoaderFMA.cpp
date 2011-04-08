// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../Sources/Utils/FTic.hpp"

#include "../Sources/Containers/FOctree.hpp"
#include "../Sources/Containers/FList.hpp"

#include "../Sources/Utils/FAssertable.hpp"
#include "../Sources/Utils/F3DPosition.hpp"

#include "../Sources/Core/FFmaParticule.hpp"
#include "../Sources/Core/FBasicCell.hpp"

#include "../Sources/Core/FSimpleLeaf.hpp"

#include "../Sources/Files/FFMALoader.hpp"

// Compile by : g++ testLoaderFMA.cpp ../Sources/Utils/FAssertable.cpp -O2 -o testLoaderFMA.exe


/**
  * In this file we show an example of FBasicLoader use
* Démarrage de /home/berenger/Dropbox/Personnel/FMB++/FMB++-build-desktop/FMB++...
* Inserting 2000000 particules ...
* Done  (5.77996).
* Deleting particules ...
* Done  (0.171918).
  */

int main(int argc, char ** argv ){
    // we store all particules to be able to dealloc
    FList<FFmaParticule*> particules;
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

    // open basic particules loader
    FFMALoader<FFmaParticule> loader(filename);
    if(!loader.isValide()){
        std::cout << "Loader Error, " << filename << "is missing\n";
        return 1;
    }

    // otree
    FOctree<FFmaParticule, FBasicCell, FSimpleLeaf, 10, 3> tree(loader.getBoxWidth(),loader.getCenterOfBox());

    // -----------------------------------------------------
    std::cout << "Inserting " << loader.getNumberOfParticules() << " particules ..." << std::endl;
    counter.tic();
    for(int idx = 0 ; idx < loader.getNumberOfParticules() ; ++idx){
        FFmaParticule* const part = new FFmaParticule();
        particules.pushFront(part);
        loader.fillParticule(part);
        tree.insert(part);
    }
    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << ")." << std::endl;

    // -----------------------------------------------------
    std::cout << "Deleting particules ..." << std::endl;
    counter.tic();
    while(particules.getSize()){
        delete particules.popFront();
    }
    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << ")." << std::endl;
    // -----------------------------------------------------

    return 0;
}


// [--LICENSE--]
