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

#include "../Sources/Core/FBasicParticule.hpp"
#include "../Sources/Core/FBasicCell.hpp"

#include "../Sources/Files/FFMALoader.hpp"

// Compile by : g++ testLoaderFMA.cpp ../Sources/Utils/FAssertable.cpp -O2 -o testLoaderFMA.exe


/**
  * In this file we show an example of BasicParticule and BasicLoader use
* Démarrage de /home/berenger/Dropbox/Personnel/FMB++/FMB++-build-desktop/FMB++...
* Inserting 2000000 particules ...
* Done  (5.77996).
* Deleting particules ...
* Done  (0.171918).
  */

int main(int , char ** ){
    // we store all particules to be able to dealloc
    FList<FBasicParticule*> particules;
    // Use testLoaderCreate.exe to create this file
    const char* const filename = "testLoaderFMA.fma";
    FTic counter;

    // open basic particules loader
    FFMALoader<FBasicParticule> loader(filename);
    if(!loader.isValide()){
        std::cout << "Loader Error, " << filename << "is missing\n";
        return 1;
    }

    // otree
    FOctree<FBasicParticule, FBasicCell, 10, 3> tree(loader.getBoxWidth(),loader.getCenterOfBox());

    // -----------------------------------------------------
    std::cout << "Inserting " << loader.getNumberOfParticules() << " particules ..." << std::endl;
    counter.tic();
    for(int idx = 0 ; idx < loader.getNumberOfParticules() ; ++idx){
        FBasicParticule* const part = new FBasicParticule();
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
