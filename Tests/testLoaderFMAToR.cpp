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

#include "../Sources/Extenssions/FExtendParticuleType.hpp"

#include "../Sources/Core/FSimpleLeaf.hpp"

#include "../Sources/Files/FFMAToRLoader.hpp"

// Compile by : g++ testLoaderFMAToR.cpp ../Sources/Utils/FAssertable.cpp -O2 -o testLoaderFMAToR.exe


class ParticuleToR : public FFmaParticule, public FExtendParticuleType {
};

/**
  * In this file we show an example of FFmaLoader use
  */

int main(int , char ** ){
    // we store all particules to be able to dealloc
    FList<ParticuleToR*> particules;
    // Use testLoaderCreate.exe to create this file
    const char* const filename = "testLoaderFMA.tor.fma";
    FTic counter;

    // open basic particules loader
    FFMAToRLoader<ParticuleToR> loader(filename);
    if(!loader.isValide()){
        std::cout << "Loader Error, " << filename << "is missing\n";
        return 1;
    }

    // otree
    FOctree<ParticuleToR, FBasicCell, FSimpleLeaf, 10, 3> tree(loader.getBoxWidth(),loader.getCenterOfBox());

    // -----------------------------------------------------
    std::cout << "Inserting " << loader.getNumberOfParticules() << " particules ..." << std::endl;
    counter.tic();
    for(int idx = 0 ; idx < loader.getNumberOfParticules() ; ++idx){
        ParticuleToR* const part = new ParticuleToR();
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
