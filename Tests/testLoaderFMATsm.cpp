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

#include "../Sources/Components/FFmaParticule.hpp"
#include "../Sources/Components/FBasicCell.hpp"

#include "../Sources/Extenssions/FExtendParticuleType.hpp"

#include "../Sources/Components/FSimpleLeaf.hpp"

#include "../Sources/Files/FFMATsmLoader.hpp"

// Compile by : g++ testLoaderFMATsm.cpp ../Sources/Utils/FAssertable.cpp -O2 -o testLoaderFMATsm.exe


class ParticuleTsm : public FFmaParticule, public FExtendParticuleType {
};

/**
  * In this file we show an example of FFmaLoader use
  */

int main(int argc, char ** argv ){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable is useless to execute.\n";
    std::cout << ">> It is only interesting to wath the code to understand\n";
    std::cout << ">> how to use the Tsm loader\n";
    //////////////////////////////////////////////////////////////

    // we store all particules to be able to dealloc
    FList<ParticuleTsm*> particules;
    // Use testLoaderCreate.exe to create this file
    FTic counter;
    const char* const defaultFilename = "testLoaderFMA.tsm.fma";
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

    // open basic particules loader
    FFMATsmLoader<ParticuleTsm> loader(filename);
    if(!loader.isValide()){
        std::cout << "Loader Error, " << filename << "is missing\n";
        return 1;
    }

    // otree
    FOctree<ParticuleTsm, FBasicCell, FSimpleLeaf, 10, 3> tree(loader.getBoxWidth(),loader.getCenterOfBox());

    // -----------------------------------------------------
    std::cout << "Inserting " << loader.getNumberOfParticules() << " particules ..." << std::endl;
    counter.tic();
    for(int idx = 0 ; idx < loader.getNumberOfParticules() ; ++idx){
        ParticuleTsm* const part = new ParticuleTsm();
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
