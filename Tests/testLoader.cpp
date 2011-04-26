// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../Src/Utils/FTic.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FList.hpp"

#include "../Src/Utils/FAssertable.hpp"
#include "../Src/Utils/F3DPosition.hpp"

#include "../Src/Components/FBasicParticle.hpp"
#include "../Src/Components/FBasicCell.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"

#include "../Src/Files/FBasicLoader.hpp"

// Compile by : g++ testLoader.cpp ../Src/Utils/FAssertable.cpp -O2 -o testLoader.exe


/**
  * In this file we show an example of FBasicLoader use
* Démarrage de /home/berenger/Dropbox/Personnel/FMB++/FMB++-build-desktop/FMB++...
* Inserting 2000000 particles ...
* Done  (5.77996).
* Deleting particles ...
* Done  (0.171918).
  */

int main(int argc, char ** argv){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable is useless to execute.\n";
    std::cout << ">> It is only interesting to wath the code to understand\n";
    std::cout << ">> how to use the basic loader.\n";
    //////////////////////////////////////////////////////////////

    // we store all particles to be able to dealloc
    FList<FBasicParticle*> particles;
    // Use testLoaderCreate.exe to create this file
    FTic counter;
    const char* const defaultFilename = "../../Data/testLoader.basic.temp";
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
    FBasicLoader<FBasicParticle> loader(filename);
    if(!loader.isValide()){
        std::cout << "Loader Error, " << filename << "is missing\n";
        return 1;
    }

    // otree
    FOctree<FBasicParticle, FBasicCell, FSimpleLeaf, 10, 3> tree(loader.getBoxWidth(),loader.getCenterOfBox());

    // -----------------------------------------------------
    std::cout << "Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    counter.tic();
    for(int idx = 0 ; idx < loader.getNumberOfParticles() ; ++idx){
        FBasicParticle* const part = new FBasicParticle();
        particles.pushFront(part);
        loader.fillParticle(part);
        tree.insert(part);
    }
    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << ")." << std::endl;

    // -----------------------------------------------------
    std::cout << "Deleting particles ..." << std::endl;
    counter.tic();
    while(particles.getSize()){
        delete particles.popFront();
    }
    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << ")." << std::endl;
    // -----------------------------------------------------

    return 0;
}


// [--LICENSE--]
