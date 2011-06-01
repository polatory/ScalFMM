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

#include "../Src/Components/FFmaParticle.hpp"
#include "../Src/Components/FBasicCell.hpp"

#include "../Src/Extensions/FExtendParticleType.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"

#include "../Src/Files/FFmaTsmLoader.hpp"

// Compile by : g++ testLoaderFMATsm.cpp ../Src/Utils/FAssertable.cpp -O2 -o testLoaderFMATsm.exe


class ParticleTsm : public FFmaParticle, public FExtendParticleType {
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

    // we store all particles to be able to dealloc
    FList<ParticleTsm*> particles;
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

    // open basic particles loader
    FFmaTsmLoader<ParticleTsm> loader(filename);
    if(!loader.hasNotFinished()){
        std::cout << "Loader Error, " << filename << "is missing\n";
        return 1;
    }
    {
        // otree
        FOctree<ParticleTsm, FBasicCell, FSimpleLeaf> tree(10, 3,loader.getBoxWidth(),loader.getCenterOfBox());

        // -----------------------------------------------------
        std::cout << "Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
        counter.tic();
        for(int idx = 0 ; idx < loader.getNumberOfParticles() ; ++idx){
            ParticleTsm part;
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
