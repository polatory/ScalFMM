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
#include "../Sources/Core/FSimpleLeaf.hpp"

// Compile by : g++ testOctree.cpp ../Sources/Utils/FAssertable.cpp -O2 -o testOctree.exe

/**
* In this file we show how to use octree
*/

int main(int , char ** ){
        const long NbPart = 2000000;
        FList<FBasicParticule*> particules;
        FTic counter;

        srand ( time(NULL) );

        // -----------------------------------------------------
        std::cout << "Creating " << NbPart << " particules ..." << std::endl;
        counter.tic();
        for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
            FBasicParticule* const particule = new FBasicParticule();
            particule->setPosition(FReal(rand())/RAND_MAX,FReal(rand())/RAND_MAX,FReal(rand())/RAND_MAX);
            particules.pushFront(particule);
        }
        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << ")." << std::endl;
        // -----------------------------------------------------

        FOctree<FBasicParticule, FBasicCell, FSimpleLeaf, 10, 3> tree(1.0,F3DPosition(0.5,0.5,0.5));
        FList<FBasicParticule*>::BasicIterator iter(particules);

        // -----------------------------------------------------
        std::cout << "Inserting particules ..." << std::endl;
        counter.tic();
        while( iter.isValide() ){
            tree.insert(iter.value());
            iter.progress();
        }
        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << ")." << std::endl;
        // -----------------------------------------------------

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
