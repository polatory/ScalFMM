// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>


#include "../Sources/Containers/FOctree.hpp"
#include "../Sources/Containers/FList.hpp"

#include "../Sources/Utils/FAssertable.hpp"
#include "../Sources/Utils/F3DPosition.hpp"

#include "../Sources/Core/FBasicParticule.hpp"
#include "../Sources/Core/FBasicCell.hpp"

#include "../Sources/Utils/FTic.hpp"

// Compile by : g++ testOctreeIter.cpp ../Sources/Utils/FAssertable.cpp -O2 -o testOctreeIter.exe

/**
* In this file we show how to use octree with iteration
* This is a good example to understand FOctree::Iterator.
*/

int main(int , char ** ){
        const int NbLevels = 9;
        const int NbSubLevels = 3;
        const long NbPart = 2E6;
        FList<FBasicParticule*> particules;
        FTic counter;

        srand ( 1 ); // volontary set seed to constant

        // -----------------------------------------------------
        std::cout << "Creating " << NbPart << " particules ..." << std::endl;
        counter.tic();
        for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
            FBasicParticule* const particule = new FBasicParticule();
            particule->setPosition(FReal(rand())/RAND_MAX,FReal(rand())/RAND_MAX,FReal(rand())/RAND_MAX);
            particules.pushFront(particule);
        }
        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;
        // -----------------------------------------------------

        FOctree<FBasicParticule, FBasicCell, NbLevels, NbSubLevels> tree(1.0,F3DPosition(0.5,0.5,0.5));
        FList<FBasicParticule*>::BasicIterator iter(particules);

        // -----------------------------------------------------
        std::cout << "Inserting particules ..." << std::endl;
        counter.tic();
        while( iter.isValide() ){
            tree.insert(iter.value());
            iter.progress();
        }
        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;
        // -----------------------------------------------------
        {
            std::cout << "Itering on particules ..." << std::endl;
            counter.tic();

            FOctree<FBasicParticule, FBasicCell, NbLevels, NbSubLevels>::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();
            for(int idx = 0 ; idx < NbLevels - 1; ++idx ){
                int counter = 0;
                do{
                    ++counter;
                    //counter += octreeIterator.getCurrentList()->getSize();
                } while(octreeIterator.moveRight());
                octreeIterator.moveUp();
                octreeIterator.gotoLeft();
                std::cout << "Cells at this level " << counter << " ...\n";
            }
            counter.tac();
            std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;
        }
        // -----------------------------------------------------
        {
            std::cout << "Itering on particules fast ..." << std::endl;
            counter.tic();

            FOctree<FBasicParticule, FBasicCell, NbLevels, NbSubLevels>::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();

            FOctree<FBasicParticule, FBasicCell, NbLevels, NbSubLevels>::Iterator avoidGoLeft(octreeIterator);

            for(int idx = 0 ; idx < NbLevels - 1; ++idx ){
                int counter = 0;
                do{
                } while(octreeIterator.moveRight());
                avoidGoLeft.moveUp();
                octreeIterator = avoidGoLeft;
            }
            counter.tac();
            std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;
        }
        // -----------------------------------------------------
        std::cout << "Deleting particules ..." << std::endl;
        counter.tic();
        while(particules.getSize()){
            delete particules.popFront();
        }
        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;
        // -----------------------------------------------------

	return 0;
}


// [--LICENSE--]
