// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>


#include "../Sources/Containers/FOctree.hpp"
#include "../Sources/Containers/FList.hpp"
#include "../Sources/Components/FSimpleLeaf.hpp"

#include "../Sources/Utils/FAssertable.hpp"
#include "../Sources/Utils/F3DPosition.hpp"

#include "../Sources/Components/FBasicParticle.hpp"
#include "../Sources/Components/FBasicCell.hpp"

#include "../Sources/Utils/FTic.hpp"

// Compile by : g++ testOctreeIter.cpp ../Sources/Utils/FAssertable.cpp -O2 -o testOctreeIter.exe

/**
* In this file we show how to use octree with iteration
* This is a good example to understand FOctree::Iterator.
*/

int main(int , char ** ){    
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable is useless to execute.\n";
    std::cout << ">> It is only interesting to wath the code to understand\n";
    std::cout << ">> how to use octree iterator.\n";
    //////////////////////////////////////////////////////////////

        const int NbLevels = 9;
        const int NbSubLevels = 3;
        const long NbPart = 2E6;
        FList<FBasicParticle*> particles;
        FTic counterTime;

        srand ( 1 ); // volontary set seed to constant

        // -----------------------------------------------------
        std::cout << "Creating " << NbPart << " particles ..." << std::endl;
        counterTime.tic();
        for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
            FBasicParticle* const particle = new FBasicParticle();
            particle->setPosition(FReal(rand())/RAND_MAX,FReal(rand())/RAND_MAX,FReal(rand())/RAND_MAX);
            particles.pushFront(particle);
        }
        counterTime.tac();
        std::cout << "Done  " << "(" << counterTime.elapsed() << "s)." << std::endl;
        // -----------------------------------------------------

        FOctree<FBasicParticle, FBasicCell, FSimpleLeaf, NbLevels, NbSubLevels> tree(1.0,F3DPosition(0.5,0.5,0.5));
        FList<FBasicParticle*>::BasicIterator iter(particles);

        // -----------------------------------------------------
        std::cout << "Inserting particles ..." << std::endl;
        counterTime.tic();
        while( iter.isValide() ){
            tree.insert(iter.value());
            iter.progress();
        }
        counterTime.tac();
        std::cout << "Done  " << "(" << counterTime.elapsed() << "s)." << std::endl;
        // -----------------------------------------------------
        {
            std::cout << "Itering on particles ..." << std::endl;
            counterTime.tic();

            FOctree<FBasicParticle, FBasicCell, FSimpleLeaf, NbLevels, NbSubLevels>::Iterator octreeIterator(&tree);
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
            counterTime.tac();
            std::cout << "Done  " << "(" << counterTime.elapsed() << "s)." << std::endl;
        }
        // -----------------------------------------------------
        {
            std::cout << "Itering on particles fast ..." << std::endl;
            counterTime.tic();

            FOctree<FBasicParticle, FBasicCell, FSimpleLeaf, NbLevels, NbSubLevels>::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();

            FOctree<FBasicParticle, FBasicCell, FSimpleLeaf, NbLevels, NbSubLevels>::Iterator avoidGoLeft(octreeIterator);

            for(int idx = 0 ; idx < NbLevels - 1; ++idx ){
                do{
                } while(octreeIterator.moveRight());
                avoidGoLeft.moveUp();
                octreeIterator = avoidGoLeft;
            }
            counterTime.tac();
            std::cout << "Done  " << "(" << counterTime.elapsed() << "s)." << std::endl;
        }
        // -----------------------------------------------------
        std::cout << "Deleting particles ..." << std::endl;
        counterTime.tic();
        while(particles.getSize()){
            delete particles.popFront();
        }
        counterTime.tac();
        std::cout << "Done  " << "(" << counterTime.elapsed() << "s)." << std::endl;
        // -----------------------------------------------------

	return 0;
}


// [--LICENSE--]
