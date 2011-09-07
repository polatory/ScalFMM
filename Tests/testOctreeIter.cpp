// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>


#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"
#include "../Src/Components/FSimpleLeaf.hpp"

#include "../Src/Utils/FAssertable.hpp"
#include "../Src/Utils/F3DPosition.hpp"

#include "../Src/Components/FBasicParticle.hpp"
#include "../Src/Components/FBasicCell.hpp"

#include "../Src/Utils/FTic.hpp"

// Compile by : g++ testOctreeIter.cpp ../Src/Utils/FAssertable.cpp -O2 -o testOctreeIter.exe

/**
* In this file we show how to use octree with iteration
* This is a good example to understand FOctree::Iterator.
*/

int main(int , char ** ){
    typedef FVector<FBasicParticle>      ContainerClass;
    typedef FSimpleLeaf<FBasicParticle, ContainerClass >                     LeafClass;
    typedef FOctree<FBasicParticle, FBasicCell, ContainerClass , LeafClass >  OctreeClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable is useless to execute.\n";
    std::cout << ">> It is only interesting to wath the code to understand\n";
    std::cout << ">> how to use octree iterator.\n";
    //////////////////////////////////////////////////////////////

        const int NbLevels = 9;
        const int NbSubLevels = 3;
        const long NbPart = 2E6;
        const FReal FRandMax = FReal(RAND_MAX);

        FTic counterTime;

        srand ( 1 ); // volontary set seed to constant
        // -----------------------------------------------------

        OctreeClass tree(NbLevels, NbSubLevels, 1.0,F3DPosition(0.5,0.5,0.5));

        // -----------------------------------------------------
        std::cout << "Creating and inserting " << NbPart << " particles ..." << std::endl;
        counterTime.tic();
        {
            FBasicParticle particle;
            for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
                particle.setPosition(FReal(rand())/FRandMax,FReal(rand())/FRandMax,FReal(rand())/FRandMax);
                tree.insert(particle);
            }
        }
        counterTime.tac();
        std::cout << "Done  " << "(" << counterTime.elapsed() << "s)." << std::endl;

        // -----------------------------------------------------
        {
            std::cout << "Itering on particles ..." << std::endl;
            counterTime.tic();

            typename OctreeClass::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();
            for(int idxLevel = NbLevels - 1 ; idxLevel >= 1 ; --idxLevel ){
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

            typename OctreeClass::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();

            typename OctreeClass::Iterator avoidGoLeft(octreeIterator);

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

	return 0;
}


// [--LICENSE--]
