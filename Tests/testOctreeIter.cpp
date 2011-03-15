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

#include "../Sources/Core/FAbstractParticule.hpp"
#include "../Sources/Core/FAbstractCell.hpp"

// We use openmp to count time (portable and easy to manage)
// Compile by : g++ testOctreeIter.cpp ../Sources/Utils/FAssertable.cpp -O2 -o testOctreeIter.exe

/**
* In this file we show how to use octree with iteration

NbLevels = 5 & NbSubLevels = 2 & NbPart = 2000000
Creating particules ...
Done  (0.263944s).
Inserting particules ...
Done  (0.349972s).
Itering on particules ...
Next level (4096)...
Next level (512)...
Next level (64)...
Next level (8)...
Done  (5.742e-05s).
Deleting particules ...
Done  (0.075429s).

NbLevels = 10 & NbSubLevels = 3 & NbPart = 2000000
Creating particules ...
Done  (0.261834s).
Inserting particules ...
Done  (2.68187s).
Itering on particules ...
Next level (1985191)...
Next level (1885276)...
Next level (1289097)...
Next level (262019)...
Next level (32768)...
Next level (4096)...
Next level (512)...
Next level (64)...
Next level (8)...
Done  (0.716064s).
Deleting particules ...
Done  (0.0830964s).
*/

// Fake particule class
class TestParticule : public FAbstractParticule {
    F3DPosition pos;
public:
    TestParticule(const F3DPosition& inPos) : pos(inPos) {
    }

    F3DPosition getPosition() const {
            return pos;
    }
};

// Fake cell class
class TestCell : public FAbstractCell {
    MortonIndex index;
public:
    MortonIndex getMortonIndex() const {
        return index;
    }

    void setMortonIndex(const MortonIndex inIndex) {
        index = inIndex;
    }
};

int main(int , char ** ){
        const int NbLevels = 10;
        const int NbSubLevels = 3;
        const long NbPart = 2E6;
        FList<TestParticule*> particules;
        FTic counter;

        srand ( 1 ); // volontary set seed to constant

        // -----------------------------------------------------
        std::cout << "Creating " << NbPart << " particules ..." << std::endl;
        counter.tic();
        for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
            particules.pushFront(new TestParticule(F3DPosition(double(rand())/RAND_MAX,double(rand())/RAND_MAX,double(rand())/RAND_MAX)));
        }
        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;
        // -----------------------------------------------------

        FOctree<TestParticule, TestCell, NbLevels, NbSubLevels> tree(1.0,F3DPosition(0.5,0.5,0.5));
        FList<TestParticule*>::BasicIterator iter(particules);

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

        std::cout << "Itering on particules ..." << std::endl;
        counter.tic();

        FOctree<TestParticule, TestCell, NbLevels, NbSubLevels>::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        for(int idx = 0 ; idx < NbLevels - 1; ++idx ){
            int counter = 0;
            do{
                ++counter;
                //counter += octreeIterator.getCurrentList()->getSize();
            } while(octreeIterator.moveRight());
            octreeIterator.moveUp();
            octreeIterator.gotoLeft();
            std::cout << "Next level (" << counter << ")...\n";
        }
        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

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
