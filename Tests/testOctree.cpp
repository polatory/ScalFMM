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

// Compile by : g++ testOctree.cpp ../Sources/Utils/FAssertable.cpp -O2 -o testOctree.exe

/**
* In this file we show how to use octree
* Démarrage de /home/berenger/Dropbox/Personnel/FMB++/FMB++-build-desktop/FMB++...
* Creating particules ...
* Done  (0.268113).
* Inserting particules ...
* Done  (2.70836).
* Deleting particules ...
* Done  (0.0791715).
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
	
    void setPosition(const F3DPosition&){}
};


int main(int , char ** ){
        const long NbPart = 2000000;
        FList<TestParticule*> particules;
        FTic counter;

        srand ( time(NULL) );


        // -----------------------------------------------------
        std::cout << "Creating " << NbPart << " particules ..." << std::endl;
        counter.tic();
        for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
            particules.pushFront(new TestParticule(F3DPosition(double(rand())/RAND_MAX,double(rand())/RAND_MAX,double(rand())/RAND_MAX)));
        }
        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << ")." << std::endl;
        // -----------------------------------------------------

        FOctree<TestParticule, TestCell, 10, 3> tree(1.0,F3DPosition(0.5,0.5,0.5));
        FList<TestParticule*>::BasicIterator iter(particules);

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
