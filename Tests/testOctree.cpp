// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <omp.h>

#include "../Sources/Containers/FOctree.hpp"
#include "../Sources/Containers/FList.hpp"

#include "../Sources/Utils/FAssertable.hpp"
#include "../Sources/Utils/F3DPosition.hpp"

#include "../Sources/Core/FAbstractParticule.hpp"
#include "../Sources/Core/FAbstractCell.hpp"

// We use openmp to count time (portable and easy to manage)
// Compile by : g++ testOctree.cpp ../Sources/Utils/FAssertable.cpp -O2 -lgomp -fopenmp -o testOctree.exe

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
};


int main(int , char ** ){
        const long NbPart = 2000000;
        FList<TestParticule*> particules;

        srand ( time(NULL) );


        // -----------------------------------------------------
        std::cout << "Creating " << NbPart << " particules ..." << std::endl;
        const double CreatingStartTime = omp_get_wtime();
        for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
            particules.pushFront(new TestParticule(F3DPosition(double(rand())/RAND_MAX,double(rand())/RAND_MAX,double(rand())/RAND_MAX)));
        }
        const double CreatingEndTime = omp_get_wtime();
        std::cout << "Done  " << "(" << (CreatingEndTime-CreatingStartTime) << ")." << std::endl;
        // -----------------------------------------------------

        FOctree<TestParticule, TestCell, 10, 3> tree(1.0,F3DPosition(0.5,0.5,0.5));
        FList<TestParticule*>::BasicIterator iter(particules);

        // -----------------------------------------------------
        std::cout << "Inserting particules ..." << std::endl;
        const double InsertingStartTime = omp_get_wtime();
        while( iter.isValide() ){
            tree.insert(iter.value());
            iter.progress();
        }
        const double InsertingEndTime = omp_get_wtime();
        std::cout << "Done  " << "(" << (InsertingEndTime-InsertingStartTime) << ")." << std::endl;
        // -----------------------------------------------------

        // -----------------------------------------------------
        std::cout << "Deleting particules ..." << std::endl;
        const double DeletingStartTime = omp_get_wtime();
        while(particules.getSize()){
            delete particules.popFront();
        }
        const double DeletingEndTime = omp_get_wtime();
        std::cout << "Done  " << "(" << (DeletingEndTime-DeletingStartTime) << ")." << std::endl;
        // -----------------------------------------------------



	return 0;
}


// [--LICENSE--]
