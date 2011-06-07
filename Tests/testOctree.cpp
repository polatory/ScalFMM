// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../Src/Utils/FTic.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Utils/FAssertable.hpp"
#include "../Src/Utils/F3DPosition.hpp"

#include "../Src/Components/FBasicParticle.hpp"
#include "../Src/Components/FBasicCell.hpp"
#include "../Src/Components/FSimpleLeaf.hpp"

// Compile by : g++ testOctree.cpp ../Src/Utils/FAssertable.cpp -O2 -o testOctree.exe

/**
* In this file we show how to use octree
*/

int main(int , char ** ){    
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable is useless to execute.\n";
    std::cout << ">> It is only interesting to wath the code to understand\n";
    std::cout << ">> how to use the Octree\n";
    //////////////////////////////////////////////////////////////

    const long NbPart = 2000000;
    FTic counter;

    srand ( time(NULL) );

    FOctree<FBasicParticle, FBasicCell, FVector, FSimpleLeaf> tree(10, 3,1.0,F3DPosition(0.5,0.5,0.5));

    // -----------------------------------------------------
    std::cout << "Creating and inserting " << NbPart << " particles ..." << std::endl;
    counter.tic();
    for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
        FBasicParticle particle;
        particle.setPosition(FReal(rand())/RAND_MAX,FReal(rand())/RAND_MAX,FReal(rand())/RAND_MAX);

        tree.insert(particle);
    }
    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << ")." << std::endl;
    // -----------------------------------------------------


    return 0;
}


// [--LICENSE--]
